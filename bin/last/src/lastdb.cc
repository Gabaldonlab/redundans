// Copyright 2008, 2009, 2010, 2011, 2013, 2014 Martin C. Frith

// Read fasta-format sequences; construct a suffix array of them; and
// write the results to files.

#include "last.hh"

#include "LastdbArguments.hh"
#include "TantanMasker.hh"
#include "zio.hh"
#include "stringify.hh"
#include "threadUtil.hh"
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <numeric>  // accumulate

#define ERR(x) throw std::runtime_error(x)
#define LOG(x) if( args.verbosity > 0 ) std::cerr << args.programName << ": " << x << '\n'

using namespace cbrc;

typedef unsigned long long countT;

static void openOrDie(std::ifstream &file, const std::string &name) {
  file.open(name.c_str());
  if (!file) ERR("can't open file: " + name);
}

// Set up an alphabet (e.g. DNA or protein), based on the user options
void makeAlphabet( Alphabet& alph, const LastdbArguments& args ){
  if( !args.userAlphabet.empty() )  alph.fromString( args.userAlphabet );
  else if( args.isAddStops )        alph.fromString( alph.proteinWithStop );
  else if( args.isProtein )         alph.fromString( alph.protein );
  else                              alph.fromString( alph.dna );
}

// Does the first sequence look like it isn't really DNA?
bool isDubiousDna( const Alphabet& alph, const MultiSequence& multi ){
  const uchar* seq = multi.seqReader() + multi.seqBeg(0);
  unsigned dnaCount = 0;

  for( unsigned i = 0; i < 100; ++i ){  // look at the first 100 letters
    uchar c = alph.numbersToUppercase[ seq[i] ];
    if( c == alph.size ) return false;  // we hit the end of the sequence early
    if( c < alph.size || c == alph.encode[ (uchar)'N' ] ) ++dnaCount;
  }

  if( dnaCount < 90 ) return true;  // more than 10% unexpected letters
  else return false;
}

// Set up the seed pattern(s)
static void makeSubsetSeeds( std::vector< CyclicSubsetSeed >& seeds,
			     const std::string& seedText,
			     const LastdbArguments& args,
			     const Alphabet& alph ){
  const std::string& a = alph.letters;
  bool isCaseSens = args.isCaseSensitive;

  if( !args.subsetSeedFile.empty() ){
    CyclicSubsetSeed::addPatterns(seeds, seedText, isCaseSens, alph.encode, a);
  }
  else if (!args.dnaSeedPatterns.empty()) {
    for (size_t x = 0; x < args.dnaSeedPatterns.size(); ++x) {
      const std::string &p = args.dnaSeedPatterns[x];
      std::string s = CyclicSubsetSeed::stringFromDnaPatterns(p);
      CyclicSubsetSeed::addPatterns(seeds, s, isCaseSens, alph.encode, a);
    }
  }
  else if( !args.seedPatterns.empty() ){
    for( unsigned x = 0; x < args.seedPatterns.size(); ++x ){
      const std::string& p = args.seedPatterns[x];
      std::string s = CyclicSubsetSeed::stringFromPatterns( p, a );
      CyclicSubsetSeed::addPatterns(seeds, s, isCaseSens, alph.encode, a);
    }
  }
  else{
    std::string s =
      (alph.letters == alph.dna) ? CyclicSubsetSeed::stringFromName("YASS") :
      args.isAddStops ? CyclicSubsetSeed::stringFromName("PSEUDO") :
      CyclicSubsetSeed::stringFromPatterns("1", a);
    CyclicSubsetSeed::addPatterns(seeds, s, isCaseSens, alph.encode, a);
  }

  if( seeds.empty() ) ERR( "no seed patterns" );
}

void writeLastalOptions( std::ostream& out, const std::string& seedText ){
  std::string trigger = "#lastal";
  std::istringstream iss( seedText );
  std::string line;
  while( getline( iss, line ) )
    if( line.compare( 0, trigger.size(), trigger ) == 0 )
      out << line << '\n';
}

void writePrjFile( const std::string& fileName, const LastdbArguments& args,
		   const Alphabet& alph, countT sequenceCount,
		   size_t maxSeqLen, const std::vector<countT>& letterCounts,
		   bool isFastq, unsigned volumes, unsigned numOfIndexes,
		   const std::string& seedText ){
  countT letterTotal = std::accumulate( letterCounts.begin(),
                                        letterCounts.end(), countT(0) );

  std::ofstream f( fileName.c_str() );
  f << "version=" <<
#include "version.hh"
    << '\n';
  f << "alphabet=" << alph << '\n';
  f << "numofsequences=" << sequenceCount << '\n';
  f << "numofletters=" << letterTotal << '\n';
  f << "maxsequenceletters=" << maxSeqLen << '\n';
  f << "letterfreqs=";
  for( unsigned i = 0; i < letterCounts.size(); ++i ){
    if( i > 0 ) f << ' ';
    f << letterCounts[i];
  }
  f << '\n';

  if( !args.isCountsOnly ){
    f << "maxunsortedinterval=" << args.minSeedLimit << '\n';
    f << "keeplowercase=" << args.isKeepLowercase << '\n';
    if( args.tantanSetting ){
      f << "tantansetting=" << args.tantanSetting << '\n';
    }
    f << "masklowercase=" << args.isCaseSensitive << '\n';
    if( isFastq ){
      f << "sequenceformat="
	<< (args.inputFormat % sequenceFormat::fastxKeep) << '\n';
    }
    if( args.minimizerWindow > 1 ){
      // Maybe this should be written (and read) by the indexes, so
      // each index can have a different window?
      f << "minimizerwindow=" << args.minimizerWindow << '\n';
    }
    if( volumes+1 > 0 ){
      f << "volumes=" << volumes << '\n';
    }
    else{
      f << "numofindexes=" << numOfIndexes << '\n';
    }
    f << "integersize=" << (posSize * CHAR_BIT) << '\n';
    writeLastalOptions( f, seedText );
  }

  f.close();
  if( !f ) ERR( "can't write file: " + fileName );
}

static void preprocessSomeSeqs(MultiSequence *multi,
			       const TantanMasker *masker,
			       const uchar *maskTable,
			       size_t numOfChunks,
			       size_t chunkNum) {
  size_t beg = firstSequenceInChunk(*multi, numOfChunks, chunkNum);
  size_t end = firstSequenceInChunk(*multi, numOfChunks, chunkNum + 1);
  uchar *w = multi->seqWriter();
  for (size_t i = beg; i < end; ++i)
    masker->mask(w + multi->seqBeg(i), w + multi->seqEnd(i), maskTable);
}

static void preprocessSeqs(MultiSequence &multi,
			   const TantanMasker &masker,
			   const uchar *maskTable,
			   size_t numOfChunks) {
#ifdef HAS_CXX_THREADS
  std::vector<std::thread> threads(numOfChunks - 1);
  for (size_t i = 1; i < numOfChunks; ++i)
    threads[i - 1] = std::thread(preprocessSomeSeqs,
				 &multi, &masker, maskTable, numOfChunks, i);
#endif
  preprocessSomeSeqs(&multi, &masker, maskTable, numOfChunks, 0);
#ifdef HAS_CXX_THREADS
  for (size_t i = 1; i < numOfChunks; ++i)
    threads[i - 1].join();
#endif
}

// Make one database volume, from one batch of sequences
void makeVolume(std::vector<CyclicSubsetSeed>& seeds,
		const DnaWordsFinder& wordsFinder, MultiSequence& multi,
		const LastdbArguments& args, const Alphabet& alph,
		std::vector<countT>& letterCountsSeen, size_t& maxSeqLenSeen,
		const TantanMasker& masker, unsigned numOfThreads,
		const std::string& seedText, const std::string& baseName) {
  size_t numOfIndexes = wordsFinder.wordLength ? 1 : seeds.size();
  size_t numOfSequences = multi.finishedSequences();
  size_t textLength = multi.seqBeg(numOfSequences);
  const uchar* seq = multi.seqReader();

  std::vector<countT> letterCounts(alph.size);
  size_t maxSeqLen = 0;
  size_t letterTotal = 0;
  for (size_t i = 0; i < numOfSequences; ++i) {
    alph.count(seq + multi.seqBeg(i), seq + multi.seqEnd(i), &letterCounts[0]);
    size_t t = accumulate(letterCounts.begin(), letterCounts.end(), countT(0));
    maxSeqLen = std::max(maxSeqLen, t - letterTotal);
    letterTotal = t;
  }

  for (unsigned c = 0; c < alph.size; ++c) {
    letterCountsSeen[c] += letterCounts[c];
  }
  maxSeqLenSeen = std::max(maxSeqLenSeen, maxSeqLen);

  if (args.isCountsOnly) return;

  if( args.tantanSetting ){
    LOG( "masking..." );
    preprocessSeqs( multi, masker, alph.numbersToLowercase, numOfThreads );
  }

  LOG( "writing..." );
  writePrjFile( baseName + ".prj", args, alph, numOfSequences,
		maxSeqLen, letterCounts,
		multi.qualsPerLetter(), -1, numOfIndexes, seedText );
  multi.toFiles( baseName );

  for( unsigned x = 0; x < numOfIndexes; ++x ){
    SubsetSuffixArray myIndex;
    std::vector<CyclicSubsetSeed> &indexSeeds = myIndex.getSeeds();
    size_t wordCounts[dnaWordsFinderNull + 1] = {0};

    if (wordsFinder.wordLength) {
      const uchar *seqEnd = seq + textLength;
      LOG("counting...");
      wordsFinder.count(seq, seqEnd, wordCounts);
      std::partial_sum(wordCounts, wordCounts + seeds.size(), wordCounts);
      LOG("gathering...");
      seeds.swap(indexSeeds);
      myIndex.setWordPositions(wordsFinder, wordCounts, seq, seqEnd);
    } else {
      indexSeeds.resize(1);
      seeds[x].swap(indexSeeds[0]);
      if (args.minimizerWindow > 1) {
	LOG("gathering...");
	for (size_t i = 0; i < numOfSequences; ++i) {
	  myIndex.addMinimizerPositions(seq, multi.seqBeg(i), multi.seqEnd(i),
					args.indexStep, args.minimizerWindow);
	}
      } else {
	const uchar *subsetMap = indexSeeds[0].firstMap();
	size_t count = 0;
	LOG("counting...");
	for (size_t i = 0; i < numOfSequences; ++i) {
	  size_t end = multi.seqEnd(i);
	  for (size_t j = multi.seqBeg(i); ; j += args.indexStep) {
	    count += (subsetMap[seq[j]] < CyclicSubsetSeed::DELIMITER);
	    if (end - j <= args.indexStep) break;
	  }
	}
	LOG("gathering...");
	PosPart *positions = myIndex.resizedPositions(count);
	for (size_t i = 0; i < numOfSequences; ++i) {
	  size_t end = multi.seqEnd(i);
	  for (size_t j = multi.seqBeg(i); ; j += args.indexStep) {
	    if (subsetMap[seq[j]] < CyclicSubsetSeed::DELIMITER) {
	      posSet(positions, j);
	      positions += posParts;
	    }
	    if (end - j <= args.indexStep) break;
	  }
	}
      }
      wordCounts[0] = myIndex.size();
    }

    LOG( "sorting..." );
    myIndex.sortIndex(seq, wordsFinder.wordLength, wordCounts,
		      args.minSeedLimit, args.childTableType, numOfThreads);

    LOG( "bucketing..." );
    myIndex.makeBuckets(seq, wordsFinder.wordLength, wordCounts,
			args.minIndexedPositionsPerBucket, args.bucketDepth,
			numOfThreads);

    LOG( "writing..." );
    if( numOfIndexes > 1 ){
      myIndex.toFiles( baseName + char('a' + x), false, textLength );
    }
    else{
      myIndex.toFiles( baseName, true, textLength );
    }

    if (wordsFinder.wordLength) {
      seeds.swap(indexSeeds);
    } else {
      seeds[x].swap(indexSeeds[0]);
    }
  }

  LOG( "done!" );
}

// The max number of sequence letters, such that the total volume size
// is likely to be less than volumeSize bytes.  (This is crude, it
// neglects memory for the sequence names, and the fact that
// lowercase-masked letters and DNA "N"s aren't indexed.)
static indexT maxLettersPerVolume( const LastdbArguments& args,
				   const DnaWordsFinder& wordsFinder,
				   size_t qualityCodesPerLetter,
				   unsigned numOfSeeds ){
  double b = args.minIndexedPositionsPerBucket;
  double x = posSize + offSize / b;  // bytes per indexed position
  if (wordsFinder.wordLength) {
    x = x * wordsFinder.numOfMatchedWords / wordsFinder.wordLookup.size();
  } else {
    x *= numOfSeeds;
    if (args.minimizerWindow > 1) {
      x = x * 2 / (args.minimizerWindow + 1);
    } else {
      x /= args.indexStep;
    }
  }
  double y = args.volumeSize / (1 + qualityCodesPerLetter + x);
  return (y < posLimit) ? y : posLimit;
}

static bool isRoomToDuplicateTheLastSequence(const MultiSequence &multi,
					     size_t maxSeqLen) {
  size_t n = multi.finishedSequences();
  size_t s = multi.seqBeg(n);
  return s <= maxSeqLen && s - multi.seqBeg(n - 1) <= maxSeqLen - s;
}

static void dump1(const std::string &dbName, const uchar *decode,
		  size_t seqCount, bool isFastq) {
  if (seqCount + 1 == 0) ERR("can't read file: " + dbName + ".prj");
  MultiSequence m;
  m.fromFiles(dbName, seqCount, isFastq);
  const uchar *s = m.seqReader();
  for (size_t i = 0; i < m.finishedSequences(); ++i) {
    std::cout << ">@"[isFastq] << m.seqName(i) << '\n';
    std::streambuf *buf = std::cout.rdbuf();
    size_t e = m.seqEnd(i);
    for (size_t j = m.seqBeg(i); j < e; ++j) {
      buf->sputc(decode[s[j]]);
    }
    if (isFastq) {
      std::cout << "\n+\n";
      std::cout.write((char *)m.qualityReader() + m.seqBeg(i), m.seqLen(i));
    }
    std::cout << '\n';
  }
}

static void dump(const std::string &dbName) {
  std::ios_base::sync_with_stdio(false);  // makes it much faster!
  Alphabet alph;
  unsigned volumes = -1;
  size_t seqCount = -1;
  size_t bitsPerInt = 4 * CHAR_BIT;
  sequenceFormat::Enum fmt = sequenceFormat::fasta;
  std::string line, word;
  std::ifstream file;
  openOrDie(file, dbName + ".prj");
  while (getline(file, line)) {
    std::istringstream iss(line);
    getline(iss, word, '=');
    if (word == "alphabet") iss >> alph;
    if (word == "numofsequences") iss >> seqCount;
    if (word == "sequenceformat") iss >> fmt;
    if (word == "volumes") iss >> volumes;
    if (word == "integersize") iss >> bitsPerInt;
  }
  if (alph.letters.empty()) ERR("can't read file: " + dbName + ".prj");
  size_t b = bitsPerInt / CHAR_BIT;
  if (posSize > 4 && b <= 4) ERR("please use lastdb for " + dbName);
  if (posSize <= 4 && b > 4) ERR("please use lastdb5 for " + dbName);
  if (volumes + 1 == 0) {
    dump1(dbName, alph.decode, seqCount, fmt != sequenceFormat::fasta);
  } else {
    for (unsigned i = 0; i < volumes; ++i) {
      std::string volName = dbName + stringify(i);
      std::ifstream f;
      openOrDie(f, volName + ".prj");
      seqCount = -1;
      while (getline(f, line)) {
	std::istringstream iss(line);
	getline(iss, word, '=');
	if (word == "numofsequences") iss >> seqCount;
      }
      dump1(volName, alph.decode, seqCount, fmt != sequenceFormat::fasta);
    }
  }
}

void lastdb( int argc, char** argv ){
  LastdbArguments args;
  args.fromArgs( argc, argv );
  if (args.isDump) return dump(args.lastdbName);

  std::string seedText;
  if( !args.subsetSeedFile.empty() ){
    seedText = CyclicSubsetSeed::stringFromName( args.subsetSeedFile );
    args.resetCumulativeOptions();
    args.fromString( seedText );  // read options from the seed file
    args.fromArgs( argc, argv );  // command line overrides seed file
  }

  args.setDefaults();

  unsigned numOfThreads =
    decideNumberOfThreads(args.numOfThreads, args.programName, args.verbosity);
  Alphabet alph;
  makeAlphabet( alph, args );
  TantanMasker tantanMasker;
  if( args.tantanSetting )
    tantanMasker.init(alph.isProtein(), args.tantanSetting == 2,
		      args.tantanSetting == 3, alph.letters, alph.encode);
  std::vector< CyclicSubsetSeed > seeds;
  makeSubsetSeeds( seeds, seedText, args, alph );

  DnaWordsFinder wordsFinder;
  makeWordsFinder(wordsFinder, &seeds[0], seeds.size(), alph.encode,
		  args.isCaseSensitive);
  if (wordsFinder.wordLength && alph.isProtein())
    err("error: word-restricted DNA seeds on protein");
  LOG("wordLength=" << wordsFinder.wordLength);

  MultiSequence multi;
  initSequences(multi, alph, false, args.isAddStops);
  unsigned volumeNumber = 0;
  countT sequenceCount = 0;
  std::vector<countT> letterCounts( alph.size );
  indexT maxLetters = 0;
  indexT maxSeqLen = posLimit;
  size_t maxSeqLenSeen = 0;

  char defaultInputName[] = "-";
  char* defaultInput[] = { defaultInputName, 0 };
  char** inputBegin = argv + args.inputStart;

  for( char** i = *inputBegin ? inputBegin : defaultInput; *i; ++i ){
    mcf::izstream inFileStream;
    std::istream& in = openIn( *i, inFileStream );
    LOG( "reading " << *i << "..." );

    while (appendSequence(multi, in, maxSeqLen, args.inputFormat, alph, 0)) {
      if (multi.isFinished()) {
	encodeSequences(multi, args.inputFormat, alph, args.isKeepLowercase,
			multi.finishedSequences() - 1);
	if (sequenceCount == 0) {
	  maxLetters = maxLettersPerVolume(args, wordsFinder,
					   multi.qualsPerLetter(),
					   seeds.size());
	  if (!args.isProtein && !args.isAddStops &&
	      args.userAlphabet.empty() && isDubiousDna(alph, multi)) {
	    std::cerr << args.programName
		      << ": that's some funny-lookin DNA\n";
	  }
	}
	maxSeqLen = maxLetters;
	if (args.strand != 1) {
	  if (args.strand == 2) {
	    ++sequenceCount;
	    if (isRoomToDuplicateTheLastSequence(multi, maxSeqLen)) {
	      size_t lastSeq = multi.finishedSequences() - 1;
	      multi.duplicateOneSequence(lastSeq);
	    } else {
	      std::string baseName =
		args.lastdbName + stringify(volumeNumber++);
	      makeVolume(seeds, wordsFinder, multi, args, alph, letterCounts,
			 maxSeqLenSeen, tantanMasker, numOfThreads, seedText,
			 baseName);
	      multi.eraseAllButTheLastSequence();
	    }
	  }
	  size_t lastSeq = multi.finishedSequences() - 1;
	  multi.reverseComplementOneSequence(lastSeq, alph.complement);
	}
        ++sequenceCount;
      } else {
	if (multi.finishedSequences() == 0) throwSeqTooBig();
	std::string baseName = args.lastdbName + stringify(volumeNumber++);
	makeVolume(seeds, wordsFinder, multi, args, alph, letterCounts,
		   maxSeqLenSeen, tantanMasker, numOfThreads, seedText,
		   baseName);
	multi.reinitForAppending();
	maxSeqLen = posLimit;
      }
    }
  }

  if( multi.finishedSequences() > 0 ){
    if( volumeNumber == 0 && !args.isCountsOnly ){
      makeVolume(seeds, wordsFinder, multi, args, alph, letterCounts,
		 maxSeqLenSeen, tantanMasker, numOfThreads, seedText,
		 args.lastdbName);
      return;
    }
    std::string baseName = args.lastdbName + stringify(volumeNumber++);
    makeVolume(seeds, wordsFinder, multi, args, alph, letterCounts,
	       maxSeqLenSeen, tantanMasker, numOfThreads, seedText, baseName);
  }

  writePrjFile( args.lastdbName + ".prj", args, alph, sequenceCount,
		maxSeqLenSeen, letterCounts, multi.qualsPerLetter(),
		volumeNumber, seeds.size(), seedText );
}

int main( int argc, char** argv )
try{
  lastdb( argc, argv );
  return EXIT_SUCCESS;
}
catch( const std::bad_alloc& e ) {  // bad_alloc::what() may be unfriendly
  std::cerr << argv[0] << ": out of memory\n";
  return EXIT_FAILURE;
}
catch( const std::exception& e ) {
  std::cerr << argv[0] << ": " << e.what() << '\n';
  return EXIT_FAILURE;
}
catch( int i ) {
  return i;
}
