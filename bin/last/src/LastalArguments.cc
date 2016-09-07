// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

#include "LastalArguments.hh"
#include "stringify.hh"
#include <unistd.h>  // getopt
#include <algorithm>  // max
#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <cctype>
#include <cmath>  // log
#include <cstring>  // strtok
#include <cstdlib>  // EXIT_SUCCESS

#define ERR(x) throw std::runtime_error(x)

static void badopt( char opt, const char* arg ){
  ERR( std::string("bad option value: -") + opt + ' ' + arg );
}

static int myGetopt( int argc, char** argv, const char* optstring ){
  if( optind < argc ){
    std::string nextarg = argv[optind];
    if( nextarg == "--help"    ) return 'h';
    if( nextarg == "--version" ) return 'V';
  }
  return getopt( argc, argv, optstring );
}

static char parseOutputFormat( const char* text ){
  std::string s = text;
  for( size_t i = 0; i < s.size(); ++i ){
    s[i] = std::tolower( s[i] );
  }
  if( s == "tab" || s == "0" ) return 't';
  if( s == "maf" || s == "1" ) return 'm';
  if( s == "blasttab" )        return 'b';
  if( s == "blasttab+" )       return 'B';
  return 0;
}

namespace cbrc{

LastalArguments::LastalArguments() :
  outputFormat('m'),
  outputType(3),
  strand(-1),  // depends on the alphabet
  isQueryStrandMatrix(false),
  globality(0),
  isKeepLowercase(true),  // depends on the option used with lastdb
  tantanSetting(-1),  // depends on the option used with lastdb
  maskLowercase(-1),  // depends on the lowercase option used with lastdb
  maxEvalue(-1),
  queryLettersPerRandomAlignment(1e6),
  minScoreGapped(-1),  // depends on the alphabet
  minScoreGapless(-1),  // depends on minScoreGapped and the outputType
  matchScore(-1),  // depends on the alphabet
  mismatchCost(-1),  // depends on the alphabet
  gapExistCost(-1),  // depends on the alphabet
  gapExtendCost(-1),  // depends on the alphabet
  insExistCost(-1),  // depends on gapExistCost
  insExtendCost(-1),  // depends on gapExtendCost
  gapPairCost(-1),  // this means: OFF
  frameshiftCost(-1),  // this means: ordinary, non-translated alignment
  matrixFile(""),
  maxDropGapped(-1),  // depends on minScoreGapped & maxDropGapless
  maxDropGapless(-1),  // depends on the score matrix
  maxDropFinal(-1),  // depends on maxDropGapped
  inputFormat(sequenceFormat::fasta),
  minHitDepth(1),
  maxHitDepth(-1),
  oneHitMultiplicity(10),
  maxGaplessAlignmentsPerQueryPosition(0),  // depends on oneHitMultiplicity
  cullingLimitForGaplessAlignments(0),
  cullingLimitForFinalAlignments(0),
  queryStep(1),
  batchSize(0),  // depends on the outputType, and voluming
  numOfThreads(1),
  maxRepeatDistance(1000),  // sufficiently conservative?
  temperature(-1),  // depends on the score matrix
  gamma(1),
  geneticCodeFile(""),
  verbosity(0){}

void LastalArguments::fromArgs( int argc, char** argv, bool optionsOnly ){
  std::string usage =
      "Usage: lastal [options] lastdb-name fasta-sequence-file(s)";

  std::string help = usage + "\n\
Find local sequence alignments.\n\
\n\
Score options (default settings):\n\
-r: match score   (DNA: 1, 0<Q<5:  6)\n\
-q: mismatch cost (DNA: 1, 0<Q<5: 18)\n\
-p: match/mismatch score matrix (protein-protein: BL62, DNA-protein: BL80)\n\
-a: gap existence cost (DNA: 7, protein: 11, 0<Q<5: 21)\n\
-b: gap extension cost (DNA: 1, protein:  2, 0<Q<5:  9)\n\
-A: insertion existence cost (a)\n\
-B: insertion extension cost (b)\n\
-c: unaligned residue pair cost (off)\n\
-F: frameshift cost (off)\n\
-x: maximum score drop for gapped alignments (max[y, e-1])\n\
-y: maximum score drop for gapless alignments (t*10)\n\
-z: maximum score drop for final gapped alignments (x)\n\
-d: minimum score for gapless alignments (min[e, t*ln(1000*refSize/n)])\n\
-e: minimum score for gapped alignments\n\
\n\
E-value options (default settings):\n\
-D: query letters per random alignment ("
    + stringify(queryLettersPerRandomAlignment) + ")\n\
-E: maximum expected alignments per square giga (1e+18/D/refSize/numOfStrands)\n\
\n\
Cosmetic options (default settings):\n\
-h, --help: show all options and their default settings, and exit\n\
-V, --version: show version information, and exit\n\
-v: be verbose: write messages about what lastal is doing\n\
-f: output format: TAB, MAF, BlastTab, BlastTab+ (MAF)\n\
\n\
Miscellaneous options (default settings):\n\
-s: strand: 0=reverse, 1=forward, 2=both (2 for DNA, 1 for protein)\n\
-S: score matrix applies to forward strand of: 0=reference, 1=query ("
    + stringify(isQueryStrandMatrix) + ")\n\
-T: type of alignment: 0=local, 1=overlap ("
    + stringify(globality) + ")\n\
-m: maximum initial matches per query position ("
    + stringify(oneHitMultiplicity) + ")\n\
-l: minimum length for initial matches ("
    + stringify(minHitDepth) + ")\n\
-L: maximum length for initial matches (infinity)\n\
-n: maximum gapless alignments per query position (infinity if m=0, else m)\n\
-C: omit gapless alignments in >= C others with > score-per-length (off)\n\
-K: omit alignments whose query range lies in >= K others with > score (off)\n\
-k: step-size along the query sequence ("
    + stringify(queryStep) + ")\n\
-i: query batch size (8 KiB, unless there is > 1 thread or lastdb volume)\n\
-P: number of parallel threads ("
    + stringify(numOfThreads) + ")\n\
-R: repeat-marking options (the same as was used for lastdb)\n\
-u: mask lowercase during extensions: 0=never, 1=gapless,\n\
    2=gapless+gapped but not final, 3=always (2 if lastdb -c and Q<5, else 0)\n\
-w: suppress repeats inside exact matches, offset by <= this distance ("
    + stringify(maxRepeatDistance) + ")\n\
-G: genetic code file\n\
-t: 'temperature' for calculating probabilities (1/lambda)\n\
-g: 'gamma' parameter for gamma-centroid and LAMA ("
    + stringify(gamma) + ")\n\
-j: output type: 0=match counts, 1=gapless, 2=redundant gapped, 3=gapped,\n\
                 4=column ambiguity estimates, 5=gamma-centroid, 6=LAMA,\n\
                 7=expected counts ("
    + stringify(outputType) + ")\n\
-Q: input format: 0=fasta, 1=fastq-sanger, 2=fastq-solexa, 3=fastq-illumina,\n\
                  4=prb, 5=PSSM ("
    + stringify(inputFormat) + ")\n\
\n\
Report bugs to: last-align (ATmark) googlegroups (dot) com\n\
LAST home page: http://last.cbrc.jp/\n\
";

  optind = 1;  // allows us to scan arguments more than once(???)
  int c;
  const char optionString[] = "hVvf:" "r:q:p:a:b:A:B:c:F:x:y:z:d:e:" "D:E:"
    "s:S:T:m:l:L:n:C:K:k:i:P:R:u:w:t:g:G:j:Q:";
  while( (c = myGetopt(argc, argv, optionString)) != -1 ){
    switch(c){
    case 'h':
      std::cout << help;
      throw EXIT_SUCCESS;
    case 'V':
      std::cout << "lastal "
#include "version.hh"
	"\n";
      throw EXIT_SUCCESS;
    case 'v':
      ++verbosity;
      break;
    case 'f':
      outputFormat = parseOutputFormat( optarg );
      if( !outputFormat ) badopt( c, optarg );
      break;

    case 'r':
      unstringify( matchScore, optarg );
      if( matchScore <= 0 ) badopt( c, optarg );
      break;
    case 'q':
      unstringify( mismatchCost, optarg );
      if( mismatchCost < 0 ) badopt( c, optarg );  // allow 0 for Fujibuchi-san
      break;
    case 'p':
      matrixFile = optarg;
      break;
    case 'a':
      unstringify( gapExistCost, optarg );
      if( gapExistCost < 0 ) badopt( c, optarg );
      break;
    case 'b':
      unstringify( gapExtendCost, optarg );
      if( gapExtendCost <= 0 ) badopt( c, optarg );
      break;
    case 'A':
      unstringify( insExistCost, optarg );
      if( insExistCost < 0 ) badopt( c, optarg );
      break;
    case 'B':
      unstringify( insExtendCost, optarg );
      if( insExtendCost <= 0 ) badopt( c, optarg );
      break;
    case 'c':
      unstringify( gapPairCost, optarg );
      if( gapPairCost <= 0 ) badopt( c, optarg );
      break;
    case 'F':
      unstringify( frameshiftCost, optarg );
      if( frameshiftCost <= 0 ) badopt( c, optarg );
      break;
    case 'x':
      unstringify( maxDropGapped, optarg );
      if( maxDropGapped < 0 ) badopt( c, optarg );
      break;
    case 'y':
      unstringify( maxDropGapless, optarg );
      if( maxDropGapless < 0 ) badopt( c, optarg );
      break;
    case 'z':
      unstringify( maxDropFinal, optarg );
      if( maxDropFinal < 0 ) badopt( c, optarg );
      break;
    case 'd':
      unstringify( minScoreGapless, optarg );
      if( minScoreGapless < 0 ) badopt( c, optarg );
      break;
    case 'e':
      unstringify( minScoreGapped, optarg );
      if( minScoreGapped < 0 ) badopt( c, optarg );
      break;

    case 'D':
      unstringify( queryLettersPerRandomAlignment, optarg );
      if( queryLettersPerRandomAlignment <= 0 ) badopt( c, optarg );
      break;
    case 'E':
      unstringify( maxEvalue, optarg );
      if( maxEvalue <= 0 ) badopt( c, optarg );
      break;

    case 's':
      unstringify( strand, optarg );
      if( strand < 0 || strand > 2 ) badopt( c, optarg );
      break;
    case 'S':
      unstringify( isQueryStrandMatrix, optarg );
      break;
    case 'T':
      unstringify( globality, optarg );
      if( globality < 0 || globality > 1 ) badopt( c, optarg );
      break;
    case 'm':
      unstringify( oneHitMultiplicity, optarg );
      break;
    case 'l':
      unstringify( minHitDepth, optarg );
      break;
    case 'L':
      unstringify( maxHitDepth, optarg );
      break;
    case 'n':
      unstringify( maxGaplessAlignmentsPerQueryPosition, optarg );
      if( maxGaplessAlignmentsPerQueryPosition <= 0 ) badopt( c, optarg );
      break;
    case 'C':
      unstringify( cullingLimitForGaplessAlignments, optarg );
      break;
    case 'K':
      unstringify( cullingLimitForFinalAlignments, optarg );
      break;
    case 'k':
      unstringify( queryStep, optarg );
      if( queryStep <= 0 ) badopt( c, optarg );
      break;
    case 'i':
      unstringifySize( batchSize, optarg );
      if( batchSize <= 0 ) badopt( c, optarg );  // 0 means "not specified"
      break;
    case 'P':
      unstringify( numOfThreads, optarg );
      break;
    case 'R':
      if( optarg[0] < '0' || optarg[0] > '1' ) badopt( c, optarg );
      if( optarg[1] < '0' || optarg[1] > '2' ) badopt( c, optarg );
      if( optarg[2] ) badopt( c, optarg );
      isKeepLowercase = optarg[0] - '0';
      tantanSetting = optarg[1] - '0';
      break;
    case 'u':
      unstringify( maskLowercase, optarg );
      if( maskLowercase < 0 || maskLowercase > 3 ) badopt( c, optarg );
      break;
    case 'w':
      unstringify( maxRepeatDistance, optarg );
      break;
    case 't':
      unstringify( temperature, optarg );
      if( temperature <= 0 ) badopt( c, optarg );
      break;
    case 'g':
      unstringify( gamma, optarg );
      if( gamma <= 0 ) badopt( c, optarg );
      break;
    case 'G':
      geneticCodeFile = optarg;
      break;
    case 'j':
      unstringify( outputType, optarg );
      if( outputType < 0 || outputType > 7 ) badopt( c, optarg );
      break;
    case 'Q':
      unstringify( inputFormat, optarg );
      break;

    case '?':
      ERR( "bad option" );
    }
  }

  if( maskLowercase == 1 && inputFormat == 5 )
    ERR( "can't combine option -u 1 with option -Q 5" );

  if( maskLowercase == 2 && inputFormat == 5 )
    ERR( "can't combine option -u 2 with option -Q 5" );

  if( isTranslated() && inputFormat == 5 )
    ERR( "can't combine option -F with option -Q 5" );

  if( isTranslated() && outputType > 3 )
    ERR( "can't combine option -F with option -j > 3" );

  if( isTranslated() && globality == 1 )
    ERR( "can't combine option -F with option -T 1" );

  if( isTranslated() && isQueryStrandMatrix )
    ERR( "can't combine option -F with option -S 1" );

  if( globality == 1 && outputType == 1 )
    ERR( "can't combine option -T 1 with option -j 1" );

  if( optionsOnly ) return;
  if( optind >= argc )
    ERR( "please give me a database name and sequence file(s)\n\n" + usage );
  lastdbName = argv[optind++];
  inputStart = optind;
}

void LastalArguments::fromLine( const std::string& line ){
  const char* delimiters = " \t";
  const char* s = line.c_str();
  std::vector<char> args( s, s + line.size() + 1 );
  std::vector<char*> argv;
  char* i = std::strtok( &args[0], delimiters );
  argv.push_back(i);
  while( i ){
    i = std::strtok( 0, delimiters );
    argv.push_back(i);
  }
  fromArgs( argv.size() - 1, &argv[0], true );
}

void LastalArguments::fromString( const std::string& s ){
  std::string trigger = "#last";
  std::istringstream iss(s);
  std::string line;
  while( std::getline( iss, line ) )
    if( line.compare( 0, trigger.size(), trigger ) == 0 )
      fromLine( line );
}

const char* LastalArguments::matrixName( bool isProtein ) const{
  if( matrixFile.empty() && matchScore < 0 && mismatchCost < 0 && isProtein )
    return isTranslated() ? "BL80" : "BL62";
  return matrixFile.c_str();
}

void LastalArguments::setDefaultsFromAlphabet( bool isDna, bool isProtein,
					       double numLettersInReference,
					       bool isKeepRefLowercase,
					       int refTantanSetting,
                                               bool isCaseSensitiveSeeds,
					       bool isVolumes,
					       unsigned realNumOfThreads ){
  if( strand < 0 ) strand = (isDna || isTranslated()) ? 2 : 1;

  if( isProtein ){
    // default match & mismatch scores: Blosum62 matrix
    if( matchScore < 0 && mismatchCost >= 0 ) matchScore   = 1;  // idiot-proof
    if( mismatchCost < 0 && matchScore >= 0 ) mismatchCost = 1;  // idiot-proof
    if( gapExistCost   < 0 ) gapExistCost   =  11;
    if( gapExtendCost  < 0 ) gapExtendCost  =   2;
  }
  else if( !isQuality( inputFormat ) ){
    if( matchScore     < 0 ) matchScore     =   1;
    if( mismatchCost   < 0 ) mismatchCost   =   1;
    if( gapExistCost   < 0 ) gapExistCost   =   7;
    if( gapExtendCost  < 0 ) gapExtendCost  =   1;
  }
  else{  // sequence quality scores will be used:
    if( matchScore     < 0 ) matchScore     =   6;
    if( mismatchCost   < 0 ) mismatchCost   =  18;
    if( gapExistCost   < 0 ) gapExistCost   =  21;
    if( gapExtendCost  < 0 ) gapExtendCost  =   9;
    // With this scoring scheme for DNA, gapless lambda ~= ln(10)/10,
    // so these scores should be comparable to PHRED scores.
    // Furthermore, since mismatchCost/matchScore = 3, the target
    // distribution of paired bases ~= 99% identity.  Because the
    // quality scores are unlikely to be perfect, it may be best to
    // use a lower target %identity than we otherwise would.
  }

  if( maxEvalue < 0 ){
    double r = std::max( numLettersInReference, 1.0 );  // avoid divide-by-0
    maxEvalue = 1e18 / (numOfStrands() * r * queryLettersPerRandomAlignment);
  }

  if( insExistCost < 0 ) insExistCost = gapExistCost;
  if( insExtendCost < 0 ) insExtendCost = gapExtendCost;

  if( tantanSetting < 0 ){
    isKeepLowercase = isKeepRefLowercase;
    tantanSetting = refTantanSetting;
  }

  if( maskLowercase < 0 ){
    if( isCaseSensitiveSeeds && inputFormat != sequenceFormat::pssm )
      maskLowercase = 2;
    else
      maskLowercase = 0;
  }

  if( batchSize == 0 ){
    // With voluming, we want the batches to be as large as will
    // comfortably fit into memory, because each volume gets read from
    // disk once per batch.
    if( !isVolumes && realNumOfThreads == 1 )
      batchSize = 0x2000;  // 8 Kbytes (?)
    else if( inputFormat == sequenceFormat::pssm )
      batchSize = 0x100000;   // 1 Mbyte
    else if( outputType == 0 )
      batchSize = 0x1000000;  // 16 Mbytes
    else if( inputFormat == sequenceFormat::prb )
      batchSize = 0x2000000;  // 32 Mbytes (?)
    else
      batchSize = 0x8000000;  // 128 Mbytes
    // (should we reduce the 128 Mbytes, for fewer out-of-memory errors?)
  }

  if( maxGaplessAlignmentsPerQueryPosition == 0 )
    maxGaplessAlignmentsPerQueryPosition =
      (oneHitMultiplicity > 0) ? oneHitMultiplicity : -1;

  if( isTranslated() && frameshiftCost < gapExtendCost )
    ERR( "the frameshift cost must not be less than the gap extension cost" );

  if( insExistCost != gapExistCost || insExtendCost != gapExtendCost ){
    if( isTranslated() )
      ERR( "can't combine option -F with option -A or -B" );
  }
}

void LastalArguments::setDefaultsFromMatrix( double lambda, int minScore ){
  if( outputType < 2 && minScoreGapped < 0 ) minScoreGapped = minScoreGapless;
  if( minScoreGapped < 0 ){
    if( outputType > 0 && minScore < 0 )
      ERR("\
can't calculate E-values: maybe the mismatch or gap costs are too weak.\n\
To proceed without E-values, set a score threshold with option -e.");
    minScoreGapped = minScore;
  }
  if( outputType < 2 ) minScoreGapless = minScoreGapped;

  if( temperature < 0 ) temperature = 1 / lambda;

  if( maxDropGapless < 0 ){  // should it depend on temperature or lambda?
    if( temperature < 0 ) maxDropGapless = 0;  // shouldn't happen
    else                  maxDropGapless = int( 10.0 * temperature + 0.5 );
  }

  if( maxDropGapped < 0 ){
    maxDropGapped = std::max( minScoreGapped - 1, maxDropGapless );
  }

  if( maxDropFinal < 0 ) maxDropFinal = maxDropGapped;
}

int LastalArguments::calcMinScoreGapless( double numLettersInReference,
					  double numOfIndexes ) const{
  if( minScoreGapless >= 0 ) return minScoreGapless;

  // ***** Default setting for minScoreGapless *****

  // This attempts to ensure that the gapped alignment phase will be
  // reasonably fast relative to the gapless alignment phase.

  // The expected number of gapped extensions per query position is:
  // kGapless * referenceSize * exp(-lambdaGapless * minScoreGapless).

  // The number of gapless extensions per query position is
  // proportional to: maxGaplessAlignmentsPerQueryPosition * numOfIndexes.

  // So we want exp(lambdaGapless * minScoreGapless) to be
  // proportional to: kGapless * referenceSize / (n * numOfIndexes).

  // But we crudely ignore kGapless.

  // The proportionality constant was guesstimated by some limited
  // trial-and-error.  It should depend on the relative speeds of
  // gapless and gapped extensions.

  if( temperature < 0 ) return minScoreGapped;  // shouldn't happen

  double n = maxGaplessAlignmentsPerQueryPosition;
  if( maxGaplessAlignmentsPerQueryPosition + 1 == 0 ) n = 10;  // ?
  double x = 1000.0 * numLettersInReference / (n * numOfIndexes);
  if( x < 1 ) x = 1;
  int s = int( temperature * std::log(x) + 0.5 );
  return std::min( s, minScoreGapped );
}

void LastalArguments::writeCommented( std::ostream& stream ) const{
  stream << '#';
  stream << " a=" << gapExistCost;
  stream << " b=" << gapExtendCost;
  stream << " A=" << insExistCost;
  stream << " B=" << insExtendCost;
  if( gapPairCost > 0 )
    stream << " c=" << gapPairCost;
  if( isTranslated() )
    stream << " F=" << frameshiftCost;
  stream << " e=" << minScoreGapped;
  stream << " d=" << minScoreGapless;
  stream << " x=" << maxDropGapped;
  stream << " y=" << maxDropGapless;
  stream << " z=" << maxDropFinal;
  stream << " D=" << queryLettersPerRandomAlignment;
  stream << " E=" << maxEvalue;
  stream << '\n';

  stream << '#';
  stream << " R=" << isKeepLowercase << tantanSetting;
  stream << " u=" << maskLowercase;
  stream << " s=" << strand;
  stream << " S=" << isQueryStrandMatrix;
  stream << " T=" << globality;
  stream << " m=" << oneHitMultiplicity;
  stream << " l=" << minHitDepth;
  if( maxHitDepth + 1 > 0 )
    stream << " L=" << maxHitDepth;
  stream << " n=" << maxGaplessAlignmentsPerQueryPosition;
  if( cullingLimitForGaplessAlignments )
    stream << " C=" << cullingLimitForGaplessAlignments;
  if( cullingLimitForFinalAlignments )
    stream << " K=" << cullingLimitForFinalAlignments;
  stream << " k=" << queryStep;
  stream << " i=" << batchSize;
  stream << " w=" << maxRepeatDistance;
  stream << " t=" << temperature;
  if( outputType > 4 && outputType < 7 )
    stream << " g=" << gamma;
  stream << " j=" << outputType;
  stream << " Q=" << inputFormat;
  stream << '\n';

  stream << "# " << lastdbName << '\n';
}

}  // end namespace cbrc
