// Copyright 2009, 2010, 2011, 2013 Martin C. Frith

#include "MultiSequence.hh"
#include "stringify.hh"
#include <algorithm>  // max_element
#include <cctype>  // toupper
#include <limits>  // numeric_limits
#include <streambuf>

#define ERR(x) throw std::runtime_error(x)

using namespace cbrc;

std::istream&
MultiSequence::appendFromFastx(std::istream& stream, indexT maxSeqLen,
			       bool isKeepQualityData) {
  if (names.empty()) {
    isReadingFastq = false;
    char c = '>';
    stream >> c;
    if (c == '@') {
      isReadingFastq = true;
    } else if (c != '>') {
      ERR("bad sequence data: missing '>' or '@'");
    }
    readFastxName(stream);
    if (!stream) return stream;
  }

  return isReadingFastq ?
    appendFromFastq(stream, maxSeqLen, isKeepQualityData) :
    appendFromFasta(stream, maxSeqLen);
}

std::istream&
MultiSequence::appendFromFastq(std::istream& stream, indexT maxSeqLen,
			       bool isKeepQualityData) {
  // initForAppending:
  qualityScoresPerLetter = isKeepQualityData;
  if( qualityScores.v.empty() ) appendQualPad();

  if( isFinished() ){
    char c = '@';
    stream >> c;
    if( c != '@' ) ERR( "bad FASTQ data: missing '@'" );
    readFastxName(stream);
    if( !stream ) return stream;
  }

  std::streambuf *buf = stream.rdbuf();
  int c = buf->sgetc();

  while (c != std::streambuf::traits_type::eof()) {
    if (c > ' ') {
      if (c == '+' || seq.v.size() >= maxSeqLen) break;
      seq.v.push_back(c);
    }
    c = buf->snextc();
  }

  if (isRoomToAppendPad(maxSeqLen)) {
    do {
      c = buf->sbumpc();
    } while (c != std::streambuf::traits_type::eof() && c != '\n');

    for (size_t i = seq.v.size() - ends.v.back(); i > 0; ) {
      c = buf->sbumpc();
      if (c == std::streambuf::traits_type::eof()) ERR("bad FASTQ data");
      if (c > ' ') {
	if (isKeepQualityData) {
	  if (c > 126) ERR("non-printable-ASCII in FASTQ quality data");
	  qualityScores.v.push_back(c);
	}
	--i;
      }
    }

    finish();
    appendQualPad();
  }

  return stream;
}

std::istream&
MultiSequence::appendFromPrb( std::istream& stream, indexT maxSeqLen,
			      unsigned alphSize, const uchar decode[] ){
  // initForAppending:
  qualityScoresPerLetter = alphSize;
  if( qualityScores.v.empty() ) appendQualPad();

  if( isFinished() ){
    std::string line;
    getline( stream, line );  // slow but simple
    if( !stream ) return stream;

    // give the sequence a boring name:
    static size_t lineCount = 0;
    std::string name = stringify( ++lineCount );
    addName(name);

    size_t oldSize = qualityScores.v.size();

    std::istringstream iss(line);
    int q;
    while( iss >> q ){
      if( q < -64 || q > 62 )
	ERR( "quality score too large: " + stringify(q) );
      qualityScores.v.push_back( q + 64 );  // ASCII-encode the quality score
    }

    size_t newSize = qualityScores.v.size();
    if (newSize % qualityScoresPerLetter != 0) ERR("bad PRB data");

    for (size_t i = oldSize; i < newSize; i += qualityScoresPerLetter) {
      const uchar *q = &qualityScores.v[i];
      unsigned maxIndex = std::max_element(q, q + qualityScoresPerLetter) - q;
      seq.v.push_back( decode[ maxIndex ] );
    }
  }

  if (isRoomToAppendPad(maxSeqLen)) {
    finish();
    appendQualPad();
  }

  return stream;
}

std::istream& MultiSequence::readPssmHeader( std::istream& stream ){
  // read the name of the sequence/PSSM:
  std::string name;
  stream >> name;
  stream.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

  // Look for a line with one letter per column of the PSSM:
  pssmColumnLetters.clear();
  std::string line, word;

  while( getline( stream, line ) ){
    std::istringstream iss(line);

    while( iss >> word ){
      if( word.size() == 1 ){
	uchar c = word[0];
        uchar letter = std::toupper(c);
        // allow for PSI-BLAST format, with repeated letters:
        if( pssmColumnLetters.size() && pssmColumnLetters[0] == letter ) break;
        pssmColumnLetters.push_back(letter);
      }else{
        pssmColumnLetters.clear();
        break;
      }
    }

    if( pssmColumnLetters.size() ) break;
  }

  if( !stream ) return stream;
  addName(name);
  return stream;
}

std::istream&
MultiSequence::appendFromPssm( std::istream& stream, indexT maxSeqLen,
                               const uchar* lettersToNumbers,
                               bool isMaskLowercase ){
  // initForAppending:
  if( pssm.empty() ) appendPssmPad();

  if( isFinished() ){
    readPssmHeader(stream);
    if( !stream ) return stream;
  }

  while( seq.v.size() < maxSeqLen ){
    unsigned position;
    uchar letter;
    int score;
    std::vector<int> scores;
    stream >> position >> letter;
    while( scores.size() < pssmColumnLetters.size() && stream >> score ){
      scores.push_back(score);
    }
    if( !stream ) break;

    seq.v.push_back(letter);

    int minScore = *std::min_element( scores.begin(), scores.end() );
    pssm.insert( pssm.end(), scoreMatrixRowSize, minScore );
    std::vector<int>::iterator row = pssm.end() - scoreMatrixRowSize;
    for( unsigned i = 0; i < scores.size(); ++i ){
      uchar columnLetter = pssmColumnLetters[i];
      unsigned column = lettersToNumbers[columnLetter];
      if( column >= scoreMatrixRowSize )
        ERR( std::string("bad column-letter in PSSM: ") + char(columnLetter) );
      row[column] = scores[i];
      unsigned maskColumn = lettersToNumbers[ std::tolower(columnLetter) ];
      if( maskColumn >= scoreMatrixRowSize ) continue;  // ?
      if( isMaskLowercase ) scores[i] = std::min(scores[i], 0);
      if( maskColumn != column ) row[maskColumn] = scores[i];
    }
    uchar delimiter = ' ';
    unsigned delimiterColumn = lettersToNumbers[delimiter];
    assert( delimiterColumn < scoreMatrixRowSize );
    row[delimiterColumn] = -INF;

    stream.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
  }

  if (isRoomToAppendPad(maxSeqLen)) {
    finish();
    appendPssmPad();
  }

  if( !stream.bad() ) stream.clear();

  return stream;
}
