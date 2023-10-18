// Copyright 2008, 2009, 2010, 2011, 2014 Martin C. Frith

// This class holds multiple sequences and their names.  The sequences
// are concatenated, with delimiters between them.

// The final sequence may be "unfinished".  This happens if the
// sequence data hits a memory limit before we finish reading it.

#ifndef MULTISEQUENCE_HH
#define MULTISEQUENCE_HH

#include "ScoreMatrixRow.hh"
#include "VectorOrMmap.hh"

#include <string>
#include <iosfwd>

namespace cbrc{

typedef unsigned char uchar;

class MultiSequence{
 public:

#if LAST_POS_BYTES > 4
  typedef size_t indexT;
#else
  typedef unsigned indexT;
#endif

  // initialize with leftmost delimiter pad, ready for appending sequences
  void initForAppending(indexT padSizeIn, bool isAppendStopSymbol = false);

  // re-initialize, but keep the last sequence if it is unfinished
  void reinitForAppending();

  void eraseAllButTheLastSequence() {
    ends.v.pop_back();
    reinitForAppending();
    ends.v.push_back(seq.v.size());
  }

  // read seqCount finished sequences, and their names, from binary files
  void fromFiles( const std::string& baseName, indexT seqCount,
                  size_t qualitiesPerLetter );

  // write all the finished sequences and their names to binary files
  void toFiles( const std::string& baseName ) const;

  // Append a sequence with delimiters.  Don't let the total size of
  // the concatenated sequences plus pads exceed maxSeqLen: thus it
  // may not finish reading the sequence.
  std::istream& appendFromFasta( std::istream& stream, indexT maxSeqLen );

  // As above, but read FASTQ format.
  std::istream& appendFromFastq(std::istream& stream, indexT maxSeqLen,
				bool isKeepQualityData);

  // As above, but read either FASTA or FASTQ format.  The first
  // sequence may have either format, but subsequent sequences must
  // have the same format.
  std::istream& appendFromFastx(std::istream& stream, indexT maxSeqLen,
				bool isKeepQualityData);

  // As above, but read quality scores too.
  std::istream& appendFromPrb( std::istream& stream, indexT maxSeqLen,
			       unsigned alphSize, const uchar decode[] );

  // As above, but read a PSSM too, in PSI-BLAST ASCII format.
  std::istream& appendFromPssm( std::istream& stream, indexT maxSeqLen,
                                const uchar* lettersToNumbers,
                                bool isMaskLowercase );

  // did we finish reading the last sequence?
  bool isFinished() const{ return ends.size() == nameEnds.size(); }

  // how many finished sequences are there?
  indexT finishedSequences() const{ return ends.size() - 1; }

  // total length of finished and unfinished sequences plus delimiters
  size_t unfinishedSize() const{ return seq.size(); }

  // which sequence is the coordinate in?
  indexT whichSequence( indexT coordinate ) const;

  indexT padBeg( indexT seqNum ) const{ return ends[seqNum] - padSize; }
  indexT seqBeg( indexT seqNum ) const{ return ends[seqNum]; }
  indexT seqEnd( indexT seqNum ) const{ return ends[seqNum+1] - padSize; }
  indexT padEnd( indexT seqNum ) const{ return ends[seqNum+1]; }
  indexT seqLen( indexT seqNum ) const{ return seqEnd(seqNum)-seqBeg(seqNum); }
  indexT padLen( indexT seqNum ) const{ return padEnd(seqNum)-padBeg(seqNum); }

  std::string seqName(indexT seqNum) const {
    const char *n = names.begin();
    indexT b = nameEnds[seqNum];
    indexT e = nameEnds[seqNum + 1];
    if (e > b && n[e-1] <= ' ') --e;
    return std::string(n + b, n + e);
  }

  char strand(indexT seqNum) const {
    const char *n = names.begin();
    indexT b = nameEnds[seqNum];
    indexT e = nameEnds[seqNum + 1];
    return (e > b && n[e-1] == '\t') ? '-' : '+';
  }

  // get a pointer to the start of the sequence data
  const uchar* seqReader() const{ return seq.begin(); }
  /***/ uchar* seqWriter()      { return &seq.v[0]; }

  // get a pointer to the start of the PSSM, or NULL if there is no PSSM
  // I am not totally sure about the reinterpret_cast...

  const ScoreMatrixRow* pssmReader() const{
    return pssm.empty() ? 0
        : reinterpret_cast< const ScoreMatrixRow* >( &pssm[0] );
  }

  // get a pointer to the start of the quality data
  const uchar* qualityReader() const{ return qualityScores.begin(); }

  // How many quality scores are there per letter?  There might be
  // none at all, one per letter, or several (e.g. 4) per letter.
  size_t qualsPerLetter() const { return qualityScoresPerLetter; }

  void reverseComplementOneSequence(indexT seqNum, const uchar *complement);

  void duplicateOneSequence(indexT seqNum);

 private:
  indexT padSize;  // number of delimiter chars between sequences
  VectorOrMmap<uchar> seq;  // concatenated sequences
  VectorOrMmap<indexT> ends;  // coordinates of ends of delimiter pads
  VectorOrMmap<char> names;  // concatenated sequence names (to save memory)
  VectorOrMmap<indexT> nameEnds;  // endpoints of the names

  std::vector<int> pssm;  // position-specific scoring matrix
  std::vector<uchar> pssmColumnLetters;  // which input column is which letter

  // The quality scores may be ASCII-coded: to get the real scores,
  // subtract e.g. 33 or 64.  The real scores might be related to
  // error probabilities in one of these ways:
  // Qphred = -10*log10(p)
  // Qsolexa = -10*log10(p/(1-p))
  VectorOrMmap<uchar> qualityScores;
  size_t qualityScoresPerLetter;
  bool isReadingFastq;
  bool isAppendingStopSymbol;

  // Read a fasta/fastq header: read the whole line but store just the
  // 1st word
  void readFastxName( std::istream& stream );

  // read the letters above PSSM columns, so we know which column is which
  std::istream& readPssmHeader( std::istream& stream );

  void finishName() {  // finish adding a sequence name: store its end coord
    nameEnds.v.push_back(names.v.size());
    if (nameEnds.v.back() < names.v.size()) {
      throw std::runtime_error("the sequence names are too long");
    }
  }

  void addName(const std::string &name) {  // add a new sequence name
    names.v.insert(names.v.end(), name.begin(), name.end());
    names.v.push_back('\n');
    finishName();
  }

  void appendQualPad() {  // add delimiter to the end of the quality scores
    uchar padQualityScore = 64;  // should never be used, but a valid value
    size_t s = padSize * qualityScoresPerLetter;
    qualityScores.v.insert(qualityScores.v.end(), s, padQualityScore);
  }

  void appendPssmPad() {  // add delimiter to the end of the PSSM
    pssm.insert(pssm.end(), padSize * scoreMatrixRowSize, -INF);
  }

  bool isRoomToAppendPad(indexT maxSeqLen) const {
    return seq.v.size() <= maxSeqLen && padSize <= maxSeqLen - seq.v.size();
  }

  // finish the last sequence: add final pad and end coordinate
  void finish() {
    seq.v.insert(seq.v.end(), padSize, ' ');
    ends.v.push_back(seq.v.size());
  }
};

// Divide the sequences into a given number of roughly-equally-sized
// chunks, and return the first sequence in the Nth chunk.
inline size_t firstSequenceInChunk(const MultiSequence &m,
				   size_t numOfChunks, size_t chunkNum) {
  size_t numOfSeqs = m.finishedSequences();
  size_t beg = m.seqBeg(0);
  size_t end = m.seqBeg(numOfSeqs) - 1;
  unsigned long long len = end - beg;  // try to avoid overflow
  size_t pos = beg + len * chunkNum / numOfChunks;
  size_t seqNum = m.whichSequence(pos);
  size_t begDistance = pos - m.seqBeg(seqNum);
  size_t endDistance = m.padEnd(seqNum) - pos;
  return (begDistance < endDistance) ? seqNum : seqNum + 1;
}

}

#endif
