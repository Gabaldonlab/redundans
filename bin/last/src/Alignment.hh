// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

// This struct holds a gapped, pair-wise alignment.

#ifndef ALIGNMENT_HH
#define ALIGNMENT_HH
#include "ScoreMatrixRow.hh"
#include "SegmentPair.hh"
#include <stddef.h>  // size_t
#include <string>
#include <vector>
#include <cstring>

namespace cbrc{

typedef unsigned char uchar;

class GappedXdropAligner;
class GeneralizedAffineGapCosts;
class LastEvaluer;
class MultiSequence;
class Alphabet;
class Centroid;
class TwoQualityScoreMatrix;

struct AlignmentText {
  // This holds the final text representation of an alignment, along
  // with data for sorting it.
  SegmentPair::indexT strandNum;
  SegmentPair::indexT queryBeg;
  SegmentPair::indexT queryEnd;
  int score;
  SegmentPair::indexT alnSize;
  SegmentPair::indexT matches;
  char *text;  // seems to be a bit faster than std::vector<char>

  AlignmentText() {}

  AlignmentText(size_t queryNumIn, size_t queryBegIn, size_t queryEndIn,
		char strandIn, int scoreIn,
		size_t alnSizeIn, size_t matchesIn, char *textIn) :
    strandNum(queryNumIn * 2 + (strandIn == '-')),
    queryBeg(queryBegIn), queryEnd(queryEndIn), score(scoreIn),
    alnSize(alnSizeIn), matches(matchesIn), text(textIn) {}

  size_t queryNum() const { return strandNum / 2; }

  bool operator<(const AlignmentText& r) const {
    // Order by query number (ascending), then score (descending):
    if (queryNum() != r.queryNum()) return queryNum() < r.queryNum();
    if (score      != r.score     ) return score      > r.score;

    // Requested by JGI:
    if (alnSize    != r.alnSize   ) return alnSize    > r.alnSize;
    if (matches    != r.matches   ) return matches    > r.matches;

    // Break ties, to make the sort order exactly reproducible:
    return std::strcmp(text, r.text) < 0;
  }
};

struct AlignmentExtras {
  // Optional (probabilistic) attributes of an alignment.
  // To save memory, these are outside the main Alignment struct.
  std::vector<uchar> columnAmbiguityCodes;  // char or uchar?
  std::vector<double> expectedCounts;  // expected emission & transition counts
  double fullScore;  // a.k.a. forward score, sum-of-paths score
  AlignmentExtras() : fullScore(0) {}
};

struct Alignment{
  // make a single-block alignment:
  void fromSegmentPair( const SegmentPair& sp );

  // Make an Alignment by doing gapped X-drop extension in both
  // directions starting from a seed SegmentPair.  The resulting
  // Alignment might not be "optimal" (see below).
  // If outputType > 3: calculates match probabilities.
  // If outputType > 4: does gamma-centroid alignment.
  void makeXdrop( GappedXdropAligner& aligner, Centroid& centroid,
		  const uchar* seq1, const uchar* seq2, int globality,
		  const ScoreMatrixRow* scoreMatrix, int smMax,
		  const GeneralizedAffineGapCosts& gap, int maxDrop,
		  int frameshiftCost, size_t frameSize,
		  const ScoreMatrixRow* pssm2,
                  const TwoQualityScoreMatrix& sm2qual,
                  const uchar* qual1, const uchar* qual2,
		  const Alphabet& alph, AlignmentExtras& extras,
		  double gamma = 0, int outputType = 0 );

  // Check that the Alignment has no prefix with score <= 0, no suffix
  // with score <= 0, and no sub-segment with score < -maxDrop.
  // Alignments that pass this test may be non-optimal in other ways.
  // If "globality" is non-zero, skip the prefix and suffix checks.
  bool isOptimal( const uchar* seq1, const uchar* seq2, int globality,
                  const ScoreMatrixRow* scoreMatrix, int maxDrop,
                  const GeneralizedAffineGapCosts& gap,
		  int frameshiftCost, size_t frameSize,
		  const ScoreMatrixRow* pssm2,
                  const TwoQualityScoreMatrix& sm2qual,
                  const uchar* qual1, const uchar* qual2 );

  AlignmentText write(const MultiSequence& seq1, const MultiSequence& seq2,
		      size_t seqNum2, char strand, const uchar* seqData2,
		      bool isTranslated, const Alphabet& alph,
		      const LastEvaluer& evaluer, int format,
		      const AlignmentExtras& extras) const;

  // data:
  std::vector<SegmentPair> blocks;  // the gapless blocks of the alignment
  int score;
  SegmentPair seed;  // the alignment remembers its seed

  size_t beg1() const{ return blocks.front().beg1(); }
  size_t beg2() const{ return blocks.front().beg2(); }
  size_t end1() const{ return blocks.back().end1(); }
  size_t end2() const{ return blocks.back().end2(); }

  void extend( std::vector< SegmentPair >& chunks,
	       std::vector< uchar >& ambiguityCodes,
	       GappedXdropAligner& aligner, Centroid& centroid,
	       const uchar* seq1, const uchar* seq2,
	       size_t start1, size_t start2,
	       bool isForward, int globality,
	       const ScoreMatrixRow* sm, int smMax, int maxDrop,
	       const GeneralizedAffineGapCosts& gap,
	       int frameshiftCost, size_t frameSize,
	       const ScoreMatrixRow* pssm2,
               const TwoQualityScoreMatrix& sm2qual,
               const uchar* qual1, const uchar* qual2,
	       const Alphabet& alph, AlignmentExtras& extras,
	       double gamma, int outputType );

  AlignmentText writeTab(const MultiSequence& seq1, const MultiSequence& seq2,
			 size_t seqNum2, char strand, bool isTranslated,
			 const LastEvaluer& evaluer,
			 const AlignmentExtras& extras) const;

  AlignmentText writeMaf(const MultiSequence& seq1, const MultiSequence& seq2,
			 size_t seqNum2, char strand, const uchar* seqData2,
			 bool isTranslated, const Alphabet& alph,
			 const LastEvaluer& evaluer,
			 const AlignmentExtras& extras) const;

  AlignmentText writeBlastTab(const MultiSequence& seq1,
			      const MultiSequence& seq2,
			      size_t seqNum2, char strand,
			      const uchar* seqData2,
			      bool isTranslated, const Alphabet& alph,
			      const LastEvaluer& evaluer,
			      bool isExtraColumns) const;

  size_t numColumns( size_t frameSize ) const;

  char* writeTopSeq( const uchar* seq, const Alphabet& alph,
		     size_t qualsPerBase, size_t frameSize, char* dest ) const;

  char* writeBotSeq( const uchar* seq, const Alphabet& alph,
		     size_t qualsPerBase, size_t frameSize, char* dest ) const;
};

}  // end namespace cbrc
#endif  // ALIGNMENT_HH
