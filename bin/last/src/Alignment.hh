// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

// This struct holds a gapped, pair-wise alignment.

#ifndef ALIGNMENT_HH
#define ALIGNMENT_HH

#include "Centroid.hh"
#include "GreedyXdropAligner.hh"
#include "SegmentPair.hh"
#include "mcf_frameshift_xdrop_aligner.hh"

#include <vector>
#include <cstring>

namespace cbrc{

class LastEvaluer;
class MultiSequence;
class Alphabet;
class TwoQualityScoreMatrix;

struct Aligners {
  Centroid centroid;
  FrameshiftXdropAligner frameshiftAligner;
  GreedyXdropAligner greedyAligner;
};

struct AlignmentText {
  // This holds the final text representation of an alignment, along
  // with data for sorting it.
  SegmentPair::indexT strandNum;
  SegmentPair::indexT queryBeg;
  SegmentPair::indexT queryEnd;
  double score;
  SegmentPair::indexT alnSize;
  SegmentPair::indexT matches;
  char *text;  // seems to be a bit faster than std::vector<char>

  AlignmentText() {}

  AlignmentText(size_t queryNumIn, size_t queryBegIn, size_t queryEndIn,
		char strandIn, double scoreIn,
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
  std::vector<char> columnAmbiguityCodes;
  std::vector<double> expectedCounts;  // expected emission & transition counts
  double fullScore;  // a.k.a. forward score, sum-of-paths score
  AlignmentExtras() : fullScore(0) {}
};

struct Alignment{
  // make a single-block alignment:
  void fromSegmentPair(const SegmentPair &sp) {
    blocks.assign(1, sp);
    score = sp.score;
  }

  // Make an Alignment by doing gapped X-drop extension in both
  // directions starting from a seed SegmentPair.  The resulting
  // Alignment might not be "optimal" (see below).
  // If outputType > 3: calculates match probabilities.
  // If outputType > 4: does gamma-centroid alignment.
  void makeXdrop( Aligners &aligners, bool isGreedy, bool isFullScore,
		  const uchar* seq1, const uchar* seq2, int globality,
		  const ScoreMatrixRow* scoreMatrix, int smMax, int smMin,
		  const const_dbl_ptr* probMatrix, double scale,
		  const GapCosts& gap, int maxDrop, size_t frameSize,
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
                  const GapCosts& gapCosts, size_t frameSize,
		  const ScoreMatrixRow* pssm2,
                  const TwoQualityScoreMatrix& sm2qual,
                  const uchar* qual1, const uchar* qual2 ) const;

  // Does the Alignment have any segment with score >= minScore?
  bool hasGoodSegment(const uchar *seq1, const uchar *seq2,
		      int minScore, const ScoreMatrixRow *scoreMatrix,
		      const GapCosts &gapCosts, size_t frameSize,
		      const ScoreMatrixRow *pssm2,
		      const TwoQualityScoreMatrix &sm2qual,
		      const uchar *qual1, const uchar *qual2) const;

  // translationType indicates that the 2nd sequence is: 0 = not
  // translated, 1 = translated into amino acids, 2 = translated into
  // codons (so codonToAmino is used to count matches & dnaAlph is
  // used to write the 2nd sequence).
  AlignmentText write(const MultiSequence& seq1, const MultiSequence& seq2,
		      size_t seqNum2, const uchar* seqData2,
		      const Alphabet& alph, const Alphabet& dnaAlph,
		      int translationType, const uchar *codonToAmino,
		      const LastEvaluer& evaluer, int format,
		      const AlignmentExtras& extras) const;

  // data:
  std::vector<SegmentPair> blocks;  // the gapless blocks of the alignment
  double score;
  SegmentPair seed;  // the alignment remembers its seed

  size_t beg1() const{ return blocks.front().beg1(); }
  size_t beg2() const{ return blocks.front().beg2(); }
  size_t end1() const{ return blocks.back().end1(); }
  size_t end2() const{ return blocks.back().end2(); }

  void extend( std::vector< SegmentPair >& chunks,
	       std::vector< char >& columnCodes,
	       Aligners &aligners, bool isGreedy, bool isFullScore,
	       const uchar* seq1, const uchar* seq2,
	       size_t start1, size_t start2,
	       bool isForward, int globality,
	       const ScoreMatrixRow* sm, int smMax, int smMin,
	       const const_dbl_ptr* probMat, double scale,
	       int maxDrop, const GapCosts& gap, size_t frameSize,
	       const ScoreMatrixRow* pssm2,
               const TwoQualityScoreMatrix& sm2qual,
               const uchar* qual1, const uchar* qual2,
	       const Alphabet& alph, AlignmentExtras& extras,
	       double gamma, int outputType );

  AlignmentText writeTab(const MultiSequence& seq1, const MultiSequence& seq2,
			 size_t seqNum2, bool isTranslated,
			 const LastEvaluer& evaluer,
			 const AlignmentExtras& extras) const;

  AlignmentText writeMaf(const MultiSequence& seq1, const MultiSequence& seq2,
			 size_t seqNum2, const uchar* seqData2,
			 const Alphabet& alph, const Alphabet& dnaAlph,
			 int translationType, const LastEvaluer& evaluer,
			 const AlignmentExtras& extras) const;

  AlignmentText writeBlastTab(const MultiSequence& seq1,
			      const MultiSequence& seq2,
			      size_t seqNum2, const uchar* seqData2,
			      const Alphabet& alph, int translationType,
			      const uchar *codonToAmino,
			      const LastEvaluer& evaluer,
			      const AlignmentExtras& extras,
			      bool isExtraColumns) const;

  size_t numColumns(size_t frameSize, bool isCodon) const;

  char *writeTopSeq(char *dest, const uchar *seq, const Alphabet &alph,
		    size_t qualsPerBase, size_t frameSize, bool isCodon) const;

  char *writeBotSeq(char *dest, const uchar *seq, const Alphabet &alph,
		    size_t qualsPerBase, size_t frameSize, bool isCodon) const;

  char *writeColumnProbs(char *dest, const char *probSymbols,
			 size_t frameSize, bool isCodon) const;
};

}  // end namespace cbrc
#endif  // ALIGNMENT_HH
