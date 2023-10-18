// Copyright 2011, 2012, 2013 Martin C. Frith

// These routines extend an alignment in a given direction (forward or
// reverse) from given start points in two sequences.

// The start points point at the first positions we'll try to align.

// To use: first call "align", which calculates the alignment but only
// returns its score.  To get the actual alignment, call
// "getNextChunk" to get each gapless chunk.

// The sequences had better end with sentinels: the score for matching
// a sentinel with anything should be -INF.

// The gap parameters correspond to "generalized affine gap costs"
// (Proteins 1998 32(1):88-96).
// gapExistenceCost = a
// gapExtensionCost = b
// gapUnalignedCost = c
// When c >= a + 2b, it reduces to standard affine gap costs:
// gap cost = gapExistenceCost + gapExtensionCost * (gap length).

// The insertion and deletion costs may differ.  Typically:
// delExistenceCost = insExistenceCost, and
// delExtensionCost = insExtensionCost.

// The algorithm proceeds antidiagonal-by-antidiagonal, similarly to
// Section 2 in J Comput Biol. 2000 7(1-2):203-14.  It does not allow
// the score to drop by more than maxScoreDrop below the highest score
// in any previous antidiagonal.

// If "globality" is 0, local alignment is performed: this finds the
// highest-scoring alignment ending anywhere.  Otherwise, overlap
// alignment is performed: this seeks the highest-scoring alignment
// ending at the end of either sequence.  If overlap alignment reaches
// the end of neither sequence (due to maxScoreDrop), -INF is
// returned.

// The parameter maxMatchScore should be the highest possible score
// for matching 2 letters.  This parameter is not actually necessary,
// but it provides some optimization opportunities.  If you give it a
// too-high value, the results will not change, but the run time may
// increase.

#ifndef GAPPED_XDROP_ALIGNER_HH
#define GAPPED_XDROP_ALIGNER_HH

#include "mcf_contiguous_queue.hh"
#include "mcf_reverse_queue.hh"
#include "mcf_gap_costs.hh"
#include "mcf_simd.hh"
#include "ScoreMatrixRow.hh"

#include <iosfwd>
#include <stddef.h>  // size_t
#include <vector>

namespace cbrc {

using namespace mcf;

typedef unsigned char uchar;
typedef const int *const_int_ptr;

typedef int Score;
typedef uchar TinyScore;

class TwoQualityScoreMatrix;

const int xdropPadLen = simdBytes;

const int droppedTinyScore = UCHAR_MAX;

class GappedXdropAligner {
 public:
  int align(const uchar *seq1,  // start point in the 1st sequence
            const uchar *seq2,  // start point in the 2nd sequence
            bool isForward,  // forward or reverse extension?
	    int globality,
            const ScoreMatrixRow *scorer,  // the substitution score matrix
	    int delExistenceCost,
	    int delExtensionCost,
	    int insExistenceCost,
	    int insExtensionCost,
            int gapUnalignedCost,
	    bool isAffine,
            int maxScoreDrop,
            int maxMatchScore);

  // Like "align", but it aligns a sequence to a PSSM.
  int alignPssm(const uchar *seq,
                const ScoreMatrixRow *pssm,
                bool isForward,
		int globality,
		int delExistenceCost,
		int delExtensionCost,
		int insExistenceCost,
		int insExtensionCost,
                int gapUnalignedCost,
		bool isAffine,
                int maxScoreDrop,
                int maxMatchScore);

  // Like "align", but both sequences have quality scores.
  int align2qual(const uchar *seq1,
                 const uchar *qual1,
                 const uchar *seq2,
                 const uchar *qual2,
                 bool isForward,
		 int globality,
                 const TwoQualityScoreMatrix &scorer,
		 int delExistenceCost,
		 int delExtensionCost,
		 int insExistenceCost,
		 int insExtensionCost,
                 int gapUnalignedCost,
		 bool isAffine,
                 int maxScoreDrop,
                 int maxMatchScore);

  // Like "align", but maybe faster for DNA.  Assumes that
  // scorer[i<4][j<4] fits in signed char.  Each sequence element is
  // first mapped through "toUnmasked".  If an unmasked sequence
  // element >= 4 appears, alignDna falls back to a slower algorithm.
  int alignDna(const uchar *seq1,
	       const uchar *seq2,
	       bool isForward,
	       const ScoreMatrixRow *scorer,
	       int delExistenceCost,
	       int delExtensionCost,
	       int insExistenceCost,
	       int insExtensionCost,
	       int maxScoreDrop,
	       int maxMatchScore,
	       const uchar *toUnmasked);

  // Call this repeatedly to get each gapless chunk of the alignment.
  // The chunks are returned in far-to-near order.  The chunk's end
  // coordinates in each sequence (relative to the start of extension)
  // and length are returned in the first 3 parameters.  If there are
  // no more chunks, the 3 parameters are unchanged and "false" is
  // returned.
  bool getNextChunk(size_t &end1,
                    size_t &end2,
                    size_t &length,
		    int delExistenceCost,
		    int delExtensionCost,
		    int insExistenceCost,
		    int insExtensionCost,
                    int gapUnalignedCost);

  // After "alignDna", must use this instead of "getNextChunk"
  bool getNextChunkDna(size_t &end1,
		       size_t &end2,
		       size_t &length,
		       int delExistenceCost,
		       int delExtensionCost,
		       int insExistenceCost,
		       int insExtensionCost);

  // Like "align", but it aligns a protein sequence to a DNA sequence.
  // The DNA should be provided as 3 protein sequences, one for each
  // reading frame.  seq2frame0 is the in-frame start point.
  // seq2frame1 is shifted by 1 in the direction of alignment
  // extension.  seq2frame2 is shifted by 1 in the opposite direction.
  // The algorithm is "3-frame alignment", as described in J Comput
  // Biol. 1997 4(3):339-49.
  int align3(const uchar *seq1,
             const uchar *seq2frame0,
             const uchar *seq2frame1,  // the +1 frame
             const uchar *seq2frame2,  // the -1 frame
             bool isForward,
             const ScoreMatrixRow *scorer,
             int gapExistenceCost,
             int gapExtensionCost,
             int gapUnalignedCost,
             int frameshiftCost,
             int maxScoreDrop,
             int maxMatchScore);

  // Like "align3", but it aligns a protein sequence to a 3-frame PSSM.
  int align3pssm(const uchar *seq,
                 const ScoreMatrixRow *pssmFrame0,
                 const ScoreMatrixRow *pssmFrame1,
                 const ScoreMatrixRow *pssmFrame2,
                 bool isForward,
                 int gapExistenceCost,
                 int gapExtensionCost,
                 int gapUnalignedCost,
                 int frameshiftCost,
                 int maxScoreDrop,
                 int maxMatchScore);

  // Use this version of getNextChunk for protein-versus-DNA
  // alignment.  The end1 parameter receives a protein coordinate, and
  // end2 receives a DNA coordinate relative to the in-frame start.
  // The length parameter receives the number of residues/codons, not
  // bases.
  bool getNextChunk3(size_t &end1,
                     size_t &end2,
                     size_t &length,
                     int gapExistenceCost,
                     int gapExtensionCost,
                     int gapUnalignedCost,
                     int frameshiftCost);

  // Like "align3", but does "new style" DNA-protein alignment.
  int alignFrame(const uchar *protein,
		 const uchar *frame0,  // translated DNA,  0 frame
		 const uchar *frame1,  // translated DNA, +1 frame
		 const uchar *frame2,  // translated DNA, +2 frame
		 bool isForward,
		 const ScoreMatrixRow *scorer,
		 const GapCosts &gapCosts,
		 int maxScoreDrop);

  // Use this after alignFrame.  The cost of the unaligned region
  // between this chunk and the next chunk is returned in the 4th
  // parameter.
  bool getNextChunkFrame(size_t &end1,
			 size_t &end2,
			 size_t &length,
			 int &costOfNearerGap,
			 const GapCosts &gapCosts);

  void writeShape(std::ostream &out) const;

  // The next few functions are for use by Centroid.  If the Centroid
  // code gets updated, it might make sense to change these functions too.

  // The number of antidiagonals, excluding dummy ones at the beginning.
  size_t numAntidiagonals() const
  { return numOfAntidiagonals; }

  size_t numCellsAndPads(size_t antidiagonal) const
  { return scoreEndsAndOrigins[2 * (antidiagonal + 3)]
      -    scoreEndsAndOrigins[2 * (antidiagonal + 2)]; }

  size_t scoreEndIndex(size_t antidiagonal) const
  { return scoreEndsAndOrigins[2 * (antidiagonal + 2)]; }

  // start of the x-drop region (i.e. number of skipped seq1 letters
  // before the x-drop region) for this antidiagonal
  size_t seq1start(size_t antidiagonal) const {
    size_t a = 2 * (antidiagonal + 2);
    return scoreEndsAndOrigins[a] + xdropPadLen - scoreEndsAndOrigins[a + 1];
  }

  size_t scoreOrigin(size_t antidiagonalIncludingDummies) const
  { return scoreEndsAndOrigins[2 * antidiagonalIncludingDummies + 1]; }

  // The index in the score vectors, of the previous "horizontal" cell.
  size_t hori(size_t antidiagonal, size_t seq1coordinate) const
  { return scoreOrigin(antidiagonal + 1) + seq1coordinate - 1; }

  // The index in the score vectors, of the previous "vertical" cell.
  size_t vert(size_t antidiagonal, size_t seq1coordinate) const
  { return scoreOrigin(antidiagonal + 1) + seq1coordinate; }

  // The index in the score vectors, of the previous "diagonal" cell.
  size_t diag(size_t antidiagonal, size_t seq1coordinate) const
  { return scoreOrigin(antidiagonal) + seq1coordinate - 1; }

  // The index in the score vectors, of the previous in-frame horizontal cell.
  size_t hori3(size_t antidiagonal, size_t seq1coordinate) const
  { return scoreOrigin(antidiagonal - 3) + seq1coordinate - 1; }

  // The index in the score vectors, of the previous in-frame vertical cell.
  size_t vert3(size_t antidiagonal, size_t seq1coordinate) const
  { return scoreOrigin(antidiagonal - 3) + seq1coordinate; }

  // The index in the score vectors, of the previous in-frame diagonal cell.
  size_t diag3(size_t antidiagonal, size_t seq1coordinate) const
  { return scoreOrigin(antidiagonal - 6) + seq1coordinate - 1; }

 private:
  std::vector<Score> xScores;  // best score ending with aligned letters
  std::vector<Score> yScores;  // best score ending with insertion in seq1
  std::vector<Score> zScores;  // best score ending with insertion in seq2

  std::vector<size_t> scoreEndsAndOrigins;  // data for each antidiagonal
  size_t numOfAntidiagonals;

  ContiguousQueue<const int *> pssmQueue;
  ContiguousQueue<uchar> seq1queue;
  ReverseQueue<uchar> seq2queue;

  // Our position during the trace-back:
  size_t bestAntidiagonal;
  size_t bestSeq1position;

  void resizeScoresIfSmaller(size_t size) {
    if (xScores.size() < size) {
      xScores.resize(size);
      yScores.resize(size);
      zScores.resize(size);
    }
  }

  void initAntidiagonal(size_t antidiagonalIncludingDummies,
			size_t seq1end, size_t thisEnd, int numCells) {
    const SimdInt mNegInf = simdFill(-INF);
    size_t nextEnd = thisEnd + xdropPadLen + numCells;

    size_t a = 2 * (antidiagonalIncludingDummies + 1);
    if (scoreEndsAndOrigins.size() <= a) {
      scoreEndsAndOrigins.resize(a + 1);
    }
    scoreEndsAndOrigins[a - 1] = nextEnd - seq1end;
    scoreEndsAndOrigins[a] = nextEnd;

    resizeScoresIfSmaller(nextEnd + (simdLen-1));
    for (int i = 0; i < xdropPadLen; i += simdLen) {
      simdStore(&xScores[thisEnd + i], mNegInf);
      simdStore(&yScores[thisEnd + i], mNegInf);
      simdStore(&zScores[thisEnd + i], mNegInf);
    }
  }

  // Puts 2 "dummy" antidiagonals at the start, so that we can safely
  // look-back from subsequent antidiagonals
  void init() {
    initAntidiagonal(0, 0, 0, 0);
    initAntidiagonal(1, 0, xdropPadLen, 0);
    xScores[xdropPadLen - 1] = 0;
    bestAntidiagonal = 0;
  }

  void initAntidiagonal3(size_t antidiagonalIncludingDummies,
			 size_t seq1end, size_t scoreEnd) {
    size_t a = 2 * (antidiagonalIncludingDummies + 1);
    if (scoreEndsAndOrigins.size() <= a) {
      scoreEndsAndOrigins.resize(a + 1);
    }
    scoreEndsAndOrigins[a - 1] = scoreEnd - seq1end;
    scoreEndsAndOrigins[a] = scoreEnd;
    resizeScoresIfSmaller(scoreEnd + (simdLen-1));
  }

  void updateBest(int &bestScore, int score, size_t antidiagonal,
                  const Score *x0, const Score *x0ori);

  void calcBestSeq1position(int bestScore, size_t numOfDummyAntidiagonals) {
    size_t seq1beg = seq1start(bestAntidiagonal + numOfDummyAntidiagonals - 2);
    const Score *x2 = &xScores[diag(bestAntidiagonal, seq1beg)];
    const Score *x2beg = x2;
    while (*x2 != bestScore) ++x2;
    bestSeq1position = x2 - x2beg + seq1beg;
  }

  void init3();
  void initFrame();

  // Everything below here is for alignDna & getNextChunkDna
#if defined __SSE4_1__ || defined __ARM_NEON
  std::vector<TinyScore> xTinyScores;
  std::vector<TinyScore> yTinyScores;
  std::vector<TinyScore> zTinyScores;

  std::vector<int> scoreRises;  // increase of best score, per antidiagonal

  void resizeTinyScoresIfSmaller(size_t size) {
    if (xTinyScores.size() < size) {
      xTinyScores.resize(size);
      yTinyScores.resize(size);
      zTinyScores.resize(size);
    }
  }

  void initAntidiagonalTiny(size_t antidiagonalIncludingDummies,
			    size_t seq1end, size_t thisEnd, int numCells) {
    const SimdUint1 mNegInf = simdOnes1();
    size_t nextEnd = thisEnd + xdropPadLen + numCells;

    size_t a = 2 * (antidiagonalIncludingDummies + 1);
    if (scoreRises.size() <= antidiagonalIncludingDummies) {
      scoreEndsAndOrigins.resize(a + 1);
      scoreRises.resize(antidiagonalIncludingDummies + 1);
    }
    scoreEndsAndOrigins[a - 1] = nextEnd - seq1end;
    scoreEndsAndOrigins[a] = nextEnd;

    resizeTinyScoresIfSmaller(nextEnd + (simdBytes-1));
    simdStore1(&xTinyScores[thisEnd], mNegInf);
    simdStore1(&yTinyScores[thisEnd], mNegInf);
    simdStore1(&zTinyScores[thisEnd], mNegInf);
  }

  void initTiny(int scoreOffset) {
    initAntidiagonalTiny(0, 0, 0, 0);
    initAntidiagonalTiny(1, 0, xdropPadLen, 0);
    xTinyScores[xdropPadLen - 1] = scoreOffset;
    bestAntidiagonal = 2;
  }

  void calcBestSeq1positionTiny(int scoreOffset) {
    size_t seq1beg = seq1start(bestAntidiagonal);
    const TinyScore *x2 = &xTinyScores[diag(bestAntidiagonal, seq1beg)];
    const TinyScore *x2beg = x2;
    int target = scoreOffset - scoreRises[bestAntidiagonal] -
      scoreRises[bestAntidiagonal + 1] - scoreRises[bestAntidiagonal + 2];
    while (*x2 != target) ++x2;
    bestSeq1position = x2 - x2beg + seq1beg;
  }
#endif
};

}

#endif
