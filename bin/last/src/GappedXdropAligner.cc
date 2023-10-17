// Copyright 2011, 2012, 2013 Martin C. Frith

// The algorithm is based on these recurrence formulas, for
// generalized affine gap costs.  For standard affine gap costs, set
// gup=infinity.
//
// gop = gapExistenceCost
// gep = gapExtensionCost
// gup = gapUnalignedCost
//
// The 1st sequence: s(1), s(2), s(3), ...
// The 2nd sequence: t(1), t(2), t(3), ...
//
// matchScore(i, j)  =  the score for aligning s(i) with t(j).
//
// Initialization:
// x(i, 0)  =  y(i, 0)  =  z(i, 0)  =  -INF  (for all i >= 0)
// x(0, j)  =  y(0, j)  =  z(0, j)  =  -INF  (for all j >= 0)
// x(0, 0)  =  0
//
// Recurrence (i > 0 and j > 0):
// X(i, j)  =  x(i-1, j-1)
// Y(i, j)  =  max[ y(i-1, j) - gep, y(i-1, j-1) - gup ]
// Z(i, j)  =  max[ z(i, j-1) - gep, z(i-1, j-1) - gup ]
// b(i, j)  =  max[ X(i, j), Y(i, j), Z(i, j) ]
// x(i, j)  =  b(i, j) + matchScore(i, j)
// y(i, j)  =  max[ b(i, j) - gop, Y(i, j) ]
// z(i, j)  =  max[ b(i, j) - gop, Z(i, j) ]
//
// Interpretations:
// X(i, j)  =  the best score for any alignment ending with s(i-1)
//             aligned to t(j-1).
// b(i, j)  =  the best score for any alignment ending at s(i-1) and
//             t(j-1).
// x(i, j)  =  the best score for any alignment ending with s(i)
//             aligned to t(j).

// The recurrences are calculated antidiagonal-by-antidiagonal, where:
// antidiagonal  =  i + j

// We store x(i, j), y(i, j), and z(i, j) in the following way.
// xScores: oxx2x33x444x5555x66666...
// yScores: xxx2x33x444x5555x66666...
// zScores: xxx2x33x444x5555x66666...
// "o" indicates a cell with score = 0.
// "x" indicates a pad cell with score = -INF.
// "2", "3", etc. indicate cells in antidiagonal 2, 3, etc.

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"

#include <ostream>
//#include <iostream>  // for debugging

namespace cbrc {

int GappedXdropAligner::align(const uchar *seq1,
                              const uchar *seq2,
                              bool isForward,
			      int globality,
                              const ScoreMatrixRow *scorer,
                              int delExistenceCost,
                              int delExtensionCost,
                              int insExistenceCost,
                              int insExtensionCost,
                              int gapUnalignedCost,
			      bool isAffine,
                              int maxScoreDrop,
                              int maxMatchScore) {
  const SimdInt mNegInf = simdFill(-INF);
  const SimdInt mDelOpenCost = simdFill(delExistenceCost);
  const SimdInt mDelGrowCost = simdFill(delExtensionCost);
  const SimdInt mInsOpenCost = simdFill(insExistenceCost);
  const SimdInt mInsGrowCost = simdFill(insExtensionCost);
  const int seqIncrement = isForward ? 1 : -1;

  int numCells = 1;
  size_t seq1end = 1;
  size_t diagPos = xdropPadLen - 1;
  size_t horiPos = xdropPadLen * 2 - 1;
  size_t thisPos = xdropPadLen * 2;

  int bestScore = 0;
  SimdInt mBestScore = simdZero();
  int bestEdgeScore = -INF;
  size_t bestEdgeAntidiagonal = 0;

  init();
  pssmQueue.clear();
  seq2queue.clear();

  bool isDelimiter1 = isDelimiter(0, scorer[*seq1]);
  bool isDelimiter2 = isDelimiter(*seq2, scorer[0]);

  for (int i = 0; i < simdLen; ++i) {
    const int *x = scorer[*seq1];
    pssmQueue.push(x, i);
    seq1 += seqIncrement * !isDelimiter(0, x);
    seq2queue.push(*seq2, i);
  }

  seq2 += seqIncrement;

  size_t antidiagonal;
  for (antidiagonal = 0; /* noop */; ++antidiagonal) {
    int n = numCells - 1;
    const const_int_ptr *s1 = &pssmQueue.fromEnd(n + simdLen);
    const uchar *s2 = seq2queue.begin();

    initAntidiagonal(antidiagonal + 2, seq1end, thisPos, numCells);
    thisPos += xdropPadLen;
    Score *x0 = &xScores[thisPos];
    Score *y0 = &yScores[thisPos];
    Score *z0 = &zScores[thisPos];
    const Score *y1 = &yScores[horiPos];
    const Score *z1 = &zScores[horiPos + 1];
    const Score *x2 = &xScores[diagPos];

    if (!globality && (isDelimiter1 || isDelimiter2)) {
      updateMaxScoreDrop(maxScoreDrop, n, maxMatchScore);
    }

    int minScore = bestScore - maxScoreDrop;
    SimdInt mMinScore = simdFill(minScore);

    if (globality && isDelimiter2) {
      const Score *z2 = &zScores[diagPos];
      int b = maxValue(x2[0], z1[0]-insExtensionCost, z2[0]-gapUnalignedCost);
      updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestSeq1position,
		  minScore, b, antidiagonal, seq1end - numCells);
    }

    if (globality && isDelimiter1) {
      const Score *y2 = &yScores[diagPos];
      int b = maxValue(x2[n], y1[n]-delExtensionCost, y2[n]-gapUnalignedCost);
      updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestSeq1position,
		  minScore, b, antidiagonal, seq1end-1);
    }

    if (isAffine) {
      for (int i = 0; i < numCells; i += simdLen) {
	SimdInt s = simdSet(
#if defined __SSE4_1__ || defined __ARM_NEON
#ifdef __AVX2__
			    s1[7][s2[7]],
			    s1[6][s2[6]],
			    s1[5][s2[5]],
			    s1[4][s2[4]],
#endif
			    s1[3][s2[3]],
			    s1[2][s2[2]],
			    s1[1][s2[1]],
#endif
			    s1[0][s2[0]]);
	SimdInt x = simdLoad(x2+i);
	SimdInt y = simdSub(simdLoad(y1+i), mDelGrowCost);
	SimdInt z = simdSub(simdLoad(z1+i), mInsGrowCost);
	SimdInt b = simdMax(simdMax(x, y), z);
	mBestScore = simdMax(b, mBestScore);
	SimdInt xNew = simdBlend(simdAdd(b, s), mNegInf, simdGt(mMinScore, b));
	simdStore(x0+i, xNew);
	simdStore(y0+i, simdMax(simdSub(b, mDelOpenCost), y));
	simdStore(z0+i, simdMax(simdSub(b, mInsOpenCost), z));
	s1 += simdLen;
	s2 += simdLen;
      }

      int newBestScore = simdHorizontalMax(mBestScore);
      if (newBestScore > bestScore) {
	bestScore = newBestScore;
	bestAntidiagonal = antidiagonal;
      }
    } else {
      const Score *y2 = &yScores[diagPos];
      const Score *z2 = &zScores[diagPos];
      for (int i = 0; i < numCells; ++i) {
        int x = x2[i];
        int y = maxValue(y1[i] - delExtensionCost, y2[i] - gapUnalignedCost);
        int z = maxValue(z1[i] - insExtensionCost, z2[i] - gapUnalignedCost);
        int b = maxValue(x, y, z);
        if (b >= minScore) {
	  if (b > bestScore) {
	    bestScore = b;
	    bestAntidiagonal = antidiagonal;
	  }
          x0[i] = b + s1[0][*s2];
          y0[i] = maxValue(b - delExistenceCost, y);
          z0[i] = maxValue(b - insExistenceCost, z);
        }
        else x0[i] = y0[i] = z0[i] = -INF;
	++s1;
	++s2;
      }
    }

    diagPos = horiPos;
    horiPos = thisPos - 1;
    thisPos += numCells;

    if (x0[n] > -INF / 2) {
      ++numCells;
      ++seq1end;
      const int *x = scorer[*seq1];
      pssmQueue.push(x, n + simdLen);
      seq1 += seqIncrement * !isDelimiter(0, x);
      isDelimiter1 = isDelimiter(0, pssmQueue.fromEnd(simdLen));
    }

    if (x0[0] > -INF / 2) {
      uchar y = *seq2;
      seq2queue.push(y, n + simdLen);
      seq2 += seqIncrement;
      isDelimiter2 = isDelimiter(y, scorer[0]);
    } else {
      --numCells;
      if (numCells == 0) break;
      ++diagPos;
      ++horiPos;
    }
  }

  if (globality) {
    bestAntidiagonal = bestEdgeAntidiagonal;
    bestScore = bestEdgeScore;
  } else {
    calcBestSeq1position(bestScore, 2);
  }
  numOfAntidiagonals = antidiagonal + 1;
  return bestScore;
}

bool GappedXdropAligner::getNextChunk(size_t &end1,
                                      size_t &end2,
                                      size_t &length,
				      int delExistenceCost,
				      int delExtensionCost,
				      int insExistenceCost,
				      int insExtensionCost,
                                      int gapUnalignedCost) {
  if (bestAntidiagonal == 0) return false;

  end1 = bestSeq1position;
  end2 = bestAntidiagonal - bestSeq1position;
  const size_t undefined = -1;
  length = undefined;

  int state = 0;

  while (1) {
    assert(bestSeq1position <= bestAntidiagonal);

    size_t h = hori(bestAntidiagonal, bestSeq1position);
    size_t v = vert(bestAntidiagonal, bestSeq1position);
    size_t d = diag(bestAntidiagonal, bestSeq1position);

    int x = xScores[d];
    int y = yScores[h] - delExtensionCost;
    int z = zScores[v] - insExtensionCost;
    int a = yScores[d] - gapUnalignedCost;
    int b = zScores[d] - gapUnalignedCost;

    if (state == 1 || state == 3) {
      y += delExistenceCost;
      a += delExistenceCost;
    }

    if (state == 2 || state == 4) {
      z += insExistenceCost;
      b += insExistenceCost;
    }

    state = maxIndex(x, y, z, a, b);

    if (length == undefined && (state > 0 || bestAntidiagonal == 0)) {
      length = end1 - bestSeq1position;
      assert(length != undefined);
    }

    if (length != undefined && state == 0) return true;

    if (state < 1 || state > 2) bestAntidiagonal -= 2;
    else                        bestAntidiagonal -= 1;

    if (state != 2) bestSeq1position -= 1;
  }
}

void GappedXdropAligner::writeShape(std::ostream &out) const {
  for (size_t i = 0; i < numAntidiagonals(); ++i) {
    size_t s = seq1start(i);
    out << s << "\t" << (s + numCellsAndPads(i) - xdropPadLen) << "\n";
  }
}

}
