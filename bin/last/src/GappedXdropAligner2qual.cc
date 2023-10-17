// Copyright 2011, 2012, 2013 Martin C. Frith

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"
#include "TwoQualityScoreMatrix.hh"

namespace cbrc {

static bool isDelimiter2qual(uchar c) {
  return c == 4;  // a bit ugly (hard-wired delimiter value)
}

int GappedXdropAligner::align2qual(const uchar *seq1,
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
                                   int maxMatchScore) {
  const int seqIncrement = isForward ? 1 : -1;

  int numCells = 1;
  size_t seq1end = 1;
  size_t diagPos = xdropPadLen - 1;
  size_t horiPos = xdropPadLen * 2 - 1;
  size_t thisPos = xdropPadLen * 2;

  int bestScore = 0;
  int bestEdgeScore = -INF;
  size_t bestEdgeAntidiagonal = 0;

  init();

  size_t antidiagonal;
  for (antidiagonal = 0; /* noop */; ++antidiagonal) {
    int n = numCells - 1;
    const uchar *s1 = seq1;
    const uchar *q1 = qual1;
    const uchar *s2 = seq2;
    const uchar *q2 = qual2;
    seq2 += seqIncrement;
    qual2 += seqIncrement;

    initAntidiagonal(antidiagonal + 2, seq1end, thisPos, numCells);
    thisPos += xdropPadLen;
    Score *x0 = &xScores[thisPos];
    Score *y0 = &yScores[thisPos];
    Score *z0 = &zScores[thisPos];
    const Score *y1 = &yScores[horiPos];
    const Score *z1 = &zScores[horiPos + 1];
    const Score *x2 = &xScores[diagPos];

    bool isDelimiter1 = isDelimiter2qual(s1[n * seqIncrement]);
    bool isDelimiter2 = isDelimiter2qual(s2[0]);

    if (!globality && (isDelimiter1 || isDelimiter2)) {
      updateMaxScoreDrop(maxScoreDrop, n, maxMatchScore);
    }

    int minScore = bestScore - maxScoreDrop;

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

    const Score *x0last = x0 + n;

    if (isAffine) {
      if (isForward)
        while (1) {
          int x = *x2;
          int y = *y1 - delExtensionCost;
          int z = *z1 - insExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
	    if (b > bestScore) {
	      bestScore = b;
	      bestAntidiagonal = antidiagonal;
	    }
            *x0 = b + scorer(*s1, *s2, *q1, *q2);
            *y0 = maxValue(b - delExistenceCost, y);
            *z0 = maxValue(b - insExistenceCost, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          ++s1;  ++q1;  --s2;  --q2;  ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;
        }
      else
        while (1) {
          int x = *x2;
          int y = *y1 - delExtensionCost;
          int z = *z1 - insExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
	    if (b > bestScore) {
	      bestScore = b;
	      bestAntidiagonal = antidiagonal;
	    }
            *x0 = b + scorer(*s1, *s2, *q1, *q2);
            *y0 = maxValue(b - delExistenceCost, y);
            *z0 = maxValue(b - insExistenceCost, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          --s1;  --q1;  ++s2;  ++q2;  ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;
        }
    } else {
      const Score *y2 = &yScores[diagPos];
      const Score *z2 = &zScores[diagPos];
      while (1) {
        int x = *x2;
        int y = maxValue(*y1 - delExtensionCost, *y2 - gapUnalignedCost);
        int z = maxValue(*z1 - insExtensionCost, *z2 - gapUnalignedCost);
        int b = maxValue(x, y, z);
        if (b >= minScore) {
	  if (b > bestScore) {
	    bestScore = b;
	    bestAntidiagonal = antidiagonal;
	  }
          *x0 = b + scorer(*s1, *s2, *q1, *q2);
          *y0 = maxValue(b - delExistenceCost, y);
          *z0 = maxValue(b - insExistenceCost, z);
        }
        else *x0 = *y0 = *z0 = -INF;
        if (x0 == x0last) break;
        ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;  ++y2;  ++z2;
	s1 += seqIncrement;  q1 += seqIncrement;
	s2 -= seqIncrement;  q2 -= seqIncrement;
      }
    }

    diagPos = horiPos;
    horiPos = thisPos - 1;
    thisPos += numCells;

    if (x0[0] > -INF / 2) {
      ++numCells;
      ++seq1end;
    }

    if (x0[-n] <= -INF / 2) {
      --numCells;
      if (numCells == 0) break;
      ++diagPos;
      ++horiPos;
      seq1 += seqIncrement;
      qual1 += seqIncrement;
      seq2 -= seqIncrement;
      qual2 -= seqIncrement;
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

}
