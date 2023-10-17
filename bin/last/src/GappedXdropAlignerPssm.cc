// Copyright 2011, 2012, 2013 Martin C. Frith

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"

namespace cbrc {

int GappedXdropAligner::alignPssm(const uchar *seq,
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
                                  int maxMatchScore) {
  const int *vectorOfMatchScores = *pssm;
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
  seq1queue.clear();
  pssmQueue.clear();

  bool isDelimiter1 = isDelimiter(*seq, vectorOfMatchScores);
  bool isDelimiter2 = isDelimiter(0, vectorOfMatchScores);

  for (int i = 0; i < simdLen; ++i) {
    uchar x = *seq;
    seq1queue.push(x, i);
    seq += seqIncrement * !isDelimiter(x, vectorOfMatchScores);
    pssmQueue.push(vectorOfMatchScores, i);
  }

  pssm += seqIncrement;

  size_t antidiagonal;
  for (antidiagonal = 0; /* noop */; ++antidiagonal) {
    int n = numCells - 1;
    const uchar *s1 = &seq1queue.fromEnd(n + simdLen);
    const const_int_ptr *s2 = &pssmQueue.fromEnd(1);

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
			    s2[-7][s1[7]],
			    s2[-6][s1[6]],
			    s2[-5][s1[5]],
			    s2[-4][s1[4]],
#endif
			    s2[-3][s1[3]],
			    s2[-2][s1[2]],
			    s2[-1][s1[1]],
#endif
			    s2[-0][s1[0]]);
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
	s2 -= simdLen;
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
          x0[i] = b + (*s2)[*s1];
          y0[i] = maxValue(b - delExistenceCost, y);
          z0[i] = maxValue(b - insExistenceCost, z);
        }
        else x0[i] = y0[i] = z0[i] = -INF;
	++s1;
	--s2;
      }
    }

    diagPos = horiPos;
    horiPos = thisPos - 1;
    thisPos += numCells;

    if (x0[n] > -INF / 2) {
      ++numCells;
      ++seq1end;
      uchar x = *seq;
      seq1queue.push(x, n + simdLen);
      seq += seqIncrement * !isDelimiter(x, vectorOfMatchScores);
      isDelimiter1 = isDelimiter(seq1queue.fromEnd(simdLen),
				 vectorOfMatchScores);
    }

    if (x0[0] > -INF / 2) {
      const int *y = *pssm;
      pssmQueue.push(y, n + simdLen);
      pssm += seqIncrement;
      isDelimiter2 = isDelimiter(0, y);
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

}
