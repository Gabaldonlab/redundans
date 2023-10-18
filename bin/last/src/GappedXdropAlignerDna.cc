// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"

#if defined __SSE4_1__ || defined __ARM_NEON

//#include <iostream>  // for debugging

namespace cbrc {

const int seqLoadLen = simdBytes;

const int delimiter = 4;

int GappedXdropAligner::alignDna(const uchar *seq1,
				 const uchar *seq2,
				 bool isForward,
				 const ScoreMatrixRow *scorer,
				 int delOpenCost,
				 int delGrowCost,
				 int insOpenCost,
				 int insGrowCost,
				 int maxScoreDrop,
				 int maxMatchScore,
				 const uchar *toUnmasked) {
  int badScoreDrop = maxScoreDrop + 1;

  delGrowCost = std::min(delGrowCost, badScoreDrop);
  delOpenCost = std::min(delOpenCost, badScoreDrop - delGrowCost);

  insGrowCost = std::min(insGrowCost, badScoreDrop);
  insOpenCost = std::min(insOpenCost, badScoreDrop - insGrowCost);

  const SimdUint1 mNegInf = simdOnes1();
  const SimdUint1 mDelOpenCost = simdFill1(delOpenCost);
  const SimdUint1 mDelGrowCost = simdFill1(delGrowCost);
  const SimdUint1 mInsOpenCost = simdFill1(insOpenCost);
  const SimdUint1 mInsGrowCost = simdFill1(insGrowCost);
  const int seqIncrement = isForward ? 1 : -1;
  const int scoreOffset = maxMatchScore * 2;

  const SimdUint1 scorer4x4 =
    simdSet1(
#ifdef __AVX2__
		 scorer[3][3], scorer[3][2], scorer[3][1], scorer[3][0],
		 scorer[2][3], scorer[2][2], scorer[2][1], scorer[2][0],
		 scorer[1][3], scorer[1][2], scorer[1][1], scorer[1][0],
		 scorer[0][3], scorer[0][2], scorer[0][1], scorer[0][0],
#endif
		 scorer[3][3], scorer[3][2], scorer[3][1], scorer[3][0],
		 scorer[2][3], scorer[2][2], scorer[2][1], scorer[2][0],
		 scorer[1][3], scorer[1][2], scorer[1][1], scorer[1][0],
		 scorer[0][3], scorer[0][2], scorer[0][1], scorer[0][0]);

  int numCells = 1;
  size_t seq1end = 1;
  size_t diagPos = xdropPadLen - 1;
  size_t horiPos = xdropPadLen * 2 - 1;
  size_t thisPos = xdropPadLen * 2;

  int bestScore = 0;
  SimdUint1 mBestScore = mNegInf;
  SimdUint1 mBadScore = simdFill1(scoreOffset + badScoreDrop);
  SimdUint1 mScoreRise1 = simdZero1();
  SimdUint1 mScoreRise2 = simdZero1();

  initTiny(scoreOffset);
  seq1queue.clear();
  seq2queue.clear();

  bool isDna = (toUnmasked[*seq1] < 4 && toUnmasked[*seq2] < 4);

  for (int i = 0; i < seqLoadLen; ++i) {
    uchar x = toUnmasked[*seq1];
    seq1queue.push(x, i);
    seq1 += seqIncrement * (x != delimiter);
    seq2queue.push(toUnmasked[*seq2], i);
  }

  seq2 += seqIncrement;

  size_t antidiagonal;
  for (antidiagonal = 2; /* noop */; ++antidiagonal) {
    int n = numCells - 1;
    const uchar *s1 = &seq1queue.fromEnd(n + seqLoadLen);
    const uchar *s2 = seq2queue.begin();

    initAntidiagonalTiny(antidiagonal, seq1end, thisPos, numCells);
    thisPos += xdropPadLen;
    TinyScore *x0 = &xTinyScores[thisPos];
    TinyScore *y0 = &yTinyScores[thisPos];
    TinyScore *z0 = &zTinyScores[thisPos];
    const TinyScore *y1 = &yTinyScores[horiPos];
    const TinyScore *z1 = &zTinyScores[horiPos + 1];
    const TinyScore *x2 = &xTinyScores[diagPos];

    const SimdUint1 mScoreRise12 = simdAdd1(mScoreRise1, mScoreRise2);
    const SimdUint1 mDelGrowCost1 = simdAdd1(mDelGrowCost, mScoreRise1);
    const SimdUint1 mInsGrowCost1 = simdAdd1(mInsGrowCost, mScoreRise1);

    if (isDna) {
      for (int i = 0; i < numCells; i += simdBytes) {
	SimdUint1 fwd1 = simdLoad1(s1+i);
	SimdUint1 rev2 = simdLoad1(s2+i);
	SimdUint1 j = simdOr1(simdQuadruple1(fwd1), rev2);
	SimdUint1 s = simdChoose1(scorer4x4, j);
	SimdUint1 x = simdAdds1(simdLoad1(x2+i), mScoreRise12);
	SimdUint1 y = simdAdds1(simdLoad1(y1+i), mDelGrowCost1);
	SimdUint1 z = simdAdds1(simdLoad1(z1+i), mInsGrowCost1);
	SimdUint1 b = simdMin1(simdMin1(x, y), z);
	SimdUint1 isDrop = simdGe1(b, mBadScore);
	mBestScore = simdMin1(b, mBestScore);
	simdStore1(x0+i, simdOr1(simdSub1(b, s), isDrop));
	simdStore1(y0+i, simdMin1(simdAdds1(b, mDelOpenCost), y));
	simdStore1(z0+i, simdMin1(simdAdds1(b, mInsOpenCost), z));
      }
    } else {
      bool isDelimiter1 = (s1[n] == delimiter);
      bool isDelimiter2 = (s2[0] == delimiter);
      if (isDelimiter1 || isDelimiter2) {
	badScoreDrop = std::min(badScoreDrop, n * maxMatchScore);
	mBadScore = simdFill1(scoreOffset + badScoreDrop);
      }

      for (int i = 0; i < numCells; i += simdBytes) {
	SimdUint1 s = simdSet1(
#ifdef __AVX2__
			     scorer[s1[31]][s2[31]],
			     scorer[s1[30]][s2[30]],
			     scorer[s1[29]][s2[29]],
			     scorer[s1[28]][s2[28]],
			     scorer[s1[27]][s2[27]],
			     scorer[s1[26]][s2[26]],
			     scorer[s1[25]][s2[25]],
			     scorer[s1[24]][s2[24]],
			     scorer[s1[23]][s2[23]],
			     scorer[s1[22]][s2[22]],
			     scorer[s1[21]][s2[21]],
			     scorer[s1[20]][s2[20]],
			     scorer[s1[19]][s2[19]],
			     scorer[s1[18]][s2[18]],
			     scorer[s1[17]][s2[17]],
			     scorer[s1[16]][s2[16]],
#endif
			     scorer[s1[15]][s2[15]],
			     scorer[s1[14]][s2[14]],
			     scorer[s1[13]][s2[13]],
			     scorer[s1[12]][s2[12]],
			     scorer[s1[11]][s2[11]],
			     scorer[s1[10]][s2[10]],
			     scorer[s1[9]][s2[9]],
			     scorer[s1[8]][s2[8]],
			     scorer[s1[7]][s2[7]],
			     scorer[s1[6]][s2[6]],
			     scorer[s1[5]][s2[5]],
			     scorer[s1[4]][s2[4]],
			     scorer[s1[3]][s2[3]],
			     scorer[s1[2]][s2[2]],
			     scorer[s1[1]][s2[1]],
			     scorer[s1[0]][s2[0]]);

	SimdUint1 x = simdAdds1(simdLoad1(x2+i), mScoreRise12);
	SimdUint1 y = simdAdds1(simdLoad1(y1+i), mDelGrowCost1);
	SimdUint1 z = simdAdds1(simdLoad1(z1+i), mInsGrowCost1);
	SimdUint1 b = simdMin1(simdMin1(x, y), z);
	SimdUint1 isDrop = simdGe1(b, mBadScore);
	mBestScore = simdMin1(b, mBestScore);
	simdStore1(x0+i, simdOr1(simdSub1(b, s), isDrop));
	simdStore1(y0+i, simdMin1(simdAdds1(b, mDelOpenCost), y));
	simdStore1(z0+i, simdMin1(simdAdds1(b, mInsOpenCost), z));
	s1 += simdBytes;
	s2 += simdBytes;
      }
      if (isDelimiter2) x0[0] = droppedTinyScore;
      if (isDelimiter1) x0[n] = droppedTinyScore;  // maybe n=0
    }

    mScoreRise2 = mScoreRise1;
    mScoreRise1 = simdZero1();
    int newBestScore = simdHorizontalMin1(mBestScore);
    int rise = 0;
    if (newBestScore < scoreOffset) {
      rise = scoreOffset - newBestScore;
      bestScore += rise;
      bestAntidiagonal = antidiagonal;
      mBestScore = mNegInf;
      mScoreRise1 = simdFill1(rise);
    }
    scoreRises[antidiagonal] = rise;

    diagPos = horiPos;
    horiPos = thisPos - 1;
    thisPos += numCells;

    if (x0[n] != droppedTinyScore) {
      ++numCells;
      ++seq1end;
      uchar x = toUnmasked[*seq1];
      seq1queue.push(x, n + seqLoadLen);
      seq1 += seqIncrement * (x != delimiter);
      uchar z = seq1queue.fromEnd(seqLoadLen);
      if (z >= 4) {
	isDna = false;
      }
    }

    if (x0[0] != droppedTinyScore) {
      uchar y = toUnmasked[*seq2];
      seq2queue.push(y, n + seqLoadLen);
      seq2 += seqIncrement;
      if (y >= 4) {
	isDna = false;
      }
    } else {
      --numCells;
      if (numCells == 0) break;
      ++diagPos;
      ++horiPos;
    }
  }

  bestAntidiagonal -= 2;
  calcBestSeq1positionTiny(scoreOffset);
  numOfAntidiagonals = antidiagonal - 1;
  return bestScore;
}

bool GappedXdropAligner::getNextChunkDna(size_t &end1,
					 size_t &end2,
					 size_t &length,
					 int delOpenCost,
					 int delGrowCost,
					 int insOpenCost,
					 int insGrowCost) {
  if (bestAntidiagonal == 0) return false;

  end1 = bestSeq1position;
  end2 = bestAntidiagonal - bestSeq1position;

  int x, y, z;
  while (1) {
    size_t h = hori(bestAntidiagonal, bestSeq1position);
    size_t v = vert(bestAntidiagonal, bestSeq1position);
    size_t d = diag(bestAntidiagonal, bestSeq1position);
    x = xTinyScores[d] + scoreRises[bestAntidiagonal];
    y = yTinyScores[h] + delGrowCost;
    z = zTinyScores[v] + insGrowCost;
    if (x > y || x > z || bestAntidiagonal == 0) break;
    bestAntidiagonal -= 2;
    bestSeq1position -= 1;
  }

  length = end1 - bestSeq1position;
  if (bestAntidiagonal == 0) return true;

  while (1) {
    bool isDel = (y <= z);
    bestAntidiagonal -= 1;
    if (isDel) bestSeq1position -= 1;
    size_t h = hori(bestAntidiagonal, bestSeq1position);
    size_t v = vert(bestAntidiagonal, bestSeq1position);
    size_t d = diag(bestAntidiagonal, bestSeq1position);
    x = xTinyScores[d] + scoreRises[bestAntidiagonal];
    y = yTinyScores[h] + delGrowCost;
    z = zTinyScores[v] + insGrowCost;
    if (isDel) {
      y -= delOpenCost;
    } else {
      z -= insOpenCost;
    }
    if (x <= y && x <= z) return true;
  }
}

}

#endif
