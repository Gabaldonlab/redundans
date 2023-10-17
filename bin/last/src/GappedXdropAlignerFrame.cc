// Author: Martin C. Frith 2020
// SPDX-License-Identifier: GPL-3.0-or-later

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"
//#include <iostream>  // for debugging

namespace cbrc {

typedef const uchar *const_uchar_ptr;

// Puts 6 "dummy" antidiagonals at the start, so that we can safely
// look-back from subsequent antidiagonals.
void GappedXdropAligner::initFrame() {
  initAntidiagonal(0, 0, 0, 0);
  initAntidiagonal(1, 0, xdropPadLen, 0);
  initAntidiagonal(2, 0, xdropPadLen * 2, 0);
  initAntidiagonal(3, 0, xdropPadLen * 3, 0);
  initAntidiagonal(4, 0, xdropPadLen * 4, 0);
  initAntidiagonal(5, 0, xdropPadLen * 5, 0);
  xScores[xdropPadLen - 1] = 0;
  bestAntidiagonal = 0;
}

int GappedXdropAligner::alignFrame(const uchar *protein,
				   const uchar *frame0,
				   const uchar *frame1,
				   const uchar *frame2,
				   bool isForward,
				   const ScoreMatrixRow *scorer,
				   const GapCosts &gapCosts,
				   int maxScoreDrop) {
  if (!isForward) {
    --protein; --frame0; --frame1; --frame2;
  }
  const_uchar_ptr frames[] = {frame0, frame1, frame2};
  const int seqIncrement = isForward ? 1 : -1;

  const int delOpenScore = -gapCosts.delPieces[0].openCost;
  const int insOpenScore = -gapCosts.insPieces[0].openCost;
  const int delScore1 = gapCosts.delScore1;
  const int delScore2 = gapCosts.delScore2;
  const int delScore3 = gapCosts.delScore3;
  const int insScore1 = gapCosts.insScore1;
  const int insScore2 = gapCosts.insScore2;
  const int insScore3 = gapCosts.insScore3;

  int runOfDrops = 2;
  int runOfEdges = 0;
  int numCells = 1;
  size_t proteinEnd = 1;
  size_t diagPos6 = xdropPadLen - 1;
  size_t horiPos5 = xdropPadLen * 2 - 1;
  size_t horiPos4 = xdropPadLen * 3 - 1;
  size_t horiPos3 = xdropPadLen * 4 - 1;
  size_t vertPos2 = xdropPadLen * 5;
  size_t vertPos1 = xdropPadLen * 6;
  size_t thisPos  = xdropPadLen * 6;

  int bestScore = 0;
  int newBestScore = 0;

  initFrame();

  for (size_t antidiagonal = 0; /* noop */; ++antidiagonal) {
    const uchar *s1 = protein;
    frames[antidiagonal % 3] += seqIncrement;
    const uchar *s2 = frames[antidiagonal % 3];

    initAntidiagonal(antidiagonal + 6, proteinEnd, thisPos, numCells);
    thisPos += xdropPadLen;
    Score *X0 = &xScores[thisPos];
    Score *Y0 = &yScores[thisPos];
    Score *Z0 = &zScores[thisPos];
    const Score *Z1 = &zScores[vertPos1];
    const Score *Z2 = &zScores[vertPos2];
    const Score *Y3 = &yScores[horiPos3];
    const Score *Z3 = &zScores[horiPos3 + 1];
    const Score *Y4 = &yScores[horiPos4];
    const Score *Y5 = &yScores[horiPos5];
    const Score *X6 = &xScores[diagPos6];

    int minScore = bestScore - maxScoreDrop;

    for (int i = 0; i < numCells; ++i) {
      s2 -= seqIncrement;
      int s = scorer[*s1][*s2];
      int y1 = Y5[i] + delScore1;
      int y2 = Y4[i] + delScore2;
      int y3 = Y3[i] + delScore3;
      int z1 = Z1[i] + insScore1;
      int z2 = Z2[i] + insScore2;
      int z3 = Z3[i] + insScore3;
      int b = maxValue(X6[i], y1, y2, y3, z1, z2, z3);
      bool isDrop = (b < minScore);
      newBestScore = maxValue(newBestScore, b);
      X0[i] = isDrop ? -INF : b + s;
      Y0[i] = maxValue(b + delOpenScore, y3);
      Z0[i] = maxValue(b + insOpenScore, z3);
      s1 += seqIncrement;
    }

    if (newBestScore > bestScore) {
      bestScore = newBestScore;
      bestAntidiagonal = antidiagonal;
    }

    diagPos6 = horiPos5;
    horiPos5 = horiPos4;
    horiPos4 = horiPos3;
    horiPos3 = vertPos2 - 1;
    vertPos2 = vertPos1;
    vertPos1 = thisPos;
    thisPos += numCells;

    int n = numCells - 1;
    if (X0[n] > -INF / 2 || runOfEdges) {
      ++runOfEdges;
      if (runOfEdges == 3) {
	++numCells;
	++proteinEnd;
	runOfEdges = 0;
      }
    }

    if (X0[0] > -INF / 2) {
      runOfDrops = 0;
    } else {
      ++runOfDrops;
      if (runOfDrops == 3) {
	--numCells;
	if (numCells == 0) break;
	protein += seqIncrement;
	frames[0] -= seqIncrement;
	frames[1] -= seqIncrement;
	frames[2] -= seqIncrement;
	++diagPos6;
	++horiPos5;
	++horiPos4;
	++horiPos3;
	++vertPos2;
	++vertPos1;
	runOfDrops = 0;
      }
    }
  }

  calcBestSeq1position(bestScore, 6);
  return bestScore;
}

bool GappedXdropAligner::getNextChunkFrame(size_t &end1,
					   size_t &end2,
					   size_t &length,
					   int &costOfNearerGap,
					   const GapCosts &gapCosts) {
  if (bestAntidiagonal == 0) return false;

  const int delOpenScore = -gapCosts.delPieces[0].openCost;
  const int insOpenScore = -gapCosts.insPieces[0].openCost;
  const int delScore1 = gapCosts.delScore1;
  const int delScore2 = gapCosts.delScore2;
  const int delScore3 = gapCosts.delScore3;
  const int insScore1 = gapCosts.insScore1;
  const int insScore2 = gapCosts.insScore2;
  const int insScore3 = gapCosts.insScore3;

  end1 = bestSeq1position;
  end2 = bestAntidiagonal - bestSeq1position * 3;

  int opt[7];
  int state = 0;

  while (1) {
    opt[0] = xScores[diag(bestAntidiagonal + 0, bestSeq1position)];
    opt[1] = yScores[hori(bestAntidiagonal + 0, bestSeq1position)] + delScore1;
    opt[2] = yScores[hori(bestAntidiagonal + 1, bestSeq1position)] + delScore2;
    opt[3] = yScores[hori(bestAntidiagonal + 2, bestSeq1position)] + delScore3;
    opt[4] = zScores[vert(bestAntidiagonal + 2, bestSeq1position)] + insScore3;
    opt[5] = zScores[vert(bestAntidiagonal + 3, bestSeq1position)] + insScore2;
    opt[6] = zScores[vert(bestAntidiagonal + 4, bestSeq1position)] + insScore1;
    state = std::max_element(opt, opt + 7) - opt;
    if (state != 0 || bestAntidiagonal == 0) break;
    bestAntidiagonal -= 6;
    bestSeq1position -= 1;
  }

  length = end1 - bestSeq1position;
  if (bestAntidiagonal == 0) {
    costOfNearerGap = 0;
    return true;
  }

  int scoreHere = opt[state];

  do {
    bool isDel = (state < 4);
    bestAntidiagonal -= 7 - state - isDel;
    if (isDel) bestSeq1position -= 1;
    opt[0] = xScores[diag(bestAntidiagonal + 0, bestSeq1position)];
    opt[1] = yScores[hori(bestAntidiagonal + 0, bestSeq1position)] + delScore1;
    opt[2] = yScores[hori(bestAntidiagonal + 1, bestSeq1position)] + delScore2;
    opt[3] = yScores[hori(bestAntidiagonal + 2, bestSeq1position)] + delScore3;
    opt[4] = zScores[vert(bestAntidiagonal + 2, bestSeq1position)] + insScore3;
    opt[5] = zScores[vert(bestAntidiagonal + 3, bestSeq1position)] + insScore2;
    opt[6] = zScores[vert(bestAntidiagonal + 4, bestSeq1position)] + insScore1;
    if (isDel) {
      opt[3] -= delOpenScore;
    } else {
      opt[4] -= insOpenScore;
    }
    state = std::max_element(opt, opt + 7) - opt;
  } while (state != 0);

  costOfNearerGap = opt[0] - scoreHere;
  return true;
}

}
