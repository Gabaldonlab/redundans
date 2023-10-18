// Author: Martin C. Frith 2020
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mcf_frameshift_xdrop_aligner.hh"

#include <assert.h>
#include <math.h>

#include <algorithm>
//#include <iostream>  // for debugging

namespace mcf {

double FrameshiftXdropAligner::forward(const uchar *protein,
				       const uchar *frame0,
				       const uchar *frame1,
				       const uchar *frame2,
				       bool isRightwardExtension,
				       const const_dbl_ptr *substitutionProbs,
				       const GapCosts &gapCosts,
				       double probDropLimit) {
  const int seqIncrement = isRightwardExtension ? 1 : -1;
  if (isRightwardExtension) {
    --protein; --frame0; --frame1; --frame2;
  }
  proteinPtr = protein;
  frames[0] = frame0 + seqIncrement;
  frames[1] = frame1;
  frames[2] = frame2;

  int numCells = 1;  // number of DynProg cells in this antidiagonal
  size_t proteinEnd = 1;
  size_t diagPos6 = 2 * padSize - 1;
  size_t horiPos5 = diagPos6 + padSize;
  size_t horiPos4 = horiPos5 + padSize;
  size_t horiPos3 = horiPos4 + padSize;
  size_t vertPos2 = horiPos3 + padSize + 1;
  size_t vertPos1 = vertPos2 + padSize;
  size_t thisPos  = vertPos1 + numCells;
  size_t antidiagonal = 1;

  if (xdropShape.empty()) {
    xdropShape.resize(3);
    xdropShape[0] = vertPos2;
    xdropShape[1] = proteinEnd;
    xdropShape[2] = thisPos;
    xFwdProbs.resize(thisPos);
    yFwdProbs.resize(thisPos);
    zFwdProbs.resize(thisPos);
    xFwdProbs[vertPos1] = 1;
    yFwdProbs[vertPos1] = 1;
    zFwdProbs[vertPos1] = 1;
  }

  double sumOfProbRatios = 1;
  double logSumOfProbRatios = 0;

  if (substitutionProbs[0][*(frames[0])] <= 0) {  // at end of DNA sequence
    numOfAntidiagonals = 1;
    rescaledSumOfProbRatios = 1;
    return 0;
  }

  const double delOpenProb = gapCosts.delProbPieces[0].openProb;
  const double insOpenProb = gapCosts.insProbPieces[0].openProb;
  const double delProb1 = gapCosts.delProb1;
  const double delProb2 = gapCosts.delProb2;
  const double delProb3 = gapCosts.delProb3;
  const double insProb1 = gapCosts.insProb1;
  const double insProb2 = gapCosts.insProb2;
  const double insProb3 = gapCosts.insProb3;

  int runOfDrops = 0;
  int runOfEdges = (substitutionProbs[*(protein + seqIncrement)][0] > 0);

  while (1) {
    const uchar *s1 = proteinPtr;
    frames[antidiagonal % 3] += seqIncrement;
    const uchar *s2 = frames[antidiagonal % 3];

    resizeFwdProbsIfSmaller(thisPos + padSize + numCells);
    for (int i = 0; i < padSize; ++i) {
      xFwdProbs[thisPos + i] = 0;
      yFwdProbs[thisPos + i] = 0;
      zFwdProbs[thisPos + i] = 0;
    }
    thisPos += padSize;
    Prob *X0 = &xFwdProbs[thisPos];
    Prob *Y0 = &yFwdProbs[thisPos];
    Prob *Z0 = &zFwdProbs[thisPos];
    const Prob *X1 = &xFwdProbs[vertPos1];
    const Prob *X2 = &xFwdProbs[vertPos2];
    const Prob *Y3 = &yFwdProbs[horiPos3];
    const Prob *Z3 = &zFwdProbs[horiPos3 + 1];
    const Prob *X4 = &xFwdProbs[horiPos4];
    const Prob *X5 = &xFwdProbs[horiPos5];
    const Prob *X6 = &xFwdProbs[diagPos6];

    double minProb = sumOfProbRatios * probDropLimit;

    bool isDelimiter2 = (substitutionProbs[0][*s2] <= 0);

    for (int i = 0; i < numCells; ++i) {
      s2 -= seqIncrement;
      double x = X6[i] * substitutionProbs[*s1][*s2];
      double y = X5[i] * delProb1 + X4[i] * delProb2 + Y3[i] * delProb3;
      double z = X1[i] * insProb1 + X2[i] * insProb2 + Z3[i] * insProb3;
      double b = x + y * delOpenProb + z * insOpenProb;
      X0[i] = b;
      Y0[i] = b + y;
      Z0[i] = b + z;
      sumOfProbRatios += b;
      s1 += seqIncrement;
    }

    bool isDelimiter1 = (substitutionProbs[*s1][0] <= 0);

    diagPos6 = horiPos5;
    horiPos5 = horiPos4;
    horiPos4 = horiPos3;
    horiPos3 = vertPos2 - 1;
    vertPos2 = vertPos1;
    vertPos1 = thisPos;
    thisPos += numCells;

    ++antidiagonal;

    size_t a2 = antidiagonal * 2;
    if (xdropShape.size() < a2 + 1) xdropShape.resize(a2 + 1);
    xdropShape[a2 - 1] = proteinEnd;
    xdropShape[a2] = thisPos;

    int n = numCells - 1;
    bool isOkBeg = (X0[0] >= minProb && !isDelimiter2);
    bool isOkEnd = (X0[n] >= minProb && !isDelimiter1);

    if (isOkEnd || runOfEdges) {
      ++runOfEdges;
      if (runOfEdges == 3) {
	++numCells;
	++proteinEnd;
	runOfEdges = 0;
      }
    }

    if (isOkBeg) {
      runOfDrops = 0;
    } else {
      ++runOfDrops;
      if (runOfDrops == 3) {
	--numCells;
	if (numCells == 0) break;
	proteinPtr += seqIncrement;
	frames[0] -= seqIncrement;
	frames[1] -= seqIncrement;
	frames[2] -= seqIncrement;
	++diagPos6; ++horiPos5; ++horiPos4; ++horiPos3; ++vertPos2; ++vertPos1;
	runOfDrops = 0;
      }
    }

    if (antidiagonal % rescaleStep == 0) {
      double scale = 1 / sumOfProbRatios;
      size_t numOfRescales = antidiagonal / rescaleStep;
      if (rescales.size() < numOfRescales) rescales.resize(numOfRescales);
      rescales[numOfRescales - 1] = scale;
      rescaleFwdProbs(begInProbs(antidiagonal - 6), thisPos, scale);
      logSumOfProbRatios += log(sumOfProbRatios);
      sumOfProbRatios = 1;
    }
  }

  numOfAntidiagonals = antidiagonal;
  rescaledSumOfProbRatios = sumOfProbRatios;
  return logSumOfProbRatios + log(sumOfProbRatios);
}

void FrameshiftXdropAligner::backward(bool isRightwardExtension,
				      const const_dbl_ptr *substitutionProbs,
				      const GapCosts &gapCosts) {
  const int seqIncrement = isRightwardExtension ? 1 : -1;

  const double delOpenProb = gapCosts.delProbPieces[0].openProb;
  const double insOpenProb = gapCosts.insProbPieces[0].openProb;
  const double delProb1 = gapCosts.delProb1;
  const double delProb2 = gapCosts.delProb2;
  const double delProb3 = gapCosts.delProb3;
  const double insProb1 = gapCosts.insProb1;
  const double insProb2 = gapCosts.insProb2;
  const double insProb3 = gapCosts.insProb3;

  size_t diagPos6 = padSize - 1;
  size_t horiPos5 = diagPos6 + padSize;
  size_t horiPos4 = horiPos5 + padSize;
  size_t horiPos3 = horiPos4 + padSize;
  size_t vertPos2 = horiPos3 + padSize + 1;
  size_t vertPos1 = vertPos2 + padSize;
  size_t thisPos  = vertPos1;
  size_t antidiagonal = numOfAntidiagonals;

  double scaledUnit = 1 / rescaledSumOfProbRatios;

  resizeBckProbsIfSmaller(begInProbs(numOfAntidiagonals));

  while (1) {
    --antidiagonal;
    const uchar *s1 = proteinPtr;
    frames[antidiagonal % 3] -= seqIncrement;
    const uchar *s2 = frames[antidiagonal % 3];

    int numCells = numOfCellsInAntidiagonal(antidiagonal);
    for (int i = 0; i < padSize; ++i) {
      xBckProbs[thisPos + i] = 0;
      yBckProbs[thisPos + i] = 0;
      zBckProbs[thisPos + i] = 0;
    }
    thisPos += padSize;
    Prob *X0 = &xBckProbs[thisPos];
    Prob *Y0 = &yBckProbs[thisPos];
    Prob *Z0 = &zBckProbs[thisPos];
    const Prob *Z1 = &zBckProbs[vertPos1];
    const Prob *Z2 = &zBckProbs[vertPos2];
    const Prob *Y3 = &yBckProbs[horiPos3];
    const Prob *Z3 = &zBckProbs[horiPos3 + 1];
    const Prob *Y4 = &yBckProbs[horiPos4];
    const Prob *Y5 = &yBckProbs[horiPos5];
    const Prob *X6 = &xBckProbs[diagPos6];

    for (int i = 0; i < numCells; ++i) {
      double x  = X6[i];
      double y1 = Y5[i] * delProb1;
      double y2 = Y4[i] * delProb2;
      double y3 = Y3[i] * delProb3;
      double z1 = Z1[i] * insProb1;
      double z2 = Z2[i] * insProb2;
      double z3 = Z3[i] * insProb3;
      double b  = x + y1 + y2 + y3 + z1 + z2 + z3 + scaledUnit;
      X0[i] = b * substitutionProbs[*s1][*s2];
      Y0[i] = b * delOpenProb + y3;
      Z0[i] = b * insOpenProb + z3;
      s1 -= seqIncrement;
      s2 += seqIncrement;
    }

    diagPos6 = horiPos5;
    horiPos5 = horiPos4;
    horiPos4 = horiPos3;
    horiPos3 = vertPos2 - 1;
    vertPos2 = vertPos1;
    vertPos1 = thisPos;
    thisPos += numCells;

    if (antidiagonal == 0) break;

    if (endInProtein(antidiagonal - 1) < endInProtein(antidiagonal)) {
      proteinPtr -= seqIncrement;
      frames[0] += seqIncrement;
      frames[1] += seqIncrement;
      frames[2] += seqIncrement;
      ++diagPos6; ++horiPos5; ++horiPos4; ++horiPos3; ++vertPos2; ++vertPos1;
    }

    size_t offset = 6;
    if ((antidiagonal + offset) % rescaleStep == 0 &&
	antidiagonal + offset < numOfAntidiagonals) {
      double scale = rescales[antidiagonal / rescaleStep];
      size_t length = begInProbs(antidiagonal + 6) - begInProbs(antidiagonal);
      rescaleBckProbs(thisPos - length, thisPos, scale);
      scaledUnit *= scale;
    }
  }

  assert(thisPos == begInProbs(numOfAntidiagonals));
}

void FrameshiftXdropAligner::count(bool isRightwardExtension,
				   const GapCosts &gapCosts,
				   const dbl_ptr *substitutionCounts,
				   double *transitionCounts) {
  const int seqIncrement = isRightwardExtension ? 1 : -1;

  size_t diagPos6 = padSize - 1;
  size_t horiPos5 = diagPos6 + padSize;
  size_t horiPos4 = horiPos5 + padSize;
  size_t horiPos3 = horiPos4 + padSize;
  size_t vertPos2 = horiPos3 + padSize + 1;
  size_t vertPos1 = vertPos2 + padSize;
  size_t thisPos  = vertPos1;
  size_t antidiagonal = numOfAntidiagonals;
  size_t totalNumOfCells = begInProbs(numOfAntidiagonals);

  size_t proteinEnd = endInProtein(antidiagonal - 1);
  if (isRightwardExtension) {
    proteinPtr += proteinEnd;
    frames[0] += (antidiagonal + 2) / 3 - (proteinEnd - 1);
    frames[1] += (antidiagonal + 1) / 3 - (proteinEnd - 1);
    frames[2] += (antidiagonal + 0) / 3 - (proteinEnd - 1);
  } else {
    proteinPtr -= proteinEnd;
    frames[0] -= (antidiagonal + 2) / 3 - (proteinEnd - 1);
    frames[1] -= (antidiagonal + 1) / 3 - (proteinEnd - 1);
    frames[2] -= (antidiagonal + 0) / 3 - (proteinEnd - 1);
  }

  double del1 = 0; double del2 = 0; double del3 = 0; double delNext = 0;
  double ins1 = 0; double ins2 = 0; double ins3 = 0; double insNext = 0;

  while (1) {
    --antidiagonal;
    const uchar *s1 = proteinPtr;
    const uchar *s2 = frames[antidiagonal % 3];
    frames[antidiagonal % 3] -= seqIncrement;

    int numCells = numOfCellsInAntidiagonal(antidiagonal);
    thisPos += padSize;
    const Prob *Z1 = &zBckProbs[vertPos1];
    const Prob *Z2 = &zBckProbs[vertPos2];
    const Prob *Y3 = &yBckProbs[horiPos3];
    const Prob *Z3 = &zBckProbs[horiPos3 + 1];
    const Prob *Y4 = &yBckProbs[horiPos4];
    const Prob *Y5 = &yBckProbs[horiPos5];
    const Prob *X6 = &xBckProbs[diagPos6];

    const Prob *fX0 = &xFwdProbs[totalNumOfCells - thisPos + padSize*7 - 1];
    const Prob *fY0 = &yFwdProbs[totalNumOfCells - thisPos + padSize*7 - 1];
    const Prob *fZ0 = &zFwdProbs[totalNumOfCells - thisPos + padSize*7 - 1];

    double d1 = 0; double d2 = 0; double d3 = 0; double dNext = 0;
    double i1 = 0; double i2 = 0; double i3 = 0; double iNext = 0;

    for (int i = 0; i < numCells; ++i) {
      double matchProb = X6[i] * fX0[-i];
      substitutionCounts[*s1][*s2] += matchProb;
      transitionCounts[0] += matchProb;
      d1 += Y5[i] * fX0[-i];
      d2 += Y4[i] * fX0[-i];
      d3 += Y3[i] * fX0[-i];
      dNext += Y3[i] * fY0[-i];
      i1 += Z1[i] * fX0[-i];
      i2 += Z2[i] * fX0[-i];
      i3 += Z3[i] * fX0[-i];
      iNext += Z3[i] * fZ0[-i];
      s1 -= seqIncrement;
      s2 += seqIncrement;
    }

    size_t mod = antidiagonal % rescaleStep;
    if (mod >= rescaleStep - 6 &&
	antidiagonal - mod + rescaleStep < numOfAntidiagonals) {
      double mul = 1 / rescales[antidiagonal / rescaleStep];
      if (mod < rescaleStep - 1) i1 *= mul;
      if (mod < rescaleStep - 2) i2 *= mul;
      if (mod < rescaleStep - 3) {
	i3 *= mul; iNext *= mul; d3 *= mul; dNext *= mul;
      }
      if (mod < rescaleStep - 4) d2 *= mul;
      if (mod < rescaleStep - 5) d1 *= mul;
    }

    ins1 += i1; ins2 += i2; ins3 += i3; insNext += iNext;
    del1 += d1; del2 += d2; del3 += d3; delNext += dNext;

    diagPos6 = horiPos5;
    horiPos5 = horiPos4;
    horiPos4 = horiPos3;
    horiPos3 = vertPos2 - 1;
    vertPos2 = vertPos1;
    vertPos1 = thisPos;
    thisPos += numCells;

    if (antidiagonal == 0) break;

    if (endInProtein(antidiagonal - 1) < endInProtein(antidiagonal)) {
      proteinPtr -= seqIncrement;
      frames[0] += seqIncrement;
      frames[1] += seqIncrement;
      frames[2] += seqIncrement;
      ++diagPos6; ++horiPos5; ++horiPos4; ++horiPos3; ++vertPos2; ++vertPos1;
    }
  }

  transitionCounts[1] += gapCosts.delProb3 * delNext;  // deleted whole codons
  transitionCounts[2] += gapCosts.insProb3 * insNext;  // inserted whole codons
  transitionCounts[3] += gapCosts.delProb3 * del3;  // in-frame opens/closes
  transitionCounts[4] += gapCosts.insProb3 * ins3;  // in-frame opens/closes

  transitionCounts[5] += gapCosts.delProb1 * del1;
  transitionCounts[6] += gapCosts.delProb2 * del2;
  transitionCounts[7] += gapCosts.insProb1 * ins1;
  transitionCounts[8] += gapCosts.insProb2 * ins2;
}

double FrameshiftXdropAligner::maxSumOfProbRatios(const uchar *protein,
						  int proteinLength,
						  const uchar *tranDna,
						  int origDnaLength,
						  const const_dbl_ptr *substitutionProbs,
						  const GapCosts &gapCosts) {
  const double delOpenProb = gapCosts.delProbPieces[0].openProb;
  const double insOpenProb = gapCosts.insProbPieces[0].openProb;
  const double delProb1 = gapCosts.delProb1;
  const double delProb2 = gapCosts.delProb2;
  const double delProb3 = gapCosts.delProb3;
  const double insProb1 = gapCosts.insProb1;
  const double insProb2 = gapCosts.insProb2;
  const double insProb3 = gapCosts.insProb3;

  const int proteinSize = proteinLength + 1;
  const int origDnaSize = origDnaLength + 1;
  xFwdProbs.resize(origDnaSize + proteinSize * origDnaSize);
  double *delRow = &xFwdProbs[0];
  double *newRow = delRow + origDnaSize;
  double Y1, Y2, Y3, Z1, Z2, Z3;
  double maxValue = 0;

  std::fill_n(delRow, origDnaSize, 0.0);

  for (int i = 0; i < proteinSize; ++i) {
    bool isIn = (i > 0);
    const double *substitutionRow = isIn ? substitutionProbs[protein[i-1]] : 0;
    const double *oldRow = newRow - origDnaSize;
    Y2 = Y3 = Z1 = Z2 = Z3 = 0;

    for (int j = 0; j < origDnaSize; ++j) {
      Y1 = Y2;
      Y2 = Y3;
      Y3 = delRow[j];
      bool isMain = (isIn && j > 2);
      double x  = isMain ? oldRow[j-3] * substitutionRow[tranDna[j-3]] : 0;
      double y1 = Y1 * delProb1;
      double y2 = Y2 * delProb2;
      double y3 = Y3 * delProb3;
      double z1 = Z1 * insProb1;
      double z2 = Z2 * insProb2;
      double z3 = Z3 * insProb3;
      double b  = x + y1 + y2 + y3 + z1 + z2 + z3 + 1;
      newRow[j] = b;
      delRow[j] = b * delOpenProb + y3;
      Z3 = Z2;
      Z2 = Z1;
      Z1 = b * insOpenProb + z3;
    }
    newRow += origDnaSize;
  }

  std::fill_n(delRow, origDnaSize, 0.0);

  for (int i = proteinSize; i-- > 0;) {
    newRow -= origDnaSize;
    bool isIn = (i < proteinLength);
    const double *substitutionRow = isIn ? substitutionProbs[protein[i]] : 0;
    const double *oldRow = newRow + origDnaSize;
    Y2 = Y3 = Z1 = Z2 = Z3 = 0;

    for (int j = origDnaSize; j-- > 0;) {
      Y1 = Y2;
      Y2 = Y3;
      Y3 = delRow[j];
      bool isMain = (isIn && j < origDnaLength - 2);
      double x  = isMain ? oldRow[j+3] * substitutionRow[tranDna[j]] : 0;
      double y1 = Y1 * delProb1;
      double y2 = Y2 * delProb2;
      double y3 = Y3 * delProb3;
      double z1 = Z1 * insProb1;
      double z2 = Z2 * insProb2;
      double z3 = Z3 * insProb3;
      double b  = x + y1 + y2 + y3 + z1 + z2 + z3 + 1;
      maxValue  = std::max(maxValue, b * newRow[j]);
      newRow[j] = b;
      delRow[j] = b * delOpenProb + y3;
      Z3 = Z2;
      Z2 = Z1;
      Z1 = b * insOpenProb + z3;
    }
  }

  return maxValue;
}

}
