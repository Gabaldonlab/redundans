// Copyright 2011 Martin C. Frith

#include "TwoQualityScoreMatrix.hh"

#include "qualityScoreUtil.hh"

#include <algorithm>  // min
#include <cassert>
#include <cmath>
//#include <iostream>  // for debugging

namespace cbrc {

void TwoQualityMatrixIndexer::init(const uchar *toUnmasked) {
  indexMap.resize(qualityCapacity * scoreMatrixRowSize);

  for (int quality = 0; quality < qualityCapacity; ++quality) {
    int normalStart = quality * numQualityLetters;
    int maskedStart = normalStart + numNormalLetters;
    int abnormalPos = 0;

    for (int letter = 0; letter < scoreMatrixRowSize; ++letter) {
      int unmasked = toUnmasked[letter];
      int i = subindex(letter, quality);

      if (letter < numNormalLetters)
        indexMap[i] = normalStart + letter;
      else if (unmasked < numNormalLetters)
        indexMap[i] = maskedStart + unmasked;
      else
        indexMap[i] = abnormalPos++;
    }

    assert(abnormalPos <= minQuality * numQualityLetters);
  }
}

static int qualityEnd(const TwoQualityMatrixIndexer &indexer, int letter) {
  if (indexer(0, letter, 0, 0) == indexer(0, letter, 0, 1))
    return indexer.minQuality + 1;
  else
    return indexer.qualityCapacity;
}

void TwoQualityScoreMatrix::init(const ScoreMatrixRow *scoreMatrix,
				 const mcf::SubstitutionMatrixStats &matStats,
                                 bool isPhred1,
                                 int qualityOffset1,
                                 bool isPhred2,
                                 int qualityOffset2,
                                 const uchar *toUnmasked,
                                 bool isMask,
				 bool isMatchMismatchMatrix) {
  typedef TwoQualityMatrixIndexer Indexer;

  indexer.init(toUnmasked);
  data.resize(Indexer::numSymbols * Indexer::numSymbols);

  const double lambda = matStats.lambda();
  const double matBias = matStats.bias();
  double expMat[Indexer::numNormalLetters][Indexer::numNormalLetters];

  for (int x1 = 0; x1 < Indexer::numNormalLetters; ++x1)
    for (int x2 = 0; x2 < Indexer::numNormalLetters; ++x2)
      expMat[x1][x2] = std::exp(lambda * scoreMatrix[x1][x2]);

  double certainties1[Indexer::qualityCapacity][Indexer::numNormalLetters];
  double certainties2[Indexer::qualityCapacity][Indexer::numNormalLetters];

  for (int q = Indexer::minQuality; q < Indexer::qualityCapacity; ++q) {
    double e1 = errorProbFromQual(q, qualityOffset1, isPhred1);
    double e2 = errorProbFromQual(q, qualityOffset2, isPhred2);
    for (int x = 0; x < Indexer::numNormalLetters; ++x) {
      certainties1[q][x] = qualityCertainty(e1, matStats.letterProbs1()[x]);
      certainties2[q][x] = qualityCertainty(e2, matStats.letterProbs2()[x]);
    }
  }

  // I tried pre-calculating 1/lambda, but there was little speed boost.

  for (int q1 = Indexer::minQuality; q1 < Indexer::qualityCapacity; ++q1) {
    int *dq1 = &data[q1 * Indexer::numQualityLetters * Indexer::numSymbols];
    const double *c1s = certainties1[q1];
    for (int q2 = Indexer::minQuality; q2 < Indexer::qualityCapacity; ++q2) {
      int *dq2 = dq1 + q2 * Indexer::numQualityLetters;
      const double *c2s = certainties2[q2];
      if (isMatchMismatchMatrix) {  // do this common special case faster
	double c = c1s[0] * c2s[0];
	int scoreSame = qualityPairScore(expMat[0][0], matBias, c, lambda);
	int scoreDiff = qualityPairScore(expMat[0][1], matBias, c, lambda);
	int scoreSameMask = isMask ? std::min(scoreSame, 0) : scoreSame;
	int scoreDiffMask = isMask ? std::min(scoreDiff, 0) : scoreDiff;
	for (int x1 = 0; x1 < Indexer::numNormalLetters; ++x1) {
	  int *dx1 = dq2 + x1 * Indexer::numSymbols;
	  int *dm1 = dx1 + Indexer::numNormalLetters * Indexer::numSymbols;
	  for (int x2 = 0; x2 < Indexer::numNormalLetters; ++x2) {
	    int m2 = x2 + Indexer::numNormalLetters;
	    if (x2 == x1) {
	      dx1[x2] = scoreSame;
	      dx1[m2] = scoreSameMask;
	      dm1[x2] = scoreSameMask;
	      dm1[m2] = scoreSameMask;
	    } else {
	      dx1[x2] = scoreDiff;
	      dx1[m2] = scoreDiffMask;
	      dm1[x2] = scoreDiffMask;
	      dm1[m2] = scoreDiffMask;
	    }
	  }
	}
      } else {
	for (int x1 = 0; x1 < Indexer::numNormalLetters; ++x1) {
	  int *dx1 = dq2 + x1 * Indexer::numSymbols;
	  int *dm1 = dx1 + Indexer::numNormalLetters * Indexer::numSymbols;
	  double c1 = c1s[x1];
	  for (int x2 = 0; x2 < Indexer::numNormalLetters; ++x2) {
	    int m2 = x2 + Indexer::numNormalLetters;
	    double c2 = c2s[x2];
	    double c = c1 * c2;
	    int score = qualityPairScore(expMat[x1][x2], matBias, c, lambda);
	    dx1[x2] = score;
	    if (isMask) score = std::min(score, 0);
	    dx1[m2] = score;
	    dm1[x2] = score;
	    dm1[m2] = score;
	  }
	}
      }
    }
  }

  for (int x1 = 0; x1 < scoreMatrixRowSize; ++x1) {
    for (int x2 = 0; x2 < scoreMatrixRowSize; ++x2) {
      int unmasked1 = toUnmasked[x1];
      int unmasked2 = toUnmasked[x2];

      bool isQuality1 = (unmasked1 < Indexer::numNormalLetters);
      bool isQuality2 = (unmasked2 < Indexer::numNormalLetters);
      if (isQuality1 && isQuality2) continue;

      bool isMasked = (unmasked1 != x1 || unmasked2 != x2);

      int score = scoreMatrix[unmasked1][unmasked2];
      if (isMasked && isMask) score = std::min(score, 0);

      int end1 = qualityEnd(indexer, x1);
      int end2 = qualityEnd(indexer, x2);

      for (int q1 = Indexer::minQuality; q1 < end1; ++q1)
        for (int q2 = Indexer::minQuality; q2 < end2; ++q2)
          data[indexer(x1, x2, q1, q2)] = score;
    }
  }
}

void TwoQualityExpMatrix::init(const TwoQualityScoreMatrix &m,
                               double temperature) {
  assert(temperature > 0);
  indexer = m.indexer;
  data.resize(indexer.numSymbols * indexer.numSymbols);

  for (int i1 = 0; i1 < scoreMatrixRowSize; ++i1) {
    for (int i2 = 0; i2 < scoreMatrixRowSize; ++i2) {
      int end1 = qualityEnd(indexer, i1);
      int end2 = qualityEnd(indexer, i2);

      for (int q1 = indexer.minQuality; q1 < end1; ++q1)
        for (int q2 = indexer.minQuality; q2 < end2; ++q2)
          data[indexer(i1, i2, q1, q2)] =
              std::exp(m(i1, i2, q1, q2) / temperature);
    }
  }
}

}
