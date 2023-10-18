// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mcf_substitution_matrix_stats.hh"
#include "LambdaCalculator.hh"
#include "cbrc_linalg.hh"
#include <algorithm>
#include <stdexcept>
#include <assert.h>
#include <math.h>
#include <string.h>  // strcmp

#define COUNTOF(a) (sizeof (a) / sizeof *(a))

const struct {
  const char *name;
  double scale;
} substitutionMatrixScales[] = {
  {"BL62", 3.0861133577701141},
  {"BL80", 2.8253295934982616},
  {"MIQS", 0.0},
  {"PAM10", 0.0},
  {"PAM30", 2.8848596855435313},
};

static double checkedExp(double lambda, int score) {
  double y = exp(lambda * score);
  if (!(y < HUGE_VAL)) throw std::overflow_error("exp overflow");
  return y;
}

static void permuteComplement(std::vector<double> &v,
			      const unsigned char *complement) {
  for (unsigned i = 0; i < v.size(); ++i) {
    unsigned j = complement[i];
    if (j < i) std::swap(v[i], v[j]);
  }
}

static double calcLetterProbs(std::vector<double> &probs, unsigned size,
			      double **expMat) {
  probs.assign(size, 1.0);
  cbrc::linalgSolve(expMat, &probs[0], size);
  double sum = 0;
  for (unsigned i = 0; i < size; ++i) {
    if (probs[i] < 0) {
      throw std::runtime_error("got a probability < 0 for the "
			       "substitution matrix");
    }
    sum += probs[i];
  }
  assert(sum > 0);
  double bias = 1 / sum;
  for (unsigned i = 0; i < size; ++i) {
    probs[i] = cbrc::roundToFewDigits(probs[i] * bias);
  }
  return bias;
}

namespace mcf {

void SubstitutionMatrixStats::calcFromScale(const const_int_ptr *scoreMatrix,
					    unsigned size, double scale) {
  assert(size > 0);
  assert(scale > 0);
  mLambda = 1 / scale;

  std::vector<double> expVec(size * size);
  std::vector<double *> expMat(size);
  for (unsigned i = 0; i < size; ++i)
    expMat[i] = &expVec[i * size];

  for (unsigned i = 0; i < size; ++i)
    for (unsigned j = 0; j < size; ++j)
      expMat[i][j] = checkedExp(mLambda, scoreMatrix[j][i]);
  double bias1 = calcLetterProbs(mLetterProbs1, size, &expMat[0]);

  for (unsigned i = 0; i < size; ++i)
    for (unsigned j = 0; j < size; ++j)
      expMat[i][j] = checkedExp(mLambda, scoreMatrix[i][j]);
  double bias2 = calcLetterProbs(mLetterProbs2, size, &expMat[0]);

  // bias1 and bias2 should be equal
  mBias = (bias1 + bias2) / 2;
}

void SubstitutionMatrixStats::calcBias(const const_int_ptr *scoreMatrix,
				       unsigned numOfRows, unsigned numOfCols,
				       double scale) {
  assert(scale > 0);
  mLambda = 1 / scale;

  mBias = 0;
  for (unsigned i = 0; i < numOfRows; ++i) {
    for (unsigned j = 0; j < numOfCols; ++j) {
      double e = checkedExp(mLambda, scoreMatrix[i][j]);
      mBias += mLetterProbs1[i] * mLetterProbs2[j] * e;
    }
  }
}

void SubstitutionMatrixStats::calcUnbiased(const char *matrixName,
					   const const_int_ptr *scoreMatrix,
					   unsigned size) {
  double scale = 0;
  unsigned realSize = size;

  for (size_t i = 0; i < COUNTOF(substitutionMatrixScales); ++i) {
    if (strcmp(matrixName, substitutionMatrixScales[i].name) == 0) {
      scale = substitutionMatrixScales[i].scale;
      if (size > 20) realSize = 20;
    }
  }

  if (scale > 0) {
    calcFromScale(scoreMatrix, realSize, scale);
  } else {
    cbrc::LambdaCalculator c;
    c.calculate(scoreMatrix, realSize);
    mLambda = c.lambda();
    if (isBad()) return;
    mBias = 1;
    mLetterProbs1.assign(c.letterProbs1(), c.letterProbs1() + realSize);
    mLetterProbs2.assign(c.letterProbs2(), c.letterProbs2() + realSize);
  }

  mLetterProbs1.insert(mLetterProbs1.end(), size - realSize, 0.0);
  mLetterProbs2.insert(mLetterProbs2.end(), size - realSize, 0.0);
}

void SubstitutionMatrixStats::flipDnaStrands(const unsigned char *complement) {
  permuteComplement(mLetterProbs1, complement);
  permuteComplement(mLetterProbs2, complement);
}

}
