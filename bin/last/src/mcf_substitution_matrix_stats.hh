// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

// A substitution score matrix S(x,y) specifies a score for aligning a
// letter of type x in the 1st sequence to a letter of type y in the
// 2nd sequence.  It is related to probabilities like this:

// S(x,y) = scale * ln{ bias * A(x,y) / [P(x) * Q(y)] }

// Terminology: scale = "temperature" = 1 / lambda

// This class calculates scale, bias, P and Q.
// The calculations may fail, putting it into a "bad" state.

// It assumes "homogeneous letter probabilities", i.e.
// P(x) = sum_y A(x,y)       Q(y) = sum_x A(x,y)

// Input: a (size x size) score matrix

#ifndef MCF_SUBSTITUTION_MATRIX_STATS_HH
#define MCF_SUBSTITUTION_MATRIX_STATS_HH

#include <vector>

namespace mcf {

typedef const int *const_int_ptr;

class SubstitutionMatrixStats {
public:
  SubstitutionMatrixStats() { mLambda = -1; }

  bool isBad() const { return mLambda < 0; }

  double lambda() const { return mLambda; }

  double bias() const { return mBias; }

  // probabilities of the letters that correspond to matrix rows
  const double *letterProbs1() const
  { return mLetterProbs1.empty() ? 0 : &mLetterProbs1[0]; }

  // probabilities of the letters that correspond to matrix columns
  const double *letterProbs2() const
  { return mLetterProbs2.empty() ? 0 : &mLetterProbs2[0]; }

  // set up an array for writing row-letter probabilities into
  double *sizedLetterProbs1(unsigned size)
  { mLetterProbs1.resize(size); return size ? &mLetterProbs1[0] : 0; }

  // set up an array for writing column-letter probabilities into
  double *sizedLetterProbs2(unsigned size)
  { mLetterProbs2.resize(size); return size ? &mLetterProbs2[0] : 0; }

  // Calculate the bias, given the scale and the letter probabilities
  // (which need not be homogeneous).  A non-square matrix is allowed.
  void calcBias(const const_int_ptr *scoreMatrix,
		unsigned numOfRows, unsigned numOfCols, double scale);

  // set scale/lambda and calculate the other values
  void calcFromScale(const const_int_ptr *scoreMatrix, unsigned size,
		     double scale);

  // Set bias = 1 and calculate the other values.
  // matrixName is used only to lookup pre-calculated cases.
  void calcUnbiased(const char *matrixName,
		    const const_int_ptr *scoreMatrix, unsigned size);

  // Set it to the parameters for aligning the reverse DNA strands.
  // j = complement[i] means the j-th letter is the complement of the
  // i-th letter.
  void flipDnaStrands(const unsigned char *complement);

private:
  double mLambda;
  double mBias;
  std::vector<double> mLetterProbs1;
  std::vector<double> mLetterProbs2;
};

}

#endif
