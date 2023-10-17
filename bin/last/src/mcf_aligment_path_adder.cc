// Author: Martin C. Frith 2021
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mcf_aligment_path_adder.hh"

#include <algorithm>

namespace mcf {

double AlignmentPathAdder::maxSum(const uchar *seq1, int len1,
				  const uchar *seq2, int len2,
				  const const_dbl_ptr *substitutionProbs,
				  double delInitProb, double delNextProb,
				  double insInitProb, double insNextProb) {
  int size1 = len1 + 1;
  int size2 = len2 + 1;
  dynamicProgrammingValues.resize(size2 + size1 * size2);
  double *delRow = &dynamicProgrammingValues[0];
  double *newRow = delRow + size2;
  double maxValue = 0;

  std::fill_n(delRow, size2, 0.0);

  for (int i = 0; i < size1; ++i) {
    bool isIn = (i > 0);
    const double *substitutionRow = isIn ? substitutionProbs[seq1[i-1]] : 0;
    const double *oldRow = newRow - size2;
    double ins = 0;

    for (int j = 0; j < size2; ++j) {
      bool isMain = (isIn && j > 0);
      double mat = isMain ? oldRow[j-1] * substitutionRow[seq2[j-1]] : 0;
      double del = delRow[j];
      double all = mat + del + ins + 1;
      newRow[j]  = all;
      delRow[j]  = all * delInitProb + del * delNextProb;
      ins        = all * insInitProb + ins * insNextProb;
    }
    newRow += size2;
  }

  std::fill_n(delRow, size2, 0.0);

  for (int i = size1; i-- > 0;) {
    newRow -= size2;
    bool isIn = (i < len1);
    const double *substitutionRow = isIn ? substitutionProbs[seq1[i]] : 0;
    const double *oldRow = newRow + size2;
    double ins = 0;

    for (int j = size2; j-- > 0;) {
      bool isMain = (isIn && j < len2);
      double mat = isMain ? oldRow[j+1] * substitutionRow[seq2[j]] : 0;
      double del = delRow[j];
      double all = mat + del + ins + 1;
      maxValue   = std::max(maxValue, all * newRow[j]);
      newRow[j]  = all;
      delRow[j]  = all * delInitProb + del * delNextProb;
      ins        = all * insInitProb + ins * insNextProb;
    }
  }

  return maxValue;
}

}
