// Copyright 2015 Martin C. Frith

#include "TantanMasker.hh"
#include "ScoreMatrix.hh"
#include "mcf_substitution_matrix_stats.hh"
#include <algorithm>
#include <cassert>
#include <cmath>

namespace cbrc {

void TantanMasker::init(bool isProtein,
			bool isATrichDna,
			bool isMaskWeakRepeats,
			const std::string &alphabet,
			const uchar *letterToIndex) {
  maxRepeatOffset = 100;
  repeatProb = 0.005;

  std::copy(probMatrix, probMatrix + scoreMatrixRowSize, probMatrixPointers);

  unsigned alphabetSize;
  ScoreMatrix s;
  const char *matrixName = "";

  if (isATrichDna) {
    repeatProb = 0.01;
    alphabetSize = 4;
    s.fromString("\
   A  C  G  T\n\
A  2 -5 -5 -5\n\
C -5  5 -5 -5\n\
G -5 -5  5 -5\n\
T -5 -5 -5  2\n\
");
  } else if (isProtein) {
    maxRepeatOffset = 50;
    alphabetSize = 20;
    matrixName = "BL62";
    s.fromString(ScoreMatrix::stringFromName(matrixName));
  } else {
    alphabetSize = alphabet.size();
    s.setMatchMismatch(1, 1, alphabet);
  }

  if (isMaskWeakRepeats) repeatProb = 0.02;

  s.init(letterToIndex);

  int *scoreMat[scoreMatrixRowSize];
  std::copy(s.caseInsensitive, s.caseInsensitive + alphabetSize, scoreMat);

  mcf::SubstitutionMatrixStats stats;
  stats.calcUnbiased(matrixName, scoreMat, alphabetSize);
  assert(!stats.isBad());

  for (int i = 0; i < scoreMatrixRowSize; ++i)
    for (int j = 0; j < scoreMatrixRowSize; ++j)
      probMatrix[i][j] = std::exp(stats.lambda() * s.caseInsensitive[i][j]);
}

}
