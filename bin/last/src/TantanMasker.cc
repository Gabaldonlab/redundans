// Copyright 2015 Martin C. Frith

#include "TantanMasker.hh"
#include "ScoreMatrix.hh"
#include "LambdaCalculator.hh"
#include <algorithm>
#include <cassert>
#include <cmath>

namespace cbrc {

void TantanMasker::init(bool isProtein,
			bool isATrichDna,
			const std::string &alphabet,
			const uchar *letterToIndex) {
  maxRepeatOffset = 100;
  repeatProb = 0.005;

  std::copy(probMatrix, probMatrix + scoreMatrixRowSize, probMatrixPointers);

  unsigned alphabetSize;
  ScoreMatrix s;

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
    s.fromString(ScoreMatrix::stringFromName("BL62"));
  } else {
    alphabetSize = alphabet.size();
    s.matchMismatch(1, 1, alphabet);
  }

  s.init(letterToIndex);

  LambdaCalculator c;
  c.calculate(s.caseInsensitive, alphabetSize);
  assert(!c.isBad());

  for (int i = 0; i < scoreMatrixRowSize; ++i)
    for (int j = 0; j < scoreMatrixRowSize; ++j)
      probMatrix[i][j] = std::exp(c.lambda() * s.caseInsensitive[i][j]);
}

}
