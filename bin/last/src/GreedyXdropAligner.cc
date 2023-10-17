// Copyright 2016 Martin C. Frith

#include "GreedyXdropAligner.hh"
#include <algorithm>
#include <cassert>
//#include <iostream>  // for debugging

static int maxValue(int a, int b, int c) {
  return std::max(std::max(a, b), c);
}

static void resizeIfSmaller(std::vector<int> &v, size_t s) {
  if (v.size() < s) v.resize(s);
}

const int undefined = -9;

static const int *definedBeg(const int *beg, const int *end) {
  while (beg < end && *beg == undefined) ++beg;
  return beg;
}

static const int *definedEnd(const int *beg, const int *end) {
  while (end > beg && *(end - 1) == undefined) --end;
  return end;
}

namespace cbrc {

// i:            number of seq1 letters aligned so far
// j:            number of seq2 letters aligned so far
// diagonal:     i - j
// antidiagonal: i + j
// distance:     number of differences (mismatches + gap characters)

// furthest:     a flattened, ragged 2D array, which holds the
//               furthest antidiagonal reachable on a given diagonal
//               with a given distance

// score0:       score + differences * (matchScore - mismatchScore),
//               or equivalently, antidiagonal * matchScore / 2

int GreedyXdropAligner::align(const uchar *seq1,
			      const uchar *seq2,
			      bool isForward,
			      const ScoreMatrixRow *scorer,
			      int maxScoreDrop,
			      uchar delimiter) {
  const int matchScore = scorer[0][0];
  const int mismatchScore = scorer[0][1];
  assert(matchScore % 2 == 0);
  const int halfMatchScore = matchScore / 2;
  const int differenceCost = matchScore - mismatchScore;
  const int lookBack = (maxScoreDrop + halfMatchScore) / differenceCost + 1;
  const int minScore0gain = lookBack * differenceCost - maxScoreDrop;

  minScore0s.assign(lookBack, 0);

  rowOrigins.resize(1);
  size_t furthestSize = 2;
  resizeIfSmaller(furthest, furthestSize);
  furthest[0] = -1;  // xxx tricky initialization
  furthest[1] = undefined;  // add one pad cell

  int lowerStop = INT_MIN;
  int upperStop = INT_MAX;

  int diagonalBeg = 0;
  int diagonalEnd = 1;

  bestDistance = -1;
  int bestScore0 = 0;

  int distance;
  for (distance = 0; ; ++distance) {
    int minScore0 = minScore0s[distance];

    size_t rowOrigin = furthestSize - diagonalBeg;
    rowOrigins.push_back(rowOrigin);
    size_t furthestSizeNew = rowOrigin + diagonalEnd + 2;  // + 2 pad cells
    resizeIfSmaller(furthest, furthestSizeNew);

    int *to = &furthest[furthestSize];
    const int *from = sources(distance, diagonalBeg);

    *to++ = undefined;  // add one pad cell
    const int *toBeg = to;

    for (int diagonal = diagonalBeg; diagonal < diagonalEnd; ++diagonal) {
      int antidiagonal = maxValue(from[0], from[1] + 1, from[2]) + 1;
      ++from;
      if (halfMatchScore * antidiagonal >= minScore0) {
	int i = (antidiagonal + diagonal) / 2;
	int j = (antidiagonal - diagonal) / 2;
	const uchar *s1;
	const uchar *s2;
	if (isForward) {
	  s1 = seq1 + i;
	  s2 = seq2 + j;
	  while (scorer[*s1][*s2] > 0) { ++s1; ++s2; }  // skip past matches
	  i = s1 - seq1;
	  j = s2 - seq2;
	} else {
	  s1 = seq1 - i;
	  s2 = seq2 - j;
	  while (scorer[*s1][*s2] > 0) { --s1; --s2; }  // skip past matches
	  i = seq1 - s1;
	  j = seq2 - s2;
	}
	if (*s2 == delimiter)                         lowerStop = diagonal;
	if (*s1 == delimiter && upperStop > diagonal) upperStop = diagonal;
	antidiagonal = i + j;
	*to = antidiagonal;
	int score0 = halfMatchScore * antidiagonal;
	if (score0 > bestScore0) {
	  bestScore0 = score0;
	  bestDistance = distance;
	  bestDiagonal = diagonal;
	}
      } else {
	*to = undefined;
      }
      ++to;
    }

    lowerStop += 2;  // unnecessary, but might improve speed
    upperStop -= 2;  // unnecessary, but might improve speed

    diagonalBeg += definedBeg(toBeg, to) - toBeg - 1;
    diagonalBeg = std::max(diagonalBeg, lowerStop + 1);

    diagonalEnd += definedEnd(toBeg, to) - to + 1;
    diagonalEnd = std::min(diagonalEnd, upperStop);

    if (diagonalBeg >= diagonalEnd) break;

    *to = undefined;  // add one pad cell
    minScore0s.push_back(bestScore0 + minScore0gain);
    bestScore0 += differenceCost;
    furthestSize = furthestSizeNew;
  }

  return bestScore0 - differenceCost * distance;
}

bool GreedyXdropAligner::getNextChunk(size_t &end1,
				      size_t &end2,
				      size_t &length) {
  if (bestDistance < 0) return false;

  int antidiagonal = furthest[rowOrigins[bestDistance + 1] + bestDiagonal + 1];
  end1 = (antidiagonal + bestDiagonal) / 2;
  end2 = (antidiagonal - bestDiagonal) / 2;

  // skip back past substitutions, until we hit an indel or the start
  // of the extension:
  int del, mis, ins;
  do {
    const int *from = sources(bestDistance, bestDiagonal);
    del = from[0];
    mis = from[1] + 1;
    ins = from[2];
    bestDistance -= 1;
  } while (mis >= del && mis >= ins);

  int newAntidiagonal;
  if (del >= ins) {
    bestDiagonal -= 1;
    newAntidiagonal = del;
  } else {
    bestDiagonal += 1;
    newAntidiagonal = ins;
  }
  length = (antidiagonal - newAntidiagonal) / 2;  // odd number / 2

  // skip back past indels, until we hit a gapless chunk or the start
  // of the extension:
  while (bestDistance >= 0) {
    const int *from = sources(bestDistance, bestDiagonal);
    del = from[0];
    mis = from[1] + 1;
    ins = from[2];
    if (mis >= del && mis >= ins) break;
    --newAntidiagonal;
    if (del >= ins) {
      if (del < newAntidiagonal) break;
      bestDiagonal -= 1;
    } else {
      if (ins < newAntidiagonal) break;
      bestDiagonal += 1;
    }
    bestDistance -= 1;
  }

  return true;
}

}
