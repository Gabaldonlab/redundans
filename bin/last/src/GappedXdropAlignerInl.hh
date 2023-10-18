// Copyright 2011, 2013, 2014 Martin C. Frith

#ifndef GAPPED_XDROP_ALIGNER_INL_HH
#define GAPPED_XDROP_ALIGNER_INL_HH

#include <algorithm>
#include <cassert>
//#include <stdexcept>

namespace cbrc {

template<typename T, int N> T arrayMax(T (&array)[N]) {
  return *std::max_element(array, array + N);
}

template<typename T, int N> T arrayMin(T (&array)[N]) {
  return *std::min_element(array, array + N);
}

template<typename T> int maxIndex(T a, T b) {
  return b > a ? 1 : 0;
}

template<typename T> int maxIndex(T a, T b, T c) {
  return c > a ? maxIndex(b, c) + 1 : maxIndex(a, b);
}

template<typename T> int maxIndex(T a, T b, T c, T d) {
  return d > a ? maxIndex(b, c, d) + 1 : maxIndex(a, b, c);
}

template<typename T> int maxIndex(T a, T b, T c, T d, T e) {
  return e > a ? maxIndex(b, c, d, e) + 1 : maxIndex(a, b, c, d);
}

template<typename T> int maxIndex(T a, T b, T c, T d, T e, T f) {
  return f > a ? maxIndex(b, c, d, e, f) + 1 : maxIndex(a, b, c, d, e);
}

template<typename T> int maxIndex(T a, T b, T c, T d, T e, T f, T g) {
  return g > a ? maxIndex(b, c, d, e, f, g) + 1 : maxIndex(a, b, c, d, e, f);
}

static inline int maxValue(int a, int b) {
  return std::max(a, b);
}

static inline int maxValue(int a, int b, int c) {
  return maxValue(maxValue(a, b), c);
}

static inline int maxValue(int a, int b, int c, int d, int e, int f, int g) {
  return maxValue(maxValue(maxValue(a, b), maxValue(c, d)),
		  maxValue(maxValue(e, f), g));
}

template<typename T>
T whichFrame(size_t antidiagonal, T frame0, T frame1, T frame2) {
  switch (antidiagonal % 3) {
    case 0: return frame1;  // the +1 frame
    case 1: return frame2;  // the -1 frame
    case 2: return frame0;
    default: assert(0);  // keeps my compiler happy
  }
}

// The next two functions will stop the alignment at delimiters.  But
// this is not guaranteed if bestScore > INF / 2.  We could avoid this
// restriction by replacing -INF / 2 with bestScore - INF.

inline const Score *finiteBeg(const Score *beg, const Score *end) {
  while (beg < end && *beg <= -INF / 2)
    ++beg;
  return beg;
}

inline const Score *finiteEnd(const Score *beg, const Score *end) {
  while (end > beg && *(end-1) <= -INF / 2)
    --end;
  return end;
}

inline bool isDelimiter(uchar c, const int *scores) {
  return scores[c] <= -INF;
}

/*
inline void checkGappedXdropScore(int bestScore) {
  // If this happens, sentinels/delimiters might not work:
  if (bestScore > INF / 2)
    throw std::overflow_error("score got too high in gapped extension");
}
*/

inline void updateBest1(int &bestScore,
			size_t &bestAntidiagonal,
			size_t &bestSeq1position,
			int minScore,
			int score,
			size_t antidiagonal,
			size_t seq1position) {
  if (score >= minScore && score > bestScore) {
    bestScore = score;
    bestAntidiagonal = antidiagonal;
    bestSeq1position = seq1position;
  }
}

inline void GappedXdropAligner::updateBest(int &bestScore, int score,
                                           size_t antidiagonal,
                                           const Score *x0,
					   const Score *x0base) {
  if (score > bestScore) {
    bestScore = score;
    bestAntidiagonal = antidiagonal;
    bestSeq1position = static_cast<size_t>(x0 - x0base);
  }
}

inline void updateMaxScoreDrop(int &maxScoreDrop,
                               int maxMatches, int maxMatchScore) {
  // If the current antidiagonal touches a sentinel/delimiter, then
  // maxMatches is the maximum possible number of matches starting
  // from the next antidiagonal.
  maxScoreDrop = std::min(maxScoreDrop, maxMatches * maxMatchScore - 1);
}

inline void updateFiniteEdges3(size_t *maxSeq1begs, size_t *minSeq1ends,
                               const Score *x0base, const Score *x0end,
                               size_t numCells) {
  const Score *x0beg = x0end - numCells;

  maxSeq1begs[0] = maxSeq1begs[1];
  maxSeq1begs[1] = maxSeq1begs[2];
  maxSeq1begs[2] = maxSeq1begs[3];
  maxSeq1begs[3] = maxSeq1begs[4] + 1;
  maxSeq1begs[4] = maxSeq1begs[5];
  maxSeq1begs[5] = maxSeq1begs[6];
  maxSeq1begs[6] = finiteBeg(x0beg, x0end) - x0base;

  minSeq1ends[0] = minSeq1ends[1];
  minSeq1ends[1] = minSeq1ends[2];
  minSeq1ends[2] = minSeq1ends[3];
  minSeq1ends[3] = minSeq1ends[4];
  minSeq1ends[4] = minSeq1ends[5] + 1;
  minSeq1ends[5] = minSeq1ends[6];
  minSeq1ends[6] = finiteEnd(x0beg, x0end) - x0base;
}

}

#endif
