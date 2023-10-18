// Copyright 2011, 2013 Martin C. Frith

#include "gaplessTwoQualityXdrop.hh"
#include "TwoQualityScoreMatrix.hh"
#include <stdexcept>

static void err(const char *s) { throw std::overflow_error(s); }

namespace cbrc {

int forwardGaplessTwoQualityXdropScore(const uchar *seq1,
                                       const uchar *qual1,
                                       const uchar *seq2,
                                       const uchar *qual2,
                                       const TwoQualityScoreMatrix &m,
                                       int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += m(*seq1++, *seq2++, *qual1++, *qual2++);  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    err("score overflow in forward gapless extension with qualities");
  return score;
}

int reverseGaplessTwoQualityXdropScore(const uchar *seq1,
                                       const uchar *qual1,
                                       const uchar *seq2,
                                       const uchar *qual2,
                                       const TwoQualityScoreMatrix &m,
                                       int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += m(*--seq1, *--seq2, *--qual1, *--qual2);  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    err("score overflow in reverse gapless extension with qualities");
  return score;
}

const uchar *forwardGaplessTwoQualityXdropEnd(const uchar *seq1,
                                              const uchar *qual1,
                                              const uchar *seq2,
                                              const uchar *qual2,
                                              const TwoQualityScoreMatrix &m,
                                              int score) {
  int s = 0;
  while (s < score) s += m(*seq1++, *seq2++, *qual1++, *qual2++);
  return seq1;
}

const uchar *reverseGaplessTwoQualityXdropEnd(const uchar *seq1,
                                              const uchar *qual1,
                                              const uchar *seq2,
                                              const uchar *qual2,
                                              const TwoQualityScoreMatrix &m,
                                              int score) {
  int s = 0;
  while (s < score) s += m(*--seq1, *--seq2, *--qual1, *--qual2);
  return seq1;
}

bool isOptimalGaplessTwoQualityXdrop(const uchar *seq1,
                                     const uchar *seq1end,
                                     const uchar *qual1,
                                     const uchar *seq2,
                                     const uchar *qual2,
                                     const TwoQualityScoreMatrix &m,
                                     int maxScoreDrop) {
  int score = 0;
  int maxScore = 0;
  while (seq1 < seq1end) {
    score += m(*seq1++, *seq2++, *qual1++, *qual2++);
    if (score > maxScore) maxScore = score;
    else if (score <= 0 ||                       // non-optimal prefix
             seq1 == seq1end ||                  // non-optimal suffix
             score < maxScore - maxScoreDrop) {  // excessive score drop
      return false;
    }
  }
  return true;
}

int gaplessTwoQualityXdropOverlap(const uchar *seq1,
				  const uchar *qual1,
				  const uchar *seq2,
				  const uchar *qual2,
				  const TwoQualityScoreMatrix &m,
				  int maxScoreDrop,
				  size_t &reverseLength,
				  size_t &forwardLength) {
  int minScore = 0;
  int maxScore = 0;
  int score = 0;

  const uchar *rs1 = seq1;
  const uchar *rq1 = qual1;
  const uchar *rs2 = seq2;
  const uchar *rq2 = qual2;
  while (true) {
    --rs1;  --rq1;  --rs2;  --rq2;
    int s = m(*rs1, *rs2, *rq1, *rq2);
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
    else if (score < minScore) minScore = score;
  }

  maxScore = score - minScore;

  const uchar *fs1 = seq1;
  const uchar *fq1 = qual1;
  const uchar *fs2 = seq2;
  const uchar *fq2 = qual2;
  while (true) {
    int s = m(*fs1, *fs2, *fq1, *fq2);
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
    ++fs1;  ++fq1;  ++fs2;  ++fq2;
  }

  reverseLength = seq1 - (rs1 + 1);
  forwardLength = fs1 - seq1;
  return score;
}

int gaplessTwoQualityAlignmentScore(const uchar *seq1,
                                    const uchar *seq1end,
                                    const uchar *qual1,
                                    const uchar *seq2,
                                    const uchar *qual2,
                                    const TwoQualityScoreMatrix &m) {
  int score = 0;
  while (seq1 < seq1end) score += m(*seq1++, *seq2++, *qual1++, *qual2++);
  return score;
}

}
