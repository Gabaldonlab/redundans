// Copyright 2011 Martin C. Frith

// Score matrix for aligning letter1 to letter2, where both letters
// have quality codes.

// Quality codes are transformed to quality scores like this:
// q  =  qualityCode - qualityOffset

// Quality scores are transformed to error probabilities like this:
// For phred qualities:        e(q)  =  10^(-q/10)
// Else for solexa qualities:  e(q)  =  1 / (1 + 10^(q/10))

// Error probabilities are transformed to "certainties" like this:
// c = 1 - e / (1 - backgroundProbability[letter])

// The score matrix Sxy is adjusted by the quality data like this:
// Rxy     =  exp(lambda * Sxy)
// R'xyij  =  c1*c2*Rxy + (1 - c1*c2) * matrixBias
// S'xyij  =  nearestInt[ ln(R'xyij) / lambda ]

// Some letters may be considered to be "masked" versions of other
// letters (e.g. indicated by lowercase in the original sequence).  If
// masking is not applied, masked letters get the same scores as their
// unmasked versions.  If masking is applied:
// maskedScore = min(score, 0).

// If either letter is ambiguous (e.g. 'N' for DNA), then quality data
// is ignored: S'xyij = Sxy.

#ifndef TWO_QUALITY_SCORE_MATRIX_HH
#define TWO_QUALITY_SCORE_MATRIX_HH

#include "mcf_substitution_matrix_stats.hh"
#include "ScoreMatrixRow.hh"

#include <vector>

namespace cbrc {

typedef unsigned char uchar;

// The "indexer" reduces memory usage.  Which hopefully makes things
// more cache-friendly.
class TwoQualityMatrixIndexer {
 public:
  static const int qualityCapacity = 128;
  static const int numNormalLetters = 4;
  static const int numQualityLetters = numNormalLetters * 2;  // normal+masked
  static const int numSymbols = numQualityLetters * qualityCapacity;
  static const int minQuality = 32;

  void init(const uchar *toUnmasked);

  int operator()(int letter1, int letter2, int quality1, int quality2) const {
    int i = indexMap[subindex(letter1, quality1)];
    int j = indexMap[subindex(letter2, quality2)];
    return i * numSymbols + j;
  }

 private:
  std::vector<int> indexMap;

  static int subindex(int letter, int quality)
  { return quality * scoreMatrixRowSize + letter; }
};

class TwoQualityScoreMatrix {
  friend class TwoQualityExpMatrix;

 public:
  void init(const ScoreMatrixRow *scoreMatrix,
	    const mcf::SubstitutionMatrixStats &matStats,
            bool isPhred1,
            int qualityOffset1,
            bool isPhred2,
            int qualityOffset2,
            const uchar *toUnmasked,  // maps letters to unmasked letters
            bool isMask,
	    bool isMatchMismatchMatrix);  // if "true", init is faster

  // Tests whether init has been called:
  operator const void *() const { return data.empty() ? 0 : this; }

  int operator()(int letter1, int letter2, int quality1, int quality2) const
  { return data[indexer(letter1, letter2, quality1, quality2)]; }

 private:
  std::vector<int> data;
  TwoQualityMatrixIndexer indexer;
};

// This class stores: exp(score(x, y, i, j) / temperature)
class TwoQualityExpMatrix {
 public:
  void init(const TwoQualityScoreMatrix &m, double temperature);

  // Tests whether init has been called:
  operator const void *() const { return data.empty() ? 0 : this; }

  double operator()(int letter1, int letter2, int quality1, int quality2) const
  { return data[indexer(letter1, letter2, quality1, quality2)]; }

 private:
  std::vector<double> data;
  TwoQualityMatrixIndexer indexer;
};

}

#endif
