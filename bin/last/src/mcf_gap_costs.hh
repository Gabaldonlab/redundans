// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

// Gap costs for pair-wise sequence alignment.
// A gap of length k costs: openCost + k * growCost.

// "Generalized affine gap cost" is supported [SF Altschul 1998
// Proteins 32(1):88-96].  In this scheme, a "gap" may consist of
// unaligned regions in both sequences.  If these unaligned regions
// have sizes j and k, where j <= k, the cost is:
// openCost + growCost*(k-j) + pairCost*j

// Different costs for deletions and insertions are supported.

// Piecewise linear gap costs are supported.

#ifndef MCF_GAP_COSTS_HH
#define MCF_GAP_COSTS_HH

#include <vector>

namespace mcf {

struct GapCosts {
  struct Piece {
    int openCost;
    int growCost;
  };

  struct ProbPiece {
    double openProb;
    double growProb;
  };

  std::vector<Piece> delPieces;
  std::vector<Piece> insPieces;
  int pairCost;

  int frameshiftCost;

  // these are for "new-style" DNA-protein alignment with frameshifts:
  int delScore1;
  int delScore2;
  int delScore3;
  int insScore1;
  int insScore2;
  int insScore3;

  bool isAffine;

  // these are for calculating alignment probabilities:
  std::vector<ProbPiece> delProbPieces;
  std::vector<ProbPiece> insProbPieces;
  double delProb1;
  double delProb2;
  double delProb3;
  double insProb1;
  double insProb2;
  double insProb3;

  // Assign piecewise linear open and grow costs, and one pairCost.
  // If unalignedPairCost <= 0, assign non-generalized costs.
  // Throw a runtime_error if any growCost is <= 0.
  // delOpenCosts.size() must equal delGrowCosts.size(), and
  // insOpenCosts.size() must equal insGrowCosts.size().
  // The del and ins vectors must not be empty.
  // The number of frameshift costs should be: 0 (no frameshifts), 1
  // (old-style frameshifts), or 4 (new-style frameshifts).
  void assign(const std::vector<int> &delOpenCosts,
	      const std::vector<int> &delGrowCosts,
	      const std::vector<int> &insOpenCosts,
	      const std::vector<int> &insGrowCosts,
	      const std::vector<int> &frameshiftCosts,
	      int unalignedPairCost,
	      double scale);  // probability ratio  =  exp(scale * score)

  // The cost of a "gap" consisting of unaligned letters in the query
  // and/or reference sequence
  int cost(int refInsertLen, int qryInsertLen) const;

  bool isNewFrameshifts() const { return frameshiftCost == -2; }
};

}

#endif
