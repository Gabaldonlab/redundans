// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mcf_gap_costs.hh"
#include <stdexcept>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>  // size_t

static void err(const char *s) {
  throw std::runtime_error(s);
}

namespace mcf {

static void assignGapCostPieces(const std::vector<int> &openCosts,
				const std::vector<int> &growCosts,
				std::vector<GapCosts::Piece> &pieces) {
  assert(openCosts.size() == growCosts.size());
  pieces.clear();
  for (size_t i = 0; i < openCosts.size(); ++i) {
    GapCosts::Piece p = {openCosts[i], growCosts[i]};
    if (p.growCost <= 0) err("gap extension cost must be > 0");
    pieces.push_back(p);
  }
  assert(!pieces.empty());
}

static void assignProbs(std::vector<GapCosts::ProbPiece> &probPieces,
			const std::vector<GapCosts::Piece> &costPieces,
			double scale, bool isMergedPathCosts) {
  probPieces.resize(costPieces.size());
  for (size_t i = 0; i < costPieces.size(); ++i) {
    double open = exp(scale * -costPieces[i].openCost);
    double grow = exp(scale * -costPieces[i].growCost);

    if (isMergedPathCosts) {
      double init = open * grow;  // probability ratio for 1st letter in a gap

      // Path parameter from alignment parameter, as in Supplementary
      // section 3.1 of "How sequence alignment scores correspond to
      // probability models", Bioinformatics 2020:
      double next = grow - init;

      open = init;
      grow = next;
    }

    probPieces[i].openProb = open;
    probPieces[i].growProb = grow;
  }
}

void GapCosts::assign(const std::vector<int> &delOpenCosts,
		      const std::vector<int> &delGrowCosts,
		      const std::vector<int> &insOpenCosts,
		      const std::vector<int> &insGrowCosts,
		      const std::vector<int> &frameshiftCosts,
		      int unalignedPairCost,
		      double scale) {
  assignGapCostPieces(delOpenCosts, delGrowCosts, delPieces);
  assignGapCostPieces(insOpenCosts, insGrowCosts, insPieces);
  if (unalignedPairCost > 0) {
    pairCost = unalignedPairCost;
  } else {
    pairCost = delPieces[0].openCost + delPieces[0].growCost +
      insPieces[0].openCost + insPieces[0].growCost;
  }
  isAffine = (delPieces.size() < 2 && insPieces.size() < 2 &&
	      pairCost >= delPieces[0].growCost + insPieces[0].growCost +
	      std::max(delPieces[0].openCost, insPieces[0].openCost));

  bool isMergedPathCosts = (frameshiftCosts.size() < 2);
  assignProbs(delProbPieces, delPieces, scale, isMergedPathCosts);
  assignProbs(insProbPieces, insPieces, scale, isMergedPathCosts);

  if (frameshiftCosts.empty()) {
    frameshiftCost = -1;
  } else if (frameshiftCosts.size() == 1) {
    frameshiftCost = frameshiftCosts[0];
    assert(frameshiftCost >= 0);
  } else {
    delScore1 = -(1 * delPieces[0].growCost + frameshiftCosts[0]);
    delScore2 = -(2 * delPieces[0].growCost + frameshiftCosts[1]);
    delScore3 = -(3 * delPieces[0].growCost);
    insScore1 = -(1 * insPieces[0].growCost + frameshiftCosts[2]);
    insScore2 = -(2 * insPieces[0].growCost + frameshiftCosts[3]);
    insScore3 = -(3 * insPieces[0].growCost);
    frameshiftCost = -2;

    delProb1 = exp(scale * delScore1);
    delProb2 = exp(scale * delScore2);
    delProb3 = exp(scale * delScore3);
    insProb1 = exp(scale * insScore1);
    insProb2 = exp(scale * insScore2);
    insProb3 = exp(scale * insScore3);
  }
}

int GapCosts::cost(int refInsertLen, int qryInsertLen) const {
  int genCost = INT_MAX;

  int delCost = 0;
  if (refInsertLen > 0) {
    delCost = INT_MAX;
    for (size_t i = 0; i < delPieces.size(); ++i) {
      int s = delPieces[i].openCost + delPieces[i].growCost * refInsertLen;
      delCost = std::min(delCost, s);
      if (refInsertLen >= qryInsertLen) {
	int t = s + (pairCost - delPieces[i].growCost) * qryInsertLen;
	genCost = std::min(genCost, t);
      }
    }
  }

  int insCost = 0;
  if (qryInsertLen > 0) {
    insCost = INT_MAX;
    for (size_t i = 0; i < insPieces.size(); ++i) {
      int s = insPieces[i].openCost + insPieces[i].growCost * qryInsertLen;
      insCost = std::min(insCost, s);
      if (qryInsertLen >= refInsertLen) {
	int t = s + (pairCost - insPieces[i].growCost) * refInsertLen;
	genCost = std::min(genCost, t);
      }
    }
  }

  return std::min(delCost + insCost, genCost);
}

}
