// Author: Martin C. Frith 2021
// SPDX-License-Identifier: GPL-3.0-or-later

// This calculates a kind of maximum local similarity score between
// two sequences.

// These inputs are parameters of an alignment probability model,
// Model A in [Fri20]:
// substitutionProbs: S' in [Fri20] section 2.1
// delInitProb: a'_D in [Fri20] section 5.1
// delNextProb: b'_D in [Fri20] section 5.1
// insInitProb: a'_I in [Fri20] section 5.1
// insNextProb: b'_I in [Fri20] section 5.1

// Any path through the model corresponds to a local alignment.  The
// model defines a probability for each path.  A path's "probability
// ratio" is its probability divided by the probability of a "null"
// path [Fri20].

// This code calculates the maximum, over coordinates i in sequence 1
// and j in sequence 2, sum of probability ratios of all paths passing
// through (i,j).  This is just like [Yao21] section 2.5, but for
// model A in [Fri20].

// [Fri20]: "How sequence alignment scores correspond to probability
// models", MC Frith 2020, Bioinformatics 36(2):408-415

// [Yao21]: "Improved DNA-versus-protein homology search for protein
// fossils", Y Yao & MC Frith

#ifndef MCF_ALIGMENT_PATH_ADDER_HH
#define MCF_ALIGMENT_PATH_ADDER_HH

#include <vector>

namespace mcf {

typedef const double *const_dbl_ptr;

class AlignmentPathAdder {
public:
  typedef unsigned char uchar;

  double maxSum(const uchar *seq1, int len1,
		const uchar *seq2, int len2,
		const const_dbl_ptr *substitutionProbs,
		double delInitProb, double delNextProb,
		double insInitProb, double insNextProb);

private:
  std::vector<double> dynamicProgrammingValues;
};

}

#endif
