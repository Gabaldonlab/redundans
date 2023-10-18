// Author: Martin C. Frith 2013
// SPDX-License-Identifier: GPL-3.0-or-later

// This routine reads alignments of query sequences to a genome, and
// estimates which alignment parts represent the source of each query.
// It allows different parts of one query to come from different parts
// of the genome.

// The input should be in MAF format
// (http://genome.ucsc.edu/FAQ/FAQformat.html#format5).  It must
// include a header of the sort written by lastal, with score
// parameters and genome size.

#ifndef LAST_SPLIT_HH
#define LAST_SPLIT_HH

#include "last_split_options.hh"

void lastSplit(LastSplitOptions &opts);

#endif
