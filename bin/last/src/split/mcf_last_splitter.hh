// Author: Martin C. Frith 2022
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef MCF_LAST_SPLITTER_HH
#define MCF_LAST_SPLITTER_HH

#include "cbrc_split_aligner.hh"
#include "last_split_options.hh"

#include <iostream>

namespace mcf {

void setLastSplitParams(cbrc::SplitAlignerParams &params,
			const LastSplitOptions &opts,
			const std::vector< std::vector<int> > &scoreMatrix,
			const char *rowNames, const char *colNames,
			int delOpenCost, int delGrowCost,
			int insOpenCost, int insGrowCost,
			double scale, double genomeSize, int sequenceFormat);

struct SliceData {
  unsigned alnBeg;
  unsigned alnEnd;
  int score;
};

class LastSplitter {
public:
  void reserve(size_t s) { mafs.reserve(s); }

  void addMaf(char **linesBeg, char **linesEnd, bool isTopSeqQuery) {
    mafs.push_back(cbrc::UnsplitAlignment(linesBeg, linesEnd, isTopSeqQuery));
  }

  // Calculate and store output, and clear MAFs
  void split(const LastSplitOptions &opts,
	     const cbrc::SplitAlignerParams &params, bool isAlreadySplit);

  bool isOutputEmpty() const { return outputText.empty(); }

  void printOutput() const
  { std::cout.write(&outputText[0], outputText.size()); }

  void clearOutput() { outputText.clear(); }

private:
  cbrc::SplitAligner sa;
  std::vector<cbrc::UnsplitAlignment> mafs;
  std::vector<SliceData> slices;
  std::vector<char> outputText;

  void doOneQuery(const LastSplitOptions &opts,
		  const cbrc::SplitAlignerParams &params, bool isAlreadySplit,
		  const cbrc::UnsplitAlignment *beg,
		  const cbrc::UnsplitAlignment *end);

  void doOneAlignmentPart(const LastSplitOptions &opts,
			  const cbrc::SplitAlignerParams &params,
			  bool isAlreadySplit, const cbrc::UnsplitAlignment &a,
			  unsigned numOfParts, unsigned partNum,
			  const SliceData &sd, unsigned alnNum,
			  unsigned qSliceBeg, unsigned qSliceEnd,
			  bool isSenseStrand, double senseStrandLogOdds);
};

}

#endif
