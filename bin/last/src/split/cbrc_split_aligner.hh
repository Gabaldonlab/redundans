// Copyright 2012 Risa Kawaguchi
// Copyright 2013, 2014 Martin C. Frith

#ifndef CBRC_SPLIT_ALIGNER_HH
#define CBRC_SPLIT_ALIGNER_HH

#include "cbrc_unsplit_alignment.hh"
#include "cbrc_int_exponentiator.hh"

#include "Alphabet.hh"
#include "MultiSequence.hh"

#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include <map>
#include <string>
#include <vector>

namespace cbrc {

struct SplitAlignerParams {
  static int scoreFromProb(double prob, double scale) {
    return floor(scale * log(prob) + 0.5);
  }

  // A   deletion of length k scores: delOpenScore + k * delGrowScore
  // An insertion of length k scores: insOpenScore + k * insGrowScore

  // "jumpScore" is the (negative) score for a trans-splice.

  // "restartScore" is the (negative) score for re-starting an
  // alignment, in the repeated matches algorithm from chapter 2 of
  // Durbin, Eddy, et al.  In that book, it is called "-T".

  // "qualityOffset" is 33 for fastq-sanger or 64 for fastq-illumina

  void setParams(int delOpenScoreIn, int delGrowScoreIn,
		 int insOpenScoreIn, int insGrowScoreIn,
		 int jumpScoreIn, int restartScoreIn, double scaleIn,
		 int qualityOffsetIn);

  void setSpliceParams(double splicePriorIn,
		       double meanLogDistIn, double sdevLogDistIn);

  void setScoreMat(const std::vector< std::vector<int> > &sm,
		   const char *rowNames, const char *colNames);

  void readGenome(const std::string &baseName);

  // XXX this should allow us to specify scores for gt-ag, at-ac, etc.
  void setSpliceSignals();

  // Outputs some algorithm parameters on lines starting with "#"
  void print() const;

  static const int numQualCodes = 64;
  static int score_mat[64][64][numQualCodes];
  int maxMatchScore;
  int qualityOffset;
  int delOpenScore;
  int delGrowScore;
  int insOpenScore;
  int insGrowScore;
  int jumpScore;
  int restartScore;
  double jumpProb;
  double restartProb;
  double scale;
  IntExponentiator scaledExp;  // for fast calculation of exp(x / scale)

  double splicePrior;
  double meanLogDist;
  double sdevLogDist;
  double spliceTerm1;
  double spliceTerm2;
  unsigned maxSpliceDist;
  int maxSpliceScore;
  int maxSpliceBegEndScore;
  std::vector<int> spliceScoreTable;  // lookup table
  std::vector<double> spliceProbTable;  // lookup table
  unsigned spliceTableSize;
  MultiSequence genome[32];
  Alphabet alphabet;
  typedef std::map<std::string, unsigned long long> StringNumMap;
  StringNumMap chromosomeIndex;
  int spliceBegScores[(4 * 4 + 1) * 2];  // donor score for any dinucleotide
  int spliceEndScores[(4 * 4 + 1) * 2];  // acceptor score for any dinucleotide
  double spliceBegProbs[(4 * 4 + 1) * 2];
  double spliceEndProbs[(4 * 4 + 1) * 2];

  bool isSpliced() const { return restartProb <= 0; }

  bool isGenome() const { return !chromosomeIndex.empty(); }

  bool isSpliceCoords() const { return splicePrior > 0 || isGenome(); }

  void dpExtensionMinScores(size_t &minScore1, size_t &minScore2) const;

  void seqEnds(const uchar *&beg, const uchar *&end,
	       const char *seqName) const;

  int spliceScore(unsigned d) const
  { return d < spliceTableSize ? spliceScoreTable[d] : calcSpliceScore(d); }

  double spliceProb(unsigned d) const
  { return d < spliceTableSize ? spliceProbTable[d] : calcSpliceProb(d); }

  void spliceBegSignal(char *out, const char *seqName, bool isForwardStrand,
		       bool isSenseStrand, unsigned coord) const;

  void spliceEndSignal(char *out, const char *seqName, bool isForwardStrand,
		       bool isSenseStrand, unsigned coord) const;

  int calcSpliceScore(double dist) const;

  double calcSpliceProb(double dist) const
  { return scaledExp(calcSpliceScore(dist)); }

  size_t maxGenomeVolumes() const { return sizeof genome / sizeof *genome; }

  void readGenomeVolume(const std::string &baseName,
			size_t seqCount, size_t volumeNumber);
};

class SplitAligner {
public:
    SplitAligner() { maxCellsPerMatrix = 0; scMemory = dpMemory = 0; }
    ~SplitAligner() { free(scMemory); free(dpMemory); }

    // Prepares to analyze some candidate alignments for one query
    // sequence: sets the number of DP matrix cells (and thus memory)
    void layout(const SplitAlignerParams &params,
		const UnsplitAlignment *beg, const UnsplitAlignment *end);

    // The number of cells in each dynamic programming matrix
    size_t cellsPerDpMatrix() const
    { return matrixRowOrigins[numAlns-1] + dpEnd(numAlns-1) + 1; }

    // Bytes of memory needed for the current query sequence (roughly)
    size_t memory(const SplitAlignerParams &params,
		  bool isBothSpliceStrands) const;

    // Call this before viterbi/forward/backward, and after layout
    void initMatricesForOneQuery(const SplitAlignerParams &params,
				 bool isBothSpliceStrands);

    // returns the optimal split-alignment score
    long viterbi(const SplitAlignerParams &params) {
      resizeVector(Vvec);
      for (unsigned i = 0; i < numAlns; ++i) {
	Vmat[matrixRowOrigins[i] + dpBegs[i]] = INT_MIN/2;
      }
      return params.isSpliced() ? viterbiSplice(params) : viterbiSplit(params);
    }

    // Gets the chunks of an optimal split alignment.
    // For each chunk, it gets:
    // 1. The index of the candidate alignment that the chunk comes from
    // 2. The chunk's start coordinate in the query sequence
    // 3. The chunk's end coordinate in the query sequence
    // It gets the chunks in reverse order, from query end to query start.
    void traceBack(const SplitAlignerParams &params, long viterbiScore,
		   std::vector<unsigned>& alnNums,
		   std::vector<unsigned>& queryBegs,
		   std::vector<unsigned>& queryEnds) const;

    // Calculates the alignment score for a segment of an alignment
    int segmentScore(unsigned alnNum,
		     unsigned queryBeg, unsigned queryEnd) const;

    void exponentiateScores(const SplitAlignerParams &params) {
      size_t s = cellsPerDpMatrix() * 2;
      for (size_t i = 0; i < s; ++i) Sexp[i] = params.scaledExp(Smat[i]);
      // if x/scale < about -745, then exp(x/scale) will be exactly 0.0
    }

    void forwardBackward(const SplitAlignerParams &params) {
      resizeVector(rescales);
      for (unsigned i = 0; i < numAlns; ++i) {
	Fmat[matrixRowOrigins[i] + dpBegs[i]] = 0;
	Bmat[matrixRowOrigins[i] + dpEnds[i]] = 0;
      }
      if (params.isSpliced()) {
	forwardSplice(params);
	backwardSplice(params);
      } else {
	forwardSplit(params);
	backwardSplit(params);
      }
    }

    // Returns one probability per column, for a segment of an alignment
    std::vector<double> marginalProbs(unsigned queryBeg, unsigned alnNum,
				      unsigned alnBeg, unsigned alnEnd) const;

    // Toggles between forward and reverse-complement splice signals
    void flipSpliceSignals(const SplitAlignerParams &params);

    // This returns log(p / (1-p)), where p is the probability that
    // the query uses splice signals in the orientation currently set
    // by flipSpliceSignals()
    double spliceSignalStrandLogOdds() const;

    // Gets the 2 genome bases immediately downstream of queryPos in
    // alnNum, and writes them into the buffer pointed to by "out"
    void spliceBegSignal(char *out, const SplitAlignerParams &params,
			 unsigned alnNum, unsigned queryPos,
			 bool isSenseStrand) const {
      const UnsplitAlignment &a = alns[alnNum];
      params.spliceBegSignal(out, a.rname, a.isForwardStrand(), isSenseStrand,
			     cell(spliceBegCoords, alnNum, queryPos));
    }

    // Gets the 2 genome bases immediately upstream of queryPos in
    // alnNum, and writes them into the buffer pointed to by "out"
    void spliceEndSignal(char *out, const SplitAlignerParams &params,
			 unsigned alnNum, unsigned queryPos,
			 bool isSenseStrand) const {
      const UnsplitAlignment &a = alns[alnNum];
      params.spliceEndSignal(out, a.rname, a.isForwardStrand(), isSenseStrand,
			     cell(spliceEndCoords, alnNum, queryPos));
    }

private:
    unsigned numAlns;  // the number of candidate alignments (for 1 query)
    const UnsplitAlignment *alns;  // the candidates
    unsigned minBeg;  // the minimum query start coordinate of any candidate
    unsigned maxEnd;  // the maximum query end coordinate of any candidate
    std::vector<unsigned> dpBegs;  // dynamic programming begin coords
    std::vector<unsigned> dpEnds;  // dynamic programming end coords
    std::vector<size_t> matrixRowOrigins;  // layout of ragged matrices

    size_t maxCellsPerMatrix;
    void *scMemory;
    void *dpMemory;

    int *Smat;
    // Smat holds position-specific substitution, insertion, and
    // deletion scores for the candidate alignments of one query
    // sequence to a genome.  These scores, for each candidate
    // alignment i, are called Aij and Dij in [Frith&Kawaguchi 2015].
    // Aij holds scores at query bases, in Smat[i][1,3,5,...,2n-1].
    // Dij holds scores between query bases, in Smat[i][0,2,4,...,2n].

    long *Vmat;  // DP matrix for Viterbi algorithm
    std::vector<long> Vvec;  // DP vector for Viterbi algorithm

    float *Sexp;
    // Sexp holds exp(Smat / t): these values are called A'ij and D'ij
    // in [Frith&Kawaguchi 2015].

    double *Fmat;  // DP matrix for Forward algorithm
    double *Bmat;  // DP matrix for Backward algorithm
    std::vector<double> rescales;  // the usual scaling for numerical stability

    long *VmatRev;
    std::vector<long> VvecRev;
    double *FmatRev;
    double *BmatRev;
    std::vector<double> rescalesRev;

    std::vector<unsigned> sortedAlnIndices;
    std::vector<unsigned> oldInplayAlnIndices;
    std::vector<unsigned> newInplayAlnIndices;

    std::vector<unsigned> spliceBegCoords;
    std::vector<unsigned> spliceEndCoords;
    std::vector<unsigned char> spliceBegSignals;
    std::vector<unsigned char> spliceEndSignals;
    std::vector<unsigned> rBegs;  // genomic beg coordinate of each candidate
    std::vector<unsigned> rEnds;  // genomic end coordinate of each candidate
    std::vector<unsigned> rnameAndStrandIds;
    const int *spliceBegScores;
    const int *spliceEndScores;
    const double *spliceBegProbs;
    const double *spliceEndProbs;
    int spliceBegScore(bool isGenome, size_t ij) const {
      return isGenome ? spliceBegScores[spliceBegSignals[ij]] : 0;
    }
    int spliceEndScore(bool isGenome, size_t ij) const {
      return isGenome ? spliceEndScores[spliceEndSignals[ij]] : 0;
    }
    double spliceBegProb(bool isGenome, size_t ij) const {
      return isGenome ? spliceBegProbs[spliceBegSignals[ij]] : 1.0;
    }
    double spliceEndProb(bool isGenome, size_t ij) const {
      return isGenome ? spliceEndProbs[spliceEndSignals[ij]] : 1.0;
    }
    void initSpliceCoords(unsigned i);
    void initSpliceSignals(const SplitAlignerParams &params, unsigned i);
    void initRnameAndStrandIds();
    void initRbegsAndEnds();

    void updateInplayAlnIndicesF(unsigned& sortedAlnPos,
				 unsigned& oldNumInplay,
				 unsigned& newNumInplay, unsigned j);

    void updateInplayAlnIndicesB(unsigned& sortedAlnPos,
				 unsigned& oldNumInplay,
				 unsigned& newNumInplay, unsigned j);

    long viterbiSplit(const SplitAlignerParams &params);
    long viterbiSplice(const SplitAlignerParams &params);

    void forwardSplit(const SplitAlignerParams &params);
    void backwardSplit(const SplitAlignerParams &params);
    void forwardSplice(const SplitAlignerParams &params);
    void backwardSplice(const SplitAlignerParams &params);

    unsigned findScore(bool isGenome, unsigned j, long score) const;
    unsigned findSpliceScore(const SplitAlignerParams &params,
			     unsigned i, unsigned j, long score) const;
    long scoreFromSplice(const SplitAlignerParams &params,
			 unsigned i, unsigned j, unsigned oldNumInplay,
			 unsigned& oldInplayPos) const;
    long endScore() const;
    unsigned findEndScore(long score) const;

    // "dp" means "dynamic programming":
    unsigned dpBeg(unsigned i) const { return dpBegs[i]; }
    unsigned dpEnd(unsigned i) const { return dpEnds[i]; }

    template<typename T> T&
    cell(std::vector<T>& v, unsigned j) const
    { return v[j - minBeg]; }

    template<typename T> const T&
    cell(const std::vector<T>& v, unsigned j) const
    { return v[j - minBeg]; }

    // cell j in row i of a ragged matrix
    template<typename T> T&
    cell(std::vector<T>& v, unsigned i, unsigned j) const
    { return v[matrixRowOrigins[i] + j]; }

    // cell j in row i of a ragged matrix
    template<typename T> const T&
    cell(const std::vector<T>& v, unsigned i, unsigned j) const
    { return v[matrixRowOrigins[i] + j]; }

    long cell(const long *v, unsigned i, unsigned j) const
    { return v[matrixRowOrigins[i] + j]; }

    template<typename T>
    void resizeVector(T& v) const
    { v.resize(maxEnd - minBeg + 1); }

    template<typename T>
    void resizeMatrix(T& m) const {
      // This reserves size for a ragged matrix, which is actually
      // stored in a flat vector.  There are numAlns rows, and row i
      // has dpEnd(i) - dpBeg(i) + 1 cells.
      size_t s = cellsPerDpMatrix();
      if (m.size() < s) m.resize(s);
    }

    double probFromSpliceF(const SplitAlignerParams &params,
			   unsigned i, unsigned j, unsigned oldNumInplay,
			   unsigned& oldInplayPos) const;

    double probFromSpliceB(const SplitAlignerParams &params,
			   unsigned i, unsigned j, unsigned oldNumInplay,
			   unsigned& oldInplayPos) const;

    void calcBaseScores(const SplitAlignerParams &params, unsigned i);
    void initDpBounds(const SplitAlignerParams &params);
};

}

#endif
