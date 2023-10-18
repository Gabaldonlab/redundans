// Copyright 2008, 2009, 2011 Michiaki Hamada
// Copyright 2012, 2013 Toshiyuki Sato

#ifndef CENTROID_HH
#define CENTROID_HH

#include "GappedXdropAligner.hh"
#include "mcf_gap_costs.hh"
#include "SegmentPair.hh"
#include "OneQualityScoreMatrix.hh"

#include <assert.h>
#include <math.h>

#include <algorithm>
#include <vector>
#include <iostream> // for debug

namespace cbrc{

  /**
   * (1) Forward and backward algorithm on the DP region given by Xdrop algorithm
   * (2) \gamma-centroid decoding
   */
  class Centroid{
    typedef const double *const_dbl_ptr;
    typedef double *dbl_ptr;

    enum { rescaleStep = 16 };

  public:
    GappedXdropAligner& aligner() { return xa; }

    void setPssm ( const ScoreMatrixRow* pssm, size_t qsize, double T,
                   const OneQualityExpMatrix& oqem,
                   const uchar* sequenceBeg, const uchar* qualityBeg );

    // For a sequence with quality data, store the probability that
    // each position is each letter (possibly scaled by a constant per
    // position).  xxx I don't think this really belongs in Centroid.
    void setLetterProbsPerPosition(unsigned alphabetSize,
				   size_t sequenceLength,
				   const uchar *sequence,
				   const uchar *qualityCodes,
				   bool isFastq,
				   const double *qualToProbCorrect,
				   const double *letterProbs,
				   const uchar *toUnmasked);

    // start2 is the index of the first position to look at in the PSSM

    double forward(const uchar *seq1, const uchar *seq2,
		   size_t start2, bool isExtendFwd,
		   const const_dbl_ptr *substitutionProbs,
		   const GapCosts &gapCosts, int globality);

    void backward(bool isExtendFwd,
		  const const_dbl_ptr *substitutionProbs,
		  const GapCosts &gapCosts, int globality);

    double dp(int outputType, double gamma) {
      bestScore = 0;
      bestAntiDiagonal = 0;
      bestPos1 = 0;
      X.resize(fM.size());
      if (outputType == 5) return dp_centroid(gamma);
      if (outputType == 6) return dp_ama(gamma);
      return 0;
    }

    void traceback(std::vector<SegmentPair> &chunks,
		   int outputType, double gamma) const {
      if (outputType==5) traceback_centroid(chunks, gamma);
      if (outputType==6) traceback_ama(chunks, gamma);
    }

    double dp_centroid( double gamma );
    void traceback_centroid( std::vector< SegmentPair >& chunks, double gamma ) const;

    double dp_ama( double gamma );
    void traceback_ama( std::vector< SegmentPair >& chunks, double gamma ) const;

    void getMatchAmbiguities(std::vector<char>& ambiguityCodes,
			     size_t seq1end, size_t seq2end,
			     size_t length) const;

    void getDeleteAmbiguities(std::vector<char>& ambiguityCodes,
			      size_t seq1end, size_t seq1beg) const;

    void getInsertAmbiguities(std::vector<char>& ambiguityCodes,
			      size_t seq2end, size_t seq2beg) const;

    // Added by MH (2008/10/10) : compute expected counts for transitions and emissions
    void addExpectedCounts(size_t start2, bool isExtendFwd,
			   const const_dbl_ptr *substitutionProbs,
			   const GapCosts &gapCosts, unsigned alphabetSize,
			   const dbl_ptr *substitutionCounts,
			   double *transitionCounts);

  private:
    typedef double ExpMatrixRow[scoreMatrixRowSize];

    GappedXdropAligner xa;
    size_t numAntidiagonals;

    std::vector<double> pssmExp;
    ExpMatrixRow* pssmExp2; // pre-computed pssm for prob align

    std::vector<double> letterProbsPerPosition;  // for uncertain sequences

    typedef std::vector< double > dvec_t;

    dvec_t fM; // f^M(i,j)
    dvec_t fD; // f^D(i,j) Ix
    dvec_t fI; // f^I(i,j) Iy

    dvec_t bM; // b^M(i,j)
    dvec_t bD; // b^D(i,j)
    dvec_t bI; // b^I(i,j)

    dvec_t mD;
    dvec_t mI;
    dvec_t mX1;
    dvec_t mX2;

    dvec_t X; // DP tables for $gamma$-decoding

    dvec_t rescales;

    double rescaledSumOfProbRatios;

    const uchar *seq1ptr;
    const uchar *seq2ptr;
    const ExpMatrixRow *pssmPtr;

    double bestScore;
    size_t bestAntiDiagonal;
    size_t bestPos1;

    void initForward() {
      numAntidiagonals = xa.numAntidiagonals();
      assert(numAntidiagonals > 0);

      size_t totalNumOfCells = xa.scoreEndIndex(numAntidiagonals);
      if (fM.size() < totalNumOfCells) {
	fM.resize(totalNumOfCells);
	fD.resize(totalNumOfCells);
	fI.resize(totalNumOfCells);
	fM[xdropPadLen - 1] = 1;
      }

      size_t numOfRescales = (numAntidiagonals - 1) / rescaleStep;
      if (rescales.size() < numOfRescales) rescales.resize(numOfRescales);
    }

    void initBackward(size_t totalNumOfCells) {
      bM.assign(totalNumOfCells, 0.0);
      bD.assign(totalNumOfCells, 0.0);
      bI.assign(totalNumOfCells, 0.0);

      mD.assign(numAntidiagonals + 2, 0.0);
      mI.assign(numAntidiagonals + 2, 0.0);
    }

    void rescaleFwdProbs(size_t beg, size_t end, double scale) {
      for (size_t i = beg; i < end; ++i) {
	fM[i] *= scale;
	fD[i] *= scale;
	fI[i] *= scale;
      }
    }

    void rescaleBckProbs(size_t beg, size_t end, double scale) {
      for (size_t i = beg; i < end; ++i) {
	bM[i] *= scale;
	bD[i] *= scale;
	bI[i] *= scale;
      }
    }

    void updateScore(double score, size_t antiDiagonal, size_t cur) {
      if (bestScore < score) {
	bestScore = score;
	bestAntiDiagonal = antiDiagonal;
	bestPos1 = cur;
      }
    }
  };

  // Return an ASCII code representing an error probability.  The
  // printable codes are 33--126, but here we use a maximum of 125, so
  // that 126 is reserved for special cases.
  inline char asciiProbability(double probCorrect) {
    assert(probCorrect >= 0);
    //assert(probCorrect <= 1);  // can fail: floating point is imperfect
    double e = 1 - probCorrect;
    double f = std::max(e, 1e-10);  // avoid overflow errors
    double g = -10 * log10(f);
    int i = static_cast<int>(g);  // round fractions down
    int j = i + 33;
    int k = std::min(j, 125);
    return static_cast<char>(k);
  }
}

#endif
