// Copyright 2015 Martin C. Frith

// This class calculates E-values for pair-wise local alignment.
// It has 2 states: "good" and "bad".  It starts in the "bad" state.

// "database" = sequence1, "query" = sequence2

// For DNA-versus-protein alignment:
// protein = "database" = sequence1, DNA = "query" = sequence2

// "deletion" means deletion in sequence 2 relative to sequence 1
// "insertion" means insertion in sequence 2 relative to sequence 1

// A length-k deletion costs delOpen + k * delEpen
// A length-k insertion costs insOpen + k * insEpen

#ifndef LAST_EVALUER_HH
#define LAST_EVALUER_HH

#include "ScoreMatrixRow.hh"
#include "mcf_frameshift_xdrop_aligner.hh"

#include "alp/sls_alignment_evaluer.hpp"

namespace cbrc {

using namespace mcf;

class GeneticCode;

class LastEvaluer {
public:
  // This routine tries to initialize the object for a given set of
  // alignment parameters.  It may fail, i.e. set the object to the
  // "bad" state and throw an Sls::error.
  // These arguments are only used to lookup pre-calculated cases:
  // matrixName, matchScore, mismatchCost, geneticCodeName.
  // DNA-versus-protein alignment is indicated by: frameshiftCost >= 0.
  // As a special case, frameshiftCost==0 means no frameshifts.
  // For DNA-versus-protein alignment, letterFreqs2 should either be
  // NULL or point to 64 codon frequencies (aaa, aac, etc).
  void init(const char *matrixName,
	    int matchScore,
	    int mismatchCost,
	    const char *alphabet,
	    const ScoreMatrixRow *scoreMatrix,  // score[sequence1][sequence2]
	    const double *letterFreqs1,
	    const double *letterFreqs2,
	    bool isGapped,
	    int delOpen,
	    int delEpen,
	    int insOpen,
	    int insEpen,
	    int frameshiftCost,
	    const GeneticCode &geneticCode,
	    const char *geneticCodeName,
	    int verbosity);

  // The next 2 routines set up for sum-of-paths local alignment scores.
  // initFrameshift implements section 2.6 of "Improved
  // DNA-versus-protein homology search for protein fossils", Y Yao &
  // MC Frith.  (It doesn't check whether Equation 3 is satisfied.)
  // initFullScores does the equivalent for "model A" in "How sequence
  // alignment scores correspond to probability models", MC Frith
  // 2020, Bioinformatics 36(2):408-415.

  // The Freqs need not sum to 1.
  // substitutionProbs is S' in [Yao & Frith 2021].
  // scale  =  lambda  =  1/t in [Yao & Frith 2021].

  void initFullScores(const const_dbl_ptr *substitutionProbs,
		      const double *letterFreqs1, int alphabetSize1,
		      const double *letterFreqs2, int alphabetSize2,
		      const GapCosts &gapCosts, double scale, int verbosity,
		      bool isFrameshift = false);

  void initFrameshift(const const_dbl_ptr *substitutionProbs,
		      const double *proteinLetterFreqs, int numProteinLetters,
		      const double *tranDnaLetterFreqs, int numTranDnaLetters,
		      const GapCosts &gapCosts, double scale, int verbosity) {
    initFullScores(substitutionProbs,
		   proteinLetterFreqs, numProteinLetters,
		   tranDnaLetterFreqs, numTranDnaLetters,
		   gapCosts, scale, verbosity, true);
  }

  void setSearchSpace(double databaseLength,  // number of database letters
		      double databaseMaxSeqLength,  // length of longest seq
		      double numOfStrands) {  // 1 or 2
    if (databaseMaxSeqLength > 0) {
      databaseSeqLen = databaseMaxSeqLength;
      databaseSeqNum = databaseLength / databaseMaxSeqLength * numOfStrands;
    } else {
      this->databaseSeqLen = 1;  // ALP doesn't like 0
      this->databaseSeqNum = 0;
    }
  }

  bool isGood() const { return evaluer.isGood(); }

  // Don't call this in the "bad" state:
  double evaluePerArea(double score) const
  { return evaluer.evaluePerArea(score); }

  // Don't call this in the "bad" state or before calling setSearchSpace:
  double area(double score, double queryLength) const
  { return databaseSeqNum * evaluer.area(score, queryLength, databaseSeqLen); }

  // Don't call this in the "bad" state:
  double bitScore(double score) const { return evaluer.bitScore(score); }

  // Returns max(0, score with E-value == "evalue").
  // Don't call this in the "bad" state.
  double minScore(double evalue, double area) const;

  // Returns max(0, score with E-value == 1 per this many query letters).
  // Don't call this in the "bad" state or before calling setSearchSpace.
  double minScore(double queryLettersPerRandomAlignment) const;

  // Writes some parameters preceded by "#".  Does nothing in the "bad" state.
  void writeCommented(std::ostream& out) const;

  // Writes all parameters in full precision.  Does nothing in the "bad" state.
  void writeParameters(std::ostream& out) const;

private:
  Sls::AlignmentEvaluer evaluer;
  double databaseSeqLen;
  double databaseSeqNum;
};

}

#endif
