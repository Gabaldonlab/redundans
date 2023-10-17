#ifndef INCLUDED_SLS_FALP_LIB
#define INCLUDED_SLS_FALP_LIB

/* $Id: $
* ===========================================================================
*
*							PUBLIC DOMAIN NOTICE
*			   National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/*****************************************************************************

File name: sls_fsa1_lib.hpp

Authors: Martin Frith, Sergey Sheetlin

Contents: library functions of main routines

******************************************************************************/

#include <string>
#include "sls_fsa1_pvalues.hpp"

namespace Sls {

	class FrameshiftAlignmentEvaluer {


	//Write the parameters:
	friend std::ostream &operator<<(std::ostream &s_,
			const FrameshiftAlignmentEvaluer &g_);

	//Read the parameters:
	friend std::istream &operator>>(std::istream &s_,
			FrameshiftAlignmentEvaluer &g_);

	public:

	//Computes gapless Gumbel parameters:
	void initGapless(long ntAlphabetSize_,//a number of letters in the DNA alphabet; usually 4
			long aaAlphabetSize_,//a number of letters in the amino acid alphabet
			const long *codonTable_,//defines the codon table; is an array of ntAlphabetSize_^3 integers

			const long *const *substitutionScoreMatrix_,//scoring matrix; aaAlphabetSize_ X aaAlphabetSize_
			const double *ntFreqs_,//background frequencies of letters in DNA sequences
			const double *aaFreqs_,//background frequencies of letters in amino acid sequences
			double max_time_=60);//maximum allowed calculation time in seconds


	void initFrameshift(long ntAlphabetSize_,//a number of letters in the DNA alphabet; usually 4
			long aaAlphabetSize_,//a number of letters in the amino acid alphabet
			const long *codonTable_,//defines the codon table; is an array of ntAlphabetSize_^3 integers
			//letters are numerated from 0; 
			//each codon defines a code (an integer number from [0,...,ntAlphabetSize_^3-1]);
			//for example, for the DNA alphabet ACGT, 
			//AAA corresponds to codonTable_[0]; AAC corresponds to codonTable_[1] and so on.
			//codonTable_[x] contains an amino acid number such that the amino acid corresponds to a codon with the code x
			//for example, for the amino acid alphabet ARNDCQEGHILKMFPSTWYVBJZX* ("*" defines the stop codons), 
			//if AAC is coded by N (amino acid with the number 2), then codonTable_[1]=2.

			const long *const *substitutionScoreMatrix_,//scoring matrix; aaAlphabetSize_ X aaAlphabetSize_
			const double *ntFreqs_,//background frequencies of letters in DNA sequences
			const double *aaFreqs_,//background frequencies of letters in amino acid sequences
			long gapOpen1_,//gap opening penalty for sequence #1
			long gapEpen1_,//gap extension penalty for sequence #1
			long gapOpen2_,//gap opening penalty for sequence #2
			long gapEpen2_,//gap extension penalty for sequence #2
			long frameshiftCost_,//frameshift penalty
			bool insertions_after_deletions_,//if true, then insertions after deletions are permitted
			double eps_lambda_,//relative error for the parameter lambda
			double eps_K_,//relative error for the parameter K
			double max_time_,//maximum allowed calculation time in seconds
			double max_mem_,//maximum allowed memory usage in MB
			long randomSeed_,//randomizaton seed

			//in the case when max_time_<=0, the following parameters are used:
			long int seq_number_=10000,//a number of alignments
			long int nalp_=5,//a number of ascending ladder points used in the calculation
			long int number_of_subsets_for_errors_calculation_=20);//a number of subsets used in the error calculation



	//Initializes Gumbel parameters using precalculated values:
	void initParameters(
	const AlignmentEvaluerParameters &parameters_);//parameters_ must be defined by the user

	void initParameters(
	const AlignmentEvaluerParametersWithErrors &parameters_);//parameters_ must be defined by the user

	//Computes P-values/E-values
	//Gumbel parameters must be defined via d_params
	void calc(double score_,//frameshift alignment score
			double seqlen1_,//length of sequence #1
			double seqlen2_,//length of sequence #2
			double &pvalue_,//resulted P-value
			double &pvalueErr_,//P-value error
			double &evalue_,//resulted E-value
			double &evalueErr_) const;//E-value error

	//Computes P-values/E-values without errors
	//Gumbel parameters must be defined via d_params
	void calc(double score_,//frameshift alignment score
			double seqlen1_,//length of sequence #1
			double seqlen2_,//length of sequence #2
			double &pvalue_,//resulted P-value
			double &evalue_) const;//resulted E-value

	double evalue(double score_,//frameshift alignment score
			double seqlen1_,//length of sequence #1
			double seqlen2_) const//length of sequence #2
	{
	  return area(score_, seqlen1_, seqlen2_) * evaluePerArea(score_);
	}

	//Computes P-values from E-values
	static double pvalue(double evalue_)
	{
		return sls_basic::one_minus_exp_function(-evalue_);
	}

	//The "area" is approximately seqlen1_*seqlen2_, but it is
	//modified by a finite size correction
	double area(double score_,//frameshift alignment score
			double seqlen1_,//length of sequence #1
			double seqlen2_) const;//length of sequence #2

	double evaluePerArea(double score_) const
	{
	  return d_params.K*exp(-d_params.lambda*score_);
	}

	double bitScore(double score_) const
	{
	  return (d_params.lambda*score_-log(d_params.K))/log(2.0);
	}


	//returns "true" if the set of parameters "d_params" is fully defined for P-value calculation
	bool isGood() const
	{
		return d_params.d_params_flag;
	}

	//provides access to the set of Gumbel parameters
	const FALP_set_of_parameters &parameters() const { return d_params; }

	private:

	//check correctness of the input parameters for gapless alignment
	void assert_Gapless_input_parameters(
			long ntAlphabetSize_,//a number of letters in the DNA alphabet; usually 4
			long aaAlphabetSize_,//a number of letters in the amino acid alphabet
			const long *codonTable_,//defines the codon table; is an array of ntAlphabetSize_^3 integers

			const double *ntFreqs_,//background frequencies of letters in DNA sequences
			const double *aaFreqs_,//background frequencies of letters in amino acid sequences
			double *&ntFreqs_normalized_,//normalized background frequencies of letters in DNA sequences
			double *&aaFreqs_normalized_,//normalized background frequencies of letters in amino acid sequences
			const std::string function_name_);//"assert_Gapless_input_parameters" is called from "function_name_" function

	private:

		FALP_set_of_parameters d_params;//set of Gumbel parameters

	};
}

#endif //! INCLUDED

