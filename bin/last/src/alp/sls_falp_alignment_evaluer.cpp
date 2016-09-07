/* $Id: $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
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

File name: sls_fsa1_lib.cpp

Author: Sergey Sheetlin

Contents: library functions of main routines

******************************************************************************/


#include "sls_falp_alignment_evaluer.hpp"

#include "sls_fsa1_parameters.hpp"
#include "sls_fsa1_utils.hpp"
#include "sls_alp_regression.hpp"
#include "njn_localmaxstatmatrix.hpp"
#include "njn_localmaxstatutil.hpp"
#include "sls_fsa1.hpp"

using namespace std;


namespace Sls {

// Write the parameters:
std::ostream &operator<<(std::ostream &s_,
const FrameshiftAlignmentEvaluer &g_)
{

	if(!FALP_pvalues::assert_Gumbel_parameters(
	g_.d_params)||!g_.isGood())
	{
		throw error("Error - the Gumbel parameters are not defined properly in the function \"std::ostream &operator<<\"\n",1);
	};

	s_<<g_.d_params;
	return s_;
}

// Read the parameters:
std::istream &operator>>(std::istream &s_,
FrameshiftAlignmentEvaluer &g_)
{
	try
	{
		g_.d_params.d_params_flag=false;
		s_>>g_.d_params;
		g_.d_params.d_params_flag=true;

		//precompute intercepts
		FALP_pvalues::compute_intercepts(g_.d_params);

		if(!FALP_pvalues::assert_Gumbel_parameters(
		g_.d_params)||!g_.isGood())
		{
			g_.d_params.d_params_flag=false;
		};

		return s_;
	}
	catch (...)
	{ 
		g_.d_params.d_params_flag=false;
		throw;
	};
}

//check correctness of the input parameters for gapless alignment
void FrameshiftAlignmentEvaluer::assert_Gapless_input_parameters(
long ntAlphabetSize_,//a number of letters in the DNA alphabet; usually 4
long aaAlphabetSize_,//a number of letters in the amino acid alphabet
const long *codonTable_,//defines the codon table; is an array of ntAlphabetSize_^3 integers

const double *ntFreqs_,//background frequencies of letters in DNA sequences
const double *aaFreqs_,//background frequencies of letters in amino acid sequences
double *&ntFreqs_normalized_,//normalized background frequencies of letters in DNA sequences
double *&aaFreqs_normalized_,//normalized background frequencies of letters in amino acid sequences
const string function_name_)//"assert_Gapless_input_parameters" is called from "function_name_" function
{
	if(!(ntAlphabetSize_>0))
	{
		d_params.d_params_flag=false;
		throw error("Error - the parameter \"ntAlphabetSize_\" in the function \""+function_name_+"\" must be positive\n",1);
	};

	if(!(aaAlphabetSize_>0))
	{
		d_params.d_params_flag=false;
		throw error("Error - the parameter \"aaAlphabetSize_\" in the function \""+function_name_+"\" must be positive\n",1);
	};

	long int i;
	for(i=0;i<ntAlphabetSize_*ntAlphabetSize_*ntAlphabetSize_;i++)
	{
		if(!(codonTable_[i]>=0&&codonTable_[i]<aaAlphabetSize_))
		{
			d_params.d_params_flag=false;
			throw error("Error - the value \"codonTable_["+FSA_utils::long_to_string(i)+"]\" in the function \""+function_name_+"\" is incorrect\n",1);
		};
	};

	double sum_nt=0;
	for(i=0;i<ntAlphabetSize_;i++)
	{
		if(ntFreqs_[i]<0)
		{
			d_params.d_params_flag=false;
			throw error("Error - the value \"ntFreqs_["+FSA_utils::long_to_string(i)+"]\" in the function \""+function_name_+"\" must be non-negative\n",1);
		};
		sum_nt+=ntFreqs_[i];
	};

	if(sum_nt<=0)
	{
		throw error("Error - sum of the frequencies \"ntFreqs_\" is non-positive in the function \""+function_name_+"\"\n",1);
	};

	ntFreqs_normalized_=new double[ntAlphabetSize_];
	FSA_utils::assert_mem(ntFreqs_normalized_);

	for(i=0;i<ntAlphabetSize_;i++)
	{
		ntFreqs_normalized_[i]=ntFreqs_[i]/sum_nt;
	};

	double sum_aa=0;
	for(i=0;i<aaAlphabetSize_;i++)
	{
		if(aaFreqs_[i]<0)
		{
			d_params.d_params_flag=false;
			throw error("Error - the value \"aaFreqs_["+FSA_utils::long_to_string(i)+"]\" in the function \""+function_name_+"\" must be non-negative\n",1);
		};
		sum_aa+=aaFreqs_[i];
	};

	if(sum_aa<=0)
	{
		throw error("Error - sum of the frequencies \"aaFreqs_\" is non-positive in the function \""+function_name_+"\"\n",1);
	};

	aaFreqs_normalized_=new double[aaAlphabetSize_];
	FSA_utils::assert_mem(aaFreqs_normalized_);

	for(i=0;i<aaAlphabetSize_;i++)
	{
		aaFreqs_normalized_[i]=aaFreqs_[i]/sum_aa;
	};
}

//Computes gapless Gumbel parameters:
void FrameshiftAlignmentEvaluer::initGapless(long ntAlphabetSize_,//a number of letters in the DNA alphabet; usually 4
long aaAlphabetSize_,//a number of letters in the amino acid alphabet
const long *codonTable_,//defines the codon table; is an array of ntAlphabetSize_^3 integers

const long *const *substitutionScoreMatrix_,//scoring matrix; aaAlphabetSize_ X aaAlphabetSize_
const double *ntFreqs_,//background frequencies of letters in DNA sequences
const double *aaFreqs_,//background frequencies of letters in amino acid sequences
double max_time_)//maximum allowed calculation time in seconds
{

	double *RR1_AA=NULL;

	try
	{


		double CurrentTime1;
		Sls::FSA_utils::get_current_time(CurrentTime1);

		//check correctness of the input parameters for gapless alignment
		string function_name="void FrameshiftAlignmentEvaluer::initGapless";
		double *ntFreqs_normalized=NULL;//normalized background frequencies of letters in DNA sequences
		double *aaFreqs_normalized=NULL;//normalized background frequencies of letters in amino acid sequences
		assert_Gapless_input_parameters(
		ntAlphabetSize_,//a number of letters in the DNA alphabet; usually 4
		aaAlphabetSize_,//a number of letters in the amino acid alphabet
		codonTable_,//defines the codon table; is an array of ntAlphabetSize_^3 integers
		ntFreqs_,//background frequencies of letters in DNA sequences
		aaFreqs_,//background frequencies of letters in amino acid sequences
		ntFreqs_normalized,//normalized background frequencies of letters in DNA sequences
		aaFreqs_normalized,//normalized background frequencies of letters in amino acid sequences
		function_name);//"assert_Gapless_input_parameters" is called from "function_name_" function

		d_params.d_params_flag=false;

		long int codon_length=3;
//------

		FSA_utils::extract_AA_frequencies_for_DNA_sequence(
		codonTable_,//<codon code,AA number>
		codon_length,//codon length 
		ntAlphabetSize_,//number of letters for the sequence 1
		aaAlphabetSize_,//number of letters for the sequence 2
		ntFreqs_normalized,//nucleotide probabilities
		RR1_AA);//the resulted frequencies

		if(max_time_<=0)
		{
			max_time_=60;
		};
		
		Njn::LocalMaxStatMatrix local_max_stat_matrix(aaAlphabetSize_,
							  substitutionScoreMatrix_,
							  RR1_AA,
							  aaFreqs_normalized,
							  aaAlphabetSize_,
							  max_time_);


		if(local_max_stat_matrix.getTerminated()) 
		{
			throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
		};


		//calculation of a and sigma
		double calculation_error=1e-6;

		d_params.gapless_alpha_J = local_max_stat_matrix.getAlpha ();
		d_params.gapless_alpha_J=FSA_utils::Tmax(d_params.gapless_alpha_J,0.0);
		d_params.gapless_alpha_J_error = calculation_error;

		d_params.gapless_alpha_I=d_params.gapless_alpha_J*9;
		d_params.gapless_alpha_I_error = calculation_error*9;

		d_params.gapless_sigma=d_params.gapless_alpha_J*3;
		d_params.gapless_sigma_error = calculation_error*3;


		d_params.gapless_a_J = local_max_stat_matrix.getA ();
		d_params.gapless_a_J=FSA_utils::Tmax(d_params.gapless_a_J,0.0);
		d_params.gapless_a_J_error = calculation_error;

		d_params.gapless_a_I = d_params.gapless_a_J*3;
		d_params.gapless_a_I_error = calculation_error*3;

		//calculation of all required parameters for a gapless case
		d_params.G=0;
		d_params.G1=0;
		d_params.G2=0;

		d_params.lambda = local_max_stat_matrix.getLambda ();
		d_params.lambda_error = calculation_error;

		d_params.K = local_max_stat_matrix.getK ();
		d_params.K_error = calculation_error;
			
		d_params.C = local_max_stat_matrix.getC ();
		d_params.C_error = calculation_error;

		if(d_params.C!=0)
		{
			d_params.K_C = d_params.K/d_params.C;
		}
		else
		{
			d_params.K_C = 0;
		};
		d_params.K_C_error = calculation_error;


		d_params.sigma = d_params.gapless_sigma;
		d_params.sigma_error = d_params.gapless_sigma_error;

		d_params.alpha_I = d_params.gapless_alpha_I;
		d_params.alpha_I_error = d_params.gapless_alpha_I_error;

		d_params.alpha_J = d_params.gapless_alpha_J;
		d_params.alpha_J_error = d_params.gapless_alpha_J_error;

		d_params.a_I = d_params.gapless_a_I;
		d_params.a_I_error = d_params.gapless_a_I_error;

		d_params.a_J = d_params.gapless_a_J;
		d_params.a_J_error = d_params.gapless_a_J_error;


		std::vector<double > sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.lambda;
		sbs_arrays[1]=d_params.lambda + calculation_error;

		d_params.m_LambdaSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.K;
		sbs_arrays[1]=d_params.K+calculation_error;

		d_params.m_KSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.C;
		sbs_arrays[1]=d_params.C+calculation_error;

		d_params.m_CSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.sigma;
		sbs_arrays[1]=d_params.sigma + calculation_error;

		d_params.m_SigmaSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.alpha_I;
		sbs_arrays[1]=d_params.alpha_I + calculation_error;

		d_params.m_AlphaISbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.alpha_J;
		sbs_arrays[1]=d_params.alpha_J + calculation_error;

		d_params.m_AlphaJSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.a_I;
		sbs_arrays[1]=d_params.a_I + calculation_error;

		d_params.m_AISbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.a_J;
		sbs_arrays[1]=d_params.a_J + calculation_error;

		d_params.m_AJSbs=sbs_arrays;

		d_params.d_params_flag=true;

		//precompute intercepts
		FALP_pvalues::compute_intercepts(d_params);

		double CurrentTime2;
		Sls::FSA_utils::get_current_time(CurrentTime2);
		d_params.m_CalcTime=CurrentTime2-CurrentTime1;

		delete[]RR1_AA;RR1_AA=NULL;

		if(!FALP_pvalues::assert_Gumbel_parameters(
		d_params)||!isGood())
		{
			d_params.d_params_flag=false;
			throw error("Error - computation of the Gumbel parameters is unsuccessful in the function \"void FrameshiftAlignmentEvaluer::initGapless\"\n",1);
		};

		delete[]ntFreqs_normalized;
		delete[]aaFreqs_normalized;

	}
	catch (...)
	{ 
		d_params.d_params_flag=false;
		delete[]RR1_AA;RR1_AA=NULL;
		throw;
	};
}


//Computes gapped Gumbel parameters:
//The NCBI convention is used for penalizing a gap:
//For example, a gap of length k is penalized as gapOpen1_+k*gapEpen1_ for sequence #1
void FrameshiftAlignmentEvaluer::initFrameshift(long ntAlphabetSize_,//a number of letters in the DNA alphabet; usually 4
long aaAlphabetSize_,//a number of letters in the amino acid alphabet; usually 4
const long *codonTable_,//defines the codon table; is an array of ntAlphabetSize_^3 integers
//letters are numerated from 0; 
//each codon defines a code (an integer number from [0,...,ntAlphabetSize_^3 integers];
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
long int seq_number_,//a number of alignments
long int nalp_,//a number of ascending ladder points used in the calculation
long int number_of_subsets_for_errors_calculation_)//a number of subsets used in the error calculation
{
	try
	{

		double CurrentTime1;
		Sls::FSA_utils::get_current_time(CurrentTime1);

		//check the input parameters
		string function_name="void FrameshiftAlignmentEvaluer::initFrameshift";
		double *ntFreqs_normalized=NULL;//normalized background frequencies of letters in DNA sequences
		double *aaFreqs_normalized=NULL;//normalized background frequencies of letters in amino acid sequences
		assert_Gapless_input_parameters(
		ntAlphabetSize_,//a number of letters in the DNA alphabet; usually 4
		aaAlphabetSize_,//a number of letters in the amino acid alphabet
		codonTable_,//defines the codon table; is an array of ntAlphabetSize_^3 integers
		ntFreqs_,//background frequencies of letters in DNA sequences
		aaFreqs_,//background frequencies of letters in amino acid sequences
		ntFreqs_normalized,//normalized background frequencies of letters in DNA sequences
		aaFreqs_normalized,//normalized background frequencies of letters in amino acid sequences
		function_name);//"assert_Gapless_input_parameters" is called from "function_name_" function


		if(!(gapEpen1_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"gapEpen1_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(gapEpen2_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"gapEpen2_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(frameshiftCost_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"frameshiftCost_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(eps_lambda_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"eps_lambda_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(eps_K_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"eps_K_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(max_mem_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"max_mem_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(seq_number_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"seq_number_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(nalp_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"nalp_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(number_of_subsets_for_errors_calculation_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"number_of_subsets_for_errors_calculation_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		d_params.d_params_flag=false;

		bool library_call_flag=true;//if true, then the additional parameters are used
		string smatr_file_name="";//scoring matrix file name
		string RR1_file_name="";//background frequencies file name for the sequence #1
		string RR2_file_name="";//background frequencies file name for the sequence #2
		string DNA_codon_table_file_name="";//a name of a file with DNA codon table

		bool gapped_flag=true;//if true, then the gapped alingment is performed


		bool forward_and_reverse_screen_output_flag=false;//determines whether the parameters are outputted for forward and reverse calculations
		double mult_for_is_lambda=1;//multiplier for lambda in the IS


		test test1;

		test1.FSA_IS(
		//additional parameters for the library code
		library_call_flag,//if true, then the additional parameters are used
		ntAlphabetSize_,//a number of letters in the DNA alphabet; usually 4
		aaAlphabetSize_,//a number of letters in the amino acid alphabet; usually 4
		codonTable_,//defines the codon table; is an array of ntAlphabetSize_^3 integers
		substitutionScoreMatrix_,//scoring matrix; aaAlphabetSize_ X aaAlphabetSize_
		ntFreqs_normalized,//background frequencies of letters in DNA sequences
		aaFreqs_normalized,//background frequencies of letters in amino acid sequences
		//additional parameters for the library code - end

		randomSeed_,//randomization number
		gapOpen1_,//gap opening penalty for the nucleotide sequence #1
		gapOpen2_,//gap opening penalty for the amino acid sequence #2

		gapEpen1_,//gap extension penalty for the nucleotide sequence #1
		gapEpen2_,//gap extension penalty for the amino acid sequence #2

		frameshiftCost_,//frameshift penalty gamma


		smatr_file_name,//scoring matrix file name
		RR1_file_name,//background frequencies file name for the sequence #1
		RR2_file_name,//background frequencies file name for the sequence #2
		DNA_codon_table_file_name,//a name of a file with DNA codon table

		eps_lambda_,//relative error for lambda calculation
		eps_K_,//relative error for K calculation

		gapped_flag,//if true, then the gapped alingment is performed
		max_time_,//maximum allowed calculation time in seconds
		max_mem_,//maximum allowed memory usage in MB

		seq_number_,//number of tested alignments
		nalp_,//number of ALPs for the calculation
		number_of_subsets_for_errors_calculation_,//number of subsets used for the splitting method

		forward_and_reverse_screen_output_flag,//determines whether the parameters are outputted for forward and reverse calculations
		insertions_after_deletions_,//if true, then insertions after deletions are allowed

		//for test
		mult_for_is_lambda,//multiplier for lambda in the IS

		//the result
		d_params);//the resulted parameters


		d_params.d_params_flag=true;

		//precompute intercepts
		//computed in test1.FSA_IS

		delete[]ntFreqs_normalized;
		delete[]aaFreqs_normalized;


		double CurrentTime2;
		Sls::FSA_utils::get_current_time(CurrentTime2);
		d_params.m_CalcTime=CurrentTime2-CurrentTime1;

		if(!FALP_pvalues::assert_Gumbel_parameters(
		d_params)||!isGood())
		{
			d_params.d_params_flag=false;
			throw error("Error - computation of the Gumbel parameters is unsuccessful in the function \"void FrameshiftAlignmentEvaluer::initFrameshift\"\n",1);
		};



	}
	catch (...)
	{ 
		d_params.d_params_flag=false;
		throw;
	};

}

//Initializes Gumbel parameters using precalculated values:
void FrameshiftAlignmentEvaluer::initParameters(
const AlignmentEvaluerParameters &parameters_)
{
	try
	{
		double CurrentTime1;
		Sls::FSA_utils::get_current_time(CurrentTime1);

		d_params.d_params_flag=false;

		double calculation_error=1e-6;

		d_params.lambda=parameters_.d_lambda;
		d_params.lambda_error=calculation_error;

		d_params.C=0;
		d_params.C_error=0;

		d_params.lambda_last_ALP_relative_error=0.0;


		d_params.K_C=0.0;
		d_params.K_C_error=0.0;

		d_params.K=parameters_.d_k;
		d_params.K_error=calculation_error;

		d_params.a_I=parameters_.d_a1;
		d_params.a_I_error=calculation_error;

		d_params.a_J=parameters_.d_a2;
		d_params.a_J_error=calculation_error;

		d_params.sigma=parameters_.d_sigma;
		d_params.sigma_error=calculation_error;

		d_params.alpha_I=parameters_.d_alpha1;
		d_params.alpha_I_error=calculation_error;

		d_params.alpha_J=parameters_.d_alpha2;
		d_params.alpha_J_error=calculation_error;


		d_params.gapless_a_I=0.0;
		d_params.gapless_a_I_error=0.0;

		d_params.gapless_a_J=0.0;
		d_params.gapless_a_J_error=0.0;



		d_params.gapless_alpha_I=0.0;
		d_params.gapless_alpha_I_error=0.0;

		d_params.gapless_alpha_J=0.0;
		d_params.gapless_alpha_J_error=0.0;

		d_params.gapless_sigma=0.0;
		d_params.gapless_sigma_error=0.0;

		d_params.G=0;
		d_params.G1=0;
		d_params.G2=0;

		d_params.realizations_number=0;

		d_params.m_CalcTime=0;

		//intercepts
		d_params.b_I=parameters_.d_b1;
		d_params.b_I_error=calculation_error;

		d_params.b_J=parameters_.d_b2;
		d_params.b_J_error=calculation_error;

		d_params.beta_I=parameters_.d_beta1;
		d_params.beta_I_error=calculation_error;

		d_params.beta_J=parameters_.d_beta2;
		d_params.beta_J_error=calculation_error;

		d_params.tau=parameters_.d_tau;
		d_params.tau_error=calculation_error;


		//arrays initialization
		std::vector<double > sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.lambda;
		sbs_arrays[1]=d_params.lambda + calculation_error;

		d_params.m_LambdaSbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.K;
		sbs_arrays[1]=d_params.K + calculation_error;

		d_params.m_KSbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.C;
		sbs_arrays[1]=d_params.C + calculation_error;

		d_params.m_CSbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.sigma;
		sbs_arrays[1]=d_params.sigma + calculation_error;

		d_params.m_SigmaSbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.alpha_I;
		sbs_arrays[1]=d_params.alpha_I + calculation_error;

		d_params.m_AlphaISbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.alpha_J;
		sbs_arrays[1]=d_params.alpha_J + calculation_error;

		d_params.m_AlphaJSbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.a_I;
		sbs_arrays[1]=d_params.a_I + calculation_error;

		d_params.m_AISbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.a_J;
		sbs_arrays[1]=d_params.a_J + calculation_error;

		d_params.m_AJSbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.b_I;
		sbs_arrays[1]=d_params.b_I + calculation_error;

		d_params.m_BISbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.b_J;
		sbs_arrays[1]=d_params.b_J + calculation_error;

		d_params.m_BJSbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.beta_I;
		sbs_arrays[1]=d_params.beta_I + calculation_error;

		d_params.m_BetaISbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.beta_J;
		sbs_arrays[1]=d_params.beta_J + calculation_error;

		d_params.m_BetaJSbs=sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.tau;
		sbs_arrays[1]=d_params.tau + calculation_error;

		d_params.m_TauSbs=sbs_arrays;

		d_params.d_params_flag=true;//if true, then the parameters are defined and P-values can be calculated

		FALP_pvalues::compute_tmp_values(d_params);

		if(!FALP_pvalues::assert_Gumbel_parameters(
		d_params)||!isGood())
		{
			d_params.d_params_flag=false;
			throw error("Error - computation of the Gumbel parameters is unsuccessful in the function \"void FrameshiftAlignmentEvaluer::initParameters\"\n",1);
		};

		double CurrentTime2;
		Sls::FSA_utils::get_current_time(CurrentTime2);
		d_params.m_CalcTime=CurrentTime2-CurrentTime1;

	}
	catch (...)
	{ 
		d_params.d_params_flag=false;
		throw;
	};

}

void FrameshiftAlignmentEvaluer::initParameters(
const AlignmentEvaluerParametersWithErrors &parameters_)
{
	try
	{
		double CurrentTime1;
		Sls::FSA_utils::get_current_time(CurrentTime1);

		d_params.d_params_flag=false;

		long int array_dim=20;

		long int seed_tmp=12345;
		srand(seed_tmp);

		d_params.lambda=parameters_.d_lambda;
		d_params.lambda_error=parameters_.d_lambda_error;

		d_params.C=0;
		d_params.C_error=0;

		d_params.lambda_last_ALP_relative_error=0.0;


		d_params.K_C=0.0;
		d_params.K_C_error=0.0;

		d_params.K=parameters_.d_k;
		d_params.K_error=parameters_.d_k_error;

		d_params.a_I=parameters_.d_a1;
		d_params.a_I_error=parameters_.d_a1_error;

		d_params.a_J=parameters_.d_a2;
		d_params.a_J_error=parameters_.d_a2_error;

		d_params.sigma=parameters_.d_sigma;
		d_params.sigma_error=parameters_.d_sigma_error;

		d_params.alpha_I=parameters_.d_alpha1;
		d_params.alpha_I_error=parameters_.d_alpha1_error;

		d_params.alpha_J=parameters_.d_alpha2;
		d_params.alpha_J_error=parameters_.d_alpha2_error;


		d_params.gapless_a_I=0.0;
		d_params.gapless_a_I_error=0.0;

		d_params.gapless_a_J=0.0;
		d_params.gapless_a_J_error=0.0;



		d_params.gapless_alpha_I=0.0;
		d_params.gapless_alpha_I_error=0.0;

		d_params.gapless_alpha_J=0.0;
		d_params.gapless_alpha_J_error=0.0;

		d_params.gapless_sigma=0.0;
		d_params.gapless_sigma_error=0.0;

		d_params.G=0;
		d_params.G1=0;
		d_params.G2=0;

		d_params.realizations_number=0;

		d_params.m_CalcTime=0;

		//intercepts
		d_params.b_I=parameters_.d_b1;
		d_params.b_I_error=parameters_.d_b1_error;

		d_params.b_J=parameters_.d_b2;
		d_params.b_J_error=parameters_.d_b2_error;

		d_params.beta_I=parameters_.d_beta1;
		d_params.beta_I_error=parameters_.d_beta1_error;

		d_params.beta_J=parameters_.d_beta2;
		d_params.beta_J_error=parameters_.d_beta2_error;

		d_params.tau=parameters_.d_tau;
		d_params.tau_error=parameters_.d_tau_error;

		double sqrt_array_dim=sqrt((double)array_dim);

		d_params.m_LambdaSbs.clear();
		d_params.m_KSbs.clear();
		d_params.m_CSbs.clear();

		d_params.m_SigmaSbs.clear();

		d_params.m_AlphaISbs.clear();
		d_params.m_AlphaJSbs.clear();

		d_params.m_AISbs.clear();
		d_params.m_AJSbs.clear();

		d_params.m_BISbs.clear();
		d_params.m_BJSbs.clear();

		d_params.m_BetaISbs.clear();
		d_params.m_BetaJSbs.clear();

		d_params.m_TauSbs.clear();

		long int i;
		for(i=0;i<array_dim;i++)
		{
			d_params.m_LambdaSbs.push_back(d_params.lambda+d_params.lambda_error*FALP_pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_KSbs.push_back(d_params.K+d_params.K_error*FALP_pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_CSbs.push_back(0);

			d_params.m_SigmaSbs.push_back(d_params.sigma+d_params.sigma_error*FALP_pvalues::standard_normal()*sqrt_array_dim);

			d_params.m_AlphaISbs.push_back(d_params.alpha_I+d_params.alpha_I_error*FALP_pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_AlphaJSbs.push_back(d_params.alpha_J+d_params.alpha_J_error*FALP_pvalues::standard_normal()*sqrt_array_dim);

			d_params.m_AISbs.push_back(d_params.a_I+d_params.a_I_error*FALP_pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_AJSbs.push_back(d_params.a_J+d_params.a_J_error*FALP_pvalues::standard_normal()*sqrt_array_dim);

			d_params.m_BISbs.push_back(d_params.b_I+d_params.b_I_error*FALP_pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_BJSbs.push_back(d_params.b_J+d_params.b_J_error*FALP_pvalues::standard_normal()*sqrt_array_dim);

			d_params.m_BetaISbs.push_back(d_params.beta_I+d_params.beta_I_error*FALP_pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_BetaJSbs.push_back(d_params.beta_J+d_params.beta_J_error*FALP_pvalues::standard_normal()*sqrt_array_dim);

			d_params.m_TauSbs.push_back(d_params.tau+d_params.tau_error*FALP_pvalues::standard_normal()*sqrt_array_dim);
		};

		d_params.d_params_flag=true;//if true, then the parameters are defined and P-values can be calculated

		FALP_pvalues::compute_tmp_values(d_params);

		if(!FALP_pvalues::assert_Gumbel_parameters(
		d_params)||!isGood())
		{
			d_params.d_params_flag=false;
			throw error("Error - computation of the Gumbel parameters is unsuccessful in the function \"void FrameshiftAlignmentEvaluer::initParameters\"\n",1);
		};

		double CurrentTime2;
		Sls::FSA_utils::get_current_time(CurrentTime2);
		d_params.m_CalcTime=CurrentTime2-CurrentTime1;

	}
	catch (...)
	{ 
		d_params.d_params_flag=false;
		throw;
	};

}


double FrameshiftAlignmentEvaluer::area(double score_,//pairwise alignment score
double seqlen1_,//length of sequence #1
double seqlen2_) const//length of sequence #2
{
	if(seqlen1_<=0||seqlen2_<=0)
	{
		throw error("Error - seqlen1_<=0 or seq2en1_<=0 in \"double FrameshiftAlignmentEvaluer::area\"\n",2);
	};

	if(!isGood())
	{
		throw error("Unexpected error - the Gumbel parameters are not defined properly in \"double FrameshiftAlignmentEvaluer::area\"\n",1);
	};

	static Sls::FALP_pvalues pvalues_obj;

	double P;
	double E;
	double area_res;
	bool area_is_1_flag=false;
	bool compute_only_area=true;

	pvalues_obj.get_appr_tail_prob_with_cov_without_errors(
	d_params,
	pvalues_obj.blast,
	score_,
	seqlen1_,
	seqlen2_,

	P,

	E,

	area_res,
	pvalues_obj.a_normal,
	pvalues_obj.b_normal,
	pvalues_obj.h_normal,
	pvalues_obj.N_normal,
	pvalues_obj.p_normal,
	area_is_1_flag,
	compute_only_area);


	return area_res;

}


void FrameshiftAlignmentEvaluer::calc(double score_,//pairwise alignment score
double seqlen1_,//length of sequence #1
double seqlen2_,//length of sequence #2
double &pvalue_,//resulted P-value
double &pvalueErr_,//P-value error
double &evalue_,//resulted E-value
double &evalueErr_) const//E-value error
{

	if(seqlen1_<=0||seqlen2_<=0)
	{
		throw error("Error - seqlen1_<=0 or seqlen2_<=0 in \"double FrameshiftAlignmentEvaluer::calc\"\n",2);
	};

	if(!isGood())
	{
		throw error("Unexpected error - the Gumbel parameters are not defined properly in \"double FrameshiftAlignmentEvaluer::calc\"\n",1);
	};

	static Sls::FALP_pvalues pvalues_obj;

	pvalues_obj.calculate_P_values(
		score_, seqlen1_, seqlen2_,
		d_params, 
		pvalue_,
		pvalueErr_,
		evalue_,
		evalueErr_);

}

void FrameshiftAlignmentEvaluer::calc(double score_,//pairwise alignment score
double seqlen1_,//length of sequence #1
double seqlen2_,//length of sequence #2
double &pvalue_,//resulted P-value
double &evalue_) const//resulted E-value
{
	if(seqlen1_<=0||seqlen2_<=0)
	{
		throw error("Error - seqlen1_<=0 or seqlen2_<=0 in \"double FrameshiftAlignmentEvaluer::calc\"\n",2);
	};

	if(!isGood())
	{
		throw error("Unexpected error - d_params is not defined in \"double FrameshiftAlignmentEvaluer::calc\"\n",1);
	};

	static Sls::FALP_pvalues pvalues_obj;

	bool area_is_1_flag=false;

	double area;


	pvalues_obj.get_appr_tail_prob_with_cov_without_errors(
	d_params,
	pvalues_obj.blast,
	score_,
	seqlen1_,
	seqlen2_,

	pvalue_,

	evalue_,

	area,
	pvalues_obj.a_normal,
	pvalues_obj.b_normal,
	pvalues_obj.h_normal,
	pvalues_obj.N_normal,
	pvalues_obj.p_normal,
	area_is_1_flag);
}

}


