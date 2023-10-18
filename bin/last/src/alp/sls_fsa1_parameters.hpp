#ifndef INCLUDED_SLS_FSA_PARAMETERS
#define INCLUDED_SLS_FSA_PARAMETERS

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

File name: sls_fsa1_parameters.hpp

Author: Sergey Sheetlin

Contents: Calculation of Gumbel parameters

******************************************************************************/

#include "sls_fsa1_utils.hpp"
#include "sls_alp_regression.hpp"
#include "sls_fsa1_pvalues.hpp"



namespace Sls { 

	struct struct_for_lambda_calculation
	{
		void **d_alp_distr;
		void **d_alp_distr_errors;
		long int d_nalp;
		double d_f_error;

		double d_last_sum;
		double d_last_sum_error;

		double d_before_last_sum;
		double d_before_last_sum_error;

		bool d_calculate_alp_number;
		long int d_alp_number;

		
	};


	class fsa_par{


		public:

		static double function_for_lambda_calculation(
		double lambda_,
		void * data_);

		static void calculate_lambda(
		bool check_the_criteria_,
		long int nalp_,
		long int &nalp_thr_,
		bool &inside_simulation_flag_,
		void **alp_distr,
		void **alp_distr_errors,
		double ungapped_lambda_,
		double &lambda_,
		double &lambda_error_);

		static void calculate_C(
		long int starting_point,
		long int nalp_,
		void **alp_distr,
		void **alp_distr_errors,
		double lambda_,
		double lambda_error_,
		double &C_,
		double &C_error_,
		bool &inside_simulation_flag_);

		static void memory_releace_for_calculate_C(
		double *&P,
		double *&P_errors,
		double *&values_P_ratio,
		double *&errors_P_ratio,

		double *&E,
		double *&E_errors,

		double *&E_T_beta,
		double *&E_T_beta_errors);

		static double get_root(
		const std::vector<double> &res_tmp_,
		double point_);

		static void memory_release_for_calculate_FSC(
		double *&exp_array,
		double *&delta_E,
		double *&delta_E_error,
		double *&delta_E_E,
		double *&delta_E_E_error,
		double *&delta_I,
		double *&delta_I_error,
		double *&delta_J,
		double *&delta_J_error,
		double *&delta_I_I,
		double *&delta_I_I_error,
		double *&delta_I_J,
		double *&delta_I_J_error,
		double *&delta_J_J,
		double *&delta_J_J_error,
		double *&cov_J_J,
		double *&cov_J_J_error,
		double *&cov_I_J,
		double *&cov_I_J_error,
		double *&cov_I_I,
		double *&cov_I_I_error,
		double *&cov_E_E,
		double *&cov_E_E_error,
		double *&C_S,
		double *&C_S_error);

		static void calculate_FSC(
		long int nalp_,
		long int ind1_,
		long int ind2_,

		double lambda_,

		long int M_max_,

		array_positive<long int> *distance_along_direction_1_,
		array_positive<long int> *distance_along_direction_2_,

		array_positive<double> *ALP_weight_,
		array_positive<long int> *ALP_edge_max_,

		double &a_I_,
		double &a_I_error_,
		double &a_J_,
		double &a_J_error_,
		double &sigma_,
		double &sigma_error_,
		double &alpha_I_,
		double &alpha_I_error_,
		double &alpha_J_,
		double &alpha_J_error_,
		bool &inside_simulation_flag_,
		double *test_FSC_vect_=NULL);//for tests

		static void sigma_calculation(
		double delta_I_aver_,
		double delta_I_aver_error_,
		double delta_J_aver_,
		double delta_J_aver_error_,
		double delta_E_aver_,
		double delta_E_aver_error_,
		double cov_E_E_aver_,
		double cov_E_E_aver_error_,
		double cov_I_J_aver_,
		double cov_I_J_aver_error_,
		double &sigma_,
		double &sigma_error_);

		static double lambda_exp(
		long int &i_,
		double *&exp_array_);

		static void Read_Params(
		Sls::FALP_set_of_parameters &gumbel_params_,
		std::string gumbelparin_file_name_);

		static void Output_Params(
		Sls::FALP_set_of_parameters &gumbel_params_,
		std::string gumbelparout_file_name_);

		friend std::ostream &operator<<(std::ostream &s_,
		const Sls::FALP_set_of_parameters &gumbel_params_);

		friend std::istream &operator>>(std::istream &s_,
		FALP_set_of_parameters &gumbel_params_);



		static void Output_Pvalues(
		char type_,//'E' or 'P'
		long int score1_,
		long int score2_,
		std::vector<double> &pv_,
		std::vector<double> &pv_err_,
		std::string pvalout_file_name_);//P-values file name






		public:

		long int d_tmp;


	};
}


#endif

