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

File name: sls_fsa1_parameters.cpp

Author: Sergey Sheetlin

Contents: Calculation of Gumbel parameters

******************************************************************************/


#include "sls_fsa1_parameters.hpp"

using namespace Sls;
using namespace std;

//parameters calculation
double fsa_par::function_for_lambda_calculation(
double lambda_,
void * data_)
{

	double *expect=NULL;
	double *expect_errors=NULL;

	try
	{

		struct_for_lambda_calculation *tmp_struct=(struct_for_lambda_calculation *)data_;
		void **alp_distr=tmp_struct->d_alp_distr;
		void **alp_distr_errors=tmp_struct->d_alp_distr_errors;
		long int nalp=tmp_struct->d_nalp;

		expect=new double[nalp];
		FSA_utils::assert_mem(expect);
		expect_errors=new double[nalp];
		FSA_utils::assert_mem(expect_errors);

		if(nalp<1)
		{
			throw error("Unexpected error\n",4);
		};



		long int k;
		for(k=1;k<=nalp;k++)
		{
			array_v<double>* tmp=((array_v<double>*)alp_distr[k]);
			array_v<double>* tmp_errors=((array_v<double>*)alp_distr_errors[k]);

			double val=0;
			double val_error=0;

			long int j;
			for(j=0;j<=tmp->d_dim;j++)
			{
				if(tmp->d_elem[j]<=0)
				{
					continue;
				};
				double exp_tmp=exp(lambda_*(j+tmp->d_ind0));
				val+=exp_tmp*tmp->d_elem[j];
				val_error+=exp_tmp*exp_tmp*tmp_errors->d_elem[j];
			};
			val_error=alp_reg::sqrt_for_errors(val_error);
			expect[k-1]=val;
			expect_errors[k-1]=val_error;

		};

		tmp_struct->d_last_sum=expect[nalp-1];
		tmp_struct->d_last_sum_error=expect_errors[nalp-1];

		if(nalp==1)
		{
			tmp_struct->d_before_last_sum=1.0;
			tmp_struct->d_before_last_sum_error=0.0;
		}
		else
		{
			tmp_struct->d_before_last_sum=expect[nalp-2];
			tmp_struct->d_before_last_sum_error=expect_errors[nalp-2];
		};

		if(tmp_struct->d_calculate_alp_number)
		{
			double tmp=0.0;
			long int k;
			for(k=0;k<nalp;k++)
			{
				if(expect_errors[k]!=0)
				{
					tmp+=1.0/(expect_errors[k]*expect_errors[k]);
				};

			};

			long int tmp_alp=nalp;
			double tmp1=0.0;
			for(k=nalp-1;k>=0;k--)
			{
				if(expect_errors[k]!=0)
				{
					tmp1+=1.0/(expect_errors[k]*expect_errors[k]);
				};
				if(tmp1>0.2*tmp)
				{
					tmp_alp=k+1;
					break;
				};
			};

			tmp_struct->d_alp_number=tmp_alp;
		};

		if(nalp==1)
		{
			double tmp=expect[0]-1.0;
			tmp_struct->d_f_error=expect_errors[0];
			delete[]expect;expect=NULL;
			delete[]expect_errors;expect_errors=NULL;
			return tmp;
		};


		long int min_length=0;
		long int number_of_elements=nalp;
		bool cut_left_tail=true;
		bool cut_right_tail=false;
		double y=2;
		double beta0;
		double beta1;
		double beta0_error;
		double beta1_error;
		long int k1_opt;
		long int k2_opt;

		bool res_was_calculated;

		alp_reg::robust_regression_sum_with_cut_LSM(
		min_length,
		number_of_elements,
		expect,
		expect_errors,
		cut_left_tail,
		cut_right_tail,
		y,
		beta0,
		beta1,
		beta0_error,
		beta1_error,
		k1_opt,
		k2_opt,
		res_was_calculated);



		if(!res_was_calculated)
		{
			throw error("",2);
		};


		delete[]expect;expect=NULL;
		delete[]expect_errors;expect_errors=NULL;

		tmp_struct->d_f_error=beta1_error;
		return beta1;

	}
	catch (...)
	{ 
		delete[]expect;expect=NULL;
		delete[]expect_errors;expect_errors=NULL;
		throw;
	};

}

void fsa_par::calculate_lambda(
bool check_the_criteria_,
long int nalp_,
long int &nalp_thr_,
bool &inside_simulation_flag_,
void **alp_distr,
void **alp_distr_errors,
double ungapped_lambda_,
double &lambda_,
double &lambda_error_)
{

	bool ee_error_flag=false;
	error ee_error("",0);

	inside_simulation_flag_=false;

	try
	{

		long int nalp=nalp_;

		if(nalp<=0)
		{
			throw error("Unexpected error\n",4);
		};

		

		struct_for_lambda_calculation tmp_struct;
		tmp_struct.d_alp_distr=alp_distr;
		tmp_struct.d_alp_distr_errors=alp_distr_errors;
		tmp_struct.d_nalp=nalp;
		tmp_struct.d_calculate_alp_number=false;


		function_type *func=function_for_lambda_calculation;
		void* func_pointer=&tmp_struct;
		double a=0;
		double b=ungapped_lambda_*2;
		//double b=ungapped_lambda_*10;
		long int n_partition=30;
		//long int n_partition=300;
		double eps=1e-10;
		std::vector<double> res;

		alp_reg::find_tetta_general(
		func,
		func_pointer,
		a,//[a,b] is the interval for search of equation solution
		b,
		n_partition,
		eps,
		res);

		


		inside_simulation_flag_=true;
		if(res.size()==0)
		{
			inside_simulation_flag_=false;
			return;
		};

		

		lambda_=get_root(res,ungapped_lambda_);
		

		tmp_struct.d_calculate_alp_number=true;
		double f1=func(lambda_,func_pointer);
		nalp_thr_=tmp_struct.d_alp_number;
		tmp_struct.d_calculate_alp_number=false;

		double slope_error=tmp_struct.d_f_error;

		double delta_lambda=lambda_/100.0;
		double f2=func(lambda_+delta_lambda,func_pointer);
		
		

		if(delta_lambda==0||f1==f2)
		{
			lambda_error_=0.0;
		}
		else
		{
			double derivative=(f2-f1)/delta_lambda;
			lambda_error_=fabs(slope_error/derivative);
		};
		
		if(!check_the_criteria_)
		{
			return;
		};


	}
	catch (error er)
	{
		ee_error_flag=true;
		ee_error=er;		
	};


	if(ee_error_flag)
	{
		if(ee_error.st!="")
		{
			throw error(ee_error.st,ee_error.error_code);
		}
		else
		{
			inside_simulation_flag_=false;
		};
	};

}

void fsa_par::memory_releace_for_calculate_C(
double *&P,
double *&P_errors,
double *&values_P_ratio,
double *&errors_P_ratio,

double *&E,
double *&E_errors,

double *&E_T_beta,
double *&E_T_beta_errors)
{
	delete[]values_P_ratio;values_P_ratio=NULL;
	delete[]errors_P_ratio;errors_P_ratio=NULL;


	delete[]P;P=NULL;
	delete[]P_errors;P_errors=NULL;

	delete[]E;E=NULL;
	delete[]E_T_beta;E_T_beta=NULL;
	delete[]E_errors;E_errors=NULL;
	delete[]E_T_beta_errors;E_T_beta_errors=NULL;
}

void fsa_par::calculate_C(
long int starting_point,
long int nalp_,
void **alp_distr,
void **alp_distr_errors,
double lambda_,
double lambda_error_,
double &C_,
double &C_error_,
bool &inside_simulation_flag_)
{
	inside_simulation_flag_=true;

	error ee_error("",0);

	double *P=NULL;
	double *P_errors=NULL;
	double *values_P_ratio=NULL;
	double *errors_P_ratio=NULL;

	double *E=NULL;
	double *E_errors=NULL;

	double *E_T_beta=NULL;
	double *E_T_beta_errors=NULL;


	try
	{
	try
	{

		long int total_number_of_ALP=nalp_;

		if(total_number_of_ALP<1)
		{
			throw error("Unexpected error\n",4);
		};


		//1)P(beta=infinity)
		long int j;


		P=new double[total_number_of_ALP+1];
		FSA_utils::assert_mem(P);
		P_errors=new double[total_number_of_ALP+1];
		FSA_utils::assert_mem(P_errors);

		P[0]=1.0;
		P_errors[0]=0.0;

		
		for(j=1;j<=total_number_of_ALP;j++)
		{
			array_v<double>* tmp=((array_v<double>*)alp_distr[j]);
			array_v<double>* tmp_errors=((array_v<double>*)alp_distr_errors[j]);

			P[j]=0;
			P_errors[j]=0;
			long int i;
			for(i=0;i<=tmp->d_dim;i++)
			{
				P[j]+=tmp->d_elem[i];
				P_errors[j]+=tmp_errors->d_elem[i];
			};

			P_errors[j]=alp_reg::sqrt_for_errors(P_errors[j]);
		};

		

		values_P_ratio=new double[total_number_of_ALP];
		FSA_utils::assert_mem(values_P_ratio);
		errors_P_ratio=new double[total_number_of_ALP];
		FSA_utils::assert_mem(errors_P_ratio);

		

		for(j=0;j<total_number_of_ALP;j++)
		{
			if(P[j]<=0)
			{
				inside_simulation_flag_=false;
				throw error("",1);
			};
			values_P_ratio[j]=P[j+1]/P[j];
			errors_P_ratio[j]=FSA_utils::error_of_the_ratio(P[j+1],P_errors[j+1],P[j],P_errors[j]);
		};



		double beta1=0;
		double beta1_error=0;

		long int number_of_elements=total_number_of_ALP;

		bool cut_left_tail=true;
		bool cut_right_tail=false;

		double y=2;

		long int k1_opt;
		long int k2_opt;


		double P_beta_inf;
		double P_beta_inf_error=0;

		bool res_was_calculated;
		
		alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
		0,
		number_of_elements-starting_point,
		values_P_ratio+starting_point,
		errors_P_ratio+starting_point,
		cut_left_tail,
		cut_right_tail,
		y,
		P_beta_inf,
		beta1,
		P_beta_inf_error,
		beta1_error,
		k1_opt,
		k2_opt,
		res_was_calculated);


		

		if(!res_was_calculated)
		{
			//throw error("The program cannot estimate the parameters; please repeat the calculation9\n",2);
			inside_simulation_flag_=false;
			throw error("",1);
		};

		P_beta_inf=1-P_beta_inf;

		

		
		//2)E(exp(lambda*T_beta)) and E(T_beta*exp(lambda*T_beta))
		E=new double[total_number_of_ALP+1];
		FSA_utils::assert_mem(E);
		E_errors=new double[total_number_of_ALP+1];
		FSA_utils::assert_mem(E_errors);

		E_T_beta=new double[total_number_of_ALP+1];
		FSA_utils::assert_mem(E_T_beta);
		E_T_beta_errors=new double[total_number_of_ALP+1];
		FSA_utils::assert_mem(E_T_beta);


		E[0]=1;
		E_T_beta[0]=0;

		E_errors[0]=0;
		E_T_beta_errors[0]=0;

		


		for(j=1;j<=total_number_of_ALP;j++)
		{
			array_v<double>* tmp=((array_v<double>*)alp_distr[j]);
			array_v<double>* tmp_errors=((array_v<double>*)alp_distr_errors[j]);

			E[j]=0;
			E_T_beta[j]=0;

			E_errors[j]=0;
			E_T_beta_errors[j]=0;

			long int i;
			for(i=0;i<=tmp->d_dim;i++)
			{
				double tmp_double=exp(lambda_*(double)i);
				E[j]+=tmp_double*tmp->d_elem[i];
				E_errors[j]+=tmp_double*tmp_double*tmp_errors->d_elem[i];

				tmp_double=(double)i*exp(lambda_*(double)i);
				E_T_beta[j]+=tmp_double*tmp->d_elem[i];
				E_T_beta_errors[j]+=tmp_double*tmp_double*tmp_errors->d_elem[i];
			};

			E_errors[j]=alp_reg::sqrt_for_errors(E_errors[j]);
			E_T_beta_errors[j]=alp_reg::sqrt_for_errors(E_T_beta_errors[j]);

		};


		double E_aver;
		double E_aver_error;

		double E_T_beta_diff_aver;
		double E_T_beta_diff_aver_error;


		if(total_number_of_ALP==1)
		{
			E_aver=E[1];
			E_aver_error=E_errors[1];

			E_T_beta_diff_aver=E_T_beta[1]-E_T_beta[0];
			E_T_beta_diff_aver_error=E_T_beta_errors[1];

		}
		else
		{
			long int number_of_elements=total_number_of_ALP;

			bool cut_left_tail=true;
			bool cut_right_tail=false;


			double beta0;
			double beta1=0;
			double beta0_error;
			double beta1_error=0;

			bool res_was_calculated;

			alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
			0,
			number_of_elements-starting_point,
			E+1+starting_point,
			E_errors+1+starting_point,
			cut_left_tail,
			cut_right_tail,
			y,
			E_aver,
			beta1,
			E_aver_error,
			beta1_error,
			k1_opt,
			k2_opt,
			res_was_calculated);


			if(!res_was_calculated)
			{
				inside_simulation_flag_=false;
				throw error("",1);
			};



			number_of_elements=total_number_of_ALP;


			alp_reg::robust_regression_sum_with_cut_LSM(
			0,
			number_of_elements-starting_point,
			E_T_beta+1+starting_point,
			E_T_beta_errors+1+starting_point,
			cut_left_tail,
			cut_right_tail,
			y,
			beta0,
			beta1,
			beta0_error,
			beta1_error,
			k1_opt,
			k2_opt,
			res_was_calculated);

			

			if(!res_was_calculated)
			{
				//throw error("The program cannot estimate the parameters; please repeat the calculation11\n",2);
				inside_simulation_flag_=false;
				throw error("",1);
			};


			E_T_beta_diff_aver=beta1;
			E_T_beta_diff_aver_error=beta1_error;


			
			
		};


		double exp_lambda_error=exp(-lambda_)*lambda_error_;
		double exp_lambda=(1-exp(-lambda_));

		
		double den_error=FSA_utils::error_of_the_product(E_T_beta_diff_aver,E_T_beta_diff_aver_error,exp_lambda,exp_lambda_error);
		double den=(1-exp(-lambda_))*E_T_beta_diff_aver;


		double nom_error=FSA_utils::error_of_the_product(P_beta_inf,P_beta_inf_error,E_aver,E_aver_error);
		double nom=P_beta_inf*E_aver;


		C_error_=FSA_utils::error_of_the_ratio(nom,nom_error,den,den_error);

		if(den<=0)
		{
			inside_simulation_flag_=false;
			throw error("",1);
		};

		C_=nom/den;

		memory_releace_for_calculate_C(
		P,
		P_errors,
		values_P_ratio,
		errors_P_ratio,

		E,
		E_errors,

		E_T_beta,
		E_T_beta_errors);


	}
	catch (error er)
	{
		if(er.st!="")
		{
			throw;
		}
		else
		{
			memory_releace_for_calculate_C(
			P,
			P_errors,
			values_P_ratio,
			errors_P_ratio,

			E,
			E_errors,

			E_T_beta,
			E_T_beta_errors);
		};

	};
	}
	catch (...)
	{ 
		//memory release
		memory_releace_for_calculate_C(
		P,
		P_errors,
		values_P_ratio,
		errors_P_ratio,

		E,
		E_errors,

		E_T_beta,
		E_T_beta_errors);
		throw;
	};



}


double fsa_par::get_root(
const std::vector<double> &res_tmp_,
double point_)
{
	if(res_tmp_.size()==0)
	{
		throw error("Error in alp_sim::get_root - the equation does not have any solutions\n",2);
	};

	long int i;
	long int p=0;
	double d1=fabs(point_-res_tmp_[0]);
	for(i=1;i<(long int)res_tmp_.size();i++)
	{
		double d2=fabs(point_-res_tmp_[i]);
		if(d2<d1)
		{
			p=i;
			d1=d2;
		};
	};

	return res_tmp_[p];
}

//------------------------------------------
//FSC
double fsa_par::lambda_exp(
long int &i_,
double *&exp_array_)
{

	if(exp_array_[i_]==-1)
	{
		throw error("The program is not able to calculate the parameters; rescaling penalties and scoring matrix might help\n",3);
	};

	return exp_array_[i_];
}

void fsa_par::memory_release_for_calculate_FSC(
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
double *&C_S_error)
{
	//memory release
	delete[]exp_array;exp_array=NULL;

	delete[]delta_E;delta_E=NULL;
	delete[]delta_E_error;delta_E_error=NULL;
	delete[]delta_E_E;delta_E_E=NULL;
	delete[]delta_E_E_error;delta_E_E_error=NULL;

	delete[]delta_I;delta_I=NULL;
	delete[]delta_I_error;delta_I_error=NULL;
	delete[]delta_J;delta_J=NULL;
	delete[]delta_J_error;delta_J_error=NULL;

	delete[]delta_I_J;delta_I_J=NULL;
	delete[]delta_I_J_error;delta_I_J_error=NULL;
	delete[]delta_J_J;delta_J_J=NULL;
	delete[]delta_J_J_error;delta_J_J_error=NULL;
	delete[]delta_I_I;delta_I_I=NULL;
	delete[]delta_I_I_error;delta_I_I_error=NULL;

	delete[]cov_I_J;cov_I_J=NULL;
	delete[]cov_I_J_error;cov_I_J_error=NULL;
	delete[]cov_J_J;cov_J_J=NULL;
	delete[]cov_J_J_error;cov_J_J_error=NULL;
	delete[]cov_I_I;cov_I_I=NULL;
	delete[]cov_I_I_error;cov_I_I_error=NULL;
	delete[]cov_E_E;cov_E_E=NULL;
	delete[]cov_E_E_error;cov_E_E_error=NULL;

	delete[]C_S;C_S=NULL;
	delete[]C_S_error;C_S_error=NULL;
}

void fsa_par::calculate_FSC(
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
double *test_FSC_vect_)//for tests
{
	inside_simulation_flag_=true;

	double *exp_array=NULL;

	double *delta_E=NULL;
	double *delta_E_error=NULL;

	double *delta_E_E=NULL;
	double *delta_E_E_error=NULL;


	double *delta_I=NULL;
	double *delta_I_error=NULL;

	double *delta_J=NULL;
	double *delta_J_error=NULL;

	double *delta_I_I=NULL;
	double *delta_I_I_error=NULL;

	double *delta_I_J=NULL;
	double *delta_I_J_error=NULL;

	double *delta_J_J=NULL;
	double *delta_J_J_error=NULL;

	double *cov_J_J=NULL;
	double *cov_J_J_error=NULL;

	double *cov_I_J=NULL;
	double *cov_I_J_error=NULL;

	double *cov_I_I=NULL;
	double *cov_I_I_error=NULL;

	double *cov_E_E=NULL;
	double *cov_E_E_error=NULL;

	double *C_S=NULL;
	double *C_S_error=NULL;


	try
	{
	try
	{


		if(nalp_<1)
		{
			throw error("Unexpected error\n",4);
		};

		double dbl_max_log=log(DBL_MAX);


		exp_array=new double[M_max_+1];
		FSA_utils::assert_mem(exp_array);


		long int i;
		for(i=0;i<=M_max_;i++)
		{
			double tmp=(double)i*lambda_;
			if(tmp<dbl_max_log)
			{
				exp_array[i]=exp(tmp);
			}
			else
			{
				exp_array[i]=-1;
			};
		};

		


		delta_E=new double[nalp_];
		FSA_utils::assert_mem(delta_E);
		delta_E_error=new double[nalp_];
		FSA_utils::assert_mem(delta_E_error);

		delta_E_E=new double[nalp_];
		FSA_utils::assert_mem(delta_E_E);
		delta_E_E_error=new double[nalp_];
		FSA_utils::assert_mem(delta_E_E_error);

		cov_E_E=new double[nalp_];
		FSA_utils::assert_mem(cov_E_E);
		cov_E_E_error=new double[nalp_];
		FSA_utils::assert_mem(cov_E_E_error);


		delta_I=new double[nalp_];
		FSA_utils::assert_mem(delta_I);
		delta_I_error=new double[nalp_];
		FSA_utils::assert_mem(delta_I_error);

		delta_J=new double[nalp_];
		FSA_utils::assert_mem(delta_J);
		delta_J_error=new double[nalp_];
		FSA_utils::assert_mem(delta_J_error);

		delta_I_I=new double[nalp_];
		FSA_utils::assert_mem(delta_I_I);
		delta_I_I_error=new double[nalp_];
		FSA_utils::assert_mem(delta_I_I_error);

		delta_I_J=new double[nalp_];
		FSA_utils::assert_mem(delta_I_J);
		delta_I_J_error=new double[nalp_];
		FSA_utils::assert_mem(delta_I_J_error);

		delta_J_J=new double[nalp_];
		FSA_utils::assert_mem(delta_J_J);
		delta_J_J_error=new double[nalp_];
		FSA_utils::assert_mem(delta_J_J_error);

		cov_J_J=new double[nalp_];
		FSA_utils::assert_mem(cov_J_J);
		cov_J_J_error=new double[nalp_];
		FSA_utils::assert_mem(cov_J_J_error);

		cov_I_J=new double[nalp_];
		FSA_utils::assert_mem(cov_I_J);
		cov_I_J_error=new double[nalp_];
		FSA_utils::assert_mem(cov_I_J_error);

		cov_I_I=new double[nalp_];
		FSA_utils::assert_mem(cov_I_I);
		cov_I_I_error=new double[nalp_];
		FSA_utils::assert_mem(cov_I_I_error);

		C_S=new double[nalp_];
		FSA_utils::assert_mem(C_S);
		C_S_error=new double[nalp_];
		FSA_utils::assert_mem(C_S_error);


		long int j;
		for(j=0;j<nalp_;j++)
		{
			delta_E[j]=0.0;
			delta_E_error[j]=0.0;

			delta_E_E[j]=0.0;
			delta_E_E_error[j]=0.0;

			delta_I[j]=0.0;
			delta_I_error[j]=0.0;
			delta_J[j]=0.0;
			delta_J_error[j]=0.0;

			delta_I_I[j]=0.0;

			delta_I_I_error[j]=0.0;
			delta_I_J[j]=0.0;
			delta_I_J_error[j]=0.0;
			delta_J_J[j]=0.0;
			delta_J_J_error[j]=0.0;

			C_S[j]=0.0;
			C_S_error[j]=0.0;

		};

		double C_S_constant=1.0;
		//calculate the constant
		bool calculate_C_S_constant_flag=true;
		if(calculate_C_S_constant_flag)
		{
			long int i;
			for(i=ind1_;i<=ind2_;i++)
			{
				long int j;
				for(j=1;j<=nalp_;j++)
				{

					long int &E_j=ALP_edge_max_[i].d_elem[j];
					double &weight_j=ALP_weight_[i].d_elem[j];


					double tmp=lambda_exp(E_j,exp_array)*weight_j;

					C_S[j-1]+=tmp;
					C_S_error[j-1]+=tmp*tmp;
				};
			};

			double ind_diff=(double)(ind2_-ind1_+1);
			for(j=0;j<nalp_;j++)
			{
				C_S[j]/=ind_diff;
				C_S_error[j]/=ind_diff;
				C_S_error[j]-=C_S[j]*C_S[j];
				C_S_error[j]/=ind_diff;
				C_S_error[j]=alp_reg::sqrt_for_errors(C_S_error[j]);

			};

			//regression

			double beta1=0;
			double beta1_error=0;

			long int number_of_elements=nalp_;

			bool cut_left_tail=true;
			bool cut_right_tail=false;

			double y=2;

			long int k1_opt;
			long int k2_opt;


			double C_S_constant_error;


			bool res_was_calculated;
			
			alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
			0,
			number_of_elements,
			C_S,
			C_S_error,
			cut_left_tail,
			cut_right_tail,
			y,
			C_S_constant,
			beta1,
			C_S_constant_error,
			beta1_error,
			k1_opt,
			k2_opt,
			res_was_calculated);

			if(!res_was_calculated||C_S_constant<=0)
			{
				inside_simulation_flag_=false;
				throw error("",1);
			};

		};

		double one_div_C_S_constant=1.0/C_S_constant;


		
		for(i=ind1_;i<=ind2_;i++)
		{

			long int j;
			for(j=1;j<=nalp_;j++)
			{
				long int j_1=j-1;

				long int &E_j_1=ALP_edge_max_[i].d_elem[j_1];
				long int &E_j=ALP_edge_max_[i].d_elem[j];
				double &weight_j=ALP_weight_[i].d_elem[j];

				long int &I_j_1=distance_along_direction_1_[i].d_elem[j_1];
				long int &I_j=distance_along_direction_1_[i].d_elem[j];

				long int &J_j_1=distance_along_direction_2_[i].d_elem[j_1];
				long int &J_j=distance_along_direction_2_[i].d_elem[j];

				double exp_tmp=lambda_exp(E_j,exp_array)*one_div_C_S_constant;

				double delta_I_tmp=(I_j-I_j_1)*exp_tmp*weight_j;
				double delta_J_tmp=(J_j-J_j_1)*exp_tmp*weight_j;
				double delta_E_tmp=(E_j-E_j_1)*exp_tmp*weight_j;
				double delta_E_E_tmp=(E_j-E_j_1)*(E_j-E_j_1)*exp_tmp*weight_j;

				

				
				double delta_I_I_tmp=delta_I_tmp*(I_j-I_j_1);
				double delta_J_J_tmp=delta_J_tmp*(J_j-J_j_1);
				double delta_I_J_tmp=delta_I_tmp*(J_j-J_j_1);

	


				delta_E[j_1]+=delta_E_tmp;
				delta_E_error[j_1]+=delta_E_tmp*delta_E_tmp;

				delta_E_E[j_1]+=delta_E_E_tmp;
				delta_E_E_error[j_1]+=delta_E_E_tmp*delta_E_E_tmp;

				delta_I[j_1]+=delta_I_tmp;
				delta_I_error[j_1]+=delta_I_tmp*delta_I_tmp;
				delta_J[j_1]+=delta_J_tmp;
				delta_J_error[j_1]+=delta_J_tmp*delta_J_tmp;

				delta_I_I[j_1]+=delta_I_I_tmp;
				delta_I_I_error[j_1]+=delta_I_I_tmp*delta_I_I_tmp;

				delta_I_J[j_1]+=delta_I_J_tmp;
				delta_I_J_error[j_1]+=delta_I_J_tmp*delta_I_J_tmp;

				delta_J_J[j_1]+=delta_J_J_tmp;
				delta_J_J_error[j_1]+=delta_J_J_tmp*delta_J_J_tmp;
				
			};
		};


		double ind_diff=(double)(ind2_-ind1_+1);
		for(j=0;j<nalp_;j++)
		{
			delta_E[j]/=ind_diff;
			delta_E_error[j]/=ind_diff;
			delta_E_error[j]-=delta_E[j]*delta_E[j];
			delta_E_error[j]/=ind_diff;
			delta_E_error[j]=alp_reg::sqrt_for_errors(delta_E_error[j]);

			

			delta_E_E[j]/=ind_diff;
			delta_E_E_error[j]/=ind_diff;
			delta_E_E_error[j]-=delta_E_E[j]*delta_E_E[j];
			delta_E_E_error[j]/=ind_diff;


			delta_I[j]/=ind_diff;
			delta_I_error[j]/=ind_diff;
			delta_I_error[j]-=delta_I[j]*delta_I[j];
			delta_I_error[j]/=ind_diff;
			delta_I_error[j]=alp_reg::sqrt_for_errors(delta_I_error[j]);

			delta_J[j]/=ind_diff;
			delta_J_error[j]/=ind_diff;
			delta_J_error[j]-=delta_J[j]*delta_J[j];
			delta_J_error[j]/=ind_diff;
			delta_J_error[j]=alp_reg::sqrt_for_errors(delta_J_error[j]);

			delta_I_J[j]/=ind_diff;
			delta_I_J_error[j]/=ind_diff;
			delta_I_J_error[j]-=delta_I_J[j]*delta_I_J[j];
			delta_I_J_error[j]/=ind_diff;


			delta_I_I[j]/=ind_diff;
			delta_I_I_error[j]/=ind_diff;
			delta_I_I_error[j]-=delta_I_I[j]*delta_I_I[j];
			delta_I_I_error[j]/=ind_diff;


			delta_J_J[j]/=ind_diff;
			delta_J_J_error[j]/=ind_diff;
			delta_J_J_error[j]-=delta_J_J[j]*delta_J_J[j];
			delta_J_J_error[j]/=ind_diff;


			cov_I_J[j]=delta_I_J[j]-delta_I[j]*delta_J[j];
			cov_I_I[j]=delta_I_I[j]-delta_I[j]*delta_I[j];
			cov_J_J[j]=delta_J_J[j]-delta_J[j]*delta_J[j];

			cov_E_E[j]=delta_E_E[j]-delta_E[j]*delta_E[j];

			cov_I_J_error[j]=FSA_utils::error_of_the_product(delta_I[j],delta_I_error[j],delta_J[j],delta_J_error[j]);
			cov_I_J_error[j]=sqrt(delta_I_J_error[j]+cov_I_J_error[j]*cov_I_J_error[j]);

			cov_I_I_error[j]=FSA_utils::error_of_the_product(delta_I[j],delta_I_error[j],delta_I[j],delta_I_error[j]);
			cov_I_I_error[j]=sqrt(delta_I_I_error[j]+cov_I_I_error[j]*cov_I_I_error[j]);

			cov_J_J_error[j]=FSA_utils::error_of_the_product(delta_J[j],delta_J_error[j],delta_J[j],delta_J_error[j]);
			cov_J_J_error[j]=sqrt(delta_J_J_error[j]+cov_J_J_error[j]*cov_J_J_error[j]);

			cov_E_E_error[j]=FSA_utils::error_of_the_product(delta_E[j],delta_E_error[j],delta_E[j],delta_E_error[j]);
			cov_E_E_error[j]=sqrt(delta_E_E_error[j]+cov_E_E_error[j]*cov_E_E_error[j]);


		};


		//regression

		double beta1=0;
		double beta1_error=0;

		long int number_of_elements=nalp_;

		bool cut_left_tail=true;
		bool cut_right_tail=false;

		double y=2;

		long int k1_opt;
		long int k2_opt;


		double delta_I_aver;
		double delta_I_aver_error;



		bool res_was_calculated;
		
		alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
		0,
		number_of_elements,
		delta_I,
		delta_I_error,
		cut_left_tail,
		cut_right_tail,
		y,
		delta_I_aver,
		beta1,
		delta_I_aver_error,
		beta1_error,
		k1_opt,
		k2_opt,
		res_was_calculated);

		if(!res_was_calculated)
		{
			inside_simulation_flag_=false;
			throw error("",1);
		};



		double delta_J_aver;
		double delta_J_aver_error;

		
		
		alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
		0,
		number_of_elements,
		delta_J,
		delta_J_error,
		cut_left_tail,
		cut_right_tail,
		y,
		delta_J_aver,
		beta1,
		delta_J_aver_error,
		beta1_error,
		k1_opt,
		k2_opt,
		res_was_calculated);

		if(!res_was_calculated)
		{
			inside_simulation_flag_=false;
			throw error("",1);
		};


		double delta_E_aver;
		double delta_E_aver_error;
		
		alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
		0,
		number_of_elements,
		delta_E,
		delta_E_error,
		cut_left_tail,
		cut_right_tail,
		y,
		delta_E_aver,
		beta1,
		delta_E_aver_error,
		beta1_error,
		k1_opt,
		k2_opt,
		res_was_calculated);

		if(!res_was_calculated)
		{
			inside_simulation_flag_=false;
			throw error("",1);
		};

		double cov_I_I_aver;
		double cov_I_I_aver_error;

		double cov_I_J_aver;
		double cov_I_J_aver_error;

		double cov_J_J_aver;
		double cov_J_J_aver_error;

		double cov_E_E_aver;
		double cov_E_E_aver_error;


		alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
		0,
		number_of_elements,
		cov_I_J,
		cov_I_J_error,
		cut_left_tail,
		cut_right_tail,
		y,
		cov_I_J_aver,
		beta1,
		cov_I_J_aver_error,
		beta1_error,
		k1_opt,
		k2_opt,
		res_was_calculated);

		if(!res_was_calculated)
		{
			inside_simulation_flag_=false;
			throw error("",1);
		};

		alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
		0,
		number_of_elements,
		cov_I_I,
		cov_I_I_error,
		cut_left_tail,
		cut_right_tail,
		y,
		cov_I_I_aver,
		beta1,
		cov_I_I_aver_error,
		beta1_error,
		k1_opt,
		k2_opt,
		res_was_calculated);


		if(!res_was_calculated)
		{
			inside_simulation_flag_=false;
			throw error("",1);
		};

		alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
		0,
		number_of_elements,
		cov_J_J,
		cov_J_J_error,
		cut_left_tail,
		cut_right_tail,
		y,
		cov_J_J_aver,
		beta1,
		cov_J_J_aver_error,
		beta1_error,
		k1_opt,
		k2_opt,
		res_was_calculated);

		if(!res_was_calculated)
		{
			inside_simulation_flag_=false;
			throw error("",1);
		};

		alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
		0,
		number_of_elements,
		cov_E_E,
		cov_E_E_error,
		cut_left_tail,
		cut_right_tail,
		y,
		cov_E_E_aver,
		beta1,
		cov_E_E_aver_error,
		beta1_error,
		k1_opt,
		k2_opt,
		res_was_calculated);

		if(!res_was_calculated)
		{
			inside_simulation_flag_=false;
			throw error("",1);
		};



		if(delta_E_aver<=0)
		{
			inside_simulation_flag_=false;
			throw error("",1);
		};



		a_I_=delta_I_aver/delta_E_aver;
		a_I_error_=FSA_utils::error_of_the_ratio(delta_I_aver,delta_I_aver_error,delta_E_aver,delta_E_aver_error);
		a_J_=delta_J_aver/delta_E_aver;
		a_J_error_=FSA_utils::error_of_the_ratio(delta_J_aver,delta_J_aver_error,delta_E_aver,delta_E_aver_error);
		
		sigma_calculation(
		delta_I_aver,
		delta_I_aver_error,
		delta_J_aver,
		delta_J_aver_error,
		delta_E_aver,
		delta_E_aver_error,
		cov_E_E_aver,
		cov_E_E_aver_error,
		cov_I_J_aver,
		cov_I_J_aver_error,
		sigma_,
		sigma_error_);


		sigma_calculation(
		delta_I_aver,
		delta_I_aver_error,
		delta_I_aver,
		delta_I_aver_error,
		delta_E_aver,
		delta_E_aver_error,
		cov_E_E_aver,
		cov_E_E_aver_error,
		cov_I_I_aver,
		cov_I_I_aver_error,
		alpha_I_,
		alpha_I_error_);

		//for test
		if(test_FSC_vect_)
		{
			test_FSC_vect_[0]=alpha_I_;
			long int j;
			for(j=0;j<number_of_elements;j++)
			{
				test_FSC_vect_[j+1]=(delta_I[j]*delta_I[j]*cov_E_E[j]+delta_E[j]*delta_E[j]*cov_I_I[j])/(delta_E[j]*delta_E[j]*delta_E[j]);
			};
		};


		if(test_FSC_vect_)
		{
			test_FSC_vect_[0]=lambda_;
			long int i,j;
			for(j=1;j<=nalp_;j++)
			{
				test_FSC_vect_[j]=0;
			};

			for(i=ind1_;i<=ind2_;i++)
			{
				

				long int j;
				for(j=1;j<=nalp_;j++)
				{

					long int &E_j=ALP_edge_max_[i].d_elem[j];
					double &weight_j=ALP_weight_[i].d_elem[j];


					double exp_tmp=lambda_exp(E_j,exp_array)*one_div_C_S_constant;

					double tmp=exp_tmp*weight_j;
					test_FSC_vect_[j]+=tmp;
				};
			};

			for(j=1;j<=nalp_;j++)
			{
				test_FSC_vect_[j]/=(double)(ind2_-ind1_)+1;
			};

		};


		sigma_calculation(
		delta_J_aver,
		delta_J_aver_error,
		delta_J_aver,
		delta_J_aver_error,
		delta_E_aver,
		delta_E_aver_error,
		cov_E_E_aver,
		cov_E_E_aver_error,
		cov_J_J_aver,
		cov_J_J_aver_error,
		alpha_J_,
		alpha_J_error_);

		memory_release_for_calculate_FSC(
		exp_array,
		delta_E,
		delta_E_error,
		delta_E_E,
		delta_E_E_error,
		delta_I,
		delta_I_error,
		delta_J,
		delta_J_error,
		delta_I_I,
		delta_I_I_error,
		delta_I_J,
		delta_I_J_error,
		delta_J_J,
		delta_J_J_error,
		cov_J_J,
		cov_J_J_error,
		cov_I_J,
		cov_I_J_error,
		cov_I_I,
		cov_I_I_error,
		cov_E_E,
		cov_E_E_error,
		C_S,
		C_S_error);


	}
	catch (error er)
	{
		if(er.st!="")
		{
			throw;
		}
		else
		{

			//memory release
			memory_release_for_calculate_FSC(
			exp_array,
			delta_E,
			delta_E_error,
			delta_E_E,
			delta_E_E_error,
			delta_I,
			delta_I_error,
			delta_J,
			delta_J_error,
			delta_I_I,
			delta_I_I_error,
			delta_I_J,
			delta_I_J_error,
			delta_J_J,
			delta_J_J_error,
			cov_J_J,
			cov_J_J_error,
			cov_I_J,
			cov_I_J_error,
			cov_I_I,
			cov_I_I_error,
			cov_E_E,
			cov_E_E_error,
			C_S,
			C_S_error);
		};

	};
	}
	catch (...)
	{ 
		//memory release
		memory_release_for_calculate_FSC(
		exp_array,
		delta_E,
		delta_E_error,
		delta_E_E,
		delta_E_E_error,
		delta_I,
		delta_I_error,
		delta_J,
		delta_J_error,
		delta_I_I,
		delta_I_I_error,
		delta_I_J,
		delta_I_J_error,
		delta_J_J,
		delta_J_J_error,
		cov_J_J,
		cov_J_J_error,
		cov_I_J,
		cov_I_J_error,
		cov_I_I,
		cov_I_I_error,
		cov_E_E,
		cov_E_E_error,
		C_S,
		C_S_error);
		throw;
	};
}

void fsa_par::sigma_calculation(
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
double &sigma_error_)
{
	double nom1_1=delta_I_aver_*delta_J_aver_;
	double nom2_2=delta_E_aver_*delta_E_aver_;

	double den=nom2_2*delta_E_aver_;

	double nom1=nom1_1*cov_E_E_aver_;
	double nom2=nom2_2*cov_I_J_aver_;

	sigma_=(nom1+nom2)/den;

	
	double nom1_sigma_error=FSA_utils::error_of_the_product(delta_I_aver_,delta_I_aver_error_,delta_J_aver_,delta_J_aver_error_);
	nom1_sigma_error=FSA_utils::error_of_the_product(nom1_1,nom1_sigma_error,cov_E_E_aver_,cov_E_E_aver_error_);

	
	double nom2_sigma_error_2=FSA_utils::error_of_the_product(delta_E_aver_,delta_E_aver_error_,delta_E_aver_,delta_E_aver_error_);
	double nom2_sigma_error=FSA_utils::error_of_the_product(nom2_2,nom2_sigma_error_2,cov_I_J_aver_,cov_I_J_aver_error_);

	
	double den_sigma_error=FSA_utils::error_of_the_product(nom2_2,nom2_sigma_error_2,delta_E_aver_,delta_E_aver_error_);

	double nom_sigma_error=FSA_utils::error_of_the_sum(nom1_sigma_error,nom2_sigma_error);

	sigma_error_=FSA_utils::error_of_the_ratio(nom1+nom2,nom_sigma_error,den,den_sigma_error);

}

void fsa_par::Output_Pvalues(
char type_,//'E' or 'P'
long int score1_,
long int score2_,
vector<double> &pv_,
vector<double> &pv_err_,
string pvalout_file_name_)//P-values file name
{
	if(!(type_=='E'||type_=='P'))
	{
		throw error("Error - the type_ parameter can take only 'E' or 'P' values in void fsa_par::Output_Pvalues\n",1);
	};

	ofstream f(pvalout_file_name_.data());
	if(!f)
	{
		throw error("Error - the file "+pvalout_file_name_+" cannot be created\n",4);
	};

	if(type_=='E')
	{
		f<<"Score\tE-value\tE-value error\n";
	};

	if(type_=='P')
	{
		f<<"Score\tP-value\tP-value error\n";
	};

	long int i;
	for(i=score1_;i<=score2_;i++)
	{
		f<<i<<"\t"<<pv_[i-score1_]<<"\t"<<pv_err_[i-score1_]<<endl;
	};


	f.close();

}

void fsa_par::Output_Params(
Sls::FALP_set_of_parameters &gumbel_params_,
string gumbelparout_file_name_)
{
	ofstream f;
	ifstream fin;


	bool append_flag=false;

	try
	{

		if(append_flag)
		{
			fin.open(gumbelparout_file_name_.data());
			if(!fin)
			{
				append_flag=false;
			}
			else
			{
				fin.close();
			};
		};
		
		if(!append_flag)
		{
			f.open(gumbelparout_file_name_.data());
		}
		else
		{
			f.open(gumbelparout_file_name_.data(),ios::app);
		};

		if(!f)
		{
			throw error("Error - file "+gumbelparout_file_name_+" cannot be created\n",4);
		};

		f<<gumbel_params_;


		f.close();
	}
	catch (...)
	{ 
		if(f.is_open())
		{
			f.close();
		};

		if(append_flag)
		{
			if(fin.is_open())
			{
				fin.close();
			};
		};

		throw;
	};

}

void fsa_par::Read_Params(
Sls::FALP_set_of_parameters &gumbel_params_,
string gumbelparin_file_name_)
{
	ifstream f;

	try
	{

		gumbel_params_.d_params_flag=false;

		f.open(gumbelparin_file_name_.data());
		if(!f)
		{
			throw error("Error - file "+gumbelparin_file_name_+" is not found\n",4);
		};

		f>>gumbel_params_;

		f.close();
		gumbel_params_.d_params_flag=true;

	}
	catch (...)
	{ 
		gumbel_params_.d_params_flag=false;
		if(f.is_open())
		{
			f.close();
		};
		throw;
	};
}

