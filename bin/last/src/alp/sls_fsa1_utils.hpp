#ifndef GUMBEL_FSA1_UTILS
#define GUMBEL_FSA1_UTILS

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

File name: sls_fsa1_utils.hpp

Author: Sergey Sheetlin, Martin Frith

Contents: Frameshift alignment algorithms utilities

******************************************************************************/
#include "sls_basic.hpp"

#include <vector>
#include <complex>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <float.h> 
#include <algorithm>
#include <sstream>
#include <ctime>
#include <limits>
#include <climits>
#include <errno.h>

#include "njn_random.hpp"
#include "njn_uniform.hpp"

namespace Sls 
{

	class FSA_utils: public sls_basic//contains general functions
	{
		public:

		static void read_RR(
		std::string RR_file_name_,
		double *&RR_,
		double *&RR_sum_,
		long int *&RR_sum_elements_,
		long int &number_of_AA_RR_,
		long int number_of_AA_RR_default_=0);

		static void check_RR_sum(
		double sum_tmp_,
		long int number_of_AA_RR_,
		std::string RR_file_name_);

		static void calculate_RR_sum(
		double *RR_,
		long int number_of_AA_RR_,
		double *&RR_sum_,
		long int *&RR_sum_elements_);

		static void read_RR(
		std::string RR_file_name_,
		double *&RR_,
		long int &number_of_AA_RR_,
		long int number_of_AA_RR_default_=0);

		static void read_smatr(
		std::string smatr_file_name_,
		long int **&smatr_,
		long int &number_of_AA_smatr_,
		long int &smatr_min_);

		static void smatr_min(
		long int **smatr_,
		long int number_of_AA_smatr_,
		long int &smatr_min_);

		static void remove_zero_probabilities(
		double *&RR1_,
		double *&RR1_sum_,
		long int *&RR1_sum_elements_,
		long int &alphabet_letters_number1_,
		double *&RR2_,
		double *&RR2_sum_,
		long int *&RR2_sum_elements_,
		long int &alphabet_letters_number2_,
		long int **&smatr_,
		long int &number_of_AA_smatr_,
		long int &smatr_min_,

		long int &number_of_letters1_,//number of letters for the sequence 1
		long int &number_of_letters2_,//number of letters for the sequence 2

		char *&alphabet1_,//alphabet letters for the sequence #1
		char *&alphabet2_,//alphabet letters for the sequence #2

		long int *&alphabet1_to_long_,//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
		long int *&alphabet2_to_long_,//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

		long int &codon_length_,//codon length 
		long int *&codon_AA_);//<codon code,AA number>





		static std::string long_to_string(//convert interer ot std::string
		long int number_);

		static char digit_to_string(//convert interer ot std::string
		long int digit_);
		

		static bool the_value_is_double(
		std::string str_,
		double &val_);

		static bool the_value_is_long(
		std::string str_,
		long int &val_);

		static double sqrt_plus(
		double x_);

		static double error_of_the_sum(//v1_+v2_
		double v1_error_,
		double v2_error_);

		static double error_of_the_product(//v1_*v2_
		double v1_,
		double v1_error_,
		double v2_,
		double v2_error_);

		static double error_of_the_lg(//lg(v1_)
		double v1_,
		double v1_error_);

		static double error_of_the_sqrt(//sqrt(v1_)
		double v1_,
		double v1_error_);

		static double error_of_the_ratio(//v1_/v2_
		double v1_,
		double v1_error_,
		double v2_,
		double v2_error_);

		static double error_of_the_sum_with_coeff(//c1_*v1_+c2_*v2_
		double c1_,
		double v1_error_,
		double c2_,
		double v2_error_);



		static void read_alphabet(
		std::string alphabet_file_name_,
		long int &number_of_AA_alphabet_,
		char* &alphabet_);

		static void calculate_composition_frequencies(
		double number_of_AA_alphabet_,
		char* alphabet_,
		std::string fasta_file_name_,
		double *comp_frequencies_);

		static long int convert_codon_into_AA(
		long int codon_length_,//codon length 
		long int *codon_AA_,//<codon code,AA number>
		long int number_of_letters1_,//number of letters for the sequence 1
		long int *codon_);

		static long int convert_codon_into_code(
		long int codon_length_,//codon length 
		long int number_of_letters1_,//number of letters for the sequence 1
		long int *codon_);//input codon


		static void convert_code_into_codon(
		long int code_,//the input code
		long int codon_length_,//codon length 
		long int number_of_letters1_,//number of letters for the sequence 1
		long int *codon_);//must be allocated


		static void read_codon_AA_file(
		std::string file_name_,
		long int &number_of_letters1_,//number of letters for the sequence 1
		long int &number_of_letters2_,//number of letters for the sequence 2

		char *&alphabet1_,//alphabet letters for the sequence #1
		char *&alphabet2_,//alphabet letters for the sequence #2

		long int *&alphabet1_to_long_,//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
		long int *&alphabet2_to_long_,//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

		long int &codon_length_,//codon length 
		long int *&codon_AA_,//<codon code,AA number>
		bool reverse_codons_flag_=false);//if true, then the codons are reversed

		static void reverse_codons(
		long int *codon_AA_,//<codon code,AA number>; original codons
		long int alphabet_letters_number1_,//number of letters for the sequence #1
		long int codon_length_,//codon length 
		long int *&codon_AA_reversed_);//<codon code,AA number>; reversed codons


		static void reverse_sequence(//reverse the letters of the sequence
			long int *seq_,
			long int seq_length_);

		static void extract_AA_frequencies_for_DNA_sequence(
			const long int *codon_AA_,//<codon code,AA number>
			long int &codon_length_,//codon length 
			long int number_of_letters1_,//number of letters for the sequence 1
			long int number_of_letters2_,//number of letters for the sequence 2
			const double *RR1_,//nucleotide probabilities
			double *&RR1_AA_);//the resulted frequencies

		static void read_sequences_for_alingment(

			std::string input_file_name_,

			long int &number_of_letters1_,//number of letters for the sequence 1
			long int &number_of_letters2_,//number of letters for the sequence 2

			char *&alphabet1_,//alphabet letters for the sequence #1
			char *&alphabet2_,//alphabet letters for the sequence #2

			long int& number_of_sequences_,

			std::string *&headers_,
			long int *&lengths1_,//lengths of the sequences #1
			long int *&lengths2_,//lengths of the sequences #2
			long int **&sequences1_,//the first index numerates sequences; the second - sequence letters
			long int **&sequences2_);

		static void read_string(
			std::ifstream &f_,
			std::string &st_,
			bool &end_of_file_flag_);



		inline static double ran2()//generates the next random value
		{
			return Njn::Uniform::variate <double> (0,1);
			//double rand_C=(double)((double)rand()/(double)RAND_MAX);
			//return rand_C;	
		}

		inline static void srand2(long int seed_)//initializes the seed
		{
			Njn::Random::seed(seed_);
			//srand(seed_);
		}

		static void srand2_old(long int seed_);//initializes the seed

		static double ran2_old();//generates the next random value

		static void process_random_factor(
			long int &random_factor_,
			bool *rand_flag_=NULL);


		static long int power_long(//returns a_^b_
			long int a_,
			long int b_);

		static void convert_distr_into_sum(//convert distr_[0], distr_[1], distr_[2],... into distr_[0], distr_[0]+distr_[1], distr_[0]+distr_[1]+distr_[2],...
			long int dim_,
			double *distr_);

		static bool Gauss(
			std::vector<std::vector<double> > A_,//matrix n*(n+1)
			std::vector<double> &x_,//solution
			double inside_eps_,
			std::vector<std::vector<double> > *inv_A_);

		static void multiply_matrices(
			const std::vector<std::vector<double> > &A_,
			const std::vector<std::vector<double> > &B_,
			std::vector<std::vector<double> > &res_);

		static void multiply_matrix_and_vector(
			const std::vector<std::vector<double> > &A_,
			const std::vector<double> &y_,
			std::vector<double> &res_);

		static void transpose_matrix(
			const std::vector<std::vector<double> > &A_,
			std::vector<std::vector<double> > &res_);

		static void print_matrix(
			const std::vector<std::vector<double> > A_);

		static double standard_deviation(//standard deviation for average of elements of vect_
			long int dim_,
			double *vect_);

		static double average(//average of elements of vect_
			long int dim_,
			double *vect_);






//----------------------------------

		template<typename T>
		static void get_memory_for_matrix(
		long int dim1_,
		long int dim2_,
		T ** &matr_)
		{
			matr_=NULL;


			try
			{

				long int i;
				matr_=new T *[dim1_];
				Sls::FSA_utils::assert_mem(matr_);

				for(i=0;i<dim1_;i++)
				{
					matr_[i]=NULL;
				};

				for(i=0;i<dim1_;i++)
				{
					matr_[i]=new T [dim2_];
					Sls::FSA_utils::assert_mem(matr_[i]);
				};

			}
			catch (...)
			{ 
				if(matr_)
				{
					long int i;
					for(i=0;i<dim1_;i++)
					{
						delete[]matr_[i];matr_[i]=NULL;
					};
					delete[]matr_;matr_=NULL;
				};
				throw;
			};

		}

		template<typename T>
		static void delete_memory_for_matrix(
		long int dim1_,
		T ** &matr_)
		{
			long int i;
			if(matr_)
			{
				for(i=0;i<dim1_;i++)
				{
					delete []matr_[i];matr_[i]=NULL;
				};
				delete []matr_;matr_=NULL;
			};

		}

		static long int random_long(
		double value_,
		long int dim_)
		{
			if(value_<0||value_>1.0||dim_<=0)
			{
				throw error("Unexpected error\n",4);
			};

			if(dim_==1)
			{
				return 0;
			};

			long int tmp=(long int)floor(value_*(double)dim_);
			tmp=Tmin(tmp,dim_-1);
			return tmp;
		}

		template<typename T>
		static T random_long(
		double value_,
		long int dim_,
		double *sum_distr_,
		T* elements_=NULL)//sum_distr_[dim_-1] must be equal to 1
		{
			if(value_<0||value_>1)	
			{
				throw error("Unexpected error in q_elem importance_sampling::get_random_pair\n",4);
			};

			long int v1=0;
			long int v2=dim_;

			while(v2-v1>1)
			{
				long int v3=(long int)(Sls::FSA_utils::round(double(v2+v1)/2.0));
				if(sum_distr_[v3-1]==value_)
				{
					v1=v3-1;
					v2=v3;
					break;
				};

				if(sum_distr_[v3-1]>value_)
				{
					v2=v3;
				}
				else
				{
					v1=v3;
				};
			};

			if(elements_)
			{
				long int v2_1=v2-1;


				long int v2_minus=-1;

				long int j;
				for(j=v2_1;j>=1;j--)
				{
					if(sum_distr_[j]!=sum_distr_[j-1])
					{
						v2_minus=j;
						break;
					};
				};

				if(v2_minus<0)
				{
					if(sum_distr_[0]>0)
					{
						v2_minus=0;
					};
				};

				if(v2_minus>=0)
				{
					return elements_[v2_minus];
				};

				long int v2_plus=-1;
				for(j=v2;j<dim_;j++)
				{
					if(sum_distr_[j]!=sum_distr_[j-1])
					{
						v2_plus=j;
						break;
					};
				};

				if(v2_minus<0)
				{
					throw error("Unexpected error in random_long\n",1);
				}
				else
				{
					return elements_[v2_plus];
				};

			}
			else
			{
				return v2-1;
			};

		}

		template<class T>
		static void print_matrix(
		std::string file_name_,
		long int dim1_,
		long int dim2_,
		T** A_)
		{
			std::ofstream f(file_name_.data());
			if(!f)
			{
				throw error("Error - the file "+file_name_+" is not found\n",3);
			};

			long int i,j;
			for(i=0;i<dim1_;i++)
			{
				for(j=0;j<dim2_;j++)
				{
					if(j<dim2_-1)
					{
						f<<A_[i][j]<<"\t";
					}
					else
					{
						f<<A_[i][j]<<"\n";
					};
				};
			};

			f.close();
		}





		public:

		long int d_tmpl;


	};





	template<typename T> class array_positive{
	public:
		array_positive(void *fsa_=NULL)// constructor
		{ 
			d_elem=NULL;
			d_fsa=fsa_; 
			if(!fsa_)
			{
				//throw error("Unexpected error\n",4);
			};
			d_dim=-1;
			d_step=10;
		}

		void increment_array(long int ind_);

		inline void set_elem(
			long int ind_,
			T elem_)
		{
			if(ind_>d_dim)
			{
				increment_array(ind_);
			};

			d_elem[ind_]=elem_;
		}

		inline T get_elem(
			long int ind_)
		{
			if(ind_>d_dim)
			{
				increment_array(ind_);
			};

			return d_elem[ind_];
		}


		inline void increase_elem_by_1(
			long int ind_)
		{
			if(ind_>d_dim)
			{
				increment_array(ind_);
			};

			d_elem[ind_]++;
		}

		inline void increase_elem_by_x(
			long int ind_,
			T x_)
		{
			if(ind_>d_dim)
			{
				increment_array(ind_);
			};

			d_elem[ind_]+=x_;
		}

		inline void divide_elem_by_x(
			long int ind_,
			T x_)
		{
			if(ind_>d_dim)
			{
				increment_array(ind_);
			};

			d_elem[ind_]/=x_;
		}

		~array_positive()
		{
			delete[]d_elem;d_elem=NULL;
		}


	public:
			
		long int d_step;
		long int d_dim;//dimension of the array is d_dim+1
		T * d_elem;
		void *d_fsa;//initial data
	};





		



//-------------------------------------------------------------------------------------------


	template<typename T> class array_v{
	public:
		array_v(void *fsa_=NULL)// constructor
		{ 
			d_elem=NULL;
			d_fsa=fsa_; 
			d_dim=-1;
			d_ind0=0;
			d_step=10;
			d_dim_plus_d_ind0=d_dim+d_ind0;
		}

		void clear()
		{
			d_dim=-1;
			d_ind0=0;
			d_dim_plus_d_ind0=d_dim+d_ind0;
			delete[]d_elem;d_elem=NULL;
		}

		void increment_array_on_the_right(long int ind_);

		void increment_array_on_the_left(long int ind_);

		inline void set_elems(const array_v<T> *a_)
		{
			long int a0=a_->d_ind0;
			long int a1=a_->d_dim_plus_d_ind0;

			if(a0>a1)return;

			while(a1>d_dim_plus_d_ind0)
			{
				d_dim_plus_d_ind0+=d_step;
			};

			while(a0<d_ind0)
			{
				d_ind0-=d_step;
			};

			d_dim=d_dim_plus_d_ind0-d_ind0;
			d_elem=new T[d_dim+1];
			sls_basic::assert_mem(d_elem);

			long int i;
			for(i=a0;i<=a1;i++)
			{
			  d_elem[i-d_ind0]=a_->d_elem[i-a0];
			}
		}

		inline void set_elem(
			long int ind_,
			T elem_)
		{
			if(ind_>d_dim_plus_d_ind0)
			{
				increment_array_on_the_right(ind_);
			};

			if(ind_<d_ind0)
			{
				increment_array_on_the_left(ind_);
			};

			d_elem[ind_-d_ind0]=elem_;
		}


		T get_elem(
			long int ind_)
		{
			if(ind_>d_dim_plus_d_ind0||ind_<d_ind0)
			{
				throw error("Error - the index ind_ in array_v::get_elem is out of range\n",1);
			};

			return d_elem[ind_-d_ind0];
		}


		inline void increase_elem_by_1(
			long int ind_)
		{
			if(ind_>d_dim_plus_d_ind0)
			{
				increment_array_on_the_right(ind_);
			};

			if(ind_<d_ind0)
			{
				increment_array_on_the_left(ind_);
			};

			d_elem[ind_-d_ind0]++;
		}


		void increase_elem_by_x(
			long int ind_,
			double x_)
		{
			if(ind_>d_dim_plus_d_ind0)
			{
				increment_array_on_the_right(ind_);
			};

			if(ind_<d_ind0)
			{
				increment_array_on_the_left(ind_);
			};

			d_elem[ind_-d_ind0]+=x_;
		}


		void divide_elem_by_x(
			long int ind_,
			double x_)
		{
			if(ind_>d_dim_plus_d_ind0)
			{
				increment_array_on_the_right(ind_);
			};

			if(ind_<d_ind0)
			{
				increment_array_on_the_left(ind_);
			};

			d_elem[ind_-d_ind0]/=x_;
		}




		~array_v()
		{
			delete[]d_elem;d_elem=NULL;
		}

		static void copy_x2_into_x1(
			array_v<T> &x1_,
			array_v<T> &x2_)
		{
			x1_=x2_;
			if(x2_.d_dim<0)
			{
				x1_.d_elem=NULL;
				return;
			};
			x1_.d_elem=new T[x2_.d_dim+1];
			FSA_utils::assert_mem(x1_.d_elem);
			long int i;
			for(i=0;i<=x2_.d_dim;i++)
			{
				x1_.d_elem[i]=x2_.d_elem[i];
			};

		}
			



	public:
			
		long int d_step;
		long int d_dim;//dimension of the array is d_dim+1
		long int d_ind0;//the leftmost index of the array
		long int d_dim_plus_d_ind0;
		T * d_elem;
		void *d_fsa;//initial data
	};


//-------------------------------------------------------------------------------------------

	template<typename T> class array_v2{


	public:
		array_v2(// constructor
			long int subarray_length_,
			bool assign_null_elem_=false,
			T *null_elem_=NULL)
		{ 
			if(subarray_length_<=0)
			{
				throw error("Error - subarray_length_<=0 in array_v2::array_v2\n",1);
			};
			d_subarray_length=subarray_length_;
			d_min_index=-1;
			d_max_index=-1;
			d_assign_null_elem=assign_null_elem_;
			if(assign_null_elem_)
			{
				d_null_elem=*null_elem_;
			};
		}

		~array_v2()
		{
			long int i;
			for(i=0;i<(long int)d_vector_of_subarrays.size();i++)
			{
				delete[]d_vector_of_subarrays[i];
			};
		}

		inline void clear()
		{
			deallocate_memory(d_max_index+1);
			d_min_index=-1;
			d_max_index=-1;
			d_vector_of_subarrays.clear();

		}

		inline void allocate_memory(
			long int ind_)//ind_ is not deleted
		{
			if(ind_>d_max_index)
			{
				long int ind_in_d_vector_of_subarrays=(long int)floor((double)ind_/(double)d_subarray_length);
				long int old_size=(long int)d_vector_of_subarrays.size();
				long int new_size=ind_in_d_vector_of_subarrays+1;
				long int j;
				for(j=old_size;j<new_size;j++)
				{
					T*pointer_tmp=new T[d_subarray_length];
					FSA_utils::assert_mem(pointer_tmp);

					d_vector_of_subarrays.push_back(pointer_tmp);
					if(d_assign_null_elem)
					{
						long int i;
						for(i=0;i<d_subarray_length;i++)
						{
							pointer_tmp[i]=d_null_elem;
						};
					};
				};

				d_max_index=d_subarray_length*new_size-1;
			};

		}


		inline void set_elem(
			long int ind_,
			T elem_)
		{
			if(ind_<d_min_index)
			{
				throw error("Unexpected error - ind_<d_min_index in set_elem\n",1);
			};

			long int ind_in_d_vector_of_subarrays=(long int)floor((double)ind_/(double)d_subarray_length);
			long int ind_in_subarray=(long int)(ind_%d_subarray_length);

			if(ind_>d_max_index)
			{
				allocate_memory(ind_);
			};

			d_vector_of_subarrays[ind_in_d_vector_of_subarrays][ind_in_subarray]=elem_;
		}


		inline T &get_elem(
			long int ind_)
		{
			if(ind_<d_min_index)
			{
				throw error("Unexpected error - ind_<d_min_index in get_elem\n",1);
			};

			if(ind_>d_max_index)
			{
				allocate_memory(ind_);
			};

			long int ind_in_d_vector_of_subarrays=(long int)floor((double)ind_/(double)d_subarray_length);
			long int ind_in_subarray=(long int)(ind_%d_subarray_length);

			return d_vector_of_subarrays[ind_in_d_vector_of_subarrays][ind_in_subarray];
		}

		inline void deallocate_memory(
			long int ind_)//ind_ is not deleted
		{
			long int ind_in_d_vector_of_subarrays=(long int)FSA_utils::Tmin(floor((double)d_max_index/(double)d_subarray_length),floor((double)ind_/(double)d_subarray_length)-1);
			long int ind_in_d_vector_of_subarrays_min=(long int)floor((double)d_min_index/(double)d_subarray_length);//always >=0
			if(ind_in_d_vector_of_subarrays_min<0)
			{
				ind_in_d_vector_of_subarrays_min=0;
			};


			long int i;

			for(i=ind_in_d_vector_of_subarrays_min;i<=ind_in_d_vector_of_subarrays;i++)
			{
				delete[]d_vector_of_subarrays[i];d_vector_of_subarrays[i]=NULL;
			};

			if(ind_in_d_vector_of_subarrays_min<ind_in_d_vector_of_subarrays+1)
			{
				d_min_index=(ind_in_d_vector_of_subarrays+1)*d_subarray_length;
			};
		}


		
	public:

		long int d_subarray_length;
		std::vector<T*> d_vector_of_subarrays;
		long int d_min_index;
		long int d_max_index;
		bool d_assign_null_elem;
		T d_null_elem;


		
	};

	template<class T>
	void array_positive<T>::increment_array(long int ind_)
	{
		T *d_elem_new=NULL;

		try
		{
			long int o_dim=d_dim;
			do{
			  d_dim+=d_step;
			}while(ind_>d_dim);

			d_elem_new=new T[d_dim+1];
			FSA_utils::assert_mem(d_elem_new);

			long int i;
			for(i=0;i<o_dim+1;i++)
			{
				d_elem_new[i]=d_elem[i];
			};

			for(i=o_dim+1;i<d_dim+1;i++)
			{
				d_elem_new[i]=0;
			};


			delete[]d_elem;d_elem=NULL;

			d_elem=d_elem_new;d_elem_new=NULL;
		}
		catch (...)
		{ 
			delete[]d_elem_new;d_elem_new=NULL;
			throw;
		};
		
	}

	template<class T>
	void array_v<T>::increment_array_on_the_right(long int ind_)
	{
		bool ee_error_flag=false;
		error ee_error("",0);
		T *d_elem_new=NULL;

		try
		{
		try
		{


			long int o_dim=d_dim;
			do{
			  d_dim+=d_step;
			  d_dim_plus_d_ind0+=d_step;
			}while(ind_>d_dim_plus_d_ind0);

			d_elem_new=new T[d_dim+1];
			FSA_utils::assert_mem(d_elem_new);

			long int i;
			for(i=0;i<o_dim+1;i++)
			{
				d_elem_new[i]=d_elem[i];
			};

			for(i=o_dim+1;i<d_dim+1;i++)
			{
				d_elem_new[i]=0;
			};

			delete[]d_elem;d_elem=NULL;
			d_elem=d_elem_new;d_elem_new=NULL;


		}
		catch (error er)
		{
			ee_error_flag=true;
			ee_error=er;		
		};
		}
		catch (...)
		{ 
			ee_error_flag=true;
			ee_error=error("Internal error in the program\n",4);
		};

		//memory release

		if(ee_error_flag)
		{
			delete[]d_elem_new;d_elem_new=NULL;
			throw error(ee_error.st,ee_error.error_code);
		};

	}

	template<class T>
		void array_v<T>::increment_array_on_the_left(long int ind_)
	{
		T *d_elem_new=NULL;

		try
		{
			long int o_dim=d_dim;
			do{
			  d_dim+=d_step;
			  d_ind0-=d_step;
			}while(ind_<d_ind0);
			long int jump=d_dim-o_dim;

			d_elem_new=new T[d_dim+1];
			FSA_utils::assert_mem(d_elem_new);

			long int i;

			for(i=0;i<jump;i++)
			{
				d_elem_new[i]=0;
			};

			for(i=0;i<o_dim+1;i++)
			{
				d_elem_new[i+jump]=d_elem[i];
			};


			delete[]d_elem;d_elem=NULL;
			d_elem=d_elem_new;d_elem_new=NULL;

		}
		catch (...)
		{
			delete[]d_elem_new;d_elem_new=NULL;
			throw;
		};


	}


}


#endif // GUMBEL_FSA1_UTILS

