#ifndef GUMBEL_FSA1
#define GUMBEL_FSA1

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

File name: sls_fsa1.hpp

Author: Sergey Sheetlin

Contents: Frameshift alignment algorithms 

******************************************************************************/

#include "sls_alp_regression.hpp"
#include "sls_fsa1_utils.hpp"
#include "sls_fsa1_parameters.hpp"
#include "sls_fsa1_pvalues.hpp"

#include "njn_localmaxstatmatrix.hpp"
#include "njn_localmaxstatutil.hpp"

#include "sls_fsa1_utils.hpp"

const double mb_bytes=1048576.0;

namespace Sls 
{
	const long int inf=LONG_MAX/10;

	class IS1_general_simulation;
	template<typename T> class two_dim_layer_alignment_algorithm;

	struct par_test1_type
	{
		std::string d_gumbelparout_file_name;
		long int d_seqlen1;
		long int d_seqlen2;
	};

	
	struct state_type
	{
		long int d_ALP_number;//number of ALPs simulated

		array_positive<long int> *d_M_array;//global maximum corresponded to the ALP
		array_positive<long int> *d_seq1_length_array;//length of sequence #1 corresponded to the ALP
		array_positive<long int> *d_seq2_length_array;//length of sequence #2 corresponded to the ALP
		array_positive<long int> *d_distance_along_direction_1_array;
		array_positive<long int> *d_distance_along_direction_2_array;
		array_positive<double> *d_ALP_weights_array;//IS weigths for the ALP points; dimension of the array is d_ALP_number+1



		IS1_general_simulation *d_IS1_general_simulation;
		two_dim_layer_alignment_algorithm<long int> *d_two_dim_layer_alignment_algorithm;

		//for matrix expanding
		bool d_expanding_finished_flag;
		long int d_seq1_length_after_expanding;
		long int d_seq2_length_after_expanding;
		long int d_M_after_expanding;
		array_v<double> d_scores_counts_after_expanding;


		state_type()
		{
			d_M_array=NULL;
			d_seq1_length_array=NULL;
			d_seq2_length_array=NULL;
			d_distance_along_direction_1_array=NULL;
			d_distance_along_direction_2_array=NULL;
			d_ALP_weights_array=NULL;

		}
	};

	struct mult_states_type
	{
		long int d_number_of_realizations;
		long int d_number_of_ALP;

		double d_total_calculation_time;//time in seconds; includes information from the old object

		double d_average_ALP_pos1_mult_ALP_pos2;//average of d_average_ALP_pos1*d_average_ALP_pos2
		double d_average_ALP_pos1;//average of ALP position for the sequence #1
		double d_average_ALP_pos2;//average of ALP position for the sequence #2
		double d_average_expanding_length1;//average length of crude sampling expansion for the sequence #1
		double d_average_expanding_length2;//average length of crude sampling expansion for the sequence #2
		double d_average_expanding_length1_mult_expanding_length2;//average of d_average_expanding_length1*d_average_expanding_length2

		//long long d_total_number_of_ALP_cells;
		//long long d_total_number_of_exp_cells;

		double d_total_number_of_ALP_cells;
		double d_total_number_of_exp_cells;

		long int d_total_number_of_ALPs;


		

		state_type *d_states;//dimension of the array is d_number_of_realizations

		mult_states_type()
		{
			d_states=NULL;
			d_number_of_realizations=-1;
		}
	};


	class FSA;
	class FSA_utils;



	class FSA
	{
		public:

		typedef void alignment_function_type(long int *seq1_,long int *seq2_,FSA *FAS_object_,void*par_);
		struct struct_for_alignment
		{
			long int d_score;//alignment score
		};




		FSA(//constructor

		bool reversed_seq_,//whether the sequences are reversed or not
		std::string alignment_type_,//type of alignment; possible values "local", "global"

		long int rand_,//randomization number
		long int open1_,//gap opening penalty for the nucleotide sequence #1
		long int open2_,//gap opening penalty for the amino acid sequence #2

		long int epen1_,//gap extension penalty for the nucleotide sequence #1
		long int epen2_,//gap extension penalty for the amino acid sequence #2

		long int gamma_,//frameshift penalty gamma

		std::string smatr_file_name_,//scoring matrix file name
		std::string RR1_file_name_,//background frequencies file name for the sequence #1
		std::string RR2_file_name_,//background frequencies file name for the sequence #2
		std::string DNA_codon_table_file_name_,//a name of a file with DNA codon table
		long int seq1_length_,//length of sequence #1
		long int seq2_length_,//length of sequence #2
		long int seq_number_);//number of tested alignments

		~FSA();



		long int random_AA1();

		long int random_AA2();

		long int convert_codon_into_AA(//returns AA corresponded to the codon
		long int *codon_);

		static void classical_global_alignment(
		long int *seq1_,//sequence #1
		long int *seq2_,//sequence #2
		FSA *FAS_object_,//a pointer to FSA object
		void*par_);//additonal parameters

		static void a1_global_alignment(
		long int *seq1_,//sequence #1
		long int *seq2_,//sequence #2
		FSA *FAS_object_,//a pointer to FSA object
		void*par_);//additonal parameters

		void DNA_AA_global(
		std::string distr_file_name_,
		alignment_function_type *alignment_function_);

		void AA_AA_global(
		std::string distr_file_name_,
		alignment_function_type *alignment_function_);

		void DNA_to_3_frames_AA_global(
		std::string distr_file_name_,
		alignment_function_type *alignment_function_);




		public:

		
		//input parameters

		bool d_reversed_seq;//whether the sequences are reversed or not
		std::string d_alignment_type;//type of alignment; possible values "local", "global"

		long int d_open1;//gap opening penalty for the nucleotide sequence #1
		long int d_open2;//gap opening penalty for the amino acid sequence #2


		long int d_epen1;//gap extension penalty for the nucleotide sequence #1
		long int d_epen2;//gap extension penalty for the amino acid sequence #2

		long int d_gamma;//frameshift penalty gamma


		double d_max_time;//maximum allowed calculation time in seconds
		double d_max_mem;//maximum allowed memory usage in MB
		double d_memory_size_in_MB;//current allocated memory size


		
		//additional parameters

		long int d_number_of_letters1;//number of letters for the sequence 1
		long int d_number_of_letters2;//number of letters for the sequence 2

		char *d_alphabet1;//alphabet letters for the sequence #1
		char *d_alphabet2;//alphabet letters for the sequence #2

		long int *d_alphabet1_to_long;//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
		long int *d_alphabet2_to_long;//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)


		long int** d_smatr;//scoring matrix

		double *d_RR1;//nucleotide probabilities
		double *d_RR1_sum;//probability distribution function for d_RR1
		long int *d_RR1_sum_elements;//numbers of nucleotide corresponded to d_RR1

		double *d_RR2;//AA probabilities
		double *d_RR2_sum;//probability distribution function for d_RR2
		long int *d_RR2_sum_elements;//numbers of AA corresponded to d_RR2

		long int d_seq1_length;//length of sequence #1
		long int d_seq2_length;//length of sequence #2
		long int d_seq_number;//number of tested alignments

		long int d_codon_length;//codon length 
		long int *d_codon_AA;//<codon code,AA number>

		bool d_rand_flag;
		long int d_random_factor;

		//temporary arrays for alignments
		long int **d_S_a1;
		long int **d_I_a1;
		long int **d_D_a1;

		long int *d_seq1;
		long int *d_seq2;

	};

	class IS1_general
	{

		public:

		IS1_general(
		long int alphabet_letters_number1_,//number of letters in the sequence #1
		long int alphabet_letters_number2_,//number of letters in the sequence #2
		double *RR1_,//background probability for the sequence #1
		double *RR2_,//background probability for the sequence #2
		long int number_of_states_,//number of states
		double **transition_probabilities_,//transition probabilities between states; matrix d_number_of_states x d_number_of_states
		std::pair<long int, long int> *states_description_,//description of the states; the index is a state number
		double ***states_distr_);//distributions of the states; the first index is a state number
								//the second and the third indexes correspond to an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second

		~IS1_general();
		
		static long int letters_to_code(//returns a unique code for the array of letters
		long int letters_number_,//total number of letters
		long int letters_dim_,//dimension of the array with letters
		long int*letters_);//array of letters
		

		static void code_to_letters(//returns a unique code for the array of letters
		long int code_,//input code
		long int letters_number_,//total number of letters
		long int letters_dim_,//dimension of the array with letters
		long int*letters_);//array of letters; the result
		
		static long int matr_indexes_to_code(
		long int code1_,//code #1
		long int code1_number_,//the range of code1_ is [0,code1_number_-1]
		long int code2_,//code #2
		long int code2_number_);//the range of code2_ is [0,code2_number_-1]

		static void code_to_matr_indexes(
		long int code_,//input code
		long int &code1_,//code #1; the result
		long int code1_number_,//the range of code1_ is [0,code1_number_-1]
		long int &code2_,//code #2; the result
		long int code2_number_);//the range of code2_ is [0,code2_number_-1]

		static double ** allocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
			long int alphabet_letters_number1_,//number of letters in the sequence #1
			long int alphabet_letters_number2_,//number of letters in the sequence #2
			std::pair<long int, long int> state_description_);//state description

		static void deallocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
			double **&state_distr_,
			long int alphabet_letters_number1_,//number of letters in the sequence #1
			std::pair<long int, long int> state_description_);//state description

		static void generate_one_state_sum_distr(//calculates sum-distributions corresponded to the input parameters
			std::pair<long int, long int> state_distr_dims_,//dimensions of the matrix state_distr_
			double **state_distr_,//state distribution in the same format as in d_states_distr[s]
			double *&state_sum_distr_);//the result; the dimention is state_distr__dim1 x state_distr__dim2

		void generate_random_letters_in_a_given_state(
			long int state_number_,//state number
			long int *letters1_,//the resulted letters for the sequence 1; the array must be allocated
			long int *letters2_);//the resulted letters for the sequence 2; the array must be allocated

		void one_step_of_the_importance_samping(//an initial state when sequences are empty must be defined outside
		long int current_state_,//the current state
		long int *seq1_add_letters_,//letters for the sequence #1; the array must be allocated
		long int *seq2_add_letters_,//letters for the sequence #2; the array must be allocated
		long int &new_state_);//a new state

		void allocate_states_distr_sums();//allocates d_states_distr_sums
		void deallocate_states_distr_sums();//deallocates d_states_distr_sums

		void calculate_inverse_matrices_for_the_infinite_sums(
			long int code1_,
			long int code2_,
			std::vector<std::vector<double> > *A1_inv_,
			std::vector<std::vector<double> > *A2_inv_,
			long int x1_,
			long int x2_);


			
		public:

		long int d_alphabet_letters_number1;//number of letters in the sequence #1
		long int d_alphabet_letters_number2;//number of letters in the sequence #2

		double *d_RR1;//background probability for the sequence #1
		double *d_RR2;//background probability for the sequence #2

		long int d_number_of_states;//number of states
		double **d_transition_probabilities;//transition probabilities between states; matrix d_number_of_states x d_number_of_states

		double **d_transition_probabilities_sum;//sum-distributions calculated from d_transition_probabilities[s]

		std::pair<long int, long int> *d_states_description;//description of the states; the index is a state number
		//.first: number of letters generated for the sequence #1
		//.second: number of letters generated for the sequence #2

		double ***d_states_distr;//distributions of the states; the first index is a state number
								//the second and the third indexes correspond to an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second

		double *****d_states_distr_sums;//distributions of the sums for d_states_distr; the first index is a state number
		//the second and the third indexes are number of letters retained for the first (x1) and the second (x2) sequence respectively
		//x1=0..state_description_type.first; x2=0..state_description_type.second
		//d_states_distr_sums[s][x1][x2] is the distribution and has the dimension d_A_letters^x1 x d_B_letters^x2

		std::vector<std::vector<double> > d_A1_inv;
		std::vector<std::vector<double> > d_A2_inv;
		//the second index is x1 and x2 respectively; the third index is the code; 
		//x1=0..state_description_type.first; x2=0..state_description_type.second
		//the dimentions of the array d_A1_inv[s][x1] is d_A_letters^x1
		//the dimentions of the array d_A2_inv[s][x2] is d_B_letters^x2


		//calculated values
		std::pair<long int, long int> *d_states_distr_dimensions;//contains dimensions of the matrix d_states_distr[s], s=0,...,d_number_of_states

		double **d_states_sum_distr;//sum-distributions calculated from d_states_distr
		long int *d_states_sum_distr_elements_for_all_states;//an array of with numbers 0,1,2,...
		//can be used for d_states_sum_distr[s] or d_transition_probabilities_sum

		std::pair<long int*, long int*> *d_tmp_letters;//an auxiliary array of length d_number_of_states; 
												//element #s contains arrays of dimensions <d_states_description[s].first,d_states_description[s].second>

		long int d_max_dim1;//maximum of d_states_description[s].first
		long int d_max_dim2;//maximum of d_states_description[s].second

	};

	template<typename T> class two_dim_layer
	{
		public:

		two_dim_layer(//constructor
		long int max_ind1_,//max of the index #1 (minimum index is 0 by default)
		long int max_ind2_,//max of the index #2 (minimum index is 0 by default)

		long int layers_number1_,//number of layers for the index #1
		long int layers_number2_,//number of layers for the index #2

		T null_)
		{
			if(max_ind1_<=0||max_ind2_<=0||layers_number1_<=0||layers_number2_<=0)
			{
				throw error("Error - parameters of two_dim_layer::two_dim_layer must be positive\n",1);
			};

			if(layers_number1_-1>max_ind2_)
			{
				throw error("Error - layers_number1_-1>max_ind2_ in two_dim_layer::two_dim_layer\n",1);
			};

			if(layers_number2_-1>max_ind1_)
			{
				throw error("Error - layers_number2_-1>max_ind1_ in two_dim_layer::two_dim_layer\n",1);
			};

			d_max_ind1=max_ind1_;//max of the index #1
			d_max_ind2=max_ind2_;//max of the index #2

			d_layers_number1=layers_number1_;//number of layers for the index #1
			d_layers_number2=layers_number2_;//number of layers for the index #2

			d_current_max_ind1=layers_number1_-1;
			d_current_max_ind2=layers_number2_-1;

			d_current_min_ind1=0;
			d_current_min_ind2=0;

			d_layers1=NULL;
			d_layers2=NULL;

			d_null=null_;

			d_layers1=new T*[d_layers_number1];
			FSA_utils::assert_mem(d_layers1);

			long int i;
			for(i=0;i<d_layers_number1;i++)
			{
				d_layers1[i]=NULL;
				d_layers1[i]=new T[d_max_ind2+1];
				FSA_utils::assert_mem(d_layers1[i]);
			};

			d_layers2=new T*[d_layers_number2];
			FSA_utils::assert_mem(d_layers2);

			for(i=0;i<d_layers_number2;i++)
			{
				d_layers2[i]=NULL;
				d_layers2[i]=new T[d_max_ind1+1];
				FSA_utils::assert_mem(d_layers2[i]);
			};

			d_layers1_tmp=NULL;
			d_layers1_tmp=new T*[d_layers_number1];
			FSA_utils::assert_mem(d_layers1_tmp);

			d_layers2_tmp=NULL;
			d_layers2_tmp=new T*[d_layers_number2];
			FSA_utils::assert_mem(d_layers2_tmp);

		}

		two_dim_layer(//constructor; copies data from two_dim_layer_ with minimally allocated memory
		long int max_ind1_,
		long int max_ind2_,
		two_dim_layer *two_dim_layer_)
		{

			d_max_ind1=max_ind1_;//max of the index #1
			d_max_ind2=max_ind2_;//max of the index #2

			d_layers_number1=two_dim_layer_->d_layers_number1;//number of layers for the index #1
			d_layers_number2=two_dim_layer_->d_layers_number2;//number of layers for the index #2

			d_current_max_ind1=two_dim_layer_->d_current_max_ind1;
			d_current_max_ind2=two_dim_layer_->d_current_max_ind2;

			d_current_min_ind1=two_dim_layer_->d_current_min_ind1;
			d_current_min_ind2=two_dim_layer_->d_current_min_ind2;

			d_layers1=NULL;
			d_layers2=NULL;

			d_null=two_dim_layer_->d_null;

			d_layers1=new T*[d_layers_number1];
			FSA_utils::assert_mem(d_layers1);

			long int i;
			for(i=0;i<d_layers_number1;i++)
			{
				d_layers1[i]=NULL;
				d_layers1[i]=new T[d_max_ind2+1];
				FSA_utils::assert_mem(d_layers1[i]);
			};

			d_layers2=new T*[d_layers_number2];
			FSA_utils::assert_mem(d_layers2);

			for(i=0;i<d_layers_number2;i++)
			{
				d_layers2[i]=NULL;
				d_layers2[i]=new T[d_max_ind1+1];
				FSA_utils::assert_mem(d_layers2[i]);
			};

			d_layers1_tmp=NULL;
			d_layers2_tmp=NULL;

			//copying data
			for(i=0;i<d_layers_number1;i++)
			{
				long int j;
				for(j=0;j<=d_max_ind2;j++)
				{
					d_layers1[i][j]=two_dim_layer_->d_layers1[i][j];
				};
			};

			for(i=0;i<d_layers_number2;i++)
			{
				long int j;
				for(j=0;j<=d_max_ind1;j++)
				{
					d_layers2[i][j]=two_dim_layer_->d_layers2[i][j];
				};
			};

		}

		static void copy_2_into_1_two_dim_layer(//copies the layers from the second argument to the first; both objects must be defined and fully functional
		two_dim_layer *two_dim_layer1_,
		two_dim_layer *two_dim_layer2_)
		{
			if(two_dim_layer1_->d_layers_number1!=two_dim_layer2_->d_layers_number1||
				two_dim_layer1_->d_layers_number2!=two_dim_layer2_->d_layers_number2)
			{
				throw error("Unexpected error - objects have different number of layers in copy_2_into_1_two_dim_layer\n",1);
			};

			if(two_dim_layer2_->d_max_ind2>two_dim_layer1_->d_max_ind2||
				two_dim_layer2_->d_max_ind1>two_dim_layer1_->d_max_ind1)
			{
				throw error("Unexpected error in copy_2_into_1_two_dim_layer\n",1);
			};

			two_dim_layer1_->d_current_max_ind1=two_dim_layer2_->d_current_max_ind1;
			two_dim_layer1_->d_current_max_ind2=two_dim_layer2_->d_current_max_ind2;

			two_dim_layer1_->d_current_min_ind1=two_dim_layer2_->d_current_min_ind1;
			two_dim_layer1_->d_current_min_ind2=two_dim_layer2_->d_current_min_ind2;


			//copying data
			long int i;
			for(i=0;i<two_dim_layer2_->d_layers_number1;i++)
			{
				long int j;
				for(j=0;j<=two_dim_layer2_->d_max_ind2;j++)
				{
					two_dim_layer1_->d_layers1[i][j]=two_dim_layer2_->d_layers1[i][j];
				};
			};

			for(i=0;i<two_dim_layer2_->d_layers_number2;i++)
			{
				long int j;
				for(j=0;j<=two_dim_layer2_->d_max_ind1;j++)
				{
					two_dim_layer1_->d_layers2[i][j]=two_dim_layer2_->d_layers2[i][j];
				};
			};

		}


		~two_dim_layer()
		{
			long int i;
			if(d_layers1)
			{
				for(i=0;i<d_layers_number1;i++)
				{
					delete[]d_layers1[i];
				};
				delete[]d_layers1;
			};

			if(d_layers2)
			{
				for(i=0;i<d_layers_number2;i++)
				{
					delete[]d_layers2[i];
				};
				delete[]d_layers2;
			};

			delete[]d_layers1_tmp;
			delete[]d_layers2_tmp;

		}

		T get_element(
			long int i1_,
			long int i2_)
		{
			if(d_current_min_ind1<=i1_&&i1_<=d_current_max_ind1&&i2_<=d_current_max_ind2)
			{
				long int layer1=i1_-d_current_min_ind1;
				return d_layers1[layer1][i2_];
			};

			if(d_current_min_ind2<=i2_&&i2_<=d_current_max_ind2&&i1_<d_current_min_ind1)
			{
				long int layer2=i2_-d_current_min_ind2;
				return d_layers2[layer2][i1_];
			};

			throw error("Error - element cannot be extracted in two_dim_layer::get_element\n",1);

		}

		void auxiliary_function_for_set_element(
		long int i1_,
		long int &current_min_ind1_,
		long int &current_max_ind1_,
		long int layers_number1_,
		T**&layers1_,
		T**&layers1_tmp_,

		long int current_min_ind2_,
		long int current_max_ind2_,
		T**layers2_)

		{
			if(i1_>current_max_ind1_)
			{
				long int current_max_ind1_new=i1_;
				long int current_min_ind1_new=i1_-layers_number1_+1;


				long int i;
				for(i=current_min_ind1_new;i<=current_max_ind1_;i++)
				{
					layers1_tmp_[i-current_min_ind1_new]=layers1_[i-current_min_ind1_];
				};

				long int ind_st=FSA_utils::Tmax(current_max_ind1_+1,current_min_ind1_new);
				for(i=ind_st;i<=current_max_ind1_new;i++)
				{
					layers1_tmp_[i-current_min_ind1_new]=layers1_[i-ind_st];
				};

				T**tmp_array=layers1_tmp_;
				layers1_tmp_=layers1_;
				layers1_=tmp_array;

				current_min_ind1_=current_min_ind1_new;
				current_max_ind1_=current_max_ind1_new;
				
				return;

			};

			if(i1_<current_min_ind1_)
			{
				long int current_min_ind1_new=i1_;
				long int current_max_ind1_new=i1_+layers_number1_-1;

				long int i;
				for(i=current_min_ind1_;i<=current_max_ind1_new;i++)
				{
					layers1_tmp_[i-current_min_ind1_new]=layers1_[i-current_min_ind1_];
				};

				long int ind_end=FSA_utils::Tmin(current_max_ind1_new,current_min_ind1_-1);
				long int ind_tmp=FSA_utils::Tmax(-current_min_ind1_new,layers_number1_-current_min_ind1_);

				for(i=current_min_ind1_new;i<=ind_end;i++)
				{
					layers1_tmp_[i-current_min_ind1_new]=layers1_[i+ind_tmp];
				};

				//transfer values
				
				for(i=current_min_ind1_new;i<=ind_end;i++)
				{
					long int j;
					for(j=current_min_ind2_;j<=current_max_ind2_;j++)
					{
						layers1_tmp_[i-current_min_ind1_new][j]=layers2_[j-current_min_ind2_][i];
					};
				};



				T**tmp_array=layers1_tmp_;
				layers1_tmp_=layers1_;
				layers1_=tmp_array;

				current_min_ind1_=current_min_ind1_new;
				current_max_ind1_=current_max_ind1_new;

			};
		}


		void set_element(
			long int i1_,
			long int i2_,
			T elem_)
		{
			if(i1_>d_max_ind1||i1_<0||i2_>d_max_ind2||i2_<0)
			{
				throw error("Error in the parameters i1_, i2_ of two_dim_layer::set_element\n",10);
			};

			bool success_flag0=false;
			if(d_current_min_ind1<=i1_&&i1_<=d_current_max_ind1&&i2_<=d_current_max_ind2)
			{
				long int layer1=i1_-d_current_min_ind1;
				d_layers1[layer1][i2_]=elem_;
				success_flag0=true;
			};

			if(d_current_min_ind2<=i2_&&i2_<=d_current_max_ind2&&i1_<=d_current_max_ind1)
			{
				long int layer2=i2_-d_current_min_ind2;
				d_layers2[layer2][i1_]=elem_;
				success_flag0=true;
			};

			if(success_flag0)
			{
				return;
			};


		
			auxiliary_function_for_set_element(
			i1_,
			d_current_min_ind1,
			d_current_max_ind1,
			d_layers_number1,
			d_layers1,
			d_layers1_tmp,
			
			d_current_min_ind2,
			d_current_max_ind2,
			d_layers2);

			

			auxiliary_function_for_set_element(
			i2_,
			d_current_min_ind2,
			d_current_max_ind2,
			d_layers_number2,
			d_layers2,
			d_layers2_tmp,
			
			d_current_min_ind1,
			d_current_max_ind1,
			d_layers1);


			
			bool success_flag=false;
			if(d_current_min_ind1<=i1_&&i1_<=d_current_max_ind1&&i2_<=d_current_max_ind2)
			{
				long int layer1=i1_-d_current_min_ind1;
				d_layers1[layer1][i2_]=elem_;
				success_flag=true;
			};

			if(d_current_min_ind2<=i2_&&i2_<=d_current_max_ind2&&i1_<=d_current_max_ind1)
			{
				long int layer2=i2_-d_current_min_ind2;
				d_layers2[layer2][i1_]=elem_;
				success_flag=true;
			};

			if(!success_flag)
			{
				throw error("Unexpected error in two_dim_layer::set_element\n",1);
			};
			

		}

		void set_max_ind(//if the element (i1_,i2_) is not available, changes the maximum and minimum indexes
			long int i1_,
			long int i2_)
		{
			bool success_flag=false;
			if(d_current_min_ind1<=i1_&&i1_<=d_current_max_ind1&&i2_<=d_current_max_ind2)
			{
				success_flag=true;
			};

			if(d_current_min_ind2<=i2_&&i2_<=d_current_max_ind2&&i1_<=d_current_max_ind1)
			{
				success_flag=true;
			};

			if(!success_flag)
			{
				set_element(i1_,i2_,d_null);
			};


		}

		public:

		long int d_max_ind1;//max of the index #1
		long int d_max_ind2;//max of the index #2

		long int d_current_max_ind1;//max of the index #1
		long int d_current_max_ind2;//max of the index #2

		long int d_current_min_ind1;//max of the index #1
		long int d_current_min_ind2;//max of the index #2


		long int d_layers_number1;//number of layers for the index #1
		long int d_layers_number2;//number of layers for the index #2

		T **d_layers1;//layers corresponded to the index #1; the dimension is d_layers_number1 x d_max_ind2
		T **d_layers2;//layers corresponded to the index #2; the dimension is d_layers_number2 x d_max_ind1

		T **d_layers1_tmp;
		T **d_layers2_tmp;

		T d_null;


	};

	template<typename T> class two_dim_layer_alignment_algorithm
	{//the class is desinged to work with two_dim_layer object

		public:

		typedef void two_dim_layer_alignment_function_type(long int i1_,long int i2_,two_dim_layer_alignment_algorithm<T> *obj_,void*par_);
			
		two_dim_layer_alignment_algorithm(
		long int number_of_variables_,//total number of variables in dynamic equations
		long int depth1_,//the maximum difference of the first index in the dynamic equations
		long int depth2_,//the maximum difference of the second index in the dynamic equations
		long int max_ind1_,//max of the index #1 (minimum index is 0 by default)
		long int max_ind2_,//max of the index #2 (minimum index is 0 by default)
		T null_,//null element of T
		long int *var_num_);//d_variable_name_to_number[c] converts a letter 'c' (of type char) to order number of variable

		two_dim_layer_alignment_algorithm(//used to save data from two_dim_layer_alignment_algorithm_ with minimal memory allocation
		long int max_ind1_,//max of the index #1
		long int max_ind2_,//max of the index #2
		two_dim_layer_alignment_algorithm<T> *two_dim_layer_alignment_algorithm_);


		~two_dim_layer_alignment_algorithm();

		void align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
		//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
		long int target_align1_length_,//target length for the sequence #1
		long int target_align2_length_);//target length for the sequence #2

		
		void align_upto_target_lengths_step(//aligns the sequences upto the target sequence lengths including; checks whether the jump exeeds d_step
		//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
		long int target_align1_length_,//target length for the sequence #1
		long int target_align2_length_);//target length for the sequence #2

		void init();//initialization for a new alignment




		public:

		long int d_number_of_variables;//total number of variables in dynamic equations
		long int d_depth1;//the maximum difference of the first index in the dynamic equations
		long int d_depth2;//the maximum difference of the second index in the dynamic equations
		long int *d_var_num;//d_variable_name_to_number[c] converts a letter 'c' (of type char) to order number of variable

		long int d_max_ind1;//max of the index #1 (minimum index is 0 by default)
		long int d_max_ind2;//max of the index #2 (minimum index is 0 by default)

		long int d_current_align_ind1;//the alignment is calculated in the rectangle [0,d_current_ind1] x [0,d_current_ind2]
		long int d_current_align_ind2;

		long int d_step1;//maximum jump in sequence #1 lengths permitted in align_upto_target_lengths_step
		long int d_step2;//maximum jump in sequence #2 lengths permitted in align_upto_target_lengths_step

		long int d_edge_maximum_calculation_depth1;//the maximum difference of the first index for the edge maximum calculation
		long int d_edge_maximum_calculation_depth2;//the maximum difference of the second index for the edge maximum calculation


		T d_null;//null element


		two_dim_layer<T> ** d_vars;//variables; the dimension of the array is d_number_of_variables

		two_dim_layer_alignment_function_type * d_two_dim_layer_alignment_function;//the function used for the alingment
		two_dim_layer_alignment_function_type * d_two_dim_layer_boundary_function;//the function used for the boundary conditions; must be defined in the rectangle [0,depth1_-1] x [0,depth2_-1]

		void *d_par;//parameters for the last argument void*par_ of functions d_two_dim_layer_boundary_function and d_two_dim_layer_alignment_function

		//statistics
		//global maximum
		bool d_M_flag;
		long int d_M;

		//edge maximum
		bool d_E_flag;
		long int d_E;

		bool d_scores_counts_flag;
		array_v<double> d_scores_counts;

		long int d_M_threshold;

		//FSC parameters
		bool d_FSC_flag;
		long int d_distance_along_direction_1;
		long int d_distance_along_direction_2;



	};

//-----------------------------------------------------

	class IS1_general_simulation
	{
		public:
		IS1_general_simulation(
			IS1_general *IS1_general_obj_,
			long int initial_state_,//initial state for the IS
			long int max_seq1_length_,//maximum sequence length
			long int max_seq2_length_,//maximum sequence length
			double *initial_distr_=NULL);//initial distribution of states; initial_state_ is ignored if d_initial_distr_ is defined

		IS1_general_simulation(//creates a copy of IS1_general_simulation_ with smallest possible memory allocation
		IS1_general_simulation *IS1_general_simulation_);//maximum sequence length


		~IS1_general_simulation();

		void init();//initialization for a new sequence generation

		void simulate_upto_target_lengths(
			long int target_seq1_length_,//target length of the sequence #1
			long int target_seq2_length_);//target length of the sequence #2


		//weights calculation
		void calculate_weight_W1_one_step(//calculates the weight assuming all the weights used in the recursive equations are defined
			long int target_length1_,//target length for sequence #1
			long int target_length2_,//target length for sequence #2
			long int i1_,//the first index (corresponded to the sequence #1)
			long int i2_);//the second index (corresponded to the sequence #2)

		void calculate_weight_W1_upto_target_lengths(
			long int target_seq1_length_,//target length of the sequence #1
			long int target_seq2_length_,//target length of the sequence #2
			long int seq1_length_tmp_,//the weights are calculated upto this length for the sequence #1
			long int seq2_length_tmp_);//the weights are calculated upto this length for the sequence #2


		void calculate_weight_W1_for_fixed_lengths_using_infinite_sums(
			long int target_seq1_length_,//target length of the sequence #1
			long int target_seq2_length_,//target length of the sequence #2
			double &weight_);//the resulted weight

		void calculate_weight_W1_for_fixed_lengths_using_recursions(
			long int target_seq1_length_,//target length of the sequence #1
			long int target_seq2_length_,//target length of the sequence #2
			double &weight_,//the resulted weight
			bool save_matrices_=false,
			std::vector<std::vector<double> > *A1_=NULL,
			std::vector<std::vector<double> > *A2_=NULL);





		public:

		long int d_seq1_current_length;//current length of the sequence #1 simulated using the IS
		long int d_seq2_current_length;//current length of the sequence #2 simulated using the IS

		long int d_max_seq1_length;//maximum sequence length
		long int d_max_seq2_length;//maximum sequence length

		long int *d_seq1;
		long int *d_seq2;

		long int d_initial_state;
		long int d_current_state;
		double *d_initial_distr_sum;

		IS1_general *d_IS1_general_obj;

		//weights calculation

		two_dim_layer<double> ** d_W1;//array of pointers to weights; dimention of the array is d_IS1_general_obj->d_number_of_states
		long int d_W1_seq1_current_length;//current length of the sequence #1 simulated using the IS
		long int d_W1_seq2_current_length;//current length of the sequence #2 simulated using the IS

		long int d_W1_step1;//increase of number of layers of d_W1 for the index #1
		long int d_W1_step2;//increase of number of layers of d_W1 for the index #2



	};

//-------------------------------------------------------------------------------------------
	//objects for tests

	class test
	{
		public:

		struct data_for_lambda_equation//struct for lambda_equation
		{
			long int d_number_of_AA;//number of AA
			long int** d_smatr;//scoring matrix
			double *d_RR1;//AA probabilities
			double *d_RR2;//AA probabilities
		};

		struct data_for_lambda_equation_FSA//struct for lambda_equation_FSA
		{
			long int d_number_of_letters1;//number of letters for the sequence #1
			long int d_number_of_letters2;//number of letters for the sequence #2
			long int** d_smatr;//scoring matrix
			double *d_RR1;//nucleotide probabilities
			double *d_RR2;//AA probabilities

			long int d_codon_length;//codon length 
			long int *d_codon_AA;//<codon code,AA number>

		};


		struct data_for_classical_alignment
		{
			long int d_open1;//gap opening penalty for the nucleotide sequence #1
			long int d_open2;//gap opening penalty for the amino acid sequence #2

			long int d_epen1;//gap extension penalty for the nucleotide sequence #1
			long int d_epen2;//gap extension penalty for the amino acid sequence #2

			long int**d_smatr;//the scoring matrix

			long int *d_seq1;//sequence #1
			long int *d_seq2;//sequence #2

		};

		struct data_for_FSA_alignment
		{
			long int d_alphabet_letters_number1;//number of letters in the sequence #1

			long int d_open1;//gap opening penalty for the nucleotide sequence #1
			long int d_open2;//gap opening penalty for the amino acid sequence #2

			long int d_epen1;//gap extension penalty for the nucleotide sequence #1
			long int d_epen2;//gap extension penalty for the amino acid sequence #2

			long int d_gamma;//frame shift penalty

			long int**d_smatr;//the scoring matrix

			long int *d_seq1;//sequence #1
			long int *d_seq2;//sequence #2

			long int *d_codon_AA;//<codon code,AA number>

			bool d_insertions_after_deletions;//if true, then insertions after deletions are allowed

			std::string d_alignment_type;//possible values are "global" or "local"

		};


		void classical_IS(

		unsigned rand_,//randomization number
		long int open1_,//gap opening penalty for the nucleotide sequence #1
		long int open2_,//gap opening penalty for the amino acid sequence #2

		long int epen1_,//gap extension penalty for the nucleotide sequence #1
		long int epen2_,//gap extension penalty for the amino acid sequence #2


		std::string smatr_file_name_,//scoring matrix file name
		std::string RR1_file_name_,//background frequencies file name for the sequence #1
		std::string RR2_file_name_);//background frequencies file name for the sequence #2

		//long int seq1_length_,//length of sequence #1
		//long int seq2_length_,//length of sequence #2

		//long int seq_number_);//number of tested alignments

		static long int realization_number_into_set_number(//numeration of sets starts from 1
		long int &realization_number_,//realization order number; the numeration starts from 0
		long int &number_of_realizations_set_);//total number of sets

		static std::pair<long int,long int> set_number_boundaries(//boundaries of realization order numbers for the set; numeration of realizations starts from 0
		long int &set_number_,//set order number; the numeration starts from 1
		long int &number_of_realizations_set_);//total number of realizations per set

		static void combine_parameters_from_forward_and_reversed_calculations_generalized(
		Sls::FALP_set_of_parameters &par_,//parameters from forward calculation
		Sls::FALP_set_of_parameters &par_reversed_,//parameters from reversed calculation
		Sls::FALP_set_of_parameters &par_result_);//the result


		void input_data_for_the_constructor(

		std::string DNA_codon_table_file_name_,//a name of a file with DNA codon table
		std::string smatr_file_name_,//scoring matrix file name
		std::string RR1_file_name_,//probabilities1 file name
		std::string RR2_file_name_,//probabilities2 file name

		long int &alphabetSize1_,
		long int &alphabetSize2_,
		long int **&substitutionScoreMatrix_,
		double *&letterFreqs1_,
		double *&letterFreqs2_,

		char *&alphabet1_,//alphabet letters for the sequence #1
		char *&alphabet2_,//alphabet letters for the sequence #2

		long int *&alphabet1_to_long_,//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
		long int *&alphabet2_to_long_,//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

		long int &codon_length_,//codon length 
		long int *&codon_AA_);//<codon code,AA number>

		void FSA_IS(
		//additional parameters for the library code
		bool library_call_flag_,//if true, then the additional parameters are used
		long ntAlphabetSize_,//a number of letters in the DNA alphabet; usually 4
		long aaAlphabetSize_,//a number of letters in the amino acid alphabet; usually 4
		const long *codonTable_,//defines the codon table; is an array of ntAlphabetSize_^3 integers
		const long *const *substitutionScoreMatrix_,//scoring matrix; aaAlphabetSize_ X aaAlphabetSize_
		const double *ntFreqs_,//background frequencies of letters in DNA sequences
		const double *aaFreqs_,//background frequencies of letters in amino acid sequences
		//additional parameters for the library code - end

		long int rand_,//randomization number
		long int open1_,//gap opening penalty for the nucleotide sequence #1
		long int open2_,//gap opening penalty for the amino acid sequence #2

		long int epen1_,//gap extension penalty for the nucleotide sequence #1
		long int epen2_,//gap extension penalty for the amino acid sequence #2

		long int gamma_,//frameshift penalty gamma


		std::string smatr_file_name_,//scoring matrix file name
		std::string RR1_file_name_,//background frequencies file name for the sequence #1
		std::string RR2_file_name_,//background frequencies file name for the sequence #2
		std::string DNA_codon_table_file_name_,//a name of a file with DNA codon table

		double eps_lambda_,//relative error for lambda calculation
		double eps_K_,//relative error for K calculation

		bool gapped_flag_,//if true, then the gapped alingment is performed
		double max_time_,//maximum allowed calculation time in seconds
		double max_mem_,//maximum allowed memory usage in MB

		long int seq_number_,//number of tested alignments
		long int nalp_,//number of ALPs for the calculation
		long int number_of_subsets_for_errors_calculation_,//number of subsets used for the splitting method
		bool forward_and_reverse_screen_output_flag_,//determines whether the parameters are outputted for forward and reverse calculations
		bool insertions_after_deletions_,//if true, then insertions after deletions are allowed

		//for test
		double mult_for_is_lambda_,//multiplier for lambda in the IS

		//the result
		Sls::FALP_set_of_parameters &par_result_,//the resulted parameters
		Sls::par_test1_type *par_test1_=NULL);//for tests




		static void FSA_IS_transition_probabilities_calculation(
		bool &FSA_flag,
		double ***&states_distr,

		bool &cs_flag,
		double ***&states_distr_cs,

		bool &sim_flag,
		double ***&states_distr_sim,

		long int &number_of_states,
		long int &alphabet_letters_number1,
		long int &alphabet_letters_number2,
		std::pair<long int, long int> *&states_description,
		long int &codon_length,
		long int *&codon_AA,
		double *&RR1,
		double *&RR2,
		std::map<std::string, long int> &state_name_into_number,
		long int **&smatr,
		double &ungappedlambda);


		static void calculate_ungapped_lambda(
		long int number_of_AA_,
		double *RR1_,
		double *RR2_,
		long int**smatr_,
		double &ungapped_lambda_);

		static void calculate_ungapped_lambda_general(
		function_type *func_,
		void* func_pointer_,
		double &ungapped_lambda_);




		static double lambda_equation(double x_,void* func_number_);
		static double lambda_equation_FSA(double x_,void* func_number_);

		static void two_dim_layer_alignment_function_classical_global(
			//function for global classical alignment
			long int i1_,
			long int i2_,
			two_dim_layer_alignment_algorithm<long int> *obj_,
			void*par_);

		static void two_dim_layer_alignment_function_FSA(//function for the FSA a1-alignment
			long int i1_,
			long int i2_,
			two_dim_layer_alignment_algorithm<long int> *obj_,
			void*par_);

		static void collect_and_output_M_distr_upto_fixed_lengths(
			long int number_of_realizations_,
			std::string M_distr_file_name_,
			two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_,
			IS1_general_simulation &IS1_general_simulation_test_,
			long int target_seq1_length_,//target length of the sequence #1
			long int target_seq2_length_);//target length of the sequence #2

		static void collect_and_output_E_distr_upto_fixed_ALP(
			long int number_of_realizations_,
			std::string E_distr_file_name_,
			two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_,
			IS1_general_simulation &IS1_general_simulation_test_,
			long int target_ALP_,//target ALP number
			double ungapped_lambda_,
			long int limit2_,
			long int number_of_sets_,//number of sets for error calculation
			Sls::FALP_set_of_parameters &par_,
			bool screen_output_flag_,

			bool futher_expanding_=false,
			IS1_general *IS1_general_cs_=NULL,
			double *eps_K_=NULL,//relative error for K
			long int *number_of_cs_steps_for_expanding_=NULL);

		static void delete_mult_states_type(
			mult_states_type *&states_old1);

		static void restore_arrays(
			long int k_,
			mult_states_type *states_old_,
			two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_,
			IS1_general_simulation &IS1_general_simulation_test_);

		static void restore_arrays(
			two_dim_layer_alignment_algorithm<long int> *two_dim_layer_alignment_algorithm_test_from_,
			IS1_general_simulation *IS1_general_simulation_test_from_,
			two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_,
			IS1_general_simulation &IS1_general_simulation_test_);

		static void compare_mult_states(
			mult_states_type *states_old_,
			mult_states_type *states_new_);



		static void collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
			bool compute_allocated_memory_flag_,//if true then total allocated memory is computed in allocated_memory_
			double &allocated_memory_,


			long int number_of_realizations_,
			std::string E_distr_file_name_,
			two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_,
			IS1_general_simulation &IS1_general_simulation_test_,
			long int target_ALP_,//target ALP number
			double ungapped_lambda_,
			long int limit2_,
			long int number_of_sets_,//number of sets for error calculation
			Sls::FALP_set_of_parameters &par_,
			bool &inside_simulation_flag_,
			bool screen_output_flag_,
			bool screen_output_k_flag_,

			double ending_time_,
			bool &stopped_by_ending_time_flag_,
			long int &number_of_realizations_with_ending_time_,

			double ending_time_to_test_logarithmic_regime_,
			
			bool save_states_flag_,
			mult_states_type *states_old_,//NULL or defined
			mult_states_type *&states_new_,//resulted value
			long int &M_thr_estimation_,//average score threshold corresponded to target_ALP_


			bool futher_expanding_=false,
			long int M_thr_=-inf,//the expanding starts from the ALP with the score >=M_thr_; used if >=0
			IS1_general *IS1_general_cs_=NULL,
			double *eps_K_=NULL,//relative error for K
			long int *number_of_cs_steps_for_expanding_=NULL,
			bool parameters_are_not_calculated_=false);

		static void FSA_Align(

			long int open1_,//gap opening penalty for the nucleotide sequence #1
			long int open2_,//gap opening penalty for the amino acid sequence #2

			long int epen1_,//gap extension penalty for the nucleotide sequence #1
			long int epen2_,//gap extension penalty for the amino acid sequence #2

			long int gamma_,//frameshift penalty gamma

			std::string smatr_file_name_,//scoring matrix file name
			std::string DNA_codon_table_file_name_,//a name of a file with DNA codon table

			bool insertions_after_deletions_,//if true, then insertions after deletions are allowed

			std::string input_sequences_,//name of file with sequences for alignment (align mode)
			std::string output_sequences_,//name of output file with alignment scores (align mode)
			std::string gumbelparin_file_name_);//Gumbel parameters input file name


		static void compare_direct_and_reverse_sampling(
			long int number_of_realizations_,
			long int seq2_length_,
			two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_,
			IS1_general_simulation &IS1_general_simulation_test_,
			two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_reverse_,
			IS1_general_simulation &IS1_general_simulation_test_reverse_,
			IS1_general *IS1_general_cs_);



		static void normalize_state_distributions(
			long int alphabet_letters_number1_,//number of letters in the sequence #1
			long int alphabet_letters_number2_,//number of letters in the sequence #2

			long int number_of_states_,//number of states

			std::pair<long int, long int> *states_description_,//description of the states; the index is a state number
			double ***states_distr_);//distributions of the states; the index is a state number





		public:


		long int **d_alignment_matr;

		IS1_general* d_IS1;

		IS1_general_simulation *d_IS1_general_simulation;


	};




}

#endif // GUMBEL_FSA1

