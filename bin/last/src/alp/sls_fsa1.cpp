/* $Id: $
* ===========================================================================
*
* PUBLIC DOMAIN NOTICE
* National Center for Biotechnology Information
*
* This software/database is a "United States Government Work" under the
* terms of the United States Copyright Act. It was written as part of
* the author's offical duties as a United States Government employee and
* thus cannot be copyrighted. This software/database is freely available
* to the public for use. The National Library of Medicine and the U.S.
* Government have not placed any restriction on its use or reproduction.
*
* Although all reasonable efforts have been taken to ensure the accuracy
* and reliability of the software and data, the NLM and the U.S.
* Government do not and cannot warrant the performance or results that
* may be obtained by using this software or data. The NLM and the U.S.
* Government disclaim all warranties, express or implied, including
* warranties of performance, merchantability or fitness for any particular
* purpose.
*
* Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/*****************************************************************************

File name: sls_repwords.cpp

Author: Sergey Sheetlin

Contents: Frameshift alignment algorithms 

******************************************************************************/

#include "sls_fsa1.hpp"

using namespace Sls;
using namespace std;

static bool test_mode=true;//of true, then the test mode is activated

static long int small_long=(long int)((double)LONG_MIN/2.0);
static long int length_max=1000;

FSA::FSA(//constructor

bool reversed_seq_,//whether the sequences are reversed or not
string alignment_type_,//type of alignment; possible values "local", "global"

long int rand_,//randomization number
long int open1_,//gap opening penalty for the nucleotide sequence #1
long int open2_,//gap opening penalty for the amino acid sequence #2

long int epen1_,//gap extension penalty for the nucleotide sequence #1
long int epen2_,//gap extension penalty for the amino acid sequence #2

long int gamma_,//frameshift penalty gamma

string smatr_file_name_,//scoring matrix file name
string RR1_file_name_,//background frequencies file name for the sequence #1
string RR2_file_name_,//background frequencies file name for the sequence #2
string DNA_codon_table_file_name_,//a name of a file with DNA codon table
long int seq1_length_,//length of sequence #1
long int seq2_length_,//length of sequence #2
long int seq_number_)//number of tested alignments
{

	d_alphabet1=NULL;
	d_alphabet2=NULL;

	d_alphabet1_to_long=NULL;
	d_alphabet2_to_long=NULL;

	d_smatr=NULL;
	d_RR1=NULL;
	d_RR1_sum=NULL;
	d_RR1_sum_elements=NULL;

	d_RR2=NULL;
	d_RR2_sum=NULL;
	d_RR2_sum_elements=NULL;

	d_S_a1=NULL;
	d_I_a1=NULL;
	d_D_a1=NULL;

	d_codon_AA=NULL;

	d_seq1=NULL;
	d_seq2=NULL;


	try
	{

		d_reversed_seq=reversed_seq_;
		d_alignment_type=alignment_type_;

		long int number_of_AA_RR2;
		long int number_of_AA_smatr;
		long int smatr_min;

		FSA_utils::read_smatr(
		smatr_file_name_,
		d_smatr,
		number_of_AA_smatr,
		smatr_min);

		

		FSA_utils::read_RR(
		RR1_file_name_,
		d_RR1,
		d_RR1_sum,
		d_RR1_sum_elements,
		d_number_of_letters1);


		FSA_utils::read_RR(
		RR2_file_name_,
		d_RR2,
		d_RR2_sum,
		d_RR2_sum_elements,
		number_of_AA_RR2);


		if(number_of_AA_RR2==number_of_AA_smatr)
		{
			d_number_of_letters2=number_of_AA_smatr;
		}
		else
		{
			throw error("The number of letters is different in the files "+smatr_file_name_+" and "+RR2_file_name_+"\n",3);
		};

		long int number_of_letters1_tmp;//number of letters for the sequence 1
		long int number_of_letters2_tmp;//number of letters for the sequence 2

		FSA_utils::read_codon_AA_file(
		DNA_codon_table_file_name_,
		number_of_letters1_tmp,//number of letters for the sequence 1
		number_of_letters2_tmp,//number of letters for the sequence 2

		d_alphabet1,//alphabet letters for the sequence #1
		d_alphabet2,//alphabet letters for the sequence #2

		d_alphabet1_to_long,//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
		d_alphabet2_to_long,//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

		d_codon_length,//codon length 
		d_codon_AA);//<codon code,AA number>


		if(d_number_of_letters1!=number_of_letters1_tmp)
		{
			throw error("The number of letters is different in the files "+DNA_codon_table_file_name_+" and "+RR1_file_name_+"\n",3);
		};
		if(d_number_of_letters2!=number_of_letters2_tmp)
		{
			throw error("The number of letters is different in the files "+DNA_codon_table_file_name_+" and "+RR2_file_name_+"\n",3);
		};

		d_open1=open1_;
		d_open2=open2_;

		d_epen1=epen1_;
		d_epen2=epen2_;

		d_gamma=gamma_;

		d_seq1_length=seq1_length_;//length of sequence #1
		d_seq2_length=seq2_length_;//length of sequence #2
		d_seq_number=seq_number_;//number of tested alignments



		//randomization
		long int random_factor=rand_;


		if(random_factor<0)
		{
			random_factor=(long int)time(NULL);
			#ifndef _MSC_VER //UNIX program
				struct timeval tv;
				struct timezone tz;
				gettimeofday(&tv, &tz);
				random_factor+=tv.tv_usec*10000000;
			#else
				struct _timeb timebuffer;
				char *timeline;
				_ftime( &timebuffer );
				timeline = ctime( & ( timebuffer.time ) );
				rand_+=timebuffer.millitm*10000000;
			#endif

			random_factor=abs(random_factor);

			d_rand_flag=false;

		};

		d_random_factor=random_factor;
		cout<<"Random seed "<<d_random_factor<<endl;

		FSA_utils::srand2(d_random_factor);

		FSA_utils::get_memory_for_matrix(d_seq1_length+1,d_seq2_length+1,d_S_a1);
		FSA_utils::get_memory_for_matrix(d_seq1_length+1,d_seq2_length+1,d_I_a1);
		FSA_utils::get_memory_for_matrix(d_seq1_length+1,d_seq2_length+1,d_D_a1);

		d_seq1=new long int[d_seq1_length];
		d_seq2=new long int[d_seq2_length];

	}
	catch (...)
	{
		this->~FSA();
		throw;
	};

}

FSA::~FSA()
{
	if(d_smatr)
	{
		FSA_utils::delete_memory_for_matrix(d_number_of_letters2,d_smatr);
	};

	delete[]d_RR1;
	delete[]d_RR2;
	delete[]d_RR1_sum;
	delete[]d_RR2_sum;
	delete[]d_RR1_sum_elements;
	delete[]d_RR2_sum_elements;
	delete[]d_alphabet1;
	delete[]d_alphabet2;

	delete[]d_alphabet1_to_long;
	delete[]d_alphabet2_to_long;
	delete[]d_codon_AA;

	if(d_S_a1)
	{
		FSA_utils::delete_memory_for_matrix(d_seq1_length+1,d_S_a1);
	};
	if(d_I_a1)
	{
		FSA_utils::delete_memory_for_matrix(d_seq1_length+1,d_I_a1);
	};
	if(d_D_a1)
	{
		FSA_utils::delete_memory_for_matrix(d_seq1_length+1,d_D_a1);
	};


}

long int FSA::random_AA1()
{
	return FSA_utils::random_long(
		FSA_utils::ran2(),
			d_number_of_letters1,
			d_RR1_sum,
			d_RR1_sum_elements);
}

long int FSA::random_AA2()
{
	return FSA_utils::random_long(
		FSA_utils::ran2(),
			d_number_of_letters2,
			d_RR2_sum,
			d_RR2_sum_elements);
}



long int FSA::convert_codon_into_AA(
long int *codon_)
{
	return FSA_utils::convert_codon_into_AA(
			d_codon_length,//codon length 
			d_codon_AA,//<codon code,AA number>
			d_number_of_letters1,//number of letters for the sequence 1
			codon_);
}

void FSA::classical_global_alignment(
long int *seq1_,//sequence #1
long int *seq2_,//sequence #2
FSA *FAS_object_,//a pointer to FSA object
void*par_)//additonal parameters
{

	bool local_flag=(FAS_object_->d_alignment_type=="local");

	long int init_for_S=-inf;
	if(local_flag)
	{
		init_for_S=0;
	};

	long int i,j;
	//initial conditions
	for(i=0;i<=FAS_object_->d_seq1_length;i++)
	{
		if(local_flag)
		{
			FAS_object_->d_D_a1[i][0]=-inf;
		}
		else
		{
			FAS_object_->d_D_a1[i][0]=-FAS_object_->d_open1-FAS_object_->d_epen1*i;
		};

		FAS_object_->d_S_a1[i][0]=init_for_S;
		FAS_object_->d_I_a1[i][0]=-inf;
	};

	for(j=0;j<=FAS_object_->d_seq2_length;j++)
	{
		FAS_object_->d_S_a1[0][j]=init_for_S;
		FAS_object_->d_D_a1[0][j]=-inf;
		if(local_flag)
		{
			FAS_object_->d_I_a1[0][j]=-inf;
		}
		else
		{
			FAS_object_->d_I_a1[0][j]=-FAS_object_->d_open2-FAS_object_->d_epen2*j;
		};

	};

	FAS_object_->d_S_a1[0][0]=0;

	for(i=1;i<=FAS_object_->d_seq1_length;i++)
	{
		for(j=1;j<=FAS_object_->d_seq2_length;j++)
		{

			FAS_object_->d_I_a1[i][j]=FSA_utils::Tmax(FAS_object_->d_S_a1[i][j-1]-FAS_object_->d_open2-FAS_object_->d_epen2,FAS_object_->d_I_a1[i][j-1]-FAS_object_->d_epen2);
			FAS_object_->d_D_a1[i][j]=FSA_utils::Tmax(FAS_object_->d_S_a1[i-1][j]-FAS_object_->d_open1-FAS_object_->d_epen1,FAS_object_->d_D_a1[i-1][j]-FAS_object_->d_epen1);

			long int smart_score=FAS_object_->d_smatr[seq1_[i-1]][seq2_[j-1]];
			FAS_object_->d_S_a1[i][j]=
				FSA_utils::Tmax(
					init_for_S,
					FSA_utils::Tmax(FAS_object_->d_S_a1[i-1][j-1],FAS_object_->d_D_a1[i-1][j-1],FAS_object_->d_I_a1[i-1][j-1])+smart_score
					);
			
		};
	};

	long int score=FAS_object_->d_S_a1[0][0];
	for(i=0;i<=FAS_object_->d_seq1_length;i++)
	{
		for(j=0;j<=FAS_object_->d_seq2_length;j++)
		{
			score=FSA_utils::Tmax(score,
				FSA_utils::Tmax(FAS_object_->d_S_a1[i][j],FAS_object_->d_D_a1[i][j],FAS_object_->d_I_a1[i][j])
				);
		};
	};

	struct_for_alignment &tmp_obj=*((struct_for_alignment*)par_);
	tmp_obj.d_score=score;

}


void FSA::a1_global_alignment(
long int *seq1_,//sequence #1
long int *seq2_,//sequence #2
FSA *FAS_object_,//a pointer to FSA object
void*par_)//additonal parameters
{
	bool local_flag=(FAS_object_->d_alignment_type=="local");

	long int init_for_S=-inf;
	if(local_flag)
	{
		init_for_S=0;
	};


	long int i,j;
	//initial conditions
	for(i=0;i<=FAS_object_->d_seq1_length;i++)
	{
		long int i_div_3;

		if(local_flag)
		{
			FAS_object_->d_D_a1[i][0]=-inf;
		}
		else
		{
			if(i%3==0)
			{
				i_div_3=i/3;
				FAS_object_->d_D_a1[i][0]=-FAS_object_->d_open1-FAS_object_->d_epen1*i_div_3;
			};

			if(i%3==1)
			{
				i_div_3=(i-1)/3;
				FAS_object_->d_D_a1[i][0]=-FAS_object_->d_open1-FAS_object_->d_epen1*i_div_3-FAS_object_->d_gamma;
			};

			if(i%3==2)
			{
				i_div_3=(i+1)/3;
				FAS_object_->d_D_a1[i][0]=-FAS_object_->d_open1-FAS_object_->d_epen1*i_div_3-FAS_object_->d_gamma;
			};

			i_div_3=(long int)floor((double)i/3.0);
			FAS_object_->d_D_a1[i][0]=-FAS_object_->d_open1-FAS_object_->d_epen1*i_div_3;
			
		};


		FAS_object_->d_S_a1[i][0]=init_for_S;
		FAS_object_->d_I_a1[i][0]=-inf;
	};

	for(j=0;j<=FAS_object_->d_seq2_length;j++)
	{
		FAS_object_->d_S_a1[0][j]=init_for_S;
		FAS_object_->d_D_a1[0][j]=-inf;
		if(local_flag)
		{
			FAS_object_->d_I_a1[0][j]=-inf;
		}
		else
		{
			FAS_object_->d_I_a1[0][j]=-FAS_object_->d_open2-FAS_object_->d_epen2*j;
		};

		FAS_object_->d_S_a1[1][j]=init_for_S;
		FAS_object_->d_D_a1[1][j]=-inf;
		if(local_flag)
		{
			FAS_object_->d_I_a1[1][j]=-inf;
		}
		else
		{
			FAS_object_->d_I_a1[1][j]=-FAS_object_->d_open2-FAS_object_->d_epen2*j;
		};

		FAS_object_->d_S_a1[2][j]=init_for_S;
		FAS_object_->d_D_a1[2][j]=-inf;
		if(local_flag)
		{
			FAS_object_->d_I_a1[2][j]=-inf;
		}
		else
		{
			FAS_object_->d_I_a1[2][j]=-FAS_object_->d_open2-FAS_object_->d_epen2*j;
		};
	};

	FAS_object_->d_S_a1[0][0]=0;
	FAS_object_->d_S_a1[1][0]=0;
	FAS_object_->d_S_a1[2][0]=0;

	for(i=3;i<=FAS_object_->d_seq1_length;i++)
	{
		for(j=1;j<=FAS_object_->d_seq2_length;j++)
		{

			FAS_object_->d_I_a1[i][j]=FSA_utils::Tmax(FAS_object_->d_S_a1[i][j-1]-FAS_object_->d_open2-FAS_object_->d_epen2,FAS_object_->d_I_a1[i][j-1]-FAS_object_->d_epen2);
			
			if(i==2)
			{
				FAS_object_->d_S_a1[i][j]=init_for_S;

				FAS_object_->d_D_a1[i][j]=FAS_object_->d_D_a1[i-1][j];

				//FAS_object_->d_D_a1[i][j]=
				//	Tmax(FAS_object_->d_S_a1[i-2][j]-FAS_object_->d_open1-FAS_object_->d_epen1,FAS_object_->d_D_a1[i-2][j]-FAS_object_->d_epen1)-FAS_object_->d_gamma;

				continue;
			};
			if(i==3)
			{
				long int smart_score=FAS_object_->d_smatr[FAS_object_->convert_codon_into_AA(seq1_+i-3)][seq2_[j-1]];
				FAS_object_->d_S_a1[i][j]=FSA_utils::Tmax(

					//for local alignment
					init_for_S,
					
					FSA_utils::Tmax(FAS_object_->d_S_a1[i-3][j-1],FAS_object_->d_D_a1[i-3][j-1],FAS_object_->d_I_a1[i-3][j-1])+smart_score,
					FSA_utils::Tmax(FAS_object_->d_S_a1[i-2][j-1],FAS_object_->d_D_a1[i-2][j-1],FAS_object_->d_I_a1[i-2][j-1])+smart_score-FAS_object_->d_gamma
					
					);

				FAS_object_->d_D_a1[i][j]=FSA_utils::Tmax(
					
					FSA_utils::Tmax(FAS_object_->d_S_a1[i-3][j]-FAS_object_->d_open1-FAS_object_->d_epen1,FAS_object_->d_D_a1[i-3][j]-FAS_object_->d_epen1),
					FSA_utils::Tmax(FAS_object_->d_S_a1[i-2][j]-FAS_object_->d_open1-FAS_object_->d_epen1,FAS_object_->d_D_a1[i-2][j]-FAS_object_->d_epen1)-FAS_object_->d_gamma

					);
				continue;
				
			};
			if(i>=4)
			{

				long int smart_score=FAS_object_->d_smatr[FAS_object_->convert_codon_into_AA(seq1_+i-3)][seq2_[j-1]];
				FAS_object_->d_S_a1[i][j]=FSA_utils::Tmax(
					
					//for local alignment
					init_for_S,
					
					FSA_utils::Tmax(FAS_object_->d_S_a1[i-3][j-1],FAS_object_->d_D_a1[i-3][j-1],FAS_object_->d_I_a1[i-3][j-1])+smart_score,
					
					FSA_utils::Tmax(
					FSA_utils::Tmax(FAS_object_->d_S_a1[i-2][j-1],FAS_object_->d_D_a1[i-2][j-1],FAS_object_->d_I_a1[i-2][j-1]),
					FSA_utils::Tmax(FAS_object_->d_S_a1[i-4][j-1],FAS_object_->d_D_a1[i-4][j-1],FAS_object_->d_I_a1[i-4][j-1])
					)+smart_score-FAS_object_->d_gamma
					
					);

				FAS_object_->d_D_a1[i][j]=FSA_utils::Tmax(
					
					FSA_utils::Tmax(FAS_object_->d_S_a1[i-3][j]-FAS_object_->d_open1-FAS_object_->d_epen1,FAS_object_->d_D_a1[i-3][j]-FAS_object_->d_epen1),
					FSA_utils::Tmax(FAS_object_->d_S_a1[i-2][j]-FAS_object_->d_open1-FAS_object_->d_epen1,FAS_object_->d_D_a1[i-2][j]-FAS_object_->d_epen1)-FAS_object_->d_gamma,
					FSA_utils::Tmax(FAS_object_->d_S_a1[i-4][j]-FAS_object_->d_open1-FAS_object_->d_epen1,FAS_object_->d_D_a1[i-4][j]-FAS_object_->d_epen1)-FAS_object_->d_gamma

					);

			};



		};
	};

	long int score=FAS_object_->d_S_a1[0][0];
	for(i=0;i<=FAS_object_->d_seq1_length;i++)
	{
		for(j=0;j<=FAS_object_->d_seq2_length;j++)
		{
			score=FSA_utils::Tmax(score,
				FSA_utils::Tmax(FAS_object_->d_S_a1[i][j],FAS_object_->d_D_a1[i][j],FAS_object_->d_I_a1[i][j])
				);
		};
	};

	struct_for_alignment &tmp_obj=*((struct_for_alignment*)par_);
	tmp_obj.d_score=score;

	bool test_flag=false;
	if(test_flag)
	{
		FSA_utils::print_matrix(
		"d_S_a1.out",
		FAS_object_->d_seq1_length+1,
		FAS_object_->d_seq2_length+1,
		FAS_object_->d_S_a1);

		FSA_utils::print_matrix(
		"d_I_a1.out",
		FAS_object_->d_seq1_length+1,
		FAS_object_->d_seq2_length+1,
		FAS_object_->d_I_a1);

		FSA_utils::print_matrix(
		"d_D_a1.out",
		FAS_object_->d_seq1_length+1,
		FAS_object_->d_seq2_length+1,
		FAS_object_->d_D_a1);


	};
}

void FSA::DNA_AA_global(
string distr_file_name_,
alignment_function_type *alignment_function_)
{
	bool test_out=false;

	ofstream ftest;
	if(test_out)
	{
		string test_out_st="sss.out";
		ftest.open(test_out_st.data());
		if(!ftest)
		{
			throw error("Error - the file "+test_out_st+" is not found\n",1);
		};

		
	};

	ofstream f(distr_file_name_.data());
	if(!f)
	{
		throw error("Error - the file "+distr_file_name_+" is not found\n",3);
	};


	array_v<double> distr(this);


	long int s;
	long int score_max=small_long,score_min=-small_long;

	for(s=0;s<d_seq_number;s++)
	{

		if(test_out)
		{
			ftest<<">seq"<<s<<"\n";
		};

		long int i,j;
		//generate sequence #1
		for(i=0;i<d_seq1_length;i++)
		{
			d_seq1[i]=random_AA1();
			if(test_out)
			{
				ftest<<this->d_alphabet1[d_seq1[i]];
			};
		};
		if(test_out)
		{
			ftest<<"\n";
		};

		//generate sequence #2
		for(j=0;j<d_seq2_length;j++)
		{
			d_seq2[j]=random_AA2();
			if(test_out)
			{
				ftest<<this->d_alphabet2[d_seq2[j]];
			};

		};
		if(test_out)
		{
			ftest<<"\n";
		};

		if(d_reversed_seq)
		{
			FSA_utils::reverse_sequence(//reverse the letters of the sequence
			d_seq1,
			d_seq1_length);

			FSA_utils::reverse_sequence(//reverse the letters of the sequence
			d_seq2,
			d_seq2_length);
		};

		struct_for_alignment tmp_str;

		alignment_function_(
		d_seq1,//sequence #1
		d_seq2,//sequence #2
		this,
		(void*)(&tmp_str));

		long int score=tmp_str.d_score;
		if(s==0)
		{
			score_max=score;
			score_min=score;
		}
		else
		{
			score_max=FSA_utils::Tmax(score_max,score);
			score_min=FSA_utils::Tmin(score_min,score);
		};

		distr.increase_elem_by_1(score);

		if(test_out)
		{
			cout<<score<<endl;
		};

		if((s+1)%1000==0)
		{
			cout<<s+1<<endl;
		};
	};

	f<<score_max+1<<endl;
	for(s=0;s<=score_max;s++)
	{
		f<<s<<"\t"<<distr.d_elem[s-distr.d_ind0]<<endl;
	};

	f.close();
}

void FSA::AA_AA_global(
string distr_file_name_,
alignment_function_type *alignment_function_)
{
	ofstream f(distr_file_name_.data());
	if(!f)
	{
		throw error("Error - the file "+distr_file_name_+" is not found\n",3);
	};

	array_v<double> distr(this);


	long int s;
	long int score_max=small_long,score_min=-small_long;

	for(s=0;s<d_seq_number;s++)
	{
		long int i,j;
		//generate sequence #1

		for(i=0;i<d_seq1_length;i++)
		{
			d_seq1[i]=random_AA2();
		};
		//generate sequence #2
		for(j=0;j<d_seq2_length;j++)
		{
			d_seq2[j]=random_AA2();
		};

		if(d_reversed_seq)
		{
			FSA_utils::reverse_sequence(//reverse the letters of the sequence
			d_seq1,
			d_seq1_length);

			FSA_utils::reverse_sequence(//reverse the letters of the sequence
			d_seq2,
			d_seq2_length);
		};

		struct_for_alignment tmp_str;

		alignment_function_(
		d_seq1,//sequence #1
		d_seq2,//sequence #2
		this,
		(void*)(&tmp_str));

		long int score=tmp_str.d_score;
		if(s==0)
		{
			score_max=score;
			score_min=score;
		}
		else
		{
			score_max=FSA_utils::Tmax(score_max,score);
			score_min=FSA_utils::Tmin(score_min,score);
		};

		distr.increase_elem_by_1(score);

		if((s+1)%1000==0)
		{
			cout<<s+1<<endl;
		};
	};

	f<<score_max+1<<endl;
	for(s=0;s<=score_max;s++)
	{
		f<<s<<"\t"<<distr.d_elem[s-distr.d_ind0]<<endl;
	};

	f.close();
}

void FSA::DNA_to_3_frames_AA_global(
string distr_file_name_,
alignment_function_type *alignment_function_)
{
	bool frame_output_flag=false;
	long int s;
	double *freqs_from_frames=new double [d_number_of_letters2];
	for(s=0;s<d_number_of_letters2;s++)
	{
		freqs_from_frames[s]=0;
	};

	ofstream f(distr_file_name_.data());
	if(!f)
	{
		throw error("Error - the file "+distr_file_name_+" is not found\n",3);
	};

	array_v<double> distr(this);


	
	long int score_max=small_long,score_min=-small_long;

	for(s=0;s<d_seq_number;s++)
	{
		long int i,j;
		//generate sequence #1
		for(i=0;i<d_seq1_length;i++)
		{
			d_seq1[i]=random_AA1();
		};
		//generate sequence #2
		for(j=0;j<d_seq2_length;j++)
		{
			d_seq2[j]=random_AA2();
		};

		if(d_reversed_seq)
		{
			FSA_utils::reverse_sequence(//reverse the letters of the sequence
			d_seq1,
			d_seq1_length);

			FSA_utils::reverse_sequence(//reverse the letters of the sequence
			d_seq2,
			d_seq2_length);
		};

		struct_for_alignment tmp_str;

		long int score=0;

		long int v;
		for(v=0;v<d_codon_length;v++)
		{
			long int length1_tmp=this->d_seq1_length;

			long int length1_v_div_3=(long int)floor((double)(length1_tmp-v)/3.0);

			long int *seq1_AA=new long int[length1_v_div_3];

			long int j;
			for(j=0;j<length1_v_div_3;j++)
			{
				long int ind=v+j*3;
				seq1_AA[j]=convert_codon_into_AA(d_seq1+ind);
				freqs_from_frames[seq1_AA[j]]++;
			};


			this->d_seq1_length=length1_v_div_3;

			alignment_function_(
			seq1_AA,//sequence #1
			d_seq2,//sequence #2
			this,
			(void*)(&tmp_str));

			this->d_seq1_length=length1_tmp;

			if(v==0)
			{
				score=tmp_str.d_score;
			}
			else
			{
				score=FSA_utils::Tmax(score,tmp_str.d_score);
			};

			delete[]seq1_AA;
		};

		if(s==0)
		{
			score_max=score;
			score_min=score;
		}
		else
		{
			score_max=FSA_utils::Tmax(score_max,score);
			score_min=FSA_utils::Tmin(score_min,score);
		};

		distr.increase_elem_by_1(score);

		if((s+1)%1000==0)
		{
			cout<<s+1<<endl;
		};
	};

	f<<score_max+1<<endl;
	for(s=0;s<=score_max;s++)
	{
		f<<s<<"\t"<<distr.d_elem[s-distr.d_ind0]<<endl;
	};

	f.close();

	string frames_file_name="RR25_from_frames.in";
	ofstream freqf;
	
	if(frame_output_flag)
	{
		freqf.open(frames_file_name.data());
		if(!freqf)
		{
			throw error("Error - the file "+frames_file_name+" is not found\n",3);
		};
	};

	double sum=0;
	for(s=0;s<d_number_of_letters2;s++)
	{
		sum+=freqs_from_frames[s];
	};

	if(frame_output_flag)
	{
		freqf<<d_number_of_letters2<<endl;
		for(s=0;s<d_number_of_letters2;s++)
		{
			freqf<<freqs_from_frames[s]/sum<<endl;
		};

		freqf.close();
	};

	delete[]freqs_from_frames;

	
}


//general important sampling

IS1_general::IS1_general(
long int alphabet_letters_number1_,//number of letters in the sequence #1
long int alphabet_letters_number2_,//number of letters in the sequence #2
double *RR1_,//background probability for the sequence #1
double *RR2_,//background probability for the sequence #2
long int number_of_states_,//number of states
double **transition_probabilities_,//transition probabilities between states; matrix d_number_of_states x d_number_of_states
pair<long int, long int> *states_description_,//description of the states; the index is a state number
double ***states_distr_)//distributions of the states; the index is a state number
{
	d_alphabet_letters_number1=alphabet_letters_number1_;
	d_alphabet_letters_number2=alphabet_letters_number2_;
	d_RR1=RR1_;
	d_RR2=RR2_;
	d_number_of_states=number_of_states_;
	d_states_description=states_description_;
	d_states_distr=states_distr_;

	if(d_number_of_states<=0)
	{
		throw error("Error - incorrect parameter number_of_states_ in IS1_general::IS1_general\n",1);
	};

	d_states_distr_dimensions=new pair<long int, long int>[d_number_of_states];
	FSA_utils::assert_mem(d_states_distr_dimensions);

	d_states_sum_distr=new double*[d_number_of_states];
	FSA_utils::assert_mem(d_states_sum_distr);

	long int max_dim=d_number_of_states;

	d_max_dim1=-1;
	d_max_dim2=-1;

	long int s;
	for(s=0;s<d_number_of_states;s++)
	{
		d_max_dim1=FSA_utils::Tmax(d_max_dim1,d_states_description[s].first);
		d_max_dim2=FSA_utils::Tmax(d_max_dim2,d_states_description[s].second);

		long int dim1=FSA_utils::power_long(alphabet_letters_number1_,d_states_description[s].first);
		long int dim2=FSA_utils::power_long(alphabet_letters_number2_,d_states_description[s].second);


		d_states_distr_dimensions[s]=make_pair(dim1,dim2);

		max_dim=FSA_utils::Tmax(max_dim,dim1*dim2);

		IS1_general::generate_one_state_sum_distr(//calculates sum-distributions corresponded to the input parameters
		d_states_distr_dimensions[s],//dimensions of the matrix state_distr_
		d_states_distr[s],//state distribution in the same format as in d_states_distr[s]
		d_states_sum_distr[s]);//the result; the dimention is state_distr__dim1 x state_distr__dim2

	};


	d_states_sum_distr_elements_for_all_states=new long int[max_dim];
	FSA_utils::assert_mem(d_states_sum_distr_elements_for_all_states);

	

	long int i;
	for(i=0;i<max_dim;i++)
	{
		d_states_sum_distr_elements_for_all_states[i]=i;
	};

	FSA_utils::get_memory_for_matrix(d_number_of_states,d_number_of_states,d_transition_probabilities_sum);
	FSA_utils::get_memory_for_matrix(d_number_of_states,d_number_of_states,d_transition_probabilities);

	double eps=1e-5;

	for(s=0;s<d_number_of_states;s++)
	{
		double sum_tmp=0;
		long int s2;
		for(s2=0;s2<d_number_of_states;s2++)
		{
			d_transition_probabilities[s][s2]=transition_probabilities_[s][s2];
			d_transition_probabilities_sum[s][s2]=transition_probabilities_[s][s2];
			sum_tmp+=transition_probabilities_[s][s2];
		};
		if(fabs(sum_tmp-1.0)>eps)
		{
			throw error("Unexpected error in the parameter transition_probabilities_ of the function IS1_general::IS1_general\n",1);
		};

		for(s2=0;s2<d_number_of_states;s2++)
		{
			d_transition_probabilities[s][s2]/=sum_tmp;
		};

		FSA_utils::convert_distr_into_sum(//calculates and allocates sum-distribution
		d_number_of_states,//dimension
		d_transition_probabilities_sum[s]);//the result; the dimention is state_distr__dim1 x state_distr__dim2

	};


	d_tmp_letters=new pair<long int*, long int*>[d_number_of_states];
	FSA_utils::assert_mem(d_tmp_letters);

	for(s=0;s<d_number_of_states;s++)
	{
		d_tmp_letters[s].first=new long int[d_states_description[s].first];
		FSA_utils::assert_mem(d_tmp_letters[s].first);
		d_tmp_letters[s].second=new long int[d_states_description[s].second];
		FSA_utils::assert_mem(d_tmp_letters[s].second);
	};

	allocate_states_distr_sums();

}

IS1_general::~IS1_general()
{
	deallocate_states_distr_sums();

	delete[]d_states_distr_dimensions;

	long int s;
	for(s=0;s<d_number_of_states;s++)
	{
		delete[]d_states_sum_distr[s];
	};

	delete[]d_states_sum_distr_elements_for_all_states;
	
	FSA_utils::delete_memory_for_matrix(d_number_of_states,d_transition_probabilities_sum);
	FSA_utils::delete_memory_for_matrix(d_number_of_states,d_transition_probabilities);


	for(s=0;s<d_number_of_states;s++)
	{
		delete[]d_tmp_letters[s].first;
		delete[]d_tmp_letters[s].second;
	};
	delete[]d_tmp_letters;

}

void IS1_general::allocate_states_distr_sums()//allocates d_states_distr_sums
{
	long int*letters1=new long int[d_max_dim1];
	FSA_utils::assert_mem(letters1);

	long int*letters2=new long int[d_max_dim2];
	FSA_utils::assert_mem(letters2);


	long int s,x1,x2;

	d_states_distr_sums=new double ****[d_number_of_states];
	FSA_utils::assert_mem(d_states_distr_sums);

	for(s=0;s<d_number_of_states;s++)
	{
		d_states_distr_sums[s]=new double ***[d_states_description[s].first+1];
		FSA_utils::assert_mem(d_states_distr_sums[s]);

		long int dim1;
		long int dim2;

		for(x1=0;x1<=d_states_description[s].first;x1++)
		{
			dim1=FSA_utils::power_long(d_alphabet_letters_number1,x1);

			d_states_distr_sums[s][x1]=new double **[d_states_description[s].second+1];
			FSA_utils::assert_mem(d_states_distr_sums[s][x1]);

			for(x2=0;x2<=d_states_description[s].second;x2++)
			{
				dim2=FSA_utils::power_long(d_alphabet_letters_number2,x2);

				FSA_utils::get_memory_for_matrix(
				dim1,
				dim2,
				d_states_distr_sums[s][x1][x2]);

				long int i1,i2;
				for(i1=0;i1<dim1;i1++)
				{

					long int v;

					code_to_letters(//returns a unique code for the array of letters
					i1,//input code
					d_alphabet_letters_number1,//total number of letters
					x1,//dimension of the array with letters
					letters1);//array of letters; the result

					for(v=x1;v<d_states_description[s].first;v++)
					{
						letters1[v]=0;
					};

					long int code1=letters_to_code(//returns a unique code for the array of letters
					d_alphabet_letters_number1,//total number of letters
					d_states_description[s].first,//dimension of the array with letters
					letters1);//array of letters

					long int len1=FSA_utils::power_long(d_alphabet_letters_number1,d_states_description[s].first-x1)-1;

					double mult1=1.0;
					for(v=0;v<x1;v++)
					{
						mult1*=d_RR1[letters1[v]];
					};


					for(i2=0;i2<dim2;i2++)
					{

						code_to_letters(//returns a unique code for the array of letters
						i2,//input code
						d_alphabet_letters_number2,//total number of letters
						x2,//dimension of the array with letters
						letters2);//array of letters; the result

						for(v=x2;v<d_states_description[s].second;v++)
						{
							letters2[v]=0;
						};

						long int code2=letters_to_code(//returns a unique code for the array of letters
						d_alphabet_letters_number2,//total number of letters
						d_states_description[s].second,//dimension of the array with letters
						letters2);//array of letters

						long int len2=FSA_utils::power_long(d_alphabet_letters_number2,d_states_description[s].second-x2)-1;

						double mult2=1.0;
						for(v=0;v<x2;v++)
						{
							mult2*=d_RR2[letters2[v]];
						};

						double sum_tmp=0;
						long int v1,v2;
						for(v1=code1;v1<=code1+len1;v1++)
						{
							for(v2=code2;v2<=code2+len2;v2++)
							{
								sum_tmp+=d_states_distr[s][v1][v2];
							};
						};

						if(sum_tmp!=0&&(mult1==0||mult2==0))
						{
							throw error("The parameters d_states_distr and d_RR1 and d_RR2 are contradictory in IS1_general::allocate_states_distr_sums()\n",1);
						};

						if(mult1==0||mult2==0)
						{
							d_states_distr_sums[s][x1][x2][i1][i2]=0;
						}
						else
						{
							d_states_distr_sums[s][x1][x2][i1][i2]=sum_tmp/(mult1*mult2);
						};





					};
				};
			};
		};
	};

	delete[]letters1;
	delete[]letters2;

	//calculate matrices for infinite sums

	
	{
		long int code1=0;
		long int code2=0;
		long int x1=0;
		long int x2=0;

		calculate_inverse_matrices_for_the_infinite_sums(
			code1,
			code2,
			&d_A1_inv,
			&d_A2_inv,
			x1,
			x2);
	};

}

void IS1_general::deallocate_states_distr_sums()//deallocates d_states_distr_sums
{
	long int s,x1,x2;

	for(s=0;s<d_number_of_states;s++)
	{

		long int dim1;

		for(x1=0;x1<=d_states_description[s].first;x1++)
		{
			dim1=FSA_utils::power_long(d_alphabet_letters_number1,x1);


			for(x2=0;x2<=d_states_description[s].second;x2++)
			{

				FSA_utils::delete_memory_for_matrix(
				dim1,
				d_states_distr_sums[s][x1][x2]);
			};

			delete[]d_states_distr_sums[s][x1];
		};

		delete[] d_states_distr_sums[s];
	};

	delete[]d_states_distr_sums;
}


long int IS1_general::letters_to_code(//returns a unique code for the array of letters
long int letters_number_,//total number of letters
long int letters_dim_,//dimension of the array with letters
long int* letters_)//array of letters
{
	if(test_mode)
	{
		if(letters_dim_<=0)
		{
			if(letters_dim_<0)
			{
				throw error("Unexpected error\n",1);
			};
			return 0;
		};

		if(letters_number_<=0)
		{
			throw error("Unexpected error\n",1);
		};

		long int i;
		for(i=0;i<letters_dim_;i++)
		{
			if(letters_[i]<0||letters_[i]>=letters_number_)
			{
				throw error("Unexpected error\n",1);
			};
		};
	};

	long int res=letters_[0];
	long int i;
	for(i=1;i<letters_dim_;i++)
	{
		res=res*letters_number_+letters_[i];
	};

	return res;
}

void IS1_general::code_to_letters(//returns a unique code for the array of letters
long int code_,//input code
long int letters_number_,//total number of letters
long int letters_dim_,//dimension of the array with letters
long int*letters_)//array of letters; the result
{
	if(test_mode)
	{
		if(letters_dim_<=0)
		{
			if(letters_dim_<0)
			{
				throw error("Unexpected error\n",1);
			};
			return;
		};

		if(code_<0)
		{
			throw error("Unexpected error\n",1);
		};

		if(letters_number_<=0)
		{
			throw error("Unexpected error\n",1);
		};

	};

	long int i;
	for(i=letters_dim_-1;i>=0;i--)
	{
		letters_[i]=code_%letters_number_;
		code_=(code_-letters_[i])/letters_number_;
	};

	if(test_mode)
	{
		if(code_!=0)
		{
			throw error("Unexpected error\n",1);
		};
	};

}

long int IS1_general::matr_indexes_to_code(
long int code1_,//code #1
long int code1_number_,//the range of code1_ is [0,code1_number_-1]
long int code2_,//code #2
long int code2_number_)//the range of code2_ is [0,code2_number_-1]
{
	if(test_mode)
	{
		if(code1_number_<=0||code2_number_<=0)
		{
			throw error("Unexpected error\n",1);
		};

		if(code2_<0||code2_>=code2_number_||code1_<0||code1_>=code1_number_)
		{
			throw error("Unexpected error\n",1);
		};
	};

	return code1_*code2_number_+code2_;
}

void IS1_general::code_to_matr_indexes(
long int code_,//input code
long int &code1_,//code #1; the result
long int code1_number_,//the range of code1_ is [0,code1_number_-1]
long int &code2_,//code #2; the result
long int code2_number_)//the range of code2_ is [0,code2_number_-1]
{
	if(test_mode)
	{
		if(code1_number_<=0||code2_number_<=0)
		{
			throw error("Unexpected error\n",1);
		};
	};

	code2_=code_%code2_number_;
	code1_=(code_-code2_)/code2_number_;

	if(test_mode)
	{

		if(code2_<0||code2_>=code2_number_||code1_<0||code1_>=code1_number_)
		{
			throw error("Unexpected error\n",1);
		};
	};
}

double ** IS1_general::allocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
long int alphabet_letters_number1_,//number of letters in the sequence #1
long int alphabet_letters_number2_,//number of letters in the sequence #2
pair<long int, long int> state_description_)//state description
{
	double ** res=NULL;

	long int dim1=FSA_utils::power_long(alphabet_letters_number1_,state_description_.first);
	long int dim2=FSA_utils::power_long(alphabet_letters_number2_,state_description_.second);

	FSA_utils::get_memory_for_matrix(dim1,dim2,res);

	return res;
}

void IS1_general::deallocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
double **&state_distr_,
long int alphabet_letters_number1_,//number of letters in the sequence #1
pair<long int, long int> state_description_)//state description
{
	long int dim1=FSA_utils::power_long(alphabet_letters_number1_,state_description_.first);

	FSA_utils::delete_memory_for_matrix(dim1,state_distr_);
	state_distr_=NULL;
}


void IS1_general::generate_one_state_sum_distr(//calculates sum-distributions corresponded to the input parameters
pair<long int, long int> state_distr_dims_,//dimensions of the matrix state_distr_
double **state_distr_,//state distribution in the same format as in d_states_distr[s]
double *&state_sum_distr_)//the result; the dimention is state_distr__dim1 x state_distr__dim2
{
	long int dim1=state_distr_dims_.first;
	long int dim2=state_distr_dims_.second;

	long int dim=dim1*dim2;

	state_sum_distr_=new double [dim];
	FSA_utils::assert_mem(state_sum_distr_);

	long int i,j;
	for(i=0;i<dim1;i++)
	{
		for(j=0;j<dim2;j++)
		{
			long int code=matr_indexes_to_code(
			i,//code #1
			dim1,//the range of code1_ is [0,code1_number_-1]
			j,//code #2
			dim2);//the range of code2_ is [0,code2_number_-1]

			state_sum_distr_[code]=state_distr_[i][j];

			if(test_mode)
			{
				if(code<0||code>=dim)
				{
					throw error("Error - code<0||code>=dim\n",1);
				};
			};

		};
	};
	
	FSA_utils::convert_distr_into_sum(//convert distr_[0], distr_[1], distr_[2],... into distr_[0], distr_[0]+distr_[1], distr_[0]+distr_[1]+distr_[2],...
	dim,
	state_sum_distr_);
}

void IS1_general::generate_random_letters_in_a_given_state(
long int state_number_,//state number
long int *letters1_,//the resulted letters for the sequence 1; the array must be allocated
long int *letters2_)//the resulted letters for the sequence 2; the array must be allocated
{
	long int dim1=d_states_distr_dimensions[state_number_].first;
	long int dim2=d_states_distr_dimensions[state_number_].second;

	long int dim=dim1*dim2;

	long int code=FSA_utils::random_long(
		FSA_utils::ran2(),
		dim,
		d_states_sum_distr[state_number_],
		d_states_sum_distr_elements_for_all_states);

		long int code1;
		long int code2;

		code_to_matr_indexes(
		code,
		code1,
		dim1,
		code2,
		dim2);

		code_to_letters(//returns a unique code for the array of letters
		code1,//input code
		//dim1,//total number of letters
		d_alphabet_letters_number1,
		d_states_description[state_number_].first,//dimension of the array with letters
		letters1_);//array of letters; the result

		code_to_letters(//returns a unique code for the array of letters
		code2,//input code
		//dim2,//total number of letters
		d_alphabet_letters_number2,
		d_states_description[state_number_].second,//dimension of the array with letters
		letters2_);//array of letters; the result


}

void IS1_general::one_step_of_the_importance_samping(//an initial state when sequences are empty must be defined outside
long int current_state_,//the current state
long int *seq1_add_letters_,//letters for the sequence #1; the array must be allocated
long int *seq2_add_letters_,//letters for the sequence #2; the array must be allocated
long int &new_state_)//a new state
{
	//generate a new state
	new_state_=FSA_utils::random_long(
		FSA_utils::ran2(),
		d_number_of_states,
		d_transition_probabilities_sum[current_state_],
		d_states_sum_distr_elements_for_all_states);

	generate_random_letters_in_a_given_state(
		new_state_,//state number
		seq1_add_letters_,//the resulted letters for the sequence 1; the array must be allocated
		seq2_add_letters_);//the resulted letters for the sequence 2; the array must be allocated


}

//-------------------------------------------------------------------
//the importance sampling simulation object

IS1_general_simulation::IS1_general_simulation(
IS1_general *IS1_general_obj_,
long int initial_state_,//initial state for the IS
long int max_seq1_length_,//maximum sequence length
long int max_seq2_length_,//maximum sequence length
double *initial_distr_)//initial distribution of states; initial_state_ is ignored if d_initial_distr_ is defined
{
	d_IS1_general_obj=IS1_general_obj_;
	d_max_seq1_length=max_seq1_length_;
	d_max_seq2_length=max_seq2_length_;

	d_seq1=new long int [d_max_seq1_length];
	FSA_utils::assert_mem(d_seq1);

	d_seq2=new long int [d_max_seq2_length];
	FSA_utils::assert_mem(d_seq2);

	d_seq1_current_length=0;
	d_seq2_current_length=0;

	d_initial_state=initial_state_;
	d_current_state=initial_state_;

	d_initial_distr_sum=NULL;
	if(initial_distr_)
	{
		double eps_tmp=1e-10;
		double sum_tmp=0;
		d_initial_distr_sum=new double[IS1_general_obj_->d_number_of_states];
		FSA_utils::assert_mem(d_initial_distr_sum);
		long int i;
		for(i=0;i<IS1_general_obj_->d_number_of_states;i++)
		{
			d_initial_distr_sum[i]=initial_distr_[i];
			sum_tmp+=initial_distr_[i];
		};
		if(fabs(sum_tmp-1.0)>=eps_tmp)
		{
			throw error("Unexpected error - the sum of probablities in initial_distr_ is not 1.0 in IS1_general_simulation::IS1_general_simulation\n",1);
		};

		FSA_utils::convert_distr_into_sum(//convert distr_[0], distr_[1], distr_[2],... into distr_[0], distr_[0]+distr_[1], distr_[0]+distr_[1]+distr_[2],...
			IS1_general_obj_->d_number_of_states,
			d_initial_distr_sum);

		d_current_state=FSA_utils::random_long<long int>(
		FSA_utils::ran2(),
		d_IS1_general_obj->d_number_of_states,
		d_initial_distr_sum);


	};


	//weights calculation
	long int depth1=-1;
	long int depth2=-1;

	long int s;
	for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
	{
		depth1=FSA_utils::Tmax(depth1,d_IS1_general_obj->d_states_description[s].first);
		depth2=FSA_utils::Tmax(depth2,d_IS1_general_obj->d_states_description[s].second);
	};

	d_W1_step1=depth1+1;
	d_W1_step2=depth2+1;

	depth1+=d_W1_step1;
	depth2+=d_W1_step2;

	d_W1=new two_dim_layer<double>* [d_IS1_general_obj->d_number_of_states];
	FSA_utils::assert_mem(d_W1);

	for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
	{
		d_W1[s]=new two_dim_layer<double> (
			max_seq1_length_,
			max_seq2_length_,
			depth1,
			depth2,
			0);

	};

}


IS1_general_simulation::IS1_general_simulation(//creates a copy of IS1_general_simulation_ with smallest possible memory allocation
IS1_general_simulation *IS1_general_simulation_)//maximum sequence length
{
	d_IS1_general_obj=IS1_general_simulation_->d_IS1_general_obj;
	d_max_seq1_length=IS1_general_simulation_->d_W1_seq1_current_length;
	d_max_seq2_length=IS1_general_simulation_->d_W1_seq2_current_length;

	d_seq1_current_length=IS1_general_simulation_->d_seq1_current_length;
	d_seq2_current_length=IS1_general_simulation_->d_seq2_current_length;

	d_seq1=new long int [d_seq1_current_length];
	FSA_utils::assert_mem(d_seq1);

	d_seq2=new long int [d_seq2_current_length];
	FSA_utils::assert_mem(d_seq2);


	d_initial_state=IS1_general_simulation_->d_initial_state;
	d_current_state=IS1_general_simulation_->d_current_state;

	d_initial_distr_sum=NULL;
	if(IS1_general_simulation_->d_initial_distr_sum)
	{
		d_initial_distr_sum=new double[IS1_general_simulation_->d_IS1_general_obj->d_number_of_states];
		FSA_utils::assert_mem(d_initial_distr_sum);
		long int i;
		for(i=0;i<IS1_general_simulation_->d_IS1_general_obj->d_number_of_states;i++)
		{
			d_initial_distr_sum[i]=IS1_general_simulation_->d_initial_distr_sum[i];
		};
	};

	//weights calculation
	long int depth1=-1;
	long int depth2=-1;

	long int s;
	for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
	{
		depth1=FSA_utils::Tmax(depth1,d_IS1_general_obj->d_states_description[s].first);
		depth2=FSA_utils::Tmax(depth2,d_IS1_general_obj->d_states_description[s].second);
	};

	d_W1_step1=depth1+1;
	d_W1_step2=depth2+1;

	depth1+=d_W1_step1;
	depth2+=d_W1_step2;

	d_W1=new two_dim_layer<double>* [d_IS1_general_obj->d_number_of_states];
	FSA_utils::assert_mem(d_W1);

	d_W1_seq1_current_length=IS1_general_simulation_->d_W1_seq1_current_length;;
	d_W1_seq2_current_length=IS1_general_simulation_->d_W1_seq2_current_length;;


	for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
	{
		d_W1[s]=new two_dim_layer<double> (//constructor; copies data from two_dim_layer_ with minimally allocated memory
		d_W1_seq1_current_length,
		d_W1_seq2_current_length,
		IS1_general_simulation_->d_W1[s]);
	};


	//copy sequences
	for(s=0;s<d_seq1_current_length;s++)
	{
		d_seq1[s]=IS1_general_simulation_->d_seq1[s];
	};

	for(s=0;s<d_seq2_current_length;s++)
	{
		d_seq2[s]=IS1_general_simulation_->d_seq2[s];
	};


}


IS1_general_simulation::~IS1_general_simulation()
{
	delete[]d_seq1;
	delete[]d_seq2;
	delete[]d_initial_distr_sum;

	//weights calculation
	long int s;
	for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
	{
		delete d_W1[s];
	};
}

void IS1_general_simulation::init()//initialization for a new sequence generation
{
	d_seq1_current_length=0;
	d_seq2_current_length=0;

	d_current_state=d_initial_state;

	if(d_initial_distr_sum)
	{
		d_current_state=FSA_utils::random_long<long int>(
		FSA_utils::ran2(),
		d_IS1_general_obj->d_number_of_states,
		d_initial_distr_sum);
	};

	//weights calculation
	d_W1_seq1_current_length=-1;
	d_W1_seq2_current_length=-1;

}

void IS1_general_simulation::simulate_upto_target_lengths(
long int target_seq1_length_,//target length of the sequence #1
long int target_seq2_length_)//target length of the sequence #2
{
	
	long int new_state;

	while((d_seq1_current_length<target_seq1_length_)||(d_seq2_current_length<target_seq2_length_))
	{
		long int expected_len1=d_seq1_current_length+d_IS1_general_obj->d_max_dim1;
		long int expected_len2=d_seq2_current_length+d_IS1_general_obj->d_max_dim2;


		if(expected_len1>d_max_seq1_length||expected_len2>d_max_seq2_length)
		{
			throw error("Error - please increase maximum allowed sequence length in IS1_general_simulation::simulate_upto_target_lengths\n",1);
		};

		d_IS1_general_obj->one_step_of_the_importance_samping(//an initial state when sequences are empty must be defined outside
		d_current_state,//the current state
		d_seq1+d_seq1_current_length,//letters for the sequence #1; the array must be allocated
		d_seq2+d_seq2_current_length,//letters for the sequence #2; the array must be allocated
		new_state);//a new state


		d_current_state=new_state;

		d_seq1_current_length+=d_IS1_general_obj->d_states_description[d_current_state].first;
		d_seq2_current_length+=d_IS1_general_obj->d_states_description[d_current_state].second;

	};
	
}

//weights calculation
void IS1_general_simulation::calculate_weight_W1_one_step(//calculates the weight assuming all the weights used in the recursive equations are defined
long int target_length1_,//target length for sequence #1
long int target_length2_,//target length for sequence #2
long int i1_,//the first index (corresponded to the sequence #1)
long int i2_)//the second index (corresponded to the sequence #2)
{
	if(i1_>d_max_seq1_length||i2_>d_max_seq2_length)
	{
		throw error("Unexpected error - i1_>d_max_seq1_length||i2_>d_max_seq2_length in IS1_general_simulation::calculate_weight_W1_one_step\n",1);
	};

	if(i1_==0&&i2_==0)
	{

		long int s;
		for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
		{
			d_W1[s]->set_element(i1_,i2_,0);
		};

		d_W1[d_initial_state]->set_element(i1_,i2_,1.0);

		if(d_initial_distr_sum)
		{
			s=0;
			d_W1[s]->set_element(i1_,i2_,d_initial_distr_sum[s]);

			for(s=1;s<d_IS1_general_obj->d_number_of_states;s++)
			{
				d_W1[s]->set_element(i1_,i2_,d_initial_distr_sum[s]-d_initial_distr_sum[s-1]);
			};
		};

		return;
	};
	

	long int s;
	for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
	{

		long int dim1_tmp=d_IS1_general_obj->d_states_description[s].first;
		long int diff1=-dim1_tmp+i1_;

		long int dim2_tmp=d_IS1_general_obj->d_states_description[s].second;
		long int diff2=-dim2_tmp+i2_;

		long int code1=-1;
		long int code2=-1;


		double w_tmp=0;

		if(diff1>=0&&diff2>=0)
		{

			
			long int x1;
			long int tmp_d=i1_-target_length1_;
			if(-tmp_d>=0)
			{
				x1=dim1_tmp;
			}
			else
			{
				if(tmp_d>=dim1_tmp)
				{
					x1=0;
				}
				else
				{
					x1=dim1_tmp-tmp_d;
				};
			};

			long int x2;
			tmp_d=i2_-target_length2_;
			if(-tmp_d>=0)
			{
				x2=dim2_tmp;
			}
			else
			{
				if(tmp_d>=dim2_tmp)
				{
					x2=0;
				}
				else
				{
					x2=dim2_tmp-tmp_d;
				};
			};

		code1=d_IS1_general_obj->letters_to_code(//returns a unique code for the array of letters
			d_IS1_general_obj->d_alphabet_letters_number1,//total number of letters
			x1,//dimension of the array with letters
			d_seq1+diff1);//array of letters

		code2=d_IS1_general_obj->letters_to_code(//returns a unique code for the array of letters
			d_IS1_general_obj->d_alphabet_letters_number2,//total number of letters
			x2,//dimension of the array with letters
			d_seq2+diff2);//array of letters


			double q2=d_IS1_general_obj->d_states_distr_sums[s][x1][x2][code1][code2];




			long int x;
			for(x=0;x<d_IS1_general_obj->d_number_of_states;x++)
			{
				w_tmp+=d_IS1_general_obj->d_transition_probabilities[x][s]*d_W1[x]->get_element(i1_-dim1_tmp,i2_-dim2_tmp);
			};
			w_tmp*=q2;
		};


		d_W1[s]->set_element(i1_,i2_,w_tmp);


	};

}

void IS1_general_simulation::calculate_weight_W1_upto_target_lengths(
long int target_seq1_length_,//target length of the sequence #1
long int target_seq2_length_,//target length of the sequence #2
long int seq1_length_tmp_,//the weights are calculated upto this length for the sequence #1
long int seq2_length_tmp_)//the weights are calculated upto this length for the sequence #2
{
	if(d_W1_seq1_current_length==-1&&d_W1_seq2_current_length==-1)
	{
		calculate_weight_W1_one_step(//calculates the weight assuming all the weights used in the recursive equations are defined
		target_seq1_length_,
		target_seq2_length_,
		0,//the first index (corresponded to the sequence #1)
		0);//the second index (corresponded to the sequence #2)

		d_W1_seq1_current_length=0;
		d_W1_seq2_current_length=0;
	};


	long int ind1,ind2;

	for(ind1=d_W1_seq1_current_length+1;ind1<=seq1_length_tmp_;ind1++)
	{


		long int s;
		for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
		{
			d_W1[s]->set_max_ind(ind1,d_W1_seq2_current_length);
		};

		long int i2;
		for(i2=0;i2<=d_W1_seq2_current_length;i2++)
		{
			calculate_weight_W1_one_step(//calculates the weight assuming all the weights used in the recursive equations are defined
			target_seq1_length_,
			target_seq2_length_,
			ind1,//the first index (corresponded to the sequence #1)
			i2);//the second index (corresponded to the sequence #2)
		};

	};

	d_W1_seq1_current_length=target_seq1_length_;

	for(ind2=d_W1_seq2_current_length+1;ind2<=seq2_length_tmp_;ind2++)
	{

		long int s;
		for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
		{
			d_W1[s]->set_max_ind(seq1_length_tmp_,ind2);
		};

		long int i1;
		for(i1=0;i1<=seq1_length_tmp_;i1++)
		{
			calculate_weight_W1_one_step(//calculates the weight assuming all the weights used in the recursive equations are defined
			target_seq1_length_,
			target_seq2_length_,
			i1,//the first index (corresponded to the sequence #1)
			ind2);//the second index (corresponded to the sequence #2)
		};

	};
	
	d_W1_seq2_current_length=target_seq2_length_;

}

void IS1_general_simulation::calculate_weight_W1_for_fixed_lengths_using_recursions(
long int target_seq1_length_,//target length of the sequence #1
long int target_seq2_length_,//target length of the sequence #2
double &weight_,//the resulted weight
bool save_matrices_,
std::vector<std::vector<double> > *A1_,
std::vector<std::vector<double> > *A2_)
{
	map<long int,bool> ind1_flag;
	map<long int,bool> ind2_flag;
	if(save_matrices_)
	{
		long int dimA1=0;
		long int dimA2=0;

		long int s;
		for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
		{
			long int dim1_tmp=d_IS1_general_obj->d_states_description[s].first;
			long int dim2_tmp=d_IS1_general_obj->d_states_description[s].second;

			if(dim1_tmp==0)
			{
				dimA2++;
				ind2_flag[s]=true;
			};

			if(dim2_tmp==0)
			{
				dimA1++;
				ind1_flag[s]=true;
			};

		};

		if(save_matrices_)
		{
			A1_->clear();
			A2_->clear();

			if(dimA1>0)
			{
				A1_->resize(dimA1,std::vector<double> (dimA1+1,1));
			};

			if(dimA2>0)
			{
				A2_->resize(dimA2,std::vector<double> (dimA2+1,1));
			};
		};

	};
	
	long int i1_ind_max=target_seq1_length_+d_IS1_general_obj->d_max_dim1-1;
	long int i2_ind_max=target_seq2_length_+d_IS1_general_obj->d_max_dim2-1;

	long int s;
	for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
	{
		d_W1[s]->set_max_ind(i1_ind_max,i2_ind_max);
	};

	long int seq1_length_tmp=target_seq1_length_+d_IS1_general_obj->d_max_dim1-1;
	long int seq2_length_tmp=target_seq2_length_+d_IS1_general_obj->d_max_dim2-1;

	calculate_weight_W1_upto_target_lengths(
		target_seq1_length_,//target length of the sequence #1
		target_seq2_length_,//target length of the sequence #2
		seq1_length_tmp,
		seq2_length_tmp);


	array_v<double> *sum_over_ind1=new array_v<double>[d_IS1_general_obj->d_number_of_states];
	array_v<double> *sum_over_ind2=new array_v<double>[d_IS1_general_obj->d_number_of_states];

	
	{
		
		//index #1
		long int s;

		long int i1;
		for(i1=-d_IS1_general_obj->d_max_dim1;i1<0;i1++)
		{
			for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
			{
				sum_over_ind2[s].set_elem(i1,0);
			};
		};


		
		for(i1=0;i1<=i1_ind_max;i1++)
		{

			if(save_matrices_)
			{
				if(i1!=target_seq1_length_)
				{
					continue;
				};
			};

			bool explicit_dependency_flag=true;

			loop1:;

			vector<double> x2_vect;

			long int ind_count=-1;
			
			

			for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
			{

				long int dim1_tmp=d_IS1_general_obj->d_states_description[s].first;
				long int dim2_tmp=d_IS1_general_obj->d_states_description[s].second;

				

				if(explicit_dependency_flag)
				{
					if(dim1_tmp==0)
					{
						continue;
					};
				}
				else
				{
					if(dim1_tmp!=0)
					{
						continue;
					};
					
				};

				double x2_vect_elem=0;


	//------------------------------------------
				long int diff1=-dim1_tmp+i1;


				double w_tmp=0;
				double q2=0;

				if(diff1>=0)
				{
					long int code1=-1;
					long int code2=0;

					
					long int x1;
					long int tmp_d=i1-target_seq1_length_;
					if(-tmp_d>=0)
					{
						x1=dim1_tmp;
					}
					else
					{
						if(tmp_d>=dim1_tmp)
						{
							x1=0;
						}
						else
						{
							x1=dim1_tmp-tmp_d;
						};
					};

					

					long int x2=0;

				code1=d_IS1_general_obj->letters_to_code(//returns a unique code for the array of letters
					d_IS1_general_obj->d_alphabet_letters_number1,//total number of letters
					x1,//dimension of the array with letters
					d_seq1+diff1);//array of letters


					q2=d_IS1_general_obj->d_states_distr_sums[s][x1][x2][code1][code2];
					


					if(save_matrices_)
					{
						if(i1==target_seq1_length_)
						{
							if(!explicit_dependency_flag)
							{
								if(dim1_tmp==0)
								{
									ind_count++;

									long int ind_count2=-1;
									long int s1;
									for(s1=0;s1<d_IS1_general_obj->d_number_of_states;s1++)
									{
										if(ind2_flag[s1])
										{
											ind_count2++;

											(*A2_)[ind_count][ind_count2]=d_IS1_general_obj->d_transition_probabilities[s1][s]*q2;

											if(s==s1)
											{
												(*A2_)[ind_count][ind_count2]-=1.0;
											};

										};
									};
								};
							};
						};
						continue;
					};


		//------------------------------------------


					long int s1;
					if(explicit_dependency_flag)
					{
						for(s1=0;s1<d_IS1_general_obj->d_number_of_states;s1++)
						{
							w_tmp+=d_IS1_general_obj->d_transition_probabilities[s1][s]*sum_over_ind2[s1].get_elem(i1-dim1_tmp);
						};

						w_tmp*=q2;
					}
					else
					{
						long int s1;
						for(s1=0;s1<d_IS1_general_obj->d_number_of_states;s1++)
						{
							if(!(d_IS1_general_obj->d_states_description[s1].first==0))
							{
								x2_vect_elem-=d_IS1_general_obj->d_transition_probabilities[s1][s]*sum_over_ind2[s1].get_elem(i1-dim1_tmp);
							};

						};

						x2_vect_elem*=q2;

					};




					
					

				};

				long int ii2;
				for(ii2=0;ii2<dim2_tmp;ii2++)
				{
					double ww=d_W1[s]->get_element(i1,target_seq2_length_+ii2);
					if(explicit_dependency_flag)
					{
						w_tmp+=ww;
					}
					else
					{
						x2_vect_elem-=ww;
					};
				};




				if(explicit_dependency_flag)
				{
					sum_over_ind2[s].set_elem(i1,w_tmp);
				}
				else
				{
					x2_vect.push_back(x2_vect_elem);
				};

			};

			if(explicit_dependency_flag)
			{
				explicit_dependency_flag=false;
				goto loop1;
			};

			if(!explicit_dependency_flag)
			{
				if(d_IS1_general_obj->d_A2_inv.size()>0)
				{
					
					std::vector<double> res2;
					
					FSA_utils::multiply_matrix_and_vector(
						d_IS1_general_obj->d_A2_inv,
						x2_vect,
						res2);
						


					long int dimA2=0;

					long int s;
					for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
					{

						long int dim1_t=d_IS1_general_obj->d_states_description[s].first;


						if(dim1_t==0)
						{
							sum_over_ind2[s].set_elem(i1,res2[dimA2]);
							dimA2++;
						};


					};

				};

			};

		};
	};


//============================================================================================

	{
		//index #2
		long int s;

		long int i2;
		for(i2=-d_IS1_general_obj->d_max_dim2;i2<0;i2++)
		{
			for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
			{
				sum_over_ind1[s].set_elem(i2,0);
			};
		};


		
		for(i2=0;i2<=i2_ind_max;i2++)
		{

			if(save_matrices_)
			{
				if(i2!=target_seq2_length_)
				{
					continue;
				};
			};

			bool explicit_dependency_flag=true;

			loop2:;

			vector<double> x1_vect;

			long int ind_count=-1;
			

			for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
			{


				long int dim1_tmp=d_IS1_general_obj->d_states_description[s].first;
				long int dim2_tmp=d_IS1_general_obj->d_states_description[s].second;

				if(explicit_dependency_flag)
				{
					if(dim2_tmp==0)
					{
						continue;
					};
				}
				else
				{
					if(dim2_tmp!=0)
					{
						continue;
					};
				};

				double x1_vect_elem=0;


	//------------------------------------------
				long int diff2=-dim2_tmp+i2;


				double w_tmp=0;

				double q2=0;

				if(diff2>=0)
				{
					long int code1=0;
					long int code2=-1;

					
					long int x2;
					long int tmp_d=i2-target_seq2_length_;
					if(-tmp_d>=0)
					{
						x2=dim2_tmp;
					}
					else
					{
						if(tmp_d>=dim2_tmp)
						{
							x2=0;
						}
						else
						{
							x2=dim2_tmp-tmp_d;
						};
					};

					

					long int x1=0;


				code2=d_IS1_general_obj->letters_to_code(//returns a unique code for the array of letters
					d_IS1_general_obj->d_alphabet_letters_number2,//total number of letters
					x2,//dimension of the array with letters
					d_seq2+diff2);//array of letters


					q2=d_IS1_general_obj->d_states_distr_sums[s][x1][x2][code1][code2];

					if(save_matrices_)
					{
						if(i2==target_seq2_length_)
						{
							if(!explicit_dependency_flag)
							{
								if(dim2_tmp==0)
								{
									ind_count++;

									long int ind_count2=-1;
									long int s1;
									for(s1=0;s1<d_IS1_general_obj->d_number_of_states;s1++)
									{
										if(ind1_flag[s1])
										{
											ind_count2++;

											(*A1_)[ind_count][ind_count2]=d_IS1_general_obj->d_transition_probabilities[s1][s]*q2;

											if(s==s1)
											{
												(*A1_)[ind_count][ind_count2]-=1.0;
											};

										};
									};
								};
							};
						};
						continue;
					};


		//------------------------------------------


					long int s1;
					if(explicit_dependency_flag)
					{
						for(s1=0;s1<d_IS1_general_obj->d_number_of_states;s1++)
						{
							w_tmp+=d_IS1_general_obj->d_transition_probabilities[s1][s]*sum_over_ind1[s1].get_elem(i2-dim2_tmp);
						};

						w_tmp*=q2;
					}
					else
					{
						long int s1;
						for(s1=0;s1<d_IS1_general_obj->d_number_of_states;s1++)
						{
							if(!(d_IS1_general_obj->d_states_description[s1].second==0))
							{
								x1_vect_elem-=d_IS1_general_obj->d_transition_probabilities[s1][s]*sum_over_ind1[s1].get_elem(i2-dim2_tmp);
							};

						};

						x1_vect_elem*=q2;

					};



				};

				long int ii1;
				for(ii1=0;ii1<dim1_tmp;ii1++)
				{
					double ww=d_W1[s]->get_element(target_seq1_length_+ii1,i2);
					if(explicit_dependency_flag)
					{
						w_tmp+=ww;
					}
					else
					{
						x1_vect_elem-=ww;
					};

				};

				if(explicit_dependency_flag)
				{
					sum_over_ind1[s].set_elem(i2,w_tmp);
				}
				if(!explicit_dependency_flag)
				{
					x1_vect.push_back(x1_vect_elem);
				};


			};

			if(explicit_dependency_flag)
			{
				explicit_dependency_flag=false;
				goto loop2;
			};

			if(!explicit_dependency_flag)
			{
				if(d_IS1_general_obj->d_A2_inv.size()>0)
				{
					
					std::vector<double> res1;
					
					FSA_utils::multiply_matrix_and_vector(
						d_IS1_general_obj->d_A1_inv,
						x1_vect,
						res1);
						


					long int dimA1=0;
					

					long int s;
					for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
					{
						
						long int dim2_t=d_IS1_general_obj->d_states_description[s].second;


						if(dim2_t==0)
						{
							sum_over_ind1[s].set_elem(i2,res1[dimA1]);
							
							dimA1++;
						};

					};

				};

			};



		};

	};

	weight_=0;

	for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
	{
		long int dim1_tmp=d_IS1_general_obj->d_states_description[s].first;
		long int dim2_tmp=d_IS1_general_obj->d_states_description[s].second;

		long int i1;
		for(i1=0;i1<dim1_tmp;i1++)
		{
			weight_+=sum_over_ind2[s].get_elem(target_seq1_length_+i1);
		};

		long int i2;
		for(i2=0;i2<dim2_tmp;i2++)
		{
			weight_+=sum_over_ind1[s].get_elem(target_seq2_length_+i2);
		};

		for(i1=0;i1<dim1_tmp;i1++)
		{
			for(i2=0;i2<dim2_tmp;i2++)
			{
				weight_-=d_W1[s]->get_element(target_seq1_length_+i1,target_seq2_length_+i2);
			};
		};

		
	};


	delete[]sum_over_ind1;
	delete[]sum_over_ind2;



}

//================================
void IS1_general::calculate_inverse_matrices_for_the_infinite_sums(
long int code1_,
long int code2_,
std::vector<std::vector<double> > *A1_inv_,
std::vector<std::vector<double> > *A2_inv_,
long int x1_,
long int x2_)
{
	std::vector<std::vector<double> > A1;
	std::vector<std::vector<double> > A2;

	map<long int,bool> ind1_flag;
	map<long int,bool> ind2_flag;

	{
		long int dimA1=0;
		long int dimA2=0;

		long int s;
		for(s=0;s<d_number_of_states;s++)
		{
			long int dim1_tmp=d_states_description[s].first;
			long int dim2_tmp=d_states_description[s].second;

			if(dim1_tmp==0)
			{
				dimA2++;
				ind2_flag[s]=true;
			};

			if(dim2_tmp==0)
			{
				dimA1++;
				ind1_flag[s]=true;
			};

		};

		
		

		if(x1_>=0)
		{
			A1.clear();
			(*A1_inv_).clear();
			if(dimA1>0)
			{
				A1.resize(dimA1,std::vector<double> (dimA1+1,1));
			};
			
		};


		if(x2_>=0)
		{
			A2.clear();
			(*A2_inv_).clear();
			if(dimA2>0)
			{
				A2.resize(dimA2,std::vector<double> (dimA2+1,1));
			};
		};


	};
	

	if(x1_>=0)
	{
		
		//index #1
		long int s;
		long int ind_count=-1;

		for(s=0;s<d_number_of_states;s++)
		{
			long int dim1_tmp=d_states_description[s].first;

			if(dim1_tmp!=0)
			{
				continue;
			};


			double q2=0;

			long int code1=code1_;
			long int code2=0;

			
			long int x1=x1_;
			long int x2=0;

			q2=d_states_distr_sums[s][x1][x2][code1][code2];
			


			ind_count++;

			long int ind_count2=-1;
			long int s1;
			for(s1=0;s1<d_number_of_states;s1++)
			{
				if(ind2_flag[s1])
				{
					ind_count2++;

					A2[ind_count][ind_count2]=d_transition_probabilities[s1][s]*q2;

					if(s==s1)
					{
						A2[ind_count][ind_count2]-=1.0;
					};

				};
			};

		};




		if(A2.size()>0)
		{
			double inside_eps=1e-12;
			std::vector<double> x2;
			FSA_utils::Gauss(
				A2,//matrix n*(n+1)
				x2,//solution
				inside_eps,
				A2_inv_);
		};

		

	};

//============================================================================================

	if(x2_>=0)
	{
		//index #2
		long int s;

		

		long int ind_count=-1;

		for(s=0;s<d_number_of_states;s++)
		{
			long int dim2_tmp=d_states_description[s].second;

			if(dim2_tmp!=0)
			{
				continue;
			};



			double q2=0;

			long int code1=0;
			long int code2=code2_;

			
			long int x2=x2_;
			long int x1=0;



			q2=d_states_distr_sums[s][x1][x2][code1][code2];

			ind_count++;

			long int ind_count2=-1;
			long int s1;
			for(s1=0;s1<d_number_of_states;s1++)
			{
				if(ind1_flag[s1])
				{
					ind_count2++;

					A1[ind_count][ind_count2]=d_transition_probabilities[s1][s]*q2;

					if(s==s1)
					{
						A1[ind_count][ind_count2]-=1.0;
					};

				};
			};



		};



		if(A1.size()>0)
		{
			double inside_eps=1e-12;
			std::vector<double> x1;
			FSA_utils::Gauss(
				A1,//matrix n*(n+1)
				x1,//solution
				inside_eps,
				A1_inv_);
		};

	};


}

//================================
void IS1_general_simulation::calculate_weight_W1_for_fixed_lengths_using_infinite_sums(
long int target_seq1_length_,//target length of the sequence #1
long int target_seq2_length_,//target length of the sequence #2
double &weight_)//the resulted weight
{
	weight_=0;

	long int margin2=FSA_utils::Tmax(d_IS1_general_obj->d_max_dim2,(long int)200);
	long int margin1=FSA_utils::Tmax(d_IS1_general_obj->d_max_dim1,3*margin2);

	margin1=FSA_utils::Tmin(margin1,d_W1[0]->d_max_ind1-target_seq1_length_);
	margin2=FSA_utils::Tmin(margin2,d_W1[0]->d_max_ind2-target_seq2_length_);

	//if(d_W1_seq1_current_length==0&&d_W1_seq2_current_length==0)
	{
		calculate_weight_W1_one_step(//calculates the weight assuming all the weights used in the recursive equations are defined
		target_seq1_length_,
		target_seq2_length_,
		0,//the first index (corresponded to the sequence #1)
		0);//the second index (corresponded to the sequence #2)
	};

	d_W1_seq1_current_length=-1;
	d_W1_seq2_current_length=-1;

	while(d_W1_seq1_current_length<target_seq1_length_+margin1)
	{

		d_W1_seq1_current_length++;

		long int i1_ind_max=d_W1_seq1_current_length;

		long int i2_ind_max=FSA_utils::Tmax(target_seq2_length_+d_IS1_general_obj->d_max_dim2-1,d_W1_seq2_current_length);

		if(d_W1_seq1_current_length>=target_seq1_length_||
			d_W1_seq1_current_length<=target_seq1_length_+d_IS1_general_obj->d_max_dim1-1)
		{
			i2_ind_max=FSA_utils::Tmax(target_seq2_length_+margin2,d_W1_seq2_current_length);
		};

		long int s;
		for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
		{
			d_W1[s]->set_max_ind(i1_ind_max,i2_ind_max);
		};

		long int i2;
		for(i2=0;i2<=i2_ind_max;i2++)
		{


			calculate_weight_W1_one_step(//calculates the weight assuming all the weights used in the recursive equations are defined
			target_seq1_length_,
			target_seq2_length_,
			d_W1_seq1_current_length,//the first index (corresponded to the sequence #1)
			i2);//the second index (corresponded to the sequence #2)

			{
				

				long int i1=d_W1_seq1_current_length;

				long int i1_diff=i1-target_seq1_length_;
				long int i2_diff=i2-target_seq2_length_;

				for(s=0;s<d_IS1_general_obj->d_number_of_states;s++)
				{

					if(i1_diff>=0&&i2_diff>=0&&
						(i1_diff<d_IS1_general_obj->d_states_description[s].first||
						i2_diff<d_IS1_general_obj->d_states_description[s].second)
						)
					{
						weight_+=d_W1[s]->get_element(i1,i2);
					};

				};
			};

		};
	};

}

void test::normalize_state_distributions(
long int alphabet_letters_number1_,//number of letters in the sequence #1
long int alphabet_letters_number2_,//number of letters in the sequence #2

long int number_of_states_,//number of states

pair<long int, long int> *states_description_,//description of the states; the index is a state number
double ***states_distr_)//distributions of the states; the index is a state number
{
	long int s;
	for(s=0;s<number_of_states_;s++)
	{
		long int dim1=FSA_utils::power_long(alphabet_letters_number1_,states_description_[s].first);
		long int dim2=FSA_utils::power_long(alphabet_letters_number2_,states_description_[s].second);

		double sum_tmp=0;
		long int i1,i2;
		for(i1=0;i1<dim1;i1++)
		{
			for(i2=0;i2<dim2;i2++)
			{
				sum_tmp+=states_distr_[s][i1][i2];
			};
		};

		if(sum_tmp<=0)
		{
			throw error("Error: the distribution defined by states_distr_ in test::normalize_state_distributions is incorrect\n",1);
		};

		for(i1=0;i1<dim1;i1++)
		{
			for(i2=0;i2<dim2;i2++)
			{
				states_distr_[s][i1][i2]/=sum_tmp;
			};
		};

	};
}

//=============================
void test::FSA_IS_transition_probabilities_calculation(
bool &FSA_flag,
double ***&states_distr,

bool &cs_flag,
double ***&states_distr_cs,

bool &sim_flag,
double ***&states_distr_sim,

long int &number_of_states,
long int &alphabet_letters_number1,
long int &alphabet_letters_number2,
pair<long int, long int> *&states_description,
long int &codon_length,
long int *&codon_AA,
double *&RR1,
double *&RR2,
map<string, long int> &state_name_into_number,
long int **&smatr,
double &ungappedlambda)
{
	long int *letters1=new long int[4];
	long int *letters2=new long int[1];
	long int *codon_tmp=new long int[3];

	long int s;
	for(s=0;s<number_of_states;s++)
	{

		long int i1,i2;

		long int max_ind1=FSA_utils::power_long(alphabet_letters_number1,states_description[s].first);
		long int max_ind2=FSA_utils::power_long(alphabet_letters_number2,states_description[s].second);

		for(i1=0;i1<max_ind1;i1++)
		{
			IS1_general::code_to_letters(//returns a unique code for the array of letters
			i1,//input code
			alphabet_letters_number1,//total number of letters
			states_description[s].first,//dimension of the array with letters
			letters1);//array of letters; the result

			double tmp1=1;
			long int t1;
			for(t1=0;t1<states_description[s].first;t1++)
			{
				tmp1*=RR1[letters1[t1]];
			};

			long int AA1=0;

			if(s==state_name_into_number["S1"])
			{
				AA1=FSA_utils::convert_codon_into_AA(
				codon_length,//codon length 
				codon_AA,//<codon code,AA number>
				alphabet_letters_number1,//number of letters for the sequence 1
				letters1);
			};

			if(s==state_name_into_number["S3"])
			{
				AA1=FSA_utils::convert_codon_into_AA(
				codon_length,//codon length 
				codon_AA,//<codon code,AA number>
				alphabet_letters_number1,//number of letters for the sequence 1
				letters1+1);
			};

			//simplified sampling
			if(sim_flag)
			{
				if(s==0)
				{
					states_distr_sim[1][i1][0]=tmp1;
				};
			};
			//simplified sampling - end



			for(i2=0;i2<max_ind2;i2++)
			{
				IS1_general::code_to_letters(//returns a unique code for the array of letters
				i2,//input code
				alphabet_letters_number2,//total number of letters
				states_description[s].second,//dimension of the array with letters
				letters2);//array of letters; the result

				double tmp2=1;
				long int t2;
				for(t2=0;t2<states_description[s].second;t2++)
				{
					tmp2*=RR2[letters2[t2]];
				};

					if(s==state_name_into_number["D1"]||s==state_name_into_number["D2"]||
						s==state_name_into_number["D3"]||s==state_name_into_number["I"]||
					s==state_name_into_number["F"])
					{
						if(FSA_flag)
						{
							states_distr[s][i1][i2]=tmp1*tmp2;
						};
						continue;
					};

				//crude sampling
				if(cs_flag)
				{
					if(s==0)
					{
						states_distr_cs[0][i1][i2]=tmp1*tmp2;
					};
				};
				//crude sampling - end

				//simplified sampling
				if(sim_flag)
				{
					if(s==0)
					{
						states_distr_sim[0][i1][i2]=tmp1*tmp2*exp(ungappedlambda*smatr[AA1][letters2[0]]);

						if(i1==0)
						{
							states_distr_sim[2][0][i2]=tmp2;
						};

					};
				};
				//simplified sampling - end



				if(s==state_name_into_number["S1"]||s==state_name_into_number["S3"])
				{
					if(FSA_flag)
					{
						states_distr[s][i1][i2]=tmp1*tmp2*exp(ungappedlambda*smatr[AA1][letters2[0]]);
					};
					continue;
				};

				if(s==state_name_into_number["S2"])
				{
					double tmp3=0;

					long int a1;
					for(a1=0;a1<alphabet_letters_number1;a1++)
					{
						codon_tmp[0]=a1;
						codon_tmp[1]=letters1[0];
						codon_tmp[2]=letters1[1];


						AA1=FSA_utils::convert_codon_into_AA(
						codon_length,//codon length 
						codon_AA,//<codon code,AA number>
						alphabet_letters_number1,//number of letters for the sequence 1
						codon_tmp);

						tmp3+=RR1[a1]*exp(ungappedlambda*smatr[AA1][letters2[0]]);


					};

					if(FSA_flag)
					{
						states_distr[s][i1][i2]=tmp1*tmp2*tmp3;
					};
					continue;
				};



			};

		};


	};

	delete[]letters1;
	delete[]letters2;
	delete[]codon_tmp;
}

void test::combine_parameters_from_forward_and_reversed_calculations_generalized(
Sls::FALP_set_of_parameters &par_,//parameters from forward calculation
Sls::FALP_set_of_parameters &par_reversed_,//parameters from reversed calculation
Sls::FALP_set_of_parameters &par_result_)//the result
{

	long int c12=par_.realizations_number+par_reversed_.realizations_number;
	if(c12<=0)
	{
		throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
	};
	double c1=(double)par_.realizations_number/(double)c12;
	double c2=(double)par_reversed_.realizations_number/(double)c12;

	par_result_.lambda=c1*par_.lambda+c2*par_reversed_.lambda;
	par_result_.lambda_error=FSA_utils::error_of_the_sum(//v1_+v2_
								c1*par_.lambda_error,
								c2*par_reversed_.lambda_error);

		

	par_result_.C=par_.C;
	par_result_.C_error=par_.C_error;

	par_result_.K_C=par_reversed_.K_C;
	par_result_.K_C_error=par_reversed_.K_C_error;

	par_result_.K=par_.C*par_reversed_.K_C;
	par_result_.K_error=FSA_utils::error_of_the_product(//v1_*v2_
								par_.C,
								par_.C_error,
								par_reversed_.K_C,
								par_reversed_.K_C_error);

	par_result_.a_I=c1*par_.a_I+c2*par_reversed_.a_I;
	par_result_.a_I=FSA_utils::Tmax(par_result_.a_I,0.0);
	par_result_.a_I_error=FSA_utils::error_of_the_sum(//v1_+v2_
								c1*par_.a_I_error,
								c2*par_reversed_.a_I_error);

	par_result_.a_J=c1*par_.a_J+c2*par_reversed_.a_J;
	par_result_.a_J=FSA_utils::Tmax(par_result_.a_J,0.0);
	par_result_.a_J_error=FSA_utils::error_of_the_sum(//v1_+v2_
								c1*par_.a_J_error,
								c2*par_reversed_.a_J_error);

	par_result_.sigma=c1*par_.sigma+c2*par_reversed_.sigma;
	par_result_.sigma=FSA_utils::Tmax(par_result_.sigma,0.0);
	par_result_.sigma_error=FSA_utils::error_of_the_sum(//v1_+v2_
								c1*par_.sigma_error,
								c2*par_reversed_.sigma_error);

	par_result_.alpha_I=c1*par_.alpha_I+c2*par_reversed_.alpha_I;
	par_result_.alpha_I=FSA_utils::Tmax(par_result_.alpha_I,0.0);
	par_result_.alpha_I_error=FSA_utils::error_of_the_sum(//v1_+v2_
								c1*par_.alpha_I_error,
								c2*par_reversed_.alpha_I_error);

	par_result_.alpha_J=c1*par_.alpha_J+c2*par_reversed_.alpha_J;
	par_result_.alpha_J=FSA_utils::Tmax(par_result_.alpha_J,0.0);
	par_result_.alpha_J_error=FSA_utils::error_of_the_sum(//v1_+v2_
								c1*par_.alpha_J_error,
								c2*par_reversed_.alpha_J_error);



	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	par_result_.G=0;
	par_result_.G1=0;
	par_result_.G2=0;

	/*
	par_result_.G=par_.G;
	par_result_.G1=par_.G1;
	par_result_.G2=par_.G2;
	*/


	size_t size_tmp=par_.m_LambdaSbs.size();

	if(size_tmp!=par_reversed_.m_LambdaSbs.size()||
		size_tmp!=par_reversed_.m_KSbs.size()||
		size_tmp!=par_reversed_.m_K_CSbs.size()||
		size_tmp!=par_.m_CSbs.size()||
		size_tmp!=par_reversed_.m_CSbs.size()||
		size_tmp!=par_.m_SigmaSbs.size()||
		size_tmp!=par_reversed_.m_SigmaSbs.size()||
		size_tmp!=par_.m_AlphaISbs.size()||
		size_tmp!=par_reversed_.m_AlphaISbs.size()||
		size_tmp!=par_.m_AlphaJSbs.size()||
		size_tmp!=par_reversed_.m_AlphaJSbs.size()||
		size_tmp!=par_.m_AISbs.size()||
		size_tmp!=par_reversed_.m_AISbs.size()||
		size_tmp!=par_.m_AJSbs.size()||
		size_tmp!=par_reversed_.m_AJSbs.size()
		)
	{
		throw error("Unexpected error in combine_parameters_from_forward_and_reversed_calculations_generalized\n",1);
	};

	par_result_.m_LambdaSbs.resize(size_tmp);
	par_result_.m_KSbs.resize(size_tmp);
	
	par_result_.m_K_CSbs.resize(size_tmp);

	par_result_.m_CSbs.resize(size_tmp);

	par_result_.m_SigmaSbs.resize(size_tmp);
	par_result_.m_AlphaISbs.resize(size_tmp);
	par_result_.m_AlphaJSbs.resize(size_tmp);

	par_result_.m_AISbs.resize(size_tmp);
	par_result_.m_AJSbs.resize(size_tmp);


	long int k;
	for(k=0;k<(long int)size_tmp;k++)
	{
		par_result_.m_LambdaSbs[k]=c1*par_.m_LambdaSbs[k]+c2*par_reversed_.m_LambdaSbs[k];
		par_result_.m_KSbs[k]=par_.m_CSbs[k]*par_reversed_.m_K_CSbs[k];
		
		par_result_.m_K_CSbs[k]=par_reversed_.m_K_CSbs[k];

		par_result_.m_CSbs[k]=par_.m_CSbs[k];

		par_result_.m_SigmaSbs[k]=c1*par_.m_SigmaSbs[k]+c2*par_reversed_.m_SigmaSbs[k];
		par_result_.m_SigmaSbs[k]=FSA_utils::Tmax(par_result_.m_SigmaSbs[k],0.0);

		par_result_.m_AlphaISbs[k]=c1*par_.m_AlphaISbs[k]+c2*par_reversed_.m_AlphaISbs[k];
		par_result_.m_AlphaISbs[k]=FSA_utils::Tmax(par_result_.m_AlphaISbs[k],0.0);

		par_result_.m_AlphaJSbs[k]=c1*par_.m_AlphaJSbs[k]+c2*par_reversed_.m_AlphaJSbs[k];
		par_result_.m_AlphaJSbs[k]=FSA_utils::Tmax(par_result_.m_AlphaJSbs[k],0.0);

		par_result_.m_AISbs[k]=c1*par_.m_AISbs[k]+c2*par_reversed_.m_AISbs[k];
		par_result_.m_AISbs[k]=FSA_utils::Tmax(par_result_.m_AISbs[k],0.0);

		par_result_.m_AJSbs[k]=c1*par_.m_AJSbs[k]+c2*par_reversed_.m_AJSbs[k];
		par_result_.m_AJSbs[k]=FSA_utils::Tmax(par_result_.m_AJSbs[k],0.0);
	};


}

void test::delete_mult_states_type(
mult_states_type *&states_old1)
{
	if(states_old1)
	{
		long int kk1;
		for(kk1=0;kk1<states_old1->d_number_of_realizations;kk1++)
		{
			

			delete states_old1->d_states[kk1].d_M_array;
			delete states_old1->d_states[kk1].d_seq1_length_array;
			delete states_old1->d_states[kk1].d_seq2_length_array;
			delete states_old1->d_states[kk1].d_distance_along_direction_1_array;
			delete states_old1->d_states[kk1].d_distance_along_direction_2_array;
			delete states_old1->d_states[kk1].d_ALP_weights_array;

			delete states_old1->d_states[kk1].d_IS1_general_simulation;
			delete states_old1->d_states[kk1].d_two_dim_layer_alignment_algorithm;

		};

		delete[]states_old1->d_states;

		delete states_old1;states_old1=NULL;
	};
}

//for test
void test::compare_mult_states(
mult_states_type *states_old_,
mult_states_type *states_new_)
{
	if(states_old_->d_average_ALP_pos1!=states_new_->d_average_ALP_pos1)
	{
		cout<<"states_old_->d_average_ALP_pos1!=states_new_->d_average_ALP_pos1\n";
	};

	if(states_old_->d_average_ALP_pos2!=states_new_->d_average_ALP_pos2)
	{
		cout<<"states_old_->d_average_ALP_pos2!=states_new_->d_average_ALP_pos2\n";
	};
	if(states_old_->d_average_ALP_pos1_mult_ALP_pos2!=states_new_->d_average_ALP_pos1_mult_ALP_pos2)
	{
		cout<<"states_old_->d_average_ALP_pos1_mult_ALP_pos2!=states_new_->d_average_ALP_pos1_mult_ALP_pos2\n";
	};
	if(states_old_->d_average_expanding_length1!=states_new_->d_average_expanding_length1)
	{
		cout<<"states_old_->d_average_expanding_length1!=states_new_->d_average_expanding_length1\n";
	};
	if(states_old_->d_average_expanding_length1_mult_expanding_length2!=states_new_->d_average_expanding_length1_mult_expanding_length2)
	{
		cout<<"states_old_->d_average_expanding_length1_mult_expanding_length2!=states_new_->d_average_expanding_length1_mult_expanding_length2\n";
	};
	if(states_old_->d_average_expanding_length2!=states_new_->d_average_expanding_length2)
	{
		cout<<"states_old_->d_average_expanding_length2!=states_new_->d_average_expanding_length2\n";
	};
	if(states_old_->d_number_of_ALP!=states_new_->d_number_of_ALP)
	{
		cout<<"states_old_->d_number_of_ALP!=states_new_->d_number_of_ALP\n";
	};
	if(states_old_->d_number_of_realizations!=states_new_->d_number_of_realizations)
	{
		cout<<"states_old_->d_number_of_realizations!=states_new_->d_number_of_realizations\n";
	};
	if(states_old_->d_total_calculation_time!=states_new_->d_total_calculation_time)
	{
		//cout<<"states_old_->d_total_calculation_time!=states_new_->d_total_calculation_time\n";
	};
	if(states_old_->d_states!=states_new_->d_states)
	{
		//cout<<"states_old_->d_states!=states_new_->d_states\n";
	};

	long int k;

	for(k=0;
		k<FSA_utils::Tmin(states_old_->d_number_of_realizations,states_new_->d_number_of_realizations);
		//k<1;
		k++)
	{
		if(states_old_->d_states[k].d_ALP_number!=states_new_->d_states[k].d_ALP_number)
		{
			//cout<<k<<"\tstates_old_->d_states[k].d_ALP_number!=states_new_->d_states[k].d_ALP_number\n";
		};

		if(states_old_->d_states[k].d_two_dim_layer_alignment_algorithm!=states_new_->d_states[k].d_two_dim_layer_alignment_algorithm)
		{
			//cout<<k<<"\tstates_old_->d_states[k].d_two_dim_layer_alignment_algorithm!=states_new_->d_states[k].d_two_dim_layer_alignment_algorithm\n";
		};

		if(states_old_->d_states[k].d_IS1_general_simulation!=states_new_->d_states[k].d_IS1_general_simulation)
		{
			//cout<<k<<"\tstates_old_->d_states[k].d_IS1_general_simulation!=states_new_->d_states[k].d_IS1_general_simulation\n";
		};


		long int i;
		for(i=0;i<=FSA_utils::Tmin(states_old_->d_states[k].d_ALP_number,states_new_->d_states[k].d_ALP_number);i++)
		{
			if(states_old_->d_states[k].d_ALP_weights_array->get_elem(i)!=states_new_->d_states[k].d_ALP_weights_array->get_elem(i))
			{
				cout<<k<<"\t"<<i<<"\tstates_old_->d_states[k].d_ALP_weights_array->get_elem(i)!=states_new_->d_states[k].d_ALP_weights_array->get_elem(i)\n";
			};

			if(states_old_->d_states[k].d_distance_along_direction_1_array->get_elem(i)!=states_new_->d_states[k].d_distance_along_direction_1_array->get_elem(i))
			{
				cout<<k<<"\t"<<i<<"\tstates_old_->d_states[k].d_distance_along_direction_1_array->get_elem(i)!=states_new_->d_states[k].d_distance_along_direction_1_array->get_elem(i)\n";
			};
			if(states_old_->d_states[k].d_distance_along_direction_2_array->get_elem(i)!=states_new_->d_states[k].d_distance_along_direction_2_array->get_elem(i))
			{
				cout<<k<<"\t"<<i<<"\tstates_old_->d_states[k].d_distance_along_direction_2_array->get_elem(i)!=states_new_->d_states[k].d_distance_along_direction_2_array->get_elem(i)\n";
			};
			if(states_old_->d_states[k].d_M_array->get_elem(i)!=states_new_->d_states[k].d_M_array->get_elem(i))
			{
				cout<<k<<"\t"<<i<<"\tstates_old_->d_states[k].d_M_array->get_elem(i)!=states_new_->d_states[k].d_M_array->get_elem(i)\n";
			};
			if(states_old_->d_states[k].d_seq1_length_array->get_elem(i)!=states_new_->d_states[k].d_seq1_length_array->get_elem(i))
			{
				cout<<k<<"\t"<<i<<"\tstates_old_->d_states[k].d_seq1_length_array->get_elem(i)!=states_new_->d_states[k].d_seq1_length_array->get_elem(i)\n";
			};
			if(states_old_->d_states[k].d_seq2_length_array->get_elem(i)!=states_new_->d_states[k].d_seq2_length_array->get_elem(i))
			{
				cout<<k<<"\t"<<i<<"\tstates_old_->d_states[k].d_seq2_length_array->get_elem(i)!=states_new_->d_states[k].d_seq2_length_array->get_elem(i)\n";
			};

		};



	};

}

//*****************************

void test::FSA_Align(

long int open1_,//gap opening penalty for the nucleotide sequence #1
long int open2_,//gap opening penalty for the amino acid sequence #2

long int epen1_,//gap extension penalty for the nucleotide sequence #1
long int epen2_,//gap extension penalty for the amino acid sequence #2

long int gamma_,//frameshift penalty gamma

string smatr_file_name_,//scoring matrix file name
string DNA_codon_table_file_name_,//a name of a file with DNA codon table

bool insertions_after_deletions_,//if true, then insertions after deletions are allowed

string input_sequences_,//name of file with sequences for alignment (align mode)
string output_sequences_,//name of output file with alignment scores (align mode)
string gumbelparin_file_name_)//Gumbel parameters input file name
{

	double time0_start;
	FSA_utils::get_current_time(time0_start);

	Sls::FALP_set_of_parameters gumbelparin;//the resulted parameters

	if(gumbelparin_file_name_!="")
	{
		fsa_par::Read_Params(
		gumbelparin,
		gumbelparin_file_name_);

		FALP_pvalues::compute_intercepts(gumbelparin);
	};

	static Sls::FALP_pvalues pvalues_obj;

	long int **smatr=NULL;


	char *alphabet1=NULL;//alphabet letters for the sequence #1
	char *alphabet2=NULL;//alphabet letters for the sequence #2

	long int *alphabet1_to_long=NULL;//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
	long int *alphabet2_to_long=NULL;//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

	long int codon_length;//codon length 
	long int *codon_AA=NULL;//<codon code,AA number>


//====================================================


	long int number_of_AA_smatr;
	long int smatr_min;

	FSA_utils::read_smatr(
	smatr_file_name_,
	smatr,
	number_of_AA_smatr,
	smatr_min);

//====================================================================================

	long int number_of_letters1_tmp;//number of letters for the sequence 1
	long int number_of_letters2_tmp;//number of letters for the sequence 2

	FSA_utils::read_codon_AA_file(
	DNA_codon_table_file_name_,
	number_of_letters1_tmp,//number of letters for the sequence 1
	number_of_letters2_tmp,//number of letters for the sequence 2

	alphabet1,//alphabet letters for the sequence #1
	alphabet2,//alphabet letters for the sequence #2

	alphabet1_to_long,//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
	alphabet2_to_long,//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

	codon_length,//codon length 
	codon_AA);//<codon code,AA number>

	if(number_of_AA_smatr!=number_of_letters2_tmp)
	{
		throw error("Error - different numbers of amino acids in the files "+smatr_file_name_+", "+DNA_codon_table_file_name_+"\n",1);
	};



	long int number_of_sequences;
	string *headers;
	long int *lengths1;//lengths of the sequences #1
	long int *lengths2;//lengths of the sequences #2
	long int **sequences1;
	long int **sequences2;

	FSA_utils::read_sequences_for_alingment(

	input_sequences_,

	number_of_letters1_tmp,//number of letters for the sequence 1
	number_of_letters2_tmp,//number of letters for the sequence 2

	alphabet1,//alphabet letters for the sequence #1
	alphabet2,//alphabet letters for the sequence #2

	number_of_sequences,

	headers,
	lengths1,
	lengths2,
	sequences1,//the first index numerates sequences; the second - sequence letters
	sequences2);

	long int k;

	long int max_ind1=-1;
	long int max_ind2=-1;
	for(k=0;k<number_of_sequences;k++)
	{
		max_ind1=FSA_utils::Tmax(max_ind1,lengths1[k]);
		max_ind2=FSA_utils::Tmax(max_ind2,lengths2[k]);
	};

	if(max_ind1==0||max_ind2==0)
	{
		throw error("Error in the file "+input_sequences_+"\n",1);
	};


//---------two dim alignment algorithm---------------

	long int var_num_dim=1000;
	long int *var_num=new long int[var_num_dim];
	long int i;
	for(i=0;i<var_num_dim;i++)
	{
		var_num[i]=-1;
	};

	var_num['S']=0;
	var_num['D']=1;
	var_num['I']=2;


	data_for_FSA_alignment data_test;

	data_test.d_alphabet_letters_number1=number_of_letters1_tmp;

	data_test.d_open1=open1_;
	data_test.d_open2=open2_;

	data_test.d_epen1=epen1_;
	data_test.d_epen2=epen2_;

	data_test.d_gamma=gamma_;

	data_test.d_smatr=smatr;

	data_test.d_codon_AA=codon_AA;

	data_test.d_insertions_after_deletions=insertions_after_deletions_;

	//data_test.d_alignment_type="global";
	data_test.d_alignment_type="local";

	long int depth1=4;//the maximum difference of the first index in the dynamic equations
	long int depth2=1;//the maximum difference of the second index in the dynamic equations



	long int number_of_variables_for_the_alignment=3;

	two_dim_layer_alignment_algorithm<long int>* two_dim_layer_alignment_algorithm_test=NULL;



	two_dim_layer_alignment_algorithm_test=
	new two_dim_layer_alignment_algorithm<long int>(
	number_of_variables_for_the_alignment,//total number of variables in dynamic equations
	depth1,//the maximum difference of the first index in the dynamic equations
	depth2,//the maximum difference of the second index in the dynamic equations
	max_ind1,//max of the index #1 (minimum index is 0 by default)
	max_ind2,//max of the index #2 (minimum index is 0 by default)
	0,//null element of T
	var_num);//d_variable_name_to_number[c] converts a letter 'c' (of type char) to order number of variable

	two_dim_layer_alignment_algorithm_test->d_two_dim_layer_alignment_function=test::two_dim_layer_alignment_function_FSA;
	two_dim_layer_alignment_algorithm_test->d_two_dim_layer_boundary_function=test::two_dim_layer_alignment_function_FSA;
	two_dim_layer_alignment_algorithm_test->d_par=&data_test;


	two_dim_layer_alignment_algorithm_test->d_M_flag=true;
	two_dim_layer_alignment_algorithm_test->d_E_flag=false;

	two_dim_layer_alignment_algorithm_test->d_scores_counts_flag=false;

	ofstream fout(output_sequences_.data());
	if(!fout)
	{
		throw error("Error - the file "+output_sequences_+" is not found\n",1);
	};

	//calculate alingment scores and P-values
	for(k=0;k<number_of_sequences;k++)
	{
		two_dim_layer_alignment_algorithm_test->init();

		data_test.d_seq1=sequences1[k];
		data_test.d_seq2=sequences2[k];

		two_dim_layer_alignment_algorithm_test->align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
		//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
		lengths1[k],//target length of the sequence #1
		lengths2[k]);//target length of the sequence #2
		fout<<headers[k]<<endl;
		fout<<two_dim_layer_alignment_algorithm_test->d_M;

		if(gumbelparin_file_name_!="")
		{
				vector<double> P_values;
				vector<double> P_values_errors;

				vector<double> E_values;
				vector<double> E_values_errors;


				pvalues_obj.calculate_P_values(
				two_dim_layer_alignment_algorithm_test->d_M,
				two_dim_layer_alignment_algorithm_test->d_M,
				lengths1[k],
				lengths2[k],
				gumbelparin,
				P_values,
				P_values_errors,
				E_values,
				E_values_errors);


				if(P_values.size()==1&&P_values_errors.size()==1)
				{
					fout<<"\t"<<P_values[0]<<"\t"<<P_values_errors[0];
				};

				if(E_values.size()==1&&E_values_errors.size()==1)
				{
					fout<<"\t"<<E_values[0]<<"\t"<<E_values_errors[0];
				};

		};
		fout<<endl;
	};

	fout.close();

	delete[]headers;
	delete[]lengths1;
	delete[]lengths2;

	
	for(k=0;k<number_of_sequences;k++)
	{
		delete[]sequences1[k];
		delete[]sequences2[k];
	};

	delete[]sequences1;
	delete[]sequences2;
	delete[]var_num;

	delete two_dim_layer_alignment_algorithm_test;

	double time0_end;
	FSA_utils::get_current_time(time0_end);

	cout<<"Total calculation time\t"<<time0_end-time0_start<<endl;

}

//*****************************
void test::input_data_for_the_constructor(

string DNA_codon_table_file_name_,//a name of a file with DNA codon table
string smatr_file_name_,//scoring matrix file name
string RR1_file_name_,//probabilities1 file name
string RR2_file_name_,//probabilities2 file name

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
long int *&codon_AA_)//<codon code,AA number>
{

	long int number_of_AA_smatr;
	long int smatr_min;


	FSA_utils::read_smatr(
	smatr_file_name_,
	substitutionScoreMatrix_,
	number_of_AA_smatr,
	smatr_min);


	FSA_utils::read_RR(
	RR1_file_name_,
	letterFreqs1_,
	alphabetSize1_,
	4);


	FSA_utils::read_RR(
	RR2_file_name_,
	letterFreqs2_,
	alphabetSize2_,
	25);



	if(number_of_AA_smatr!=alphabetSize2_)
	{
		throw error("Error - different number of letters in the files "+smatr_file_name_+", "+RR2_file_name_+"\n",1);
	};

//====================================================================================

	long int number_of_letters1_tmp;//number of letters for the sequence 1
	long int number_of_letters2_tmp;//number of letters for the sequence 2


	FSA_utils::read_codon_AA_file(
	DNA_codon_table_file_name_,
	number_of_letters1_tmp,//number of letters for the sequence 1
	number_of_letters2_tmp,//number of letters for the sequence 2

	alphabet1_,//alphabet letters for the sequence #1
	alphabet2_,//alphabet letters for the sequence #2

	alphabet1_to_long_,//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
	alphabet2_to_long_,//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

	codon_length_,//codon length 
	codon_AA_);//<codon code,AA number>



	if(alphabetSize1_!=number_of_letters1_tmp)
	{
		throw error("The numbers of letters is different in the files "+DNA_codon_table_file_name_+" and "+RR1_file_name_+"\n",3);
	};
	if(alphabetSize2_!=number_of_letters2_tmp)
	{
		throw error("The numbers of letters is different in the files "+DNA_codon_table_file_name_+" and "+RR2_file_name_+"\n",3);
	};

}

void test::FSA_IS(
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

string smatr_file_name_,//scoring matrix file name
string RR1_file_name_,//background frequencies file name for the sequence #1
string RR2_file_name_,//background frequencies file name for the sequence #2
string DNA_codon_table_file_name_,//a name of a file with DNA codon table

double eps_lambda_,//relative error for lambda calculation
double eps_K_,//relative error for K calculation

bool gapped_flag_,//if true, then the gapped alingment is performed
double max_time_,//maximum allowed calculation time in seconds
double max_mem_,//maximum allowed memory usage in MB

//additional parameters
long int seq_number_,//number of tested alignments
long int nalp_,//number of ALPs for the calculation
long int number_of_subsets_for_errors_calculation_,//number of subsets used for the splitting method
bool forward_and_reverse_screen_output_flag_,//determines whether the parameters are outputted for forward and reverse calculations
bool insertions_after_deletions_,//if true, then insertions after deletions are allowed

//for test
double mult_for_is_lambda_,//multiplier for lambda in the IS

//the result
Sls::FALP_set_of_parameters &par_result_,//the resulted parameters
Sls::par_test1_type *par_test1_)//for tests
{
	double seconds1;
	FSA_utils::get_current_time(seconds1);

	FSA_utils::process_random_factor(
	rand_);

	double empirical_mem_coeff=2.5;

	double max_mem_for_states_alloc=max_mem_*1048576/empirical_mem_coeff;//maximum allocatin size for states in bytes
	double max_number_of_states=FSA_utils::Tmin((double)inf,max_mem_for_states_alloc/sizeof(state_type));

	//the parameters

	bool FSA_flag=true;//if true, then FSA IS is performed
	bool cs_flag=true;//if true, then crude sampling simulation is also performed
	bool sim_flag=false;//flag for the simplified IS technique

	double eps_K=eps_K_;//relative error for K

	long int depth1=4;//the maximum difference of the first index in the dynamic equations
	long int depth2=1;//the maximum difference of the second index in the dynamic equations

	long int max_ind2=10000;//max of the index #2
	long int max_ind1=max_ind2*3;//max of the index #1


	long int target_ALP=nalp_;
	

	long int number_of_realizations=seq_number_;
	long int number_of_sets=number_of_subsets_for_errors_calculation_;//number of sets for error calculation

	bool futher_expanding=true;
	long int number_of_cs_steps_for_expanding=max_ind2;

	long int mult_margin=5;
	long int add_margin=100;

	bool lambda_output_flag=false;

#ifdef _MSDOS_
	string lambda_file_name="lambda_"+FSA_utils::long_to_string(rand_)+".out";
#else
	string lambda_file_name="/net/frosty/vol/export1/sls/ff1/res1/lambda_"+FSA_utils::long_to_string(rand_)+".out";
#endif

//-------------------------------------------------

	//bool add_fake_state_flag=true;
	bool add_fake_state_flag=false;

	FSA_utils::srand2(rand_);
	if(!library_call_flag_)
	{
		cout<<"Randomization factor\t"<<rand_<<endl;
	};

	long int alphabet_letters_number1;//number of letters in the sequence #1
	long int alphabet_letters_number2;//number of letters in the sequence #2
	double *RR1=NULL;//background probability for the sequence #1
	double *RR2=NULL;//background probability for the sequence #2
	long int number_of_states=7;//number of states
	if(add_fake_state_flag)
	{
		number_of_states++;
	};
	double **transition_probabilities=NULL;//transition probabilities between states; matrix d_number_of_states x d_number_of_states
	pair<long int, long int> *states_description=NULL;//description of the states; the index is a state number
	double ***states_distr=NULL;//distributions of the states; the index is a state number
								//the second and the third indexes correspond to an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second

	double ***states_distr_reversed=NULL;

//crude sampling
	pair<long int, long int> *states_description_cs=NULL;
	double ***states_distr_cs=NULL;
	double **transition_probabilities_cs=NULL;
//crude sampling - end

//simplified sampling
	double ***states_distr_sim=NULL;
//simplified sampling - end

	double *RR1_sum=NULL;
	long int *RR1_sum_elements=NULL;
	double *RR2_sum=NULL;
	long int *RR2_sum_elements=NULL;

	long int **smatr=NULL;

	double ungappedlambda;

	d_IS1=NULL;
	d_IS1_general_simulation=NULL;

	double *RR1_AA=NULL;


//===============================================

	char *alphabet1=NULL;//alphabet letters for the sequence #1
	char *alphabet2=NULL;//alphabet letters for the sequence #2

	long int *alphabet1_to_long=NULL;//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
	long int *alphabet2_to_long=NULL;//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

	long int codon_length;//codon length 
	long int *codon_AA=NULL;//<codon code,AA number>
	long int *codon_reversed_AA=NULL;//<codon code,AA number>


//====================================================

	{
		
		if(!library_call_flag_)
		{
			input_data_for_the_constructor(

			DNA_codon_table_file_name_,//a name of a file with DNA codon table
			smatr_file_name_,//scoring matrix file name
			RR1_file_name_,//probabilities1 file name
			RR2_file_name_,//probabilities2 file name

			alphabet_letters_number1,
			alphabet_letters_number2,
			smatr,
			RR1,
			RR2,

			alphabet1,//alphabet letters for the sequence #1
			alphabet2,//alphabet letters for the sequence #2

			alphabet1_to_long,//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
			alphabet2_to_long,//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

			codon_length,//codon length 
			codon_AA);//<codon code,AA number>
		}
		else
		{
			alphabet_letters_number1=ntAlphabetSize_;
			alphabet_letters_number2=aaAlphabetSize_;

			RR1=new double[ntAlphabetSize_];
			RR2=new double[aaAlphabetSize_];

			alphabet1=new char[ntAlphabetSize_];//alphabet letters for the sequence #1
			alphabet2=new char[aaAlphabetSize_];//alphabet letters for the sequence #2


			long int i,j;
			for(i=0;i<ntAlphabetSize_;i++)
			{
				RR1[i]=ntFreqs_[i];
				alphabet1[i]=(char)(i+1);
			};

			for(i=0;i<aaAlphabetSize_;i++)
			{
				RR2[i]=aaFreqs_[i];
				alphabet2[i]=(char)(i+1);
			};

			FSA_utils::get_memory_for_matrix(aaAlphabetSize_,aaAlphabetSize_,smatr);
			for(i=0;i<aaAlphabetSize_;i++)
			{
				for(j=0;j<aaAlphabetSize_;j++)
				{
					smatr[i][j]=substitutionScoreMatrix_[i][j];
				};
			};
			

			alphabet1_to_long=new long int [length_max];
			FSA_utils::assert_mem(alphabet1_to_long);
			alphabet2_to_long=new long int [length_max];
			FSA_utils::assert_mem(alphabet2_to_long);

			long int k;
			for(k=0;k<length_max;k++)
			{
				alphabet1_to_long[k]=-1;
				alphabet2_to_long[k]=-1;
			};

			for(k=0;k<alphabet_letters_number1;k++)
			{
				alphabet1_to_long[(size_t)alphabet1[k]]=k;
			};
			for(k=0;k<alphabet_letters_number2;k++)
			{
				alphabet2_to_long[(size_t)alphabet2[k]]=k;
			};


			codon_length=3;//codon length 

			long int number_of_codons=1;
			for(k=0;k<codon_length;k++)
			{
				number_of_codons*=alphabet_letters_number1;
			};

			codon_AA=new long int [number_of_codons];
			FSA_utils::assert_mem(codon_AA);

			for(k=0;k<number_of_codons;k++)
			{
				codon_AA[k]=codonTable_[k];//<codon code,AA number>
			};


		};

		long int number_of_AA_smatr=alphabet_letters_number2;
		long int smatr_min;

		long int number_of_letters1_tmp=alphabet_letters_number1;
		long int number_of_letters2_tmp=alphabet_letters_number2;

		FSA_utils::calculate_RR_sum(
		RR1,
		number_of_letters1_tmp,
		RR1_sum,
		RR1_sum_elements);

		FSA_utils::calculate_RR_sum(
		RR2,
		number_of_letters2_tmp,
		RR2_sum,
		RR2_sum_elements);

		FSA_utils::smatr_min(
		smatr,
		number_of_AA_smatr,
		smatr_min);

		{
			FSA_utils::remove_zero_probabilities(
			RR1,
			RR1_sum,
			RR1_sum_elements,
			alphabet_letters_number1,
			RR2,
			RR2_sum,
			RR2_sum_elements,
			alphabet_letters_number2,
			smatr,
			number_of_AA_smatr,
			smatr_min,

			number_of_letters1_tmp,//number of letters for the sequence 1
			number_of_letters2_tmp,//number of letters for the sequence 2

			alphabet1,//alphabet letters for the sequence #1
			alphabet2,//alphabet letters for the sequence #2

			alphabet1_to_long,//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
			alphabet2_to_long,//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

			codon_length,//codon length 
			codon_AA);//<codon code,AA number>

		};


//====================================================================================

		data_for_lambda_equation_FSA func_pointer;

		func_pointer.d_number_of_letters1=alphabet_letters_number1;
		func_pointer.d_number_of_letters2=alphabet_letters_number2;
		func_pointer.d_smatr=smatr;
		func_pointer.d_RR1=RR1;
		func_pointer.d_RR2=RR2;

		func_pointer.d_codon_length=codon_length;


		if(futher_expanding)
		{
			FSA_utils::reverse_codons(
			codon_AA,//<codon code,AA number>; original codons
			alphabet_letters_number1,//number of letters for the sequence #1
			codon_length,//codon length 
			codon_reversed_AA);//<codon code,AA number>; reversed codons

		};

		func_pointer.d_codon_AA=codon_AA;

		calculate_ungapped_lambda_general(
		lambda_equation_FSA,
		(void*) &func_pointer,
		ungappedlambda);

		ungappedlambda*=mult_for_is_lambda_;
		//cout<<"Multiplier for IS lambda\t"<<mult_for_is_lambda_<<endl;

		//gapless calculaton begin
		double GaplessTimePortion=0.5;

		double GaplessCalculationTime=max_time_;

		if(max_time_<=0)
		{
			GaplessCalculationTime=120;
		};

		if(gapped_flag_)
		{
			//Gapless calculation may take only a portion of maximum allowed calculation time in the case of gapped calculation 
			GaplessCalculationTime*=GaplessTimePortion;
		};



		

		FSA_utils::extract_AA_frequencies_for_DNA_sequence(
		codon_AA,//<codon code,AA number>
		codon_length,//codon length 
		number_of_letters1_tmp,//number of letters for the sequence 1
		number_of_letters2_tmp,//number of letters for the sequence 2
		RR1,//nucleotide probabilities
		RR1_AA);//the resulted frequencies

		
		Njn::LocalMaxStatMatrix local_max_stat_matrix(number_of_letters2_tmp,
							  smatr,
							  RR1_AA,
							  RR2,
							  number_of_letters2_tmp,
							  GaplessCalculationTime);

		if(local_max_stat_matrix.getTerminated()) 
		{
			throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
		};

		//calculation of a and sigma
		double calculation_error=1e-6;

		par_result_.gapless_alpha_J = local_max_stat_matrix.getAlpha ();
		par_result_.gapless_alpha_J=FSA_utils::Tmax(par_result_.gapless_alpha_J,0.0);
		par_result_.gapless_alpha_J_error = calculation_error;

		par_result_.gapless_alpha_I=par_result_.gapless_alpha_J*9;
		par_result_.gapless_alpha_I_error = calculation_error*9;

		par_result_.gapless_sigma=par_result_.gapless_alpha_J*3;
		par_result_.gapless_sigma_error = calculation_error*3;


		par_result_.gapless_a_J = local_max_stat_matrix.getA ();
		par_result_.gapless_a_J=FSA_utils::Tmax(par_result_.gapless_a_J,0.0);
		par_result_.gapless_a_J_error = calculation_error;

		par_result_.gapless_a_I = par_result_.gapless_a_J*3;
		par_result_.gapless_a_I_error = calculation_error*3;



		if(!gapped_flag_) {


			//calculation of all required parameters for a gapless case
			par_result_.G=0;
			par_result_.G1=0;
			par_result_.G2=0;

			par_result_.lambda = local_max_stat_matrix.getLambda ();
			par_result_.lambda_error = calculation_error;

			par_result_.K = local_max_stat_matrix.getK ();
			par_result_.K_error = calculation_error;
				
			par_result_.C = local_max_stat_matrix.getC ();;
			par_result_.C_error = calculation_error;

			if(par_result_.C!=0)
			{
				par_result_.K_C = par_result_.K/par_result_.C;
			}
			else
			{
				par_result_.K_C = 0;
			};
			par_result_.K_C_error = calculation_error;


			par_result_.sigma = par_result_.gapless_sigma;
			par_result_.sigma_error = par_result_.gapless_sigma_error;

			par_result_.alpha_I = par_result_.gapless_alpha_I;
			par_result_.alpha_I_error = par_result_.gapless_alpha_I_error;

			par_result_.alpha_J = par_result_.gapless_alpha_J;
			par_result_.alpha_J_error = par_result_.gapless_alpha_J_error;

			par_result_.a_I = par_result_.gapless_a_I;
			par_result_.a_I_error = par_result_.gapless_a_I_error;

			par_result_.a_J = par_result_.gapless_a_J;
			par_result_.a_J_error = par_result_.gapless_a_J_error;


			std::vector<double > sbs_arrays;

			sbs_arrays.resize(2);
			sbs_arrays[0]=par_result_.lambda;
			sbs_arrays[1]=par_result_.lambda + calculation_error;

			par_result_.m_LambdaSbs=sbs_arrays;


			sbs_arrays.resize(2);
			sbs_arrays[0]=par_result_.K;
			sbs_arrays[1]=par_result_.K+calculation_error;

			par_result_.m_KSbs=sbs_arrays;


			sbs_arrays.resize(2);
			sbs_arrays[0]=par_result_.C;
			sbs_arrays[1]=par_result_.C+calculation_error;

			par_result_.m_CSbs=sbs_arrays;


			sbs_arrays.resize(2);
			sbs_arrays[0]=par_result_.sigma;
			sbs_arrays[1]=par_result_.sigma + calculation_error;

			par_result_.m_SigmaSbs=sbs_arrays;


			sbs_arrays.resize(2);
			sbs_arrays[0]=par_result_.alpha_I;
			sbs_arrays[1]=par_result_.alpha_I + calculation_error;

			par_result_.m_AlphaISbs=sbs_arrays;


			sbs_arrays.resize(2);
			sbs_arrays[0]=par_result_.alpha_J;
			sbs_arrays[1]=par_result_.alpha_J + calculation_error;

			par_result_.m_AlphaJSbs=sbs_arrays;


			sbs_arrays.resize(2);
			sbs_arrays[0]=par_result_.a_I;
			sbs_arrays[1]=par_result_.a_I + calculation_error;

			par_result_.m_AISbs=sbs_arrays;


			sbs_arrays.resize(2);
			sbs_arrays[0]=par_result_.a_J;
			sbs_arrays[1]=par_result_.a_J + calculation_error;

			par_result_.m_AJSbs=sbs_arrays;


			FSA_flag=false;//if true, then FSA IS is performed
			cs_flag=false;//if true, then crude sampling simulation is also performed
			sim_flag=false;//flag for the simplified IS technique
			futher_expanding=false;

		};


		//gapless calculation end

	};

	//possible choice #1
	//long int minopen=FSA_utils::Tmin(open1_,open2_);
	//long int minepen=FSA_utils::Tmin(epen1_,epen2_);

	//possible choice #2
	//long int minopen=FSA_utils::Tmax(open1_,open2_);
	//long int minepen=FSA_utils::Tmax(epen1_,epen2_);

	//possible choice #3
	long int minopen=(long int)FSA_utils::round((double)(open1_+open2_)*0.5);
	long int minepen=(long int)FSA_utils::round((double)(epen1_+epen2_)*0.5);

	long int choice_of_IS_weights=1;


	long int gamma=gamma_;


	map<string, long int> state_name_into_number;
	map<long int,string> state_number_into_name;

	
	
	state_name_into_number["S1"]=0;
	state_name_into_number["S2"]=1;
	state_name_into_number["S3"]=2;
	state_name_into_number["D1"]=3;
	state_name_into_number["D2"]=4;
	state_name_into_number["D3"]=5;
	state_name_into_number["I"]=6;
	if(add_fake_state_flag)
	{
		state_name_into_number["F"]=7;
	}
	else
	{
		state_name_into_number["F"]=-1;
	};



	state_number_into_name[0]="S1";
	state_number_into_name[1]="S2";
	state_number_into_name[2]="S3";
	state_number_into_name[3]="D1";
	state_number_into_name[4]="D2";
	state_number_into_name[5]="D3";
	state_number_into_name[6]="I";
	if(add_fake_state_flag)
	{
		state_number_into_name[7]="F";
	};

	states_description=new pair<long int, long int>[number_of_states];
	FSA_utils::assert_mem(states_description);

	states_description[state_name_into_number["S1"]]=make_pair(3,1);
	states_description[state_name_into_number["S2"]]=make_pair(2,1);
	states_description[state_name_into_number["S3"]]=make_pair(4,1);
	states_description[state_name_into_number["D1"]]=make_pair(3,0);
	states_description[state_name_into_number["D2"]]=make_pair(2,0);
	states_description[state_name_into_number["D3"]]=make_pair(4,0);
	states_description[state_name_into_number["I"]]=make_pair(0,1);
	if(add_fake_state_flag)
	{
		states_description[state_name_into_number["F"]]=make_pair(0,0);
	};

	if(FSA_flag)
	{
		FSA_utils::get_memory_for_matrix(number_of_states,number_of_states,transition_probabilities);

		if(choice_of_IS_weights==1)
		{
			transition_probabilities[0][0]=((-1 + exp(gamma*ungappedlambda))*(-1 + exp(minepen*ungappedlambda)))/(1 - exp(gamma*ungappedlambda) + exp(minepen*ungappedlambda) + exp((gamma + minepen)*ungappedlambda) + 4*exp((minepen - minopen)*ungappedlambda) + 2*exp((gamma + minepen - minopen)*ungappedlambda));
			transition_probabilities[0][1]=exp((minepen + minopen)*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[0][2]=exp((minepen + minopen)*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[0][3]=exp((gamma + minepen)*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[0][4]=exp(minepen*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[0][5]=exp(minepen*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[0][6]=((2 + exp(gamma*ungappedlambda))*exp(minepen*ungappedlambda))/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[1][0]=((-1 + exp(gamma*ungappedlambda))*(-1 + exp(minepen*ungappedlambda)))/(1 - exp(gamma*ungappedlambda) + exp(minepen*ungappedlambda) + exp((gamma + minepen)*ungappedlambda) + 4*exp((minepen - minopen)*ungappedlambda) + 2*exp((gamma + minepen - minopen)*ungappedlambda));
			transition_probabilities[1][1]=exp((minepen + minopen)*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[1][2]=exp((minepen + minopen)*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[1][3]=exp((gamma + minepen)*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[1][4]=exp(minepen*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[1][5]=exp(minepen*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[1][6]=((2 + exp(gamma*ungappedlambda))*exp(minepen*ungappedlambda))/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[2][0]=((-1 + exp(gamma*ungappedlambda))*(-1 + exp(minepen*ungappedlambda)))/(1 - exp(gamma*ungappedlambda) + exp(minepen*ungappedlambda) + exp((gamma + minepen)*ungappedlambda) + 4*exp((minepen - minopen)*ungappedlambda) + 2*exp((gamma + minepen - minopen)*ungappedlambda));
			transition_probabilities[2][1]=exp((minepen + minopen)*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[2][2]=exp((minepen + minopen)*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[2][3]=exp((gamma + minepen)*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[2][4]=exp(minepen*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[2][5]=exp(minepen*ungappedlambda)/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[2][6]=((2 + exp(gamma*ungappedlambda))*exp(minepen*ungappedlambda))/(4*exp(minepen*ungappedlambda) + 2*exp((gamma + minepen)*ungappedlambda) + exp(minopen*ungappedlambda) - exp((gamma + minopen)*ungappedlambda) + exp((minepen + minopen)*ungappedlambda) + exp((gamma + minepen + minopen)*ungappedlambda));
			transition_probabilities[3][0]=(-1 + exp(gamma*ungappedlambda))*(-1 + exp(minepen*ungappedlambda))*exp(-((gamma + minepen)*ungappedlambda));
			transition_probabilities[3][1]=(exp(-(gamma*ungappedlambda)) - exp(-((gamma + minepen)*ungappedlambda)))/2.;
			transition_probabilities[3][2]=(exp(-(gamma*ungappedlambda)) - exp(-((gamma + minepen)*ungappedlambda)))/2.;
			transition_probabilities[3][3]=(-1 + exp(gamma*ungappedlambda))*exp(-((gamma + minepen)*ungappedlambda));
			transition_probabilities[3][4]=exp(-((gamma + minepen)*ungappedlambda))/2.;
			transition_probabilities[3][5]=exp(-((gamma + minepen)*ungappedlambda))/2.;
			transition_probabilities[3][6]=0;
			transition_probabilities[4][0]=(-1 + exp(gamma*ungappedlambda))*(-1 + exp(minepen*ungappedlambda))*exp(-((gamma + minepen)*ungappedlambda));
			transition_probabilities[4][1]=(exp(-(gamma*ungappedlambda)) - exp(-((gamma + minepen)*ungappedlambda)))/2.;
			transition_probabilities[4][2]=(exp(-(gamma*ungappedlambda)) - exp(-((gamma + minepen)*ungappedlambda)))/2.;
			transition_probabilities[4][3]=(-1 + exp(gamma*ungappedlambda))*exp(-((gamma + minepen)*ungappedlambda));
			transition_probabilities[4][4]=exp(-((gamma + minepen)*ungappedlambda))/2.;
			transition_probabilities[4][5]=exp(-((gamma + minepen)*ungappedlambda))/2.;
			transition_probabilities[4][6]=0;
			transition_probabilities[5][0]=(-1 + exp(gamma*ungappedlambda))*(-1 + exp(minepen*ungappedlambda))*exp(-((gamma + minepen)*ungappedlambda));
			transition_probabilities[5][1]=(exp(-(gamma*ungappedlambda)) - exp(-((gamma + minepen)*ungappedlambda)))/2.;
			transition_probabilities[5][2]=(exp(-(gamma*ungappedlambda)) - exp(-((gamma + minepen)*ungappedlambda)))/2.;
			transition_probabilities[5][3]=(-1 + exp(gamma*ungappedlambda))*exp(-((gamma + minepen)*ungappedlambda));
			transition_probabilities[5][4]=exp(-((gamma + minepen)*ungappedlambda))/2.;
			transition_probabilities[5][5]=exp(-((gamma + minepen)*ungappedlambda))/2.;
			transition_probabilities[5][6]=0;
			transition_probabilities[6][0]=1 - exp(-(gamma*ungappedlambda)) - exp(-(minepen*ungappedlambda)) + 2*exp(-((gamma + minepen)*ungappedlambda));
			transition_probabilities[6][1]=((-1 + exp(minepen*ungappedlambda))*exp(-((gamma + minepen)*ungappedlambda)))/2.;
			transition_probabilities[6][2]=(exp(-(gamma*ungappedlambda)) - exp(-((gamma + minepen)*ungappedlambda)))/2.;
			transition_probabilities[6][3]=0;
			transition_probabilities[6][4]=0;
			transition_probabilities[6][5]=0;
			transition_probabilities[6][6]=exp(-(minepen*ungappedlambda)) - exp(-((gamma + minepen)*ungappedlambda));
		};

		//for test
		if(choice_of_IS_weights==2)
		{
			transition_probabilities[0][0]=1/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[0][1]=exp(-(gamma*ungappedlambda))/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[0][2]=exp(-(gamma*ungappedlambda))/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[0][3]=exp((-epen1_ - open1_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[0][4]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[0][5]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[0][6]=exp((-epen2_ - open2_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[1][0]=1/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[1][1]=exp(-(gamma*ungappedlambda))/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[1][2]=exp(-(gamma*ungappedlambda))/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[1][3]=exp((-epen1_ - open1_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[1][4]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[1][5]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[1][6]=exp((-epen2_ - open2_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[2][0]=1/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[2][1]=exp(-(gamma*ungappedlambda))/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[2][2]=exp(-(gamma*ungappedlambda))/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[2][3]=exp((-epen1_ - open1_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[2][4]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[2][5]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[2][6]=exp((-epen2_ - open2_)*ungappedlambda)/(1 + 2*exp(-(gamma*ungappedlambda)) + exp((-epen1_ - open1_)*ungappedlambda) + 2*exp((-epen1_ - gamma - open1_)*ungappedlambda) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[3][0]=1/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[3][1]=exp(-(gamma*ungappedlambda))/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[3][2]=exp(-(gamma*ungappedlambda))/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[3][3]=exp(-(epen1_*ungappedlambda))/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[3][4]=exp((-epen1_ - gamma)*ungappedlambda)/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[3][5]=exp((-epen1_ - gamma)*ungappedlambda)/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[3][6]=exp((-epen2_ - open2_)*ungappedlambda)/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[4][0]=1/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[4][1]=exp(-(gamma*ungappedlambda))/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[4][2]=exp(-(gamma*ungappedlambda))/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[4][3]=exp(-(epen1_*ungappedlambda))/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[4][4]=exp((-epen1_ - gamma)*ungappedlambda)/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[4][5]=exp((-epen1_ - gamma)*ungappedlambda)/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[4][6]=exp((-epen2_ - open2_)*ungappedlambda)/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[5][0]=1/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[5][1]=exp(-(gamma*ungappedlambda))/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[5][2]=exp(-(gamma*ungappedlambda))/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[5][3]=exp(-(epen1_*ungappedlambda))/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[5][4]=exp((-epen1_ - gamma)*ungappedlambda)/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[5][5]=exp((-epen1_ - gamma)*ungappedlambda)/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[5][6]=exp((-epen2_ - open2_)*ungappedlambda)/(1 + exp(-(epen1_*ungappedlambda)) + 2*exp((-epen1_ - gamma)*ungappedlambda) + 2*exp(-(gamma*ungappedlambda)) + exp((-epen2_ - open2_)*ungappedlambda));
			transition_probabilities[6][0]=1/(1 + exp(-(epen2_*ungappedlambda)) + 2*exp(-(gamma*ungappedlambda)));
			transition_probabilities[6][1]=exp(-(gamma*ungappedlambda))/(1 + exp(-(epen2_*ungappedlambda)) + 2*exp(-(gamma*ungappedlambda)));
			transition_probabilities[6][2]=exp(-(gamma*ungappedlambda))/(1 + exp(-(epen2_*ungappedlambda)) + 2*exp(-(gamma*ungappedlambda)));
			transition_probabilities[6][3]=0;
			transition_probabilities[6][4]=0;
			transition_probabilities[6][5]=0;
			transition_probabilities[6][6]=exp(-(epen2_*ungappedlambda))/(1 + exp(-(epen2_*ungappedlambda)) + 2*exp(-(gamma*ungappedlambda)));
		};

		if(choice_of_IS_weights==3)
		{
			transition_probabilities[0][0]=1/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[0][1]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[0][2]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[0][3]=exp((-epen1_ - open1_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[0][4]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[0][5]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[0][6]=exp((-epen2_ - open2_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[1][0]=1/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[1][1]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[1][2]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[1][3]=exp((-epen1_ - open1_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[1][4]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[1][5]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[1][6]=exp((-epen2_ - open2_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[2][0]=1/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[2][1]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[2][2]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[2][3]=exp((-epen1_ - open1_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[2][4]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[2][5]=exp((-epen1_ - gamma - open1_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[2][6]=exp((-epen2_ - open2_)*ungappedlambda)/(3.*(0.3333333333333333 + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen1_ - open1_)*ungappedlambda)/3. + (2*exp((-epen1_ - gamma - open1_)*ungappedlambda))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[3][0]=1/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[3][1]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[3][2]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[3][3]=exp(-(epen1_*ungappedlambda))/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[3][4]=exp((-epen1_ - gamma)*ungappedlambda)/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[3][5]=exp((-epen1_ - gamma)*ungappedlambda)/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[3][6]=exp((-epen2_ - open2_)*ungappedlambda)/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[4][0]=1/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[4][1]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[4][2]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[4][3]=exp(-(epen1_*ungappedlambda))/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[4][4]=exp((-epen1_ - gamma)*ungappedlambda)/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[4][5]=exp((-epen1_ - gamma)*ungappedlambda)/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[4][6]=exp((-epen2_ - open2_)*ungappedlambda)/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[5][0]=1/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[5][1]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[5][2]=exp(-(gamma*ungappedlambda))/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[5][3]=exp(-(epen1_*ungappedlambda))/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[5][4]=exp((-epen1_ - gamma)*ungappedlambda)/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[5][5]=exp((-epen1_ - gamma)*ungappedlambda)/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[5][6]=exp((-epen2_ - open2_)*ungappedlambda)/(3.*(0.3333333333333333 + exp(-(epen1_*ungappedlambda))/3. + (2*exp((-epen1_ - gamma)*ungappedlambda))/3. + (2*exp(-(gamma*ungappedlambda)))/3. + exp((-epen2_ - open2_)*ungappedlambda)/3.));
			transition_probabilities[6][0]=1/(1 + exp(-(epen2_*ungappedlambda)) + 2*exp(-(gamma*ungappedlambda)));
			transition_probabilities[6][1]=exp(-(gamma*ungappedlambda))/(1 + exp(-(epen2_*ungappedlambda)) + 2*exp(-(gamma*ungappedlambda)));
			transition_probabilities[6][2]=exp(-(gamma*ungappedlambda))/(1 + exp(-(epen2_*ungappedlambda)) + 2*exp(-(gamma*ungappedlambda)));
			transition_probabilities[6][3]=0;
			transition_probabilities[6][4]=0;
			transition_probabilities[6][5]=0;
			transition_probabilities[6][6]=exp(-(epen2_*ungappedlambda))/(1 + exp(-(epen2_*ungappedlambda)) + 2*exp(-(gamma*ungappedlambda)));
		};
		//for test

		//introducing a fake state
		if(add_fake_state_flag)
		{
			long int i;
			for(i=0;i<number_of_states;i++)
			{
				transition_probabilities[i][state_name_into_number["F"]]=0;
				transition_probabilities[state_name_into_number["F"]][i]=0;
			};

			transition_probabilities[state_name_into_number["F"]][state_name_into_number["S1"]]=1.0/3.0;
			transition_probabilities[state_name_into_number["F"]][state_name_into_number["S2"]]=1.0/3.0;
			transition_probabilities[state_name_into_number["F"]][state_name_into_number["S3"]]=1.0-1.0/3.0-1.0/3.0;

		};


		//test transition probablities
		long int i,j;
		for(i=0;i<number_of_states;i++)
		{
			for(j=0;j<number_of_states;j++)
			{
				if(transition_probabilities[i][j]<0)
				{
					throw error("Error - transition probablities of the importance sampling are negative;\nthe method is not applicable for the input scoring scheme\n",1);
				};
			};
		};
	
	};

//crude sampling

	long int number_of_states_cs=1;

	if(cs_flag)
	{
		states_description_cs=new pair<long int, long int>[number_of_states_cs];
		FSA_utils::assert_mem(states_description_cs);

		states_description_cs[0]=make_pair(3,1);

		FSA_utils::get_memory_for_matrix(number_of_states_cs,number_of_states_cs,transition_probabilities_cs);
		transition_probabilities_cs[0][0]=1;
	};

//crude sampling - end


//=========================

	long int s;

	if(FSA_flag)
	{
		states_distr=new double **[number_of_states];

		long int s;
		for(s=0;s<number_of_states;s++)
		{
			states_distr[s]=IS1_general::allocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
				alphabet_letters_number1,//number of letters in the sequence #1
				alphabet_letters_number2,//number of letters in the sequence #2
				states_description[s]);//state description

		};
	};

	if(futher_expanding&&FSA_flag)
	{
		states_distr_reversed=new double **[number_of_states];

		long int s;
		for(s=0;s<number_of_states;s++)
		{
			states_distr_reversed[s]=IS1_general::allocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
				alphabet_letters_number1,//number of letters in the sequence #1
				alphabet_letters_number2,//number of letters in the sequence #2
				states_description[s]);//state description

		};

	};


	if(cs_flag)
	{
		states_distr_cs=new double **[number_of_states_cs];

		long int s;
		for(s=0;s<number_of_states_cs;s++)
		{
			states_distr_cs[s]=IS1_general::allocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
				alphabet_letters_number1,//number of letters in the sequence #1
				alphabet_letters_number2,//number of letters in the sequence #2
				states_description_cs[s]);//state description

		};

	};

//=========================
	FSA_IS_transition_probabilities_calculation(
	FSA_flag,
	states_distr,

	cs_flag,
	states_distr_cs,

	sim_flag,
	states_distr_sim,

	number_of_states,
	alphabet_letters_number1,
	alphabet_letters_number2,
	states_description,
	codon_length,
	codon_AA,
	RR1,
	RR2,
	state_name_into_number,
	smatr,
	ungappedlambda);

	if(futher_expanding&&FSA_flag)
	{
		bool cs_flag_tmp=false;
		bool sim_flag_tmp=false;

		FSA_IS_transition_probabilities_calculation(
		FSA_flag,
		states_distr_reversed,

		cs_flag_tmp,
		states_distr_cs,

		sim_flag_tmp,
		states_distr_sim,

		number_of_states,
		alphabet_letters_number1,
		alphabet_letters_number2,
		states_description,
		codon_length,
		codon_reversed_AA,
		RR1,
		RR2,
		state_name_into_number,
		smatr,
		ungappedlambda);

	};

//=========================

	if(futher_expanding&&FSA_flag)
	{
		normalize_state_distributions(
		alphabet_letters_number1,//number of letters in the sequence #1
		alphabet_letters_number2,//number of letters in the sequence #2
		number_of_states,//number of states
		states_description,//description of the states; the index is a state number
		states_distr_reversed);//distributions of the states; the index is a state number

	};

	if(FSA_flag)
	{
		normalize_state_distributions(
		alphabet_letters_number1,//number of letters in the sequence #1
		alphabet_letters_number2,//number of letters in the sequence #2
		number_of_states,//number of states
		states_description,//description of the states; the index is a state number
		states_distr);//distributions of the states; the index is a state number

	};

	if(cs_flag)
	{
		normalize_state_distributions(
		alphabet_letters_number1,//number of letters in the sequence #1
		alphabet_letters_number2,//number of letters in the sequence #2
		number_of_states_cs,//number of states
		states_description_cs,//description of the states; the index is a state number
		states_distr_cs);//distributions of the states; the index is a state number

	};

//=========================

	IS1_general* IS1_reversed=NULL;
	if(FSA_flag&&futher_expanding)
	{
		IS1_reversed=new IS1_general(
		alphabet_letters_number1,//number of letters in the sequence #1
		alphabet_letters_number2,//number of letters in the sequence #2
		RR1,//background probability for the sequence #1
		RR2,//background probability for the sequence #2
		number_of_states,//number of states
		transition_probabilities,//transition probabilities between states; matrix d_number_of_states x d_number_of_states
		states_description,//description of the states; the index is a state number
		states_distr_reversed);//distributions of the states; the index is a state number

	};


	if(FSA_flag)
	{
		d_IS1=new IS1_general(
		alphabet_letters_number1,//number of letters in the sequence #1
		alphabet_letters_number2,//number of letters in the sequence #2
		RR1,//background probability for the sequence #1
		RR2,//background probability for the sequence #2
		number_of_states,//number of states
		transition_probabilities,//transition probabilities between states; matrix d_number_of_states x d_number_of_states
		states_description,//description of the states; the index is a state number
		states_distr);//distributions of the states; the index is a state number

	};


	//crude sampling 


	IS1_general* IS1_cs=NULL;
	if(cs_flag)
	{
		IS1_cs=new IS1_general(
		alphabet_letters_number1,//number of letters in the sequence #1
		alphabet_letters_number2,//number of letters in the sequence #2
		RR1,//background probability for the sequence #1
		RR2,//background probability for the sequence #2
		number_of_states_cs,//number of states
		transition_probabilities_cs,//transition probabilities between states; matrix d_number_of_states x d_number_of_states
		states_description_cs,//description of the states; the index is a state number
		states_distr_cs);//distributions of the states; the index is a state number
	};
	//crude sampling - end


	ofstream lambda_out;

	if(lambda_output_flag)
	{
		lambda_out.open(lambda_file_name.data());
		if(!lambda_out)
		{
			throw error("Error - the file "+lambda_file_name+" is not found\n",1);
		};
	};


	long int var_num_dim=1000;
	long int *var_num=new long int[var_num_dim];
	long int i;
	for(i=0;i<var_num_dim;i++)
	{
		var_num[i]=-1;
	};

	var_num['S']=0;
	var_num['D']=1;
	var_num['I']=2;


	data_for_FSA_alignment data_test;
	

	data_test.d_alphabet_letters_number1=alphabet_letters_number1;

	data_test.d_open1=open1_;
	data_test.d_open2=open2_;

	data_test.d_epen1=epen1_;
	data_test.d_epen2=epen2_;

	data_test.d_gamma=gamma_;

	data_test.d_smatr=smatr;

	data_test.d_codon_AA=codon_AA;

	data_test.d_insertions_after_deletions=insertions_after_deletions_;

	data_test.d_alignment_type="global";
	//data_test.d_alignment_type="local";

	if(lambda_output_flag)
	{
		lambda_out<<open1_<<"\t"<<epen1_<<"\t"<<open2_<<"\t"<<epen2_<<"\t"<<gamma_<<"\t"<<number_of_realizations<<"\t"<<target_ALP<<endl;
		lambda_out<<"lambda\t";
		lambda_out<<"C\t";
		lambda_out<<"a_I\t";
		lambda_out<<"a_J\t";
		lambda_out<<"sigma\t";
		lambda_out<<"alpha_I\t";
		lambda_out<<"alpha_J\t";
		lambda_out<<"K_C\t";
		lambda_out<<endl;

	};


	long int number_of_variables_for_the_alignment=3;


	two_dim_layer_alignment_algorithm<long int>* two_dim_layer_alignment_algorithm_test=NULL;

	Sls::FALP_set_of_parameters par;

	double time0_start;
	FSA_utils::get_current_time(time0_start);

	bool accuracy_is_achieved_flag=false;

//--------------------------
		//crude sampling
		string M_distr_file_name_cs=lambda_file_name+"_M_cs.out";
		data_for_FSA_alignment data_test_cs;

		IS1_general_simulation *IS1_general_simulation_cs=NULL;
		two_dim_layer_alignment_algorithm<long int> *two_dim_layer_alignment_algorithm_test_cs=NULL;

		//crude sampling - end
//--------------------------

	long int M_thr_estimation_direct=0;//average score threshold corresponded to target_ALP_


	if(FSA_flag)
	{
		long int initial_state=state_name_into_number["S1"];

		if(add_fake_state_flag)
		{
			initial_state=state_name_into_number["F"];
		};



		d_IS1_general_simulation=new IS1_general_simulation(
		d_IS1,
		initial_state,//initial state for the IS
		max_ind1*mult_margin+add_margin,//maximum sequence length
		max_ind2*mult_margin+3*add_margin);//maximum sequence length


		data_test.d_seq1=d_IS1_general_simulation->d_seq1;//sequence #1
		data_test.d_seq2=d_IS1_general_simulation->d_seq2;//sequence #2
	

		two_dim_layer_alignment_algorithm_test=
		new two_dim_layer_alignment_algorithm<long int>(
		number_of_variables_for_the_alignment,//total number of variables in dynamic equations
		depth1,//the maximum difference of the first index in the dynamic equations
		depth2,//the maximum difference of the second index in the dynamic equations
		max_ind1,//max of the index #1 (minimum index is 0 by default)
		max_ind2,//max of the index #2 (minimum index is 0 by default)
		0,//null element of T
		var_num);//d_variable_name_to_number[c] converts a letter 'c' (of type char) to order number of variable




		//string M_distr_file_name="M_distr_FSA.out";
		string M_distr_file_name=lambda_file_name+"_M.out";

		if(par_test1_)
		{
			M_distr_file_name=par_test1_->d_gumbelparout_file_name;
		};

		
		two_dim_layer_alignment_algorithm_test->d_two_dim_layer_alignment_function=test::two_dim_layer_alignment_function_FSA;
		two_dim_layer_alignment_algorithm_test->d_two_dim_layer_boundary_function=test::two_dim_layer_alignment_function_FSA;
		two_dim_layer_alignment_algorithm_test->d_par=&data_test;


		two_dim_layer_alignment_algorithm_test->d_M_flag=true;
		two_dim_layer_alignment_algorithm_test->d_E_flag=false;

		two_dim_layer_alignment_algorithm_test->d_scores_counts_flag=false;


		double fraction_of_ALP_opt=0.5;


		{

			//bool futher_expanding_tmp=true;
			bool futher_expanding_tmp=false;
			

			long int limit2=max_ind2;

			two_dim_layer_alignment_algorithm_test->d_scores_counts_flag=futher_expanding_tmp;
			two_dim_layer_alignment_algorithm_test->d_FSC_flag=true;



			//additional test whether the scoring scheme is logarithmic

			double time_for_one_realization_without_killing=0;
			double time_for_one_realization_with_killing=0;
			long int M_thr_estimation_direct_from_the_minimal_test;

			{
				double time_tmp00;
				double time_tmp11;
				

				double ending_time_tmp2=-1;
				bool stopped_by_ending_time_flag_tmp2;
				long int number_of_realizations_with_ending_time_tmp2;
				double ending_time_to_test_logarithmic_regime2=-1;

				bool inside_simulation_flag=true;
				bool save_states_flag=false;
				mult_states_type *states_old1=NULL;
				mult_states_type *states_new1=NULL;
				bool futher_expanding_test=true;
				long int M_thr=-inf;

				//number of test realizations for the test
				long int number_of_realizations_test=1;
				long int number_of_sets_test=1;
				


				try
				{

					two_dim_layer_alignment_algorithm_test->d_scores_counts_flag=true;
					bool parameters_are_not_calculated=true;

					

					FSA_utils::get_current_time(time_tmp00);

					bool compute_allocated_memory_flag_tmp=false;
					double allocated_memory_tmp;


					//without expanding with calculation of M-threshold
					collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
					compute_allocated_memory_flag_tmp,
					allocated_memory_tmp,

					number_of_realizations_test,
					M_distr_file_name,
					*two_dim_layer_alignment_algorithm_test,
					*d_IS1_general_simulation,
					target_ALP,//target ALP number
					ungappedlambda,
					limit2,

					number_of_sets_test,
					par,
					inside_simulation_flag,
					false,
					false,

					ending_time_tmp2,
					stopped_by_ending_time_flag_tmp2,
					number_of_realizations_with_ending_time_tmp2,

					ending_time_to_test_logarithmic_regime2,

					save_states_flag,
					states_old1,
					states_new1,
					M_thr_estimation_direct,

					//futher_expanding_test,
					false,
					M_thr,
					IS1_cs,
					&eps_K,
					&number_of_cs_steps_for_expanding,
					parameters_are_not_calculated);

					FSA_utils::get_current_time(time_tmp11);

					time_for_one_realization_without_killing=(time_tmp11-time_tmp00)/(double)number_of_realizations_test;
					if(time_for_one_realization_without_killing<0)
					{
						time_for_one_realization_without_killing=0;
					};

					{
						

						bool compute_allocated_memory_flag_tmp=false;
						double allocated_memory_tmp;

						//test with expanding
						collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
						compute_allocated_memory_flag_tmp,
						allocated_memory_tmp,


						number_of_realizations_test,
						M_distr_file_name,
						*two_dim_layer_alignment_algorithm_test,
						*d_IS1_general_simulation,
						target_ALP,//target ALP number
						ungappedlambda,
						limit2,

						number_of_sets_test,
						par,
						inside_simulation_flag,
						false,
						false,

						ending_time_tmp2,
						stopped_by_ending_time_flag_tmp2,
						number_of_realizations_with_ending_time_tmp2,

						ending_time_to_test_logarithmic_regime2,

						save_states_flag,
						states_old1,
						states_new1,
						M_thr_estimation_direct_from_the_minimal_test,

						futher_expanding_test,
						M_thr_estimation_direct,
						IS1_cs,
						&eps_K,
						&number_of_cs_steps_for_expanding,
						parameters_are_not_calculated);

						M_thr_estimation_direct_from_the_minimal_test=(long int)FSA_utils::round(0.5*(M_thr_estimation_direct_from_the_minimal_test+
							M_thr_estimation_direct));
					};


					two_dim_layer_alignment_algorithm_test->d_scores_counts_flag=false;


				}
				catch (error er)
				{ 
					if(er.error_code==10)
					{
						throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
					}
					else
					{
						throw error(er.st,er.error_code);
					};

				};

				double time_tmp22;
				FSA_utils::get_current_time(time_tmp22);

				time_for_one_realization_with_killing=(time_tmp22-time_tmp11)/(double)number_of_realizations_test;
				if(time_for_one_realization_with_killing<0)
				{
					time_for_one_realization_with_killing=0;
				};


			};




			//automatic adjusment of ALP number
			if(max_time_>0)
			{
				if(!library_call_flag_)
				{
					cout<<"\nPreliminary stage...\n";
				};

				double mult_for_prediction=1.0;

				double mult0=0.5;
				double mult1=0.5*mult0;

				long int max_number_of_realizations0=100;
				



				long int number_of_realizations0=max_number_of_realizations0;

				//long int number_of_realizations1_without_extending=number_of_realizations0;
				long int number_of_realizations1_with_extending=100;

				long int number_of_sets1_without_extending=5;
				long int number_of_sets1_with_extending=1;
				

				//double max_calculation_time_for_preliminary_stage=60;//in seconds
				double max_calculation_time_for_preliminary_stage=inf;//in seconds

				double time_preliminary_stage_expected=
					//time_for_one_realization_with_killing*number_of_realizations1_with_extending+
					//number_of_realizations1_without_extending*time_for_one_realization_without_killing;
					time_for_one_realization_with_killing*number_of_realizations1_with_extending;

				if(max_calculation_time_for_preliminary_stage<time_preliminary_stage_expected)
				{
					max_calculation_time_for_preliminary_stage=1.5*time_preliminary_stage_expected;
				};

				if(max_calculation_time_for_preliminary_stage>max_time_*mult0)
				{
					max_calculation_time_for_preliminary_stage=max_time_*mult0;
				};

				if(!library_call_flag_)
				{
					cout<<"Maximum time for the preliminary stage: "<<max_calculation_time_for_preliminary_stage<<" sec\n\n";
				};



				bool test_output1=false;

				
				ofstream ff;

				if(test_output1)
				{
					string fn="ALP_stat_test1.out";

					ff.open(fn.data());
					if(!ff)
					{
						throw error("Error - file "+fn+" is not found\n",1);
					};


					ff<<"number of realizations\t\
number of ALP\t\
time per one cell (ALP)\t\
ALP position for seq #1\t\
ALP position for seq #2\t\
product of positions\t\
expected relative error of lambda\t\
time per one cell (expanding)\t\
expanding length for seq #1\t\
expanding length for seq #2\t\
product of lengths\n";
				};


				bool save_states_flag=true;
				mult_states_type *states_old1=NULL;
				mult_states_type *states_new1=new mult_states_type;
				Sls::FALP_set_of_parameters par1;

				bool further_expanding_tmp1=true;

				long int ALP_max_number=10;

				long int ALP1_start=3;

				vector<bool> ALP_flag(ALP_max_number+1,false);


				long int ALP1=ALP1_start;
				
				

				long int count_tmp1=0;
				vector<double> lambda_relative_errors;
				
				long int number_of_ALP_for_errors=1;

				double L1=0;
				double L2=0;
				double L1_L2=0;

				double X1;
				double X2;
				double X1_X2;

				double time_cell_ALP=0;
				double time_cell_exp=0;
				double time_exp=0;

				long int ALP_opt=3;
				double error_opt=inf;
				mult_states_type *states_opt=NULL;
				long int expected_number_of_realizations1_opt=0;
				

				long int delete_old_object_flag=0;

				bool inside_simulation_flag;

				long int count_of_failed_relative_errors_max=3;
				long int count_of_failed_relative_errors=0;
				bool flag1=true;
				double total_time_without_extending=0;
				double total_time_without_extending_opt=0;

				long int number_of_realizations1=0;
				long int number_of_sets1=0;

				while(flag1)
				{
					
					Sls::FALP_set_of_parameters par1_old;

					if(states_old1)
					{
						par1_old=par1;
					};

					
					double ending_time_tmp1=-1.0;
					bool stopped_by_ending_time_flag_tmp1;
					long int number_of_realizations_with_ending_time_tmp1;
					double ending_time_to_test_logarithmic_regime1;
					long int M_thr=-inf;

						bool compute_allocated_memory_flag_tmp;
						double allocated_memory_tmp;

					
					if(count_tmp1==0)
					{
						compute_allocated_memory_flag_tmp=true;
						allocated_memory_tmp=0;

						number_of_realizations1=number_of_realizations1_with_extending;
						number_of_sets1=number_of_sets1_with_extending;
						ending_time_to_test_logarithmic_regime1=time0_start+FSA_utils::Tmin(max_calculation_time_for_preliminary_stage,mult1*max_time_);
						M_thr=M_thr_estimation_direct_from_the_minimal_test;
					}
					else
					{
						compute_allocated_memory_flag_tmp=false;

						//number_of_realizations1=number_of_realizations1_without_extending;
						number_of_sets1=number_of_sets1_without_extending;
						ending_time_to_test_logarithmic_regime1=-1;
					};

					double time_without_extending_tmp1;
					FSA_utils::get_current_time(time_without_extending_tmp1);

					collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
					compute_allocated_memory_flag_tmp,
					allocated_memory_tmp,

					number_of_realizations1,
					M_distr_file_name,
					*two_dim_layer_alignment_algorithm_test,
					*d_IS1_general_simulation,
					ALP1,//target ALP number
					ungappedlambda,
					limit2,
					number_of_sets1,
					par1,
					inside_simulation_flag,
					false,
					false,

					ending_time_tmp1,
					stopped_by_ending_time_flag_tmp1,
					number_of_realizations_with_ending_time_tmp1,

					ending_time_to_test_logarithmic_regime1,

					save_states_flag,
					states_old1,
					states_new1,
					M_thr_estimation_direct,

					false,
					M_thr,
					IS1_cs,
					&eps_K,
					&number_of_cs_steps_for_expanding);

					if(compute_allocated_memory_flag_tmp)
					{
						if(allocated_memory_tmp>0)
						{
							max_number_of_states=FSA_utils::round(FSA_utils::Tmin((double)inf,max_mem_for_states_alloc/allocated_memory_tmp/((double)((ALP1+1)*(ALP1+1))/(double)(ALP1*ALP1))/2.0));
						};
					}
					else
					{
						max_number_of_states/=(double)((ALP1+1)*(ALP1+1))/(double)(ALP1*ALP1);
					};



					double time_without_extending_tmp2;
					FSA_utils::get_current_time(time_without_extending_tmp2);


					total_time_without_extending+=time_without_extending_tmp2-time_without_extending_tmp1;


					if(further_expanding_tmp1)
					{
						delete_mult_states_type(states_old1);
						
						states_old1=states_new1;
						states_new1=new mult_states_type;

						bool compute_allocated_memory_flag_tmp=false;
						double allocated_memory_tmp;


						collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
						compute_allocated_memory_flag_tmp,
						allocated_memory_tmp,

						number_of_realizations1,
						M_distr_file_name,
						*two_dim_layer_alignment_algorithm_test,
						*d_IS1_general_simulation,
						ALP1,//target ALP number
						ungappedlambda,
						limit2,
						number_of_sets1,
						par1,
						inside_simulation_flag,
						false,
						false,

						ending_time_tmp1,
						stopped_by_ending_time_flag_tmp1,
						number_of_realizations_with_ending_time_tmp1,

						ending_time_to_test_logarithmic_regime1,

						save_states_flag,
						states_old1,
						states_new1,
						M_thr_estimation_direct,

						further_expanding_tmp1,
						M_thr,
						IS1_cs,
						&eps_K,
						&number_of_cs_steps_for_expanding);


						time_cell_exp=states_new1->d_total_calculation_time;
						time_exp=time_cell_exp;

						states_new1->d_total_calculation_time=states_old1->d_total_calculation_time;


						delete_mult_states_type(states_old1);
						states_old1=NULL;


					};

					count_tmp1++;

					if(inside_simulation_flag)
					{
						if(!library_call_flag_)
						{
							cout<<"ALP#="<<ALP1<<"; realizations#="<<number_of_realizations1<<"; time="<<states_new1->d_total_calculation_time<<"\n";
						};
						ALP_flag[ALP1]=true;
					}
					else
					{
						if(!library_call_flag_)
						{
							cout<<"ALP#="<<ALP1<<"; realizations#="<<number_of_realizations1<<"; time="<<states_new1->d_total_calculation_time<<": the calculation is not successful\n";
						};
					};

					X1=states_new1->d_average_ALP_pos1/(double)states_new1->d_number_of_ALP;
					X2=states_new1->d_average_ALP_pos2/(double)states_new1->d_number_of_ALP;
					X1_X2=states_new1->d_average_ALP_pos1_mult_ALP_pos2/((double)(states_new1->d_number_of_ALP*states_new1->d_number_of_ALP));

					if(further_expanding_tmp1)
					{
						L1=states_new1->d_average_expanding_length1;
						L2=states_new1->d_average_expanding_length2;
						L1_L2=states_new1->d_average_expanding_length1_mult_expanding_length2;


						if(states_new1->d_total_number_of_exp_cells>0)
						{
							time_cell_exp=time_exp/(double)states_new1->d_total_number_of_exp_cells;
						}
						else
						{
							throw error("Unexpected error\n",1);
						};
					};

					double time_cell_ALP_tmp=states_new1->d_total_calculation_time/((double)states_new1->d_total_number_of_ALP_cells);

					time_cell_ALP=time_cell_ALP_tmp;

					if(test_output1)
					{
						ff<<number_of_realizations1<<"\t"<<states_new1->d_number_of_ALP<<"\t"
							<<time_cell_ALP<<"\t"
							<<X1<<"\t"
							<<X2<<"\t"
							<<X1_X2;

					};






					//check relative errors
					if(ALP_flag[ALP1]&&ALP_flag[ALP1-1])
					{


						
						double T_per_realization=(2*ALP1*ALP1*X1_X2*time_cell_ALP+(ALP1*(X1*L2+X2*L1)+L1_L2)*time_cell_exp);

						if(T_per_realization<=0)
						{
							throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
						};

						long int expected_number_of_realizations1=(long int)FSA_utils::round(
							FSA_utils::Tmin(max_number_of_states,
							mult_for_prediction*(max_time_-time_exp)
							/T_per_realization)
							);

						expected_number_of_realizations1=FSA_utils::Tmax(expected_number_of_realizations1,(long int)1);


						double lambda_error1=par1.lambda_last_ALP_relative_error;


						lambda_error1*=sqrt((double)number_of_realizations1/(2.0*(double)expected_number_of_realizations1));

						lambda_relative_errors.push_back(lambda_error1);
						if(!library_call_flag_)
						{
							cout<<"Expected lambda relative error for ALP#="<<ALP1<<" is "<<lambda_error1<<endl<<endl;
						};

						if(test_output1)
						{
							ff<<"\t"<<lambda_error1;
						};

						long int vect_size=(long int)lambda_relative_errors.size();
						if(vect_size>=number_of_ALP_for_errors)
						{
							double average_relative_error=0;

							long int j;
							for(j=0;j<number_of_ALP_for_errors;j++)
							{
								average_relative_error+=lambda_relative_errors[vect_size-1-j];
							};

							average_relative_error/=(double)number_of_ALP_for_errors;


							if(average_relative_error<eps_lambda_)
							{
								flag1=false;
							};

							delete_old_object_flag++;


							if(average_relative_error<error_opt)
							{
								total_time_without_extending_opt=total_time_without_extending;

								count_of_failed_relative_errors=0;
								ALP_opt=ALP1;
								error_opt=average_relative_error;

								expected_number_of_realizations1_opt=expected_number_of_realizations1;

								//calculaton of the fraction of ALP
								double nom=ALP1*ALP1*X1_X2*time_cell_ALP;
								fraction_of_ALP_opt=nom/T_per_realization;


								if(states_opt==states_old1)
								{
									states_old1=NULL;
								};

								delete_mult_states_type(states_opt);
								states_opt=states_new1;

								delete_old_object_flag=-2;

							}
							else
							{
								if(count_of_failed_relative_errors_max>0)
								{
									count_of_failed_relative_errors++;
									if(count_of_failed_relative_errors>count_of_failed_relative_errors_max)
									{
										flag1=false;
									};
								};
							};
						};
					}
					else
					{
						if(test_output1)
						{
							ff<<"\t#N/A";
						};
						if(!library_call_flag_)
						{
							cout<<endl;
						};
					};


					if(test_output1)
					{

						if(further_expanding_tmp1)
						{
							ff<<"\t"<<time_cell_exp<<"\t"
							<<L1<<"\t"
							<<L2<<"\t"
							<<L1_L2;
						}
						else
						{
							ff<<"\t#N/A\t"
							<<"#N/A\t"
							<<"#N/A\t"
							<<"#N/A";

						};

						
					};



					if(ALP1>=ALP_max_number)
					{
						flag1=false;
					};
					//cout<<ALP1<<endl;


					double time0_current;
					FSA_utils::get_current_time(time0_current);

					//time limits
					if(time0_current-time0_start>=max_time_*mult1||time0_current-time0_start>=max_calculation_time_for_preliminary_stage)
					{
						flag1=false;
					};

					if(test_output1)
					{
						ff<<"\t"<<time0_current-time0_start-time_exp<<"\t"<<states_new1->d_total_calculation_time;
						ff<<"\n";
					};

					//deallocate and allocate memory
					if(states_opt!=states_old1)
					{
						delete_mult_states_type(states_old1);
					};


					states_old1=states_new1;
					states_new1=new mult_states_type;


					if(!flag1)
					{
						break;
					};

					//increase number of ALP
					ALP1++;

					//prediction of the time required for the calculation with the new number of ALP=ALP1
					long int ALP1_corrected=FSA_utils::Tmax(ALP1,ALP1_start+number_of_ALP_for_errors+1);
					double T_per_realization=(2*ALP1_corrected*ALP1_corrected*X1_X2*time_cell_ALP+(ALP1_corrected*(X1*L2+X2*L1)+L1_L2)*time_cell_exp);

					if(T_per_realization<=0)
					{
						throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
					};

					
					number_of_realizations1=(long int)FSA_utils::round(
						FSA_utils::Tmin(max_number_of_states,
						FSA_utils::Tmin(mult0*max_calculation_time_for_preliminary_stage,mult1*(max_time_-time_exp))
						/T_per_realization)
						);

					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					//number_of_realizations1=FSA_utils::Tmin(max_number_of_realizations0,number_of_realizations1);


					if(further_expanding_tmp1)
					{
						further_expanding_tmp1=false;
					};					
				};


				if(test_output1)
				{
					ff.close();
				};

				if(states_opt!=states_old1)
				{
					delete_mult_states_type(states_old1);
				};
				if(states_opt!=states_new1)
				{
					delete_mult_states_type(states_new1);
				};

				if(!states_opt)
				{
					throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
				};


				if(eps_lambda_>0)
				{
					if(error_opt<eps_lambda_)
					{
						accuracy_is_achieved_flag=true;
						//number of realizations is too large for the requested accuracy
						//decreasing number of realizations
						double coeff=error_opt/eps_lambda_;
						coeff=coeff*coeff;
						expected_number_of_realizations1_opt=(long int)FSA_utils::round(
							FSA_utils::Tmin(max_number_of_states,
							coeff*expected_number_of_realizations1_opt));
						expected_number_of_realizations1_opt=FSA_utils::Tmax(number_of_realizations0,expected_number_of_realizations1_opt);
					}
					else
					{
						expected_number_of_realizations1_opt*=8;
					};
				}
				else
				{
					expected_number_of_realizations1_opt*=8;
				};

				expected_number_of_realizations1_opt=(long int)FSA_utils::round(
					FSA_utils::Tmin((double)max_number_of_states,
					(double)expected_number_of_realizations1_opt));

				if(!library_call_flag_)
				{
					cout<<"\n--------------\n";
					cout<<"\nOptimal ALP# for the forward and reverse stages is\t"<<ALP_opt<<endl;
					cout<<"\n--------------\n";

					if(accuracy_is_achieved_flag)
					{
						cout<<"The input relative error for lambda can be acheived with this number of realizations\n";
					};
					cout<<endl;

				};


				double time_tmp2;
				FSA_utils::get_current_time(time_tmp2);

				double tmp2=(max_time_-(time_tmp2-time0_start)+states_opt->d_total_calculation_time)*fraction_of_ALP_opt
				//-states_opt->d_total_calculation_time;
				-total_time_without_extending_opt;

				if(tmp2<0)
				{
					tmp2=0;
				};

				double ending_time_tmp2=time_tmp2+tmp2;

				bool stopped_by_ending_time_flag_tmp2;
				long int number_of_realizations_with_ending_time_tmp2;
				double ending_time_to_test_logarithmic_regime2=-1;
				long int M_thr=-inf;

				//final calculation
				if(!library_call_flag_)
				{
					cout<<"Forward stage...\n\n";
				};
				save_states_flag=true;
				states_new1=NULL;
				states_new1=new mult_states_type;
				
				bool compute_allocated_memory_flag_tmp=false;
				double allocated_memory_tmp;

				collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
				compute_allocated_memory_flag_tmp,
				allocated_memory_tmp,


				expected_number_of_realizations1_opt,
				M_distr_file_name,
				*two_dim_layer_alignment_algorithm_test,
				*d_IS1_general_simulation,
				ALP_opt,//target ALP number
				ungappedlambda,
				limit2,
				number_of_sets,
				par,
				inside_simulation_flag,
				false,
				true&&forward_and_reverse_screen_output_flag_,

				ending_time_tmp2,
				stopped_by_ending_time_flag_tmp2,
				number_of_realizations_with_ending_time_tmp2,

				ending_time_to_test_logarithmic_regime2,

				save_states_flag,
				states_opt,
				states_new1,
				M_thr_estimation_direct,

				false,
				M_thr,
				IS1_cs,
				&eps_K,
				&number_of_cs_steps_for_expanding);

				number_of_realizations=expected_number_of_realizations1_opt;
				target_ALP=ALP_opt;

				if(stopped_by_ending_time_flag_tmp2)
				{
					if(!library_call_flag_)
					{
						cout<<"Number of realizations of the forward stage is\t"<<number_of_realizations_with_ending_time_tmp2<<endl<<endl;
						cout<<"Parameters estimation from the forward stage\n";
					};
					number_of_realizations=number_of_realizations_with_ending_time_tmp2;


					save_states_flag=false;
					ending_time_tmp2=-1;
					mult_states_type *states_new1_tmp=NULL;

					bool compute_allocated_memory_flag_tmp=false;
					double allocated_memory_tmp;


					collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
					compute_allocated_memory_flag_tmp,
					allocated_memory_tmp,

					number_of_realizations,
					M_distr_file_name,
					*two_dim_layer_alignment_algorithm_test,
					*d_IS1_general_simulation,
					ALP_opt,//target ALP number
					ungappedlambda,
					limit2,

					number_of_sets,
					par,
					inside_simulation_flag,
					forward_and_reverse_screen_output_flag_,
					false,

					ending_time_tmp2,
					stopped_by_ending_time_flag_tmp2,
					number_of_realizations_with_ending_time_tmp2,

					ending_time_to_test_logarithmic_regime2,

					save_states_flag,
					states_new1,
					states_new1_tmp,
					M_thr_estimation_direct,

					false,
					M_thr,
					IS1_cs,
					&eps_K,
					&number_of_cs_steps_for_expanding);
				};

				if(!inside_simulation_flag)
				{
					throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
				};


				delete_mult_states_type(states_new1);
				delete_mult_states_type(states_opt);



			};


			if(max_time_<=0)
			{
				if(!library_call_flag_)
				{
					cout<<"\nForward stage...\n\n";
				};

				double ending_time_tmp2=-1;
				bool stopped_by_ending_time_flag_tmp2;
				long int number_of_realizations_with_ending_time_tmp2;
				double ending_time_to_test_logarithmic_regime2=-1;
				long int M_thr=-inf;

				bool inside_simulation_flag=true;
				bool save_states_flag=false;
				mult_states_type *states_old1=NULL;
				mult_states_type *states_new1=NULL;

				bool compute_allocated_memory_flag_tmp=false;
				double allocated_memory_tmp;

				collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
				compute_allocated_memory_flag_tmp,
				allocated_memory_tmp,

				number_of_realizations,
				M_distr_file_name,
				*two_dim_layer_alignment_algorithm_test,
				*d_IS1_general_simulation,
				target_ALP,//target ALP number
				ungappedlambda,
				limit2,

				number_of_sets,
				par,
				inside_simulation_flag,
				forward_and_reverse_screen_output_flag_,
				true&&forward_and_reverse_screen_output_flag_,

				ending_time_tmp2,
				stopped_by_ending_time_flag_tmp2,
				number_of_realizations_with_ending_time_tmp2,

				ending_time_to_test_logarithmic_regime2,

				save_states_flag,
				states_old1,
				states_new1,
				M_thr_estimation_direct,

				futher_expanding_tmp,
				M_thr,
				IS1_cs,
				&eps_K,
				&number_of_cs_steps_for_expanding);


				if(!inside_simulation_flag)
				{
					throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
				};


			};

			double time0_end;
			FSA_utils::get_current_time(time0_end);

			if(!library_call_flag_)
			{
				cout<<"Time of the forward stage\t"<<time0_end-time0_start<<endl;
			};


			

			if(lambda_output_flag)
			{
				lambda_out<<par.lambda<<"\t";lambda_out<<par.lambda_error<<"\t";
				lambda_out<<par.C<<"\t";lambda_out<<par.C_error<<"\t";
				lambda_out<<par.a_I<<"\t";lambda_out<<par.a_I_error<<"\t";
				lambda_out<<par.a_J<<"\t";lambda_out<<par.a_J_error<<"\t";
				lambda_out<<par.sigma<<"\t";lambda_out<<par.sigma_error<<"\t";
				lambda_out<<par.alpha_I<<"\t";lambda_out<<par.alpha_I_error<<"\t";
				lambda_out<<par.alpha_J<<"\t";lambda_out<<par.alpha_J_error<<"\t";

				if(futher_expanding_tmp)
				{
					lambda_out<<par.K_C<<"\t";lambda_out<<par.K_C_error<<"\t";
				}
				else
				{
					lambda_out<<"#N/A\t";
				};

			};



		};
	};




//reversed FSA sampling
	string M_distr_file_name_reversed=lambda_file_name+"_M_reversed.out";
	data_for_FSA_alignment data_test_reversed;

	IS1_general_simulation *IS1_general_simulation_reversed=NULL;
	two_dim_layer_alignment_algorithm<long int> *two_dim_layer_alignment_algorithm_test_reversed=NULL;

	Sls::FALP_set_of_parameters par_reversed;


	if(FSA_flag&&futher_expanding)
	{
		long int initial_state=state_name_into_number["S1"];

		if(add_fake_state_flag)
		{
			initial_state=state_name_into_number["F"];
		};


		IS1_general_simulation_reversed=new IS1_general_simulation(
		IS1_reversed,
		initial_state,//initial state for the IS
		max_ind1*mult_margin+add_margin,//maximum sequence length
		max_ind2*mult_margin+3*add_margin);//maximum sequence length




		data_test_reversed=data_test;
		data_test_reversed.d_codon_AA=codon_reversed_AA;


		data_test_reversed.d_seq1=IS1_general_simulation_reversed->d_seq1;//sequence #1
		data_test_reversed.d_seq2=IS1_general_simulation_reversed->d_seq2;//sequence #2
	


		



		two_dim_layer_alignment_algorithm_test_reversed=
		new two_dim_layer_alignment_algorithm<long int>(
		number_of_variables_for_the_alignment,//total number of variables in dynamic equations
		depth1,//the maximum difference of the first index in the dynamic equations
		depth2,//the maximum difference of the second index in the dynamic equations
		max_ind1,//max of the index #1 (minimum index is 0 by default)
		max_ind2,//max of the index #2 (minimum index is 0 by default)
		0,//null element of T
		var_num);//d_variable_name_to_number[c] converts a letter 'c' (of type char) to order number of variable




		
		two_dim_layer_alignment_algorithm_test_reversed->d_two_dim_layer_alignment_function=test::two_dim_layer_alignment_function_FSA;
		two_dim_layer_alignment_algorithm_test_reversed->d_two_dim_layer_boundary_function=test::two_dim_layer_alignment_function_FSA;
		two_dim_layer_alignment_algorithm_test_reversed->d_par=&data_test_reversed;


		two_dim_layer_alignment_algorithm_test_reversed->d_M_flag=true;
		two_dim_layer_alignment_algorithm_test_reversed->d_E_flag=false;
		two_dim_layer_alignment_algorithm_test_reversed->d_scores_counts_flag=false;



		
		{
			two_dim_layer_alignment_algorithm_test_reversed->d_scores_counts_flag=true;
			two_dim_layer_alignment_algorithm_test_reversed->d_FSC_flag=true;

			if(!library_call_flag_)
			{
				cout<<"\n--------------\n";
				cout<<"\nReverse stage...\n\n";
			};

			double time1_start;
			FSA_utils::get_current_time(time1_start);

			long int limit2=max_ind2;
			long int M_thr_estimation_reverse;

			

			bool inside_simulation_flag=true;
			long int M_thr=M_thr_estimation_direct;
			if(!library_call_flag_)
			{
				cout<<"Score threshold for the reverse stage is \t"<<M_thr<<endl;
			};
			

			if(max_time_>0)
			{


				if(!accuracy_is_achieved_flag)
				{
					number_of_realizations*=8;
				};



				double ending_time_tmp3=time0_start+max_time_;
				bool stopped_by_ending_time_flag_tmp3;
				long int number_of_realizations_with_ending_time_tmp3;
				double ending_time_to_test_logarithmic_regime3=-1;
				


				bool save_states_flag=true;
				mult_states_type *states_old1=NULL;
				mult_states_type *states_new1=NULL;
				states_new1=new mult_states_type;

				bool compute_allocated_memory_flag_tmp=false;
				double allocated_memory_tmp;

				collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
				compute_allocated_memory_flag_tmp,
				allocated_memory_tmp,

				number_of_realizations,
				M_distr_file_name_reversed,
				*two_dim_layer_alignment_algorithm_test_reversed,
				*IS1_general_simulation_reversed,
				target_ALP,//target ALP number
				ungappedlambda,
				limit2,

				number_of_sets,
				par_reversed,
				inside_simulation_flag,
				//forward_and_reverse_screen_output_flag_,
				false,
				true&&forward_and_reverse_screen_output_flag_,

				ending_time_tmp3,
				stopped_by_ending_time_flag_tmp3,
				number_of_realizations_with_ending_time_tmp3,

				ending_time_to_test_logarithmic_regime3,

				save_states_flag,
				states_old1,
				states_new1,
				M_thr_estimation_reverse,

				futher_expanding,
				M_thr,
				IS1_cs,
				&eps_K,
				&number_of_cs_steps_for_expanding);

				


				if(stopped_by_ending_time_flag_tmp3)
				{
					if(!library_call_flag_)
					{
						cout<<"Number of realizations of the reverse stage is\t"<<number_of_realizations_with_ending_time_tmp3<<endl<<endl;
						cout<<"Parameters estimation from the reverse stage\n";
					};
					number_of_realizations=number_of_realizations_with_ending_time_tmp3;
					//number_of_realizations_reverse_final=number_of_realizations_with_ending_time_tmp3;

					save_states_flag=false;
					ending_time_tmp3=-1;
					mult_states_type *states_new1_tmp=NULL;

					bool compute_allocated_memory_flag_tmp=false;
					double allocated_memory_tmp;


					collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
					compute_allocated_memory_flag_tmp,
					allocated_memory_tmp,

					number_of_realizations,
					M_distr_file_name_reversed,
					*two_dim_layer_alignment_algorithm_test_reversed,
					*IS1_general_simulation_reversed,
					target_ALP,//target ALP number
					ungappedlambda,
					limit2,

					number_of_sets,
					par_reversed,
					inside_simulation_flag,
					forward_and_reverse_screen_output_flag_,
					false,

					ending_time_tmp3,
					stopped_by_ending_time_flag_tmp3,
					number_of_realizations_with_ending_time_tmp3,

					ending_time_to_test_logarithmic_regime3,

					save_states_flag,
					states_new1,
					states_new1_tmp,
					M_thr_estimation_reverse,

					futher_expanding,
					M_thr,
					IS1_cs,
					&eps_K,
					&number_of_cs_steps_for_expanding);
				};

				delete_mult_states_type(states_new1);

			}
			else
			{
				double ending_time_tmp3=-1;
				bool stopped_by_ending_time_flag_tmp3;
				long int number_of_realizations_with_ending_time_tmp3;
				double ending_time_to_test_logarithmic_regime3=-1;

				bool save_states_flag=false;
				mult_states_type *states_old1=NULL;
				mult_states_type *states_new1=NULL;

				bool compute_allocated_memory_flag_tmp=false;
				double allocated_memory_tmp;

				collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
				compute_allocated_memory_flag_tmp,
				allocated_memory_tmp,

				number_of_realizations,
				M_distr_file_name_reversed,
				*two_dim_layer_alignment_algorithm_test_reversed,
				*IS1_general_simulation_reversed,
				target_ALP,//target ALP number
				ungappedlambda,
				limit2,

				number_of_sets,
				par_reversed,
				inside_simulation_flag,
				forward_and_reverse_screen_output_flag_,
				true&&forward_and_reverse_screen_output_flag_,

				ending_time_tmp3,
				stopped_by_ending_time_flag_tmp3,
				number_of_realizations_with_ending_time_tmp3,

				ending_time_to_test_logarithmic_regime3,

				save_states_flag,
				states_old1,
				states_new1,
				M_thr_estimation_reverse,

				futher_expanding,
				M_thr,
				IS1_cs,
				&eps_K,
				&number_of_cs_steps_for_expanding);

				

			};

			if(!inside_simulation_flag)
			{
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			};

			double time1_end;
			FSA_utils::get_current_time(time1_end);

			if(!library_call_flag_)
			{
				cout<<"Time of the reverse stage\t"<<time1_end-time1_start<<endl;
			};


			if(lambda_output_flag)
			{
				lambda_out<<par_reversed.lambda<<"\t";lambda_out<<par_reversed.lambda_error<<"\t";
				lambda_out<<par_reversed.C<<"\t";lambda_out<<par_reversed.C_error<<"\t";
				lambda_out<<par_reversed.a_I<<"\t";lambda_out<<par_reversed.a_I_error<<"\t";
				lambda_out<<par_reversed.a_J<<"\t";lambda_out<<par_reversed.a_J_error<<"\t";
				lambda_out<<par_reversed.sigma<<"\t";lambda_out<<par_reversed.sigma_error<<"\t";
				lambda_out<<par_reversed.alpha_I<<"\t";lambda_out<<par_reversed.alpha_I_error<<"\t";
				lambda_out<<par_reversed.alpha_J<<"\t";lambda_out<<par_reversed.alpha_J_error<<"\t";

				if(futher_expanding)
				{
					lambda_out<<par_reversed.K_C<<"\t";lambda_out<<par_reversed.K_C_error<<"\t";
				}
				else
				{
					lambda_out<<"#N/A\t";
				};


			};


		};
	};








	if(FSA_flag&&futher_expanding)
	{
		combine_parameters_from_forward_and_reversed_calculations_generalized(
		par,//parameters from forward calculation
		par_reversed,//parameters from reversed calculation
		par_result_);

		//assign gap penalties for the FSC
		par_result_.G=FSA_utils::Tmin(data_test.d_open1+data_test.d_epen1,
			data_test.d_open2+data_test.d_epen2,
			data_test.d_gamma);

		par_result_.G1=par_result_.G;
		par_result_.G2=par_result_.G;
	};
//reversed FSA sampling - end





	delete[]var_num;


	delete[]RR1;
	delete[]RR2;

	delete[]RR1_sum;
	delete[]RR2_sum;

	delete[]RR1_sum_elements;
	delete[]RR2_sum_elements;

	FSA_utils::delete_memory_for_matrix(alphabet_letters_number1,smatr);
	FSA_utils::delete_memory_for_matrix(number_of_states,transition_probabilities);

	if(futher_expanding&&FSA_flag)
	{
		
		delete two_dim_layer_alignment_algorithm_test_reversed;
		delete IS1_general_simulation_reversed;
		delete IS1_reversed;

		for(s=0;s<number_of_states;s++)
		{
			IS1_general::deallocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
				states_distr_reversed[s],
				alphabet_letters_number1,//number of letters in the sequence #1
				states_description[s]);//state description
		};

		delete[]states_distr_reversed;

		delete[]codon_reversed_AA;

	};

	delete[]codon_AA;

	//crude sampling
	if(cs_flag)
	{

		FSA_utils::delete_memory_for_matrix(number_of_states_cs,transition_probabilities_cs);
		
		delete two_dim_layer_alignment_algorithm_test_cs;
		delete IS1_general_simulation_cs;
		delete IS1_cs;

		long int s;
		for(s=0;s<number_of_states_cs;s++)
		{
			IS1_general::deallocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
				states_distr_cs[s],
				alphabet_letters_number1,//number of letters in the sequence #1
				states_description_cs[s]);//state description
		};

		delete[]states_description_cs;

	};
	//crude sampling - end


	if(FSA_flag)
	{
		delete two_dim_layer_alignment_algorithm_test;
		delete d_IS1_general_simulation;
		delete d_IS1;


		long int s;
		for(s=0;s<number_of_states;s++)
		{
			IS1_general::deallocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
				states_distr[s],
				alphabet_letters_number1,//number of letters in the sequence #1
				states_description[s]);//state description
		};

		delete[]states_description;
	};

	delete[]RR1_AA;



	if(lambda_output_flag)
	{
		lambda_out<<endl;
		lambda_out.close();
	};

	//precompute intercepts
	par_result_.d_params_flag=true;
	FALP_pvalues::compute_intercepts(par_result_);

	double seconds2;
	FSA_utils::get_current_time(seconds2);

	par_result_.m_CalcTime=seconds2-seconds1;


}

//-------------------------------------------------------------------
void test::collect_and_output_M_distr_upto_fixed_lengths(
long int number_of_realizations_,
string M_distr_file_name_,
two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_,
IS1_general_simulation &IS1_general_simulation_test_,
long int target_seq1_length_,//target length of the sequence #1
long int target_seq2_length_)//target length of the sequence #2
{
	//for test - begin
	bool test_flag=false;

		
	string input_file_name="sss.out";

	long int number_of_letters1;//number of letters for the sequence 1
	long int number_of_letters2;//number of letters for the sequence 2

	char *alphabet1;//alphabet letters for the sequence #1
	char *alphabet2;//alphabet letters for the sequence #2

	long int number_of_sequences;

	string *headers;
	long int *lengths1;//lengths of the sequences #1
	long int *lengths2;//lengths of the sequences #2
	long int **sequences1=NULL;//the first index numerates sequences; the second - sequence letters
	long int **sequences2=NULL;

	if(test_flag)
	{

		long int *alphabet1_to_long=NULL;//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
		long int *alphabet2_to_long=NULL;//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

		long int codon_length;//codon length 
		long int *codon_AA=NULL;//<codon code,AA number>

		string DNA_codon_table_file_name="test_input_files/codon_AA.in";//a name of a file with DNA codon table




		FSA_utils::read_codon_AA_file(
		DNA_codon_table_file_name,
		number_of_letters1,//number of letters for the sequence 1
		number_of_letters2,//number of letters for the sequence 2

		alphabet1,//alphabet letters for the sequence #1
		alphabet2,//alphabet letters for the sequence #2

		alphabet1_to_long,//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
		alphabet2_to_long,//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

		codon_length,//codon length 
		codon_AA);//<codon code,AA number>

		FSA_utils::read_sequences_for_alingment(

			input_file_name,

			number_of_letters1,//number of letters for the sequence 1
			number_of_letters2,//number of letters for the sequence 2

			alphabet1,//alphabet letters for the sequence #1
			alphabet2,//alphabet letters for the sequence #2

			number_of_sequences,

			headers,
			lengths1,//lengths of the sequences #1
			lengths2,//lengths of the sequences #2
			sequences1,//the first index numerates sequences; the second - sequence letters
			sequences2);
	};

	//for test - end

	if(!two_dim_layer_alignment_algorithm_test_.d_M_flag)
	{
		throw error("Error - two_dim_layer_alignment_algorithm_test_.d_M_flag must be true in test::collect_and_output_M_distr_upto_fixed_lengths\n",1);
	};

	ofstream fM;

	array_v<double> *distrM=NULL;
	array_v<double> *distrM_error=NULL;
	long int score_max=-inf;

	//distribution of M

	fM.open(M_distr_file_name_.data());
	if(!fM)
	{
		throw error("Error - the file "+M_distr_file_name_+" is not found\n",3);
	};


	distrM=new array_v<double>(NULL);
	distrM_error=new array_v<double>(NULL);



	//--------------------

	long int k;
	for(k=1;k<=number_of_realizations_;k++)
	{

		IS1_general_simulation_test_.init();
		IS1_general_simulation_test_.simulate_upto_target_lengths(
		target_seq1_length_,
		target_seq2_length_);


		/*
		double weight;
		if(0)
		{
			IS1_general_simulation_test_.calculate_weight_W1_for_fixed_lengths_using_infinite_sums(
			target_seq1_length_,//target length of the sequence #1
			target_seq2_length_,//target length of the sequence #2
			weight);//the resulted weight
		};
		*/

		//weights calculation
		IS1_general_simulation_test_.d_W1_seq1_current_length=-1;
		IS1_general_simulation_test_.d_W1_seq2_current_length=-1;

		double weight2;

		{
			IS1_general_simulation_test_.calculate_weight_W1_for_fixed_lengths_using_recursions(
			target_seq1_length_,//target length of the sequence #1
			target_seq2_length_,//target length of the sequence #2
			weight2);//the resulted weight
		};




		two_dim_layer_alignment_algorithm_test_.init();

		if(test_flag)
		{
			data_for_FSA_alignment &data=*(data_for_FSA_alignment*)two_dim_layer_alignment_algorithm_test_.d_par;

			data.d_seq1=sequences1[k-1];
			data.d_seq2=sequences2[k-1];
		};

		two_dim_layer_alignment_algorithm_test_.align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
		//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
		target_seq1_length_,//target length of the sequence #1
		target_seq2_length_);//target length of the sequence #2

		score_max=FSA_utils::Tmax(score_max,two_dim_layer_alignment_algorithm_test_.d_M);
		double w_tmp=1/weight2;
		distrM->increase_elem_by_x(two_dim_layer_alignment_algorithm_test_.d_M,w_tmp);
		distrM_error->increase_elem_by_x(two_dim_layer_alignment_algorithm_test_.d_M,w_tmp*w_tmp);

		if(test_flag)
		{
			cout<<two_dim_layer_alignment_algorithm_test_.d_M<<endl;
		};
		

		if(k%100==0)
		{
			cout<<k<<endl;
		};

	};

	long int s;

	//convert into tail
	for(s=score_max-1;s>=0;s--)
	{
		distrM->d_elem[s-distrM->d_ind0]+=distrM->d_elem[s+1-distrM->d_ind0];
		distrM_error->d_elem[s-distrM->d_ind0]+=distrM_error->d_elem[s+1-distrM->d_ind0];
	};

	for(s=0;s<=score_max;s++)
	{
		distrM->d_elem[s-distrM->d_ind0]/=(double)number_of_realizations_;
		distrM_error->d_elem[s-distrM->d_ind0]/=(double)number_of_realizations_;
		distrM_error->d_elem[s-distrM->d_ind0]-=distrM->d_elem[s-distrM->d_ind0]*distrM->d_elem[s-distrM->d_ind0];
		distrM_error->d_elem[s-distrM->d_ind0]/=(double)number_of_realizations_;
		distrM_error->d_elem[s-distrM->d_ind0]=FSA_utils::sqrt_plus(distrM_error->d_elem[s-distrM->d_ind0]);

	};



	fM<<score_max+1<<endl;
	for(s=0;s<=score_max;s++)
	{
		fM<<s<<"\t"<<distrM->d_elem[s-distrM->d_ind0]<<"\t"<<distrM_error->d_elem[s-distrM->d_ind0]<<endl;
	};

	fM.close();


	delete distrM;
	delete distrM_error;


}

//===============================
void test::compare_direct_and_reverse_sampling(
long int number_of_realizations_,
long int seq2_length_,
two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_,
IS1_general_simulation &IS1_general_simulation_test_,
two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_reverse_,
IS1_general_simulation &IS1_general_simulation_test_reverse_,
IS1_general *IS1_general_cs_)
{

	cout<<endl;

	cout<<"Starting...\n";

	long int k;
	for(k=1;k<=number_of_realizations_;k++)
	{

		
		

		IS1_general_simulation_test_.init();

		IS1_general_simulation_test_.d_W1_seq1_current_length=-1;
		IS1_general_simulation_test_.d_W1_seq2_current_length=-1;

		two_dim_layer_alignment_algorithm_test_.init();



		long int seq1_length=0;
		long int seq2_length=0;

		long int current_state=0;//the current state
		long int new_state;

		long int M_direct;
		long int M_reverse;

		long int j;
		for(j=1;j<=seq2_length_;j++)
		{
			IS1_general_cs_->one_step_of_the_importance_samping(//an initial state when sequences are empty must be defined outside
			current_state,//the current state
			IS1_general_simulation_test_.d_seq1+seq1_length,//letters for the sequence #1; the array must be allocated
			IS1_general_simulation_test_.d_seq2+seq2_length,//letters for the sequence #2; the array must be allocated
			new_state);//a new state

			seq1_length+=two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth1;
			seq2_length+=two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth2;

			two_dim_layer_alignment_algorithm_test_.align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
			//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
			seq1_length,//target length of the sequence #1
			seq2_length);//target length of the sequence #2
		};

		M_direct=two_dim_layer_alignment_algorithm_test_.d_M;

//--------
		IS1_general_simulation_test_reverse_.init();

		IS1_general_simulation_test_reverse_.d_W1_seq1_current_length=-1;
		IS1_general_simulation_test_reverse_.d_W1_seq2_current_length=-1;

		two_dim_layer_alignment_algorithm_test_reverse_.init();
//--------

		long int t1=seq1_length;
		long int t2=seq2_length;

		for(j=1;j<=seq2_length_;j++)
		{
			
			long int ii;
			for(ii=0;ii<two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth1;ii++)
			{
				seq1_length--;
				IS1_general_simulation_test_reverse_.d_seq1[t1-1-seq1_length]=
					IS1_general_simulation_test_.d_seq1[seq1_length];

			};


			for(ii=0;ii<two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth2;ii++)
			{
				seq2_length--;
				IS1_general_simulation_test_reverse_.d_seq2[t2-1-seq2_length]=
					IS1_general_simulation_test_.d_seq2[seq2_length];

			};
			



			two_dim_layer_alignment_algorithm_test_reverse_.align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
			//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
			j*two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth1,//target length of the sequence #1
			j*two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth2);//target length of the sequence #2

			
		};

		M_reverse=two_dim_layer_alignment_algorithm_test_reverse_.d_M;

		

		if(M_direct!=M_reverse)
		{
			cout<<M_direct<<"\t"<<M_reverse<<endl;
		};


		if(k%100==0)
		{
			cout<<k<<endl;
		};

	};

	cout<<"Finished!\n";

}

//-------------------------------------------------------------------

long int test::realization_number_into_set_number(//numeration of sets starts from 1
long int &realization_number_,//realization order number; the numeration starts from 0
long int &number_of_realizations_set_)//total number of sets
{
	return (realization_number_/number_of_realizations_set_)+1;
}

pair<long int,long int> test::set_number_boundaries(//boundaries of realization order numbers for the set; numeration of realizations starts from 0
long int &set_number_,//set order number; the numeration starts from 1
long int &number_of_realizations_set_)//total number of realizations per set
{
	pair<long int,long int> res;
	res.first=(set_number_-1)*number_of_realizations_set_;
	res.second=res.first+number_of_realizations_set_-1;
	return res;
}

void test::collect_and_output_E_distr_upto_fixed_ALP(
long int number_of_realizations_,
string E_distr_file_name_,
two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_,
IS1_general_simulation &IS1_general_simulation_test_,
long int target_ALP_,//target ALP number
double ungapped_lambda_,
long int limit2_,
long int number_of_sets_,//number of sets for error calculation
Sls::FALP_set_of_parameters &par_,
bool screen_output_flag_,

bool futher_expanding_,
IS1_general *IS1_general_cs_,
double *eps_K_,//relative error for K
long int *number_of_cs_steps_for_expanding_)
{

	

	double &lambda_out_=par_.lambda;
	double &lambda_out_error_=par_.lambda_error;
	double &C_=par_.C;
	double &C_error_=par_.C_error;

	double &a_I_=par_.a_I;
	double &a_I_error_=par_.a_I_error;
	double &a_J_=par_.a_J;
	double &a_J_error_=par_.a_J_error;
	double &sigma_=par_.sigma;
	double &sigma_error_=par_.sigma_error;
	double &alpha_I_=par_.alpha_I;
	double &alpha_I_error_=par_.alpha_I_error;
	double &alpha_J_=par_.alpha_J;
	double &alpha_J_error_=par_.alpha_J_error;

	double &K_C_=par_.K_C;
	double &K_C_error_=par_.K_C_error;



	bool average_set_results_flag=false;

	{
		long int number_of_realizations_tmp=(number_of_realizations_/number_of_sets_)*number_of_sets_;
		if(number_of_realizations_tmp<number_of_realizations_)
		{
			number_of_realizations_=number_of_realizations_tmp+number_of_sets_;
		};
	};

	long int number_of_realizations_set=number_of_realizations_/number_of_sets_;

	//------------------------------------------------------

	bool output_flag=false;

	if(!two_dim_layer_alignment_algorithm_test_.d_M_flag)
	{
		throw error("Error - two_dim_layer_alignment_algorithm_test_.d_E_flag must be true in test::collect_and_output_E_distr_upto_fixed_ALP\n",1);
	};

	ofstream fE;
	ofstream fE_summary;



	array_v<double> ***distrE=new array_v<double> **[number_of_sets_+1];
	FSA_utils::assert_mem(distrE);
	array_v<double> ***distrE_errors=new array_v<double> **[number_of_sets_+1];
	FSA_utils::assert_mem(distrE_errors);

	long int **score_max=new long int*[number_of_sets_+1];//statistic
	FSA_utils::assert_mem(score_max);

	long int k1;
	for(k1=0;k1<=number_of_sets_;k1++)
	{
		distrE[k1]=new array_v<double>*[target_ALP_+1];
		FSA_utils::assert_mem(distrE[k1]);
		distrE_errors[k1]=new array_v<double>*[target_ALP_+1];
		FSA_utils::assert_mem(distrE_errors[k1]);
		score_max[k1]=new long int[target_ALP_+1];
		FSA_utils::assert_mem(score_max[k1]);

		long int k;
		for(k=0;k<=target_ALP_;k++)
		{
			distrE[k1][k]=new array_v<double>(NULL);
			FSA_utils::assert_mem(distrE[k1][k]);

			distrE_errors[k1][k]=new array_v<double>(NULL);
			FSA_utils::assert_mem(distrE_errors[k1][k]);

			score_max[k1][k]=-inf;
		};

	};

	//distribution of E
	if(output_flag)
	{
		fE.open(E_distr_file_name_.data());
		if(!fE)
		{
			throw error("Error - the file "+E_distr_file_name_+" is not found\n",3);
		};

		string E_distr_file_name_summary=E_distr_file_name_+"_summary";
		fE_summary.open(E_distr_file_name_summary.data());
		if(!fE_summary)
		{
			throw error("Error - the file "+E_distr_file_name_summary+" is not found\n",3);
		};
	};




	long int k;


	array_positive<double> *diff=new array_positive<double>[number_of_sets_+1];//statistic
	array_positive<double> *diff_errors=new array_positive<double>[number_of_sets_+1];//statistic
	


	array_positive<long int> *distance_along_direction_1=new array_positive<long int>[number_of_realizations_];//statistic
	array_positive<long int> *distance_along_direction_2=new array_positive<long int>[number_of_realizations_];//statistic

	array_positive<double> *ALP_weight=new array_positive<double>[number_of_realizations_];//statistic
	array_positive<long int> *ALP_edge_max=new array_positive<long int>[number_of_realizations_];//statistic

	//--------------------

	

	double *sum_of_weights=new double[number_of_sets_+1];//statistic;
	double *sum_of_weights_error=new double[number_of_sets_+1];//statistic;

	for(k1=0;k1<=number_of_sets_;k1++)
	{
		sum_of_weights[k1]=0;//statistic
		sum_of_weights_error[k1]=0;//statistic
	};

	long int M_threshold=(long int)ceil(log((*eps_K_)*(1.0-ungapped_lambda_))/ungapped_lambda_);
	M_threshold*=2;

	if(M_threshold>=0)
	{
		M_threshold=-1;
	};

	two_dim_layer_alignment_algorithm_test_.d_M_threshold=M_threshold;
	
	long int M_max=0;

	for(k=0;k<number_of_realizations_;k++)
	{
		long int sn=realization_number_into_set_number(k,number_of_realizations_set);

		
		

		IS1_general_simulation_test_.init();

		IS1_general_simulation_test_.d_W1_seq1_current_length=-1;
		IS1_general_simulation_test_.d_W1_seq2_current_length=-1;

		two_dim_layer_alignment_algorithm_test_.init();

		long int current_ALP_number=0;
		long int current_M=0;
		distrE[0][current_ALP_number]->increase_elem_by_x(current_M,1.0);

		distrE[sn][current_ALP_number]->increase_elem_by_x(current_M,1.0);
		
		distance_along_direction_1[k].set_elem(current_ALP_number,0);
		distance_along_direction_2[k].set_elem(current_ALP_number,0);

		ALP_weight[k].set_elem(current_ALP_number,1.0);
		ALP_edge_max[k].set_elem(current_ALP_number,0);

		score_max[0][current_ALP_number]=current_M;
		score_max[sn][current_ALP_number]=current_M;

		long int seq1_length=0;
		long int seq2_length=0;


		//weights calculation
		double weight2;
		double tmp_w=1;
		double tmp_w2=1;


		while(current_ALP_number<target_ALP_)
		{

			seq1_length+=two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth1;
			seq2_length+=two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth2;

			if(seq1_length>two_dim_layer_alignment_algorithm_test_.d_max_ind1||
				seq2_length>two_dim_layer_alignment_algorithm_test_.d_max_ind2)
			{
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			};


			IS1_general_simulation_test_.simulate_upto_target_lengths(
			seq1_length,
			seq2_length);


			IS1_general_simulation_test_.calculate_weight_W1_upto_target_lengths(
			seq1_length,//target length of the sequence #1
			seq2_length,//target length of the sequence #2
			seq1_length,//the weights are calculated upto this length for the sequence #1
			seq2_length);//the weights are calculated upto this length for the sequence #2


			two_dim_layer_alignment_algorithm_test_.align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
			//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
			seq1_length,//target length of the sequence #1
			seq2_length);//target length of the sequence #2


			if(two_dim_layer_alignment_algorithm_test_.d_M>current_M)
			{//new ALP found

				current_M=two_dim_layer_alignment_algorithm_test_.d_M;
				current_ALP_number++;



				IS1_general_simulation_test_.calculate_weight_W1_for_fixed_lengths_using_recursions(
				seq1_length,//target length of the sequence #1
				seq2_length,//target length of the sequence #2
				weight2);//the resulted weight



				score_max[0][current_ALP_number]=FSA_utils::Tmax(score_max[0][current_ALP_number],current_M);
				score_max[sn][current_ALP_number]=FSA_utils::Tmax(score_max[sn][current_ALP_number],current_M);

				tmp_w=1/weight2;
				tmp_w2=tmp_w*tmp_w;

				

				distrE[0][current_ALP_number]->increase_elem_by_x(current_M,tmp_w);
				distrE[sn][current_ALP_number]->increase_elem_by_x(current_M,tmp_w);
				
				distrE_errors[0][current_ALP_number]->increase_elem_by_x(current_M,tmp_w2);
				distrE_errors[sn][current_ALP_number]->increase_elem_by_x(current_M,tmp_w2);

				distance_along_direction_1[k].set_elem(current_ALP_number,two_dim_layer_alignment_algorithm_test_.d_distance_along_direction_1);
				distance_along_direction_2[k].set_elem(current_ALP_number,two_dim_layer_alignment_algorithm_test_.d_distance_along_direction_2);

				ALP_weight[k].set_elem(current_ALP_number,tmp_w);
				ALP_edge_max[k].set_elem(current_ALP_number,current_M);

				M_max=FSA_utils::Tmax(current_M,M_max);

		
			};

			if(seq2_length>limit2_)
			{
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			};





		};

	

		if((k+1)%100==0)
		{
			cout<<k+1<<endl;
		};

		//expanding
		if(futher_expanding_)
		{
			long int current_state=0;//the current state
			long int new_state;

			long int j;
			for(j=1;j<=*number_of_cs_steps_for_expanding_;j++)
			{
				IS1_general_cs_->one_step_of_the_importance_samping(//an initial state when sequences are empty must be defined outside
				current_state,//the current state
				IS1_general_simulation_test_.d_seq1+seq1_length,//letters for the sequence #1; the array must be allocated
				IS1_general_simulation_test_.d_seq2+seq2_length,//letters for the sequence #2; the array must be allocated
				new_state);//a new state

				seq1_length+=two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth1;
				seq2_length+=two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth2;

				two_dim_layer_alignment_algorithm_test_.align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
				//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
				seq1_length,//target length of the sequence #1
				seq2_length);//target length of the sequence #2

				long int &M=two_dim_layer_alignment_algorithm_test_.d_M;
				long int &E=two_dim_layer_alignment_algorithm_test_.d_E;

				if(E-M<M_threshold)
				{
					break;
				};


			};

			IS1_general_simulation_test_.d_seq1_current_length=seq1_length;
			IS1_general_simulation_test_.d_seq2_current_length=seq2_length;


		};
		


		//collect differences
		if(two_dim_layer_alignment_algorithm_test_.d_scores_counts_flag&&futher_expanding_)
		{
			//sum of weigths
			sum_of_weights[0]+=tmp_w;
			sum_of_weights[sn]+=tmp_w;
			sum_of_weights_error[0]+=tmp_w2;
			sum_of_weights_error[sn]+=tmp_w2;

			long int &M=two_dim_layer_alignment_algorithm_test_.d_M;
			long int ii;
			for(ii=FSA_utils::Tmax(M_threshold,two_dim_layer_alignment_algorithm_test_.d_scores_counts.d_ind0);ii<=FSA_utils::Tmin(M,two_dim_layer_alignment_algorithm_test_.d_scores_counts.d_dim_plus_d_ind0);ii++)
			{
				long int diff_tmp=(M-ii);
				double diff_tmp_weight=two_dim_layer_alignment_algorithm_test_.d_scores_counts.get_elem(ii)*tmp_w;
				double diff_tmp_weight2=diff_tmp_weight*diff_tmp_weight;

				diff[0].increase_elem_by_x(diff_tmp,diff_tmp_weight);
				diff[sn].increase_elem_by_x(diff_tmp,diff_tmp_weight);

				diff_errors[0].increase_elem_by_x(diff_tmp,diff_tmp_weight2);
				diff_errors[sn].increase_elem_by_x(diff_tmp,diff_tmp_weight2);
			};

		};


	};


	
	for(k1=0;k1<=number_of_sets_;k1++)
	{
		long int number_of_realizations_tmp=number_of_realizations_set;
		if(k1==0)
		{
			number_of_realizations_tmp=number_of_realizations_;
		};

		for(k=0;k<=target_ALP_;k++)
		{
			long int j;
			for(j=distrE[k1][k]->d_ind0;j<=score_max[k1][k];j++)
			{
				distrE[k1][k]->divide_elem_by_x(j,(double)number_of_realizations_tmp);
				distrE_errors[k1][k]->divide_elem_by_x(j,(double)number_of_realizations_tmp);
				distrE_errors[k1][k]->increase_elem_by_x(j,-distrE[k1][k]->get_elem(j)*distrE[k1][k]->get_elem(j));
				distrE_errors[k1][k]->divide_elem_by_x(j,(double)number_of_realizations_tmp);
			};
		};
	};

	double *lambda_out_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(lambda_out_vect);

	double *lambda_out_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(lambda_out_error_vect);


	double *C_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(C_vect);

	double *C_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(C_error_vect);


	{

		long int nalp=target_ALP_;

		for(k1=0;k1<=number_of_sets_;k1++)
		{
			bool check_the_criteria=false;
			long int nalp_thr;
			bool inside_simulation_flag;
			void **alp_distr=(void**)distrE[k1];
			void **alp_distr_errors=(void**)distrE_errors[k1];
			double ungapped_lambda=ungapped_lambda_;
			double lambda;
			double lambda_error;

			fsa_par::calculate_lambda(
			check_the_criteria,
			nalp,
			nalp_thr,
			inside_simulation_flag,
			alp_distr,
			alp_distr_errors,
			ungapped_lambda,
			lambda,
			lambda_error);

			if(!inside_simulation_flag)
			{
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			};


			lambda_out_vect[k1]=lambda;
			lambda_out_error_vect[k1]=lambda_error;
		};

		
		lambda_out_=lambda_out_vect[0];
		lambda_out_error_=lambda_out_error_vect[0];

		if(screen_output_flag_)
		{
			cout<<endl;
		};


		if(screen_output_flag_)
		{
			cout<<"Lambda\t"<<lambda_out_;
		};

		if(number_of_sets_>1)
		{
			lambda_out_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			lambda_out_vect+1);

			if(average_set_results_flag)
			{
				lambda_out_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				lambda_out_vect+1);


				if(screen_output_flag_)
				{
					cout<<"\t"<<lambda_out_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<lambda_out_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_LambdaSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_LambdaSbs.push_back(lambda_out_vect[k1]);
		};


		for(k1=0;k1<=number_of_sets_;k1++)
		{

			void **alp_distr=(void**)distrE[k1];
			void **alp_distr_errors=(void**)distrE_errors[k1];

			long int starting_point=0;
			bool inside_simulation_flag;

			fsa_par::calculate_C(
			starting_point,
			nalp,
			alp_distr,
			alp_distr_errors,
			lambda_out_vect[k1],
			lambda_out_error_vect[k1],
			C_vect[k1],
			C_error_vect[k1],
			inside_simulation_flag);

			if(!inside_simulation_flag)
			{
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			};
		};

		C_=C_vect[0];
		C_error_=C_error_vect[0];

		if(screen_output_flag_)
		{
			cout<<"C\t"<<C_;
		};

		if(number_of_sets_>1)
		{
			C_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			C_vect+1);

			if(average_set_results_flag)
			{
				C_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				C_vect+1);

				if(screen_output_flag_)
				{
					cout<<"\t"<<C_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<C_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_CSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_CSbs.push_back(C_vect[k1]);
		};

	};



	for(k=0;k<=target_ALP_;k++)
	{

		if(output_flag)
		{
			long int s;

			double sum_tmp=0;
			double sum_tmp_error=0;

			fE<<k<<endl;
			fE<<score_max[0][k]+1<<endl;

			for(s=0;s<=score_max[0][k];s++)
			{
				fE<<s<<"\t"<<distrE[0][k]->d_elem[s-distrE[0][k]->d_ind0]<<"\t"<<FSA_utils::sqrt_plus(distrE_errors[0][k]->d_elem[s-distrE_errors[0][k]->d_ind0])<<endl;
				sum_tmp+=distrE[0][k]->d_elem[s-distrE[0][k]->d_ind0]*exp((double)s*ungapped_lambda_);
				sum_tmp_error+=FSA_utils::sqrt_plus(distrE_errors[0][k]->d_elem[s-distrE_errors[0][k]->d_ind0])*exp((double)s*ungapped_lambda_);
			};
			fE<<endl;

			fE_summary<<k<<"\t"<<sum_tmp<<"\t"<<sum_tmp_error<<endl;
		};
	};


	double *K_C_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(K_C_vect);

	double *K_C_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(K_C_error_vect);


	//calculation of the ratio K/C
	if(two_dim_layer_alignment_algorithm_test_.d_scores_counts_flag&&futher_expanding_)
	{
		for(k1=0;k1<=number_of_sets_;k1++)
		{
			long int number_of_realizations_tmp=number_of_realizations_set;
			if(k1==0)
			{
				number_of_realizations_tmp=number_of_realizations_;
			};

			long int ii;
			for(ii=0;ii<=diff[k1].d_dim;ii++)
			{
				diff[k1].divide_elem_by_x(ii,(double)number_of_realizations_tmp);
				diff_errors[k1].divide_elem_by_x(ii,(double)number_of_realizations_tmp);
				diff_errors[k1].increase_elem_by_x(ii,-diff[k1].get_elem(ii)*diff[k1].get_elem(ii));
				diff_errors[k1].divide_elem_by_x(ii,(double)number_of_realizations_tmp);
				double tmp=FSA_utils::sqrt_plus(diff_errors[k1].get_elem(ii));
				diff_errors[k1].set_elem(ii,tmp);
			};

			

			sum_of_weights[k1]/=(double)number_of_realizations_tmp;
			sum_of_weights_error[k1]/=(double)number_of_realizations_tmp;
			sum_of_weights_error[k1]-=sum_of_weights[k1]*sum_of_weights[k1];
			sum_of_weights_error[k1]/=(double)number_of_realizations_tmp;
			sum_of_weights_error[k1]=FSA_utils::sqrt_plus(sum_of_weights_error[k1]);
		

			double den=0;
			double den_error=0;
			for(ii=0;ii<=diff[k1].d_dim;ii++)
			{
				double tmp=exp(-lambda_out_vect[k1]*(double)ii);
				den+=tmp*diff[k1].get_elem(ii);
				den_error+=tmp*tmp*diff_errors[k1].get_elem(ii)*diff_errors[k1].get_elem(ii);

			};


			den_error=FSA_utils::sqrt_plus(den_error);

			K_C_vect[k1]=sum_of_weights[k1]/den;
			K_C_error_vect[k1]=FSA_utils::error_of_the_ratio(sum_of_weights[k1],sum_of_weights_error[k1],den,den_error);

		};

		K_C_=K_C_vect[0];
		K_C_error_=K_C_error_vect[0];


		if(output_flag)
		{
			fE<<endl;
			fE<<"K/C distributions\n";
			fE<<diff[0].d_dim+1<<endl;
			long int ii;
			for(ii=0;ii<=diff[0].d_dim;ii++)
			{
				//double tmp=exp(-lambda_out_*(double)ii);
				fE<<ii<<"\t"<<diff[0].get_elem(ii)/sum_of_weights[0]<<"\t"<<diff_errors[0].get_elem(ii)/sum_of_weights[0]<<endl;
			};

		};

		if(screen_output_flag_)
		{
			cout<<"K/C\t"<<K_C_;
		};

		if(number_of_sets_>1)
		{
			K_C_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			K_C_vect+1);

			if(average_set_results_flag)
			{
				K_C_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				K_C_vect+1);

				if(screen_output_flag_)
				{
					cout<<"\t"<<K_C_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<K_C_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_K_CSbs.clear();
		par_.m_KSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_K_CSbs.push_back(K_C_vect[k1]);
			par_.m_KSbs.push_back(K_C_vect[k1]*C_vect[k1]);
		};




	};


	double *a_I_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(a_I_vect);
	double *a_I_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(a_I_error_vect);
	double *a_J_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(a_J_vect);
	double *a_J_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(a_J_error_vect);
	double *sigma_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(sigma_vect);
	double *sigma_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(sigma_error_vect);
	double *alpha_I_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(alpha_I_vect);
	double *alpha_I_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(alpha_I_error_vect);
	double *alpha_J_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(alpha_J_vect);
	double *alpha_J_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(alpha_J_error_vect);


	//FSC
	{
		long int k1;
		for(k1=0;k1<=number_of_sets_;k1++)
		{

			pair<long int,long int> boundaries;
			if(k1==0)
			{
				boundaries.first=0;
				boundaries.second=number_of_realizations_-1;
			}
			else
			{
				boundaries=set_number_boundaries(//boundaries of realization order numbers for the set; numeration of realizations starts from 0
				k1,//set order number; the numeration starts from 1
				number_of_realizations_set);//total number of sets
			};

			bool inside_simulation_flag;

			fsa_par::calculate_FSC(
			target_ALP_,
			boundaries.first,
			boundaries.second,
			

			lambda_out_vect[k1],

			M_max,

			distance_along_direction_1,
			distance_along_direction_2,

			ALP_weight,
			ALP_edge_max,

			a_I_vect[k1],
			a_I_error_vect[k1],
			a_J_vect[k1],
			a_J_error_vect[k1],
			sigma_vect[k1],
			sigma_error_vect[k1],
			alpha_I_vect[k1],
			alpha_I_error_vect[k1],
			alpha_J_vect[k1],
			alpha_J_error_vect[k1],
			inside_simulation_flag);

			if(!inside_simulation_flag)
			{
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			};
		};

		a_I_=a_I_vect[0];
		a_I_error_=a_I_error_vect[0];
		a_J_=a_J_vect[0];
		a_J_error_=a_J_error_vect[0];
		sigma_=sigma_vect[0];
		sigma_error_=sigma_error_vect[0];
		alpha_I_=alpha_I_vect[0];
		alpha_I_error_=alpha_I_error_vect[0];
		alpha_J_=alpha_J_vect[0];
		alpha_J_error_=alpha_J_error_vect[0];


		if(screen_output_flag_)
		{
			cout<<"a_I\t"<<a_I_;
		};

		if(number_of_sets_>1)
		{
			a_I_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			a_I_vect+1);

			if(average_set_results_flag)
			{
				a_I_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				a_I_vect+1);


				if(screen_output_flag_)
				{
					cout<<"\t"<<a_I_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<a_I_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_AISbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_AISbs.push_back(a_I_vect[k1]);
		};



		if(screen_output_flag_)
		{
			cout<<"a_J\t"<<a_J_;
		};

		if(number_of_sets_>1)
		{
			a_J_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			a_J_vect+1);

			if(average_set_results_flag)
			{
				a_J_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				a_J_vect+1);


				if(screen_output_flag_)
				{
					cout<<"\t"<<a_J_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<a_J_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_AJSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_AJSbs.push_back(a_J_vect[k1]);
		};


		if(screen_output_flag_)
		{
			cout<<"sigma\t"<<sigma_;
		};

		if(number_of_sets_>1)
		{
			sigma_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			sigma_vect+1);

			if(average_set_results_flag)
			{
				sigma_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				sigma_vect+1);


				if(screen_output_flag_)
				{
					cout<<"\t"<<sigma_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<sigma_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_SigmaSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_SigmaSbs.push_back(sigma_vect[k1]);
		};



		if(screen_output_flag_)
		{
			cout<<"alpha_I\t"<<alpha_I_;
		};

		if(number_of_sets_>1)
		{
			alpha_I_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			alpha_I_vect+1);

			if(average_set_results_flag)
			{
				alpha_I_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				alpha_I_vect+1);

				if(screen_output_flag_)
				{
					cout<<"\t"<<alpha_I_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<alpha_I_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_AlphaISbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_AlphaISbs.push_back(alpha_I_vect[k1]);
		};


		if(screen_output_flag_)
		{
			cout<<"alpha_J\t"<<alpha_J_;
		};

		if(number_of_sets_>1)
		{
			alpha_J_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			alpha_J_vect+1);

			if(average_set_results_flag)
			{
				alpha_J_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				alpha_J_vect+1);

				cout<<"\t"<<alpha_J_;
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<alpha_J_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_AlphaJSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_AlphaJSbs.push_back(alpha_J_vect[k1]);
		};

	};



	if(output_flag)
	{
		fE.close();
		fE_summary.close();
	};


	if(score_max)
	{
		for(k1=0;k1<=number_of_sets_;k1++)
		{
			delete[]score_max[k1];
		};
		delete []score_max;
	};

	if(distrE)
	{
		for(k1=0;k1<=number_of_sets_;k1++)
		{
			for(k=0;k<=target_ALP_;k++)
			{
				delete distrE[k1][k];
			};
			delete []distrE[k1];
		};
		delete []distrE;
	};

	if(distrE_errors)
	{
		for(k1=0;k1<=number_of_sets_;k1++)
		{
			for(k=0;k<=target_ALP_;k++)
			{
				delete distrE_errors[k1][k];
			};
			delete []distrE_errors[k1];
		};
		delete []distrE_errors;
	};
	


	//FSC
	delete[]distance_along_direction_1;
	delete[]distance_along_direction_2;
	delete[]ALP_weight;
	delete[]ALP_edge_max;



	delete[]lambda_out_vect;
	delete[]lambda_out_error_vect;
	delete[]C_vect;
	delete[]C_error_vect;
	delete[]K_C_vect;
	delete[]K_C_error_vect;
	delete[]a_I_vect;
	delete[]a_I_error_vect;
	delete[]a_J_vect;
	delete[]a_J_error_vect;
	delete[]sigma_vect;
	delete[]sigma_error_vect;
	delete[]alpha_I_vect;
	delete[]alpha_I_error_vect;
	delete[]alpha_J_vect;
	delete[]alpha_J_error_vect;

}

void test::restore_arrays(
long int k_,
mult_states_type *states_old_,
two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_,
IS1_general_simulation &IS1_general_simulation_test_)
{

	if(!states_old_)
	{
		throw error("Unexpected error in test::restore_arrays\n",1);
	};

	test::restore_arrays(
	states_old_->d_states[k_].d_two_dim_layer_alignment_algorithm,
	states_old_->d_states[k_].d_IS1_general_simulation,
	two_dim_layer_alignment_algorithm_test_,
	IS1_general_simulation_test_);
	
}


void test::restore_arrays(
two_dim_layer_alignment_algorithm<long int> *two_dim_layer_alignment_algorithm_test_from_,
IS1_general_simulation *IS1_general_simulation_test_from_,
two_dim_layer_alignment_algorithm<long int> &two_dim_layer_alignment_algorithm_test_,
IS1_general_simulation &IS1_general_simulation_test_)
{

	if(!two_dim_layer_alignment_algorithm_test_from_||!IS1_general_simulation_test_from_)
	{
		throw error("Unexpected error in test::restore_arrays\n",1);
	};
	
	//alignment object - begin
	long int i;
	for(i=0;i<two_dim_layer_alignment_algorithm_test_.d_number_of_variables;i++)
	{
		two_dim_layer<long int>::copy_2_into_1_two_dim_layer(
		two_dim_layer_alignment_algorithm_test_.d_vars[i],
		two_dim_layer_alignment_algorithm_test_from_->d_vars[i]);
	};
	two_dim_layer_alignment_algorithm_test_.d_current_align_ind1=two_dim_layer_alignment_algorithm_test_from_->d_current_align_ind1;
	two_dim_layer_alignment_algorithm_test_.d_current_align_ind2=two_dim_layer_alignment_algorithm_test_from_->d_current_align_ind2;
	//alignment object - end

	//IS1_general_simulation object - begin
	//restore the sequences
	IS1_general_simulation_test_.d_seq1_current_length=IS1_general_simulation_test_from_->d_seq1_current_length;
	IS1_general_simulation_test_.d_seq2_current_length=IS1_general_simulation_test_from_->d_seq2_current_length;

	for(i=0;i<IS1_general_simulation_test_from_->d_seq1_current_length;i++)
	{
		IS1_general_simulation_test_.d_seq1[i]=IS1_general_simulation_test_from_->d_seq1[i];
	};

	for(i=0;i<IS1_general_simulation_test_from_->d_seq2_current_length;i++)
	{
		IS1_general_simulation_test_.d_seq2[i]=IS1_general_simulation_test_from_->d_seq2[i];
	};

	//restore the layers for weights
	IS1_general_simulation_test_.d_current_state=IS1_general_simulation_test_from_->d_current_state;
	IS1_general_simulation_test_.d_W1_seq1_current_length=IS1_general_simulation_test_from_->d_W1_seq1_current_length;
	IS1_general_simulation_test_.d_W1_seq2_current_length=IS1_general_simulation_test_from_->d_W1_seq2_current_length;

	for(i=0;i<IS1_general_simulation_test_.d_IS1_general_obj->d_number_of_states;i++)
	{
		two_dim_layer<double>::copy_2_into_1_two_dim_layer(
		IS1_general_simulation_test_.d_W1[i],
		IS1_general_simulation_test_from_->d_W1[i]);
	};



	//IS1_general_simulation object - end

}

void test::collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states(
bool compute_allocated_memory_flag_,//if true then total allocated memory is computed in allocated_memory_
double &allocated_memory_,

long int number_of_realizations_,
string E_distr_file_name_,
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

bool futher_expanding_,
long int M_thr_,//the expanding starts from the ALP with the score >=M_thr_; used if >=0
IS1_general *IS1_general_cs_,
double *eps_K_,//relative error for K
long int *number_of_cs_steps_for_expanding_,
bool parameters_are_not_calculated_)
{

	if(!save_states_flag_&&ending_time_>=0)
	{
		throw error("Error - !save_states_flag_&&ending_time_>=0 in test::collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states\n",1);
	};

	if(compute_allocated_memory_flag_)
	{
		allocated_memory_=sizeof(state_type);
	};

	stopped_by_ending_time_flag_=false;

	inside_simulation_flag_=true;

	double time1=0;
	if(save_states_flag_||ending_time_>=0)
	{
		FSA_utils::get_current_time(time1);
	};
	

	double &lambda_out_=par_.lambda;
	double &lambda_out_error_=par_.lambda_error;
	double &C_=par_.C;
	double &C_error_=par_.C_error;

	double &a_I_=par_.a_I;
	double &a_I_error_=par_.a_I_error;
	double &a_J_=par_.a_J;
	double &a_J_error_=par_.a_J_error;
	double &sigma_=par_.sigma;
	double &sigma_error_=par_.sigma_error;
	double &alpha_I_=par_.alpha_I;
	double &alpha_I_error_=par_.alpha_I_error;
	double &alpha_J_=par_.alpha_J;
	double &alpha_J_error_=par_.alpha_J_error;

	double &K_C_=par_.K_C;
	double &K_C_error_=par_.K_C_error;



	bool average_set_results_flag=false;

	if(save_states_flag_)
	{
		states_new_->d_total_number_of_exp_cells=0;

		if(!states_old_)
		{
			states_new_->d_total_number_of_ALP_cells=0;
			states_new_->d_total_number_of_ALPs=0;
		}
		else
		{
			states_new_->d_total_number_of_ALP_cells=states_old_->d_total_number_of_ALP_cells;
			states_new_->d_total_number_of_ALPs=states_old_->d_total_number_of_ALPs;
		};
	};


	{
		long int number_of_realizations_tmp=(number_of_realizations_/number_of_sets_)*number_of_sets_;
		if(number_of_realizations_tmp<number_of_realizations_)
		{
			number_of_realizations_=number_of_realizations_tmp+number_of_sets_;
		};
	};

	long int number_of_realizations_set=number_of_realizations_/number_of_sets_;

	long int number_of_realizations_old=0;
	long int number_of_realizations_max=number_of_realizations_;

	if(states_old_)
	{
		//number_of_realizations_min=FSA_utils::Tmin(number_of_realizations_,states_old_->d_number_of_realizations);
		number_of_realizations_max=FSA_utils::Tmax(number_of_realizations_,states_old_->d_number_of_realizations);
		number_of_realizations_old=states_old_->d_number_of_realizations;


	};

	long int number_of_realizations_before_loop=number_of_realizations_;
	if(save_states_flag_)
	{
		states_new_->d_number_of_realizations=number_of_realizations_max;
		//allocate states
		states_new_->d_states=new state_type[states_new_->d_number_of_realizations];
		if(compute_allocated_memory_flag_)
		{
			allocated_memory_+=states_new_->d_number_of_realizations*sizeof(state_type);
		};

		long int ii;
		for(ii=0;ii<states_new_->d_number_of_realizations;ii++)
		{
			long int target_ALP_tmp;
			if(number_of_realizations_>number_of_realizations_old)
			{
				if(ii<number_of_realizations_old)
				{
					target_ALP_tmp=FSA_utils::Tmax(target_ALP_,states_old_->d_states[ii].d_ALP_number);
				}
				else
				{
					target_ALP_tmp=target_ALP_;
				};
			}
			else
			{
				if(ii<number_of_realizations_)
				{
					target_ALP_tmp=FSA_utils::Tmax(target_ALP_,states_old_->d_states[ii].d_ALP_number);
				}
				else
				{
					target_ALP_tmp=states_old_->d_states[ii].d_ALP_number;
				};
			};

			states_new_->d_states[ii].d_ALP_number=target_ALP_tmp;
			

			states_new_->d_states[ii].d_IS1_general_simulation=NULL;
			states_new_->d_states[ii].d_two_dim_layer_alignment_algorithm=NULL;
			states_new_->d_states[ii].d_expanding_finished_flag=false;


			states_new_->d_states[ii].d_ALP_weights_array=NULL;
			states_new_->d_states[ii].d_M_array=NULL;
			states_new_->d_states[ii].d_seq1_length_array=NULL;
			states_new_->d_states[ii].d_seq2_length_array=NULL;
			states_new_->d_states[ii].d_distance_along_direction_1_array=NULL;
			states_new_->d_states[ii].d_distance_along_direction_2_array=NULL;


			
		};

		

		if(states_old_)
		{


			long int ii;
			for(ii=0;ii<states_old_->d_number_of_realizations;ii++)
			{
				//allocate memory
				
				states_new_->d_states[ii].d_ALP_weights_array=new array_positive<double>;
				states_new_->d_states[ii].d_M_array=new array_positive<long int>;
				states_new_->d_states[ii].d_seq1_length_array=new array_positive<long int>;
				states_new_->d_states[ii].d_seq2_length_array=new array_positive<long int>;
				states_new_->d_states[ii].d_distance_along_direction_1_array=new array_positive<long int>;
				states_new_->d_states[ii].d_distance_along_direction_2_array=new array_positive<long int>;



				long int k;
				for(k=0;k<=states_old_->d_states[ii].d_ALP_number;k++)
				{
					states_new_->d_states[ii].d_ALP_weights_array->set_elem(k,states_old_->d_states[ii].d_ALP_weights_array->get_elem(k));
					states_new_->d_states[ii].d_M_array->set_elem(k,states_old_->d_states[ii].d_M_array->get_elem(k));
					states_new_->d_states[ii].d_seq1_length_array->set_elem(k,states_old_->d_states[ii].d_seq1_length_array->get_elem(k));
					states_new_->d_states[ii].d_seq2_length_array->set_elem(k,states_old_->d_states[ii].d_seq2_length_array->get_elem(k));
					states_new_->d_states[ii].d_distance_along_direction_1_array->set_elem(k,states_old_->d_states[ii].d_distance_along_direction_1_array->get_elem(k));
					states_new_->d_states[ii].d_distance_along_direction_2_array->set_elem(k,states_old_->d_states[ii].d_distance_along_direction_2_array->get_elem(k));

				};
			};

			for(ii=number_of_realizations_;ii<states_old_->d_number_of_realizations;ii++)
			{


				states_new_->d_states[ii].d_IS1_general_simulation=
				new IS1_general_simulation(//creates a copy of IS1_general_simulation_ with smallest possible memory allocation
				states_old_->d_states[ii].d_IS1_general_simulation);//maximum sequence length

				states_new_->d_states[ii].d_two_dim_layer_alignment_algorithm=
				new two_dim_layer_alignment_algorithm<long int>(
				states_old_->d_states[ii].d_seq1_length_array->get_elem(states_old_->d_states[ii].d_ALP_number),//target length of the sequence #1
				states_old_->d_states[ii].d_seq2_length_array->get_elem(states_old_->d_states[ii].d_ALP_number),//target length of the sequence #2
				states_old_->d_states[ii].d_two_dim_layer_alignment_algorithm);
			};
		};



	};

	//------------------------------------------------------

	bool output_flag=false;

	if(!two_dim_layer_alignment_algorithm_test_.d_M_flag)
	{
		throw error("Error - two_dim_layer_alignment_algorithm_test_.d_E_flag must be true in test::collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states\n",1);
	};

	ofstream fE;
	ofstream fE_summary;




	array_v<double> ***distrE=new array_v<double> **[number_of_sets_+1];
	FSA_utils::assert_mem(distrE);
	array_v<double> ***distrE_errors=new array_v<double> **[number_of_sets_+1];
	FSA_utils::assert_mem(distrE_errors);

	long int **score_max=new long int*[number_of_sets_+1];//statistic
	FSA_utils::assert_mem(score_max);

	long int k1;
	for(k1=0;k1<=number_of_sets_;k1++)
	{
		distrE[k1]=new array_v<double>*[target_ALP_+1];
		FSA_utils::assert_mem(distrE[k1]);
		distrE_errors[k1]=new array_v<double>*[target_ALP_+1];
		FSA_utils::assert_mem(distrE_errors[k1]);
		score_max[k1]=new long int[target_ALP_+1];
		FSA_utils::assert_mem(score_max[k1]);

		long int k;
		for(k=0;k<=target_ALP_;k++)
		{
			distrE[k1][k]=new array_v<double>(NULL);
			FSA_utils::assert_mem(distrE[k1][k]);

			distrE_errors[k1][k]=new array_v<double>(NULL);
			FSA_utils::assert_mem(distrE_errors[k1][k]);

			score_max[k1][k]=-inf;
		};

	};

	//distribution of E

	if(output_flag)
	{
		fE.open(E_distr_file_name_.data());
		if(!fE)
		{
			throw error("Error - the file "+E_distr_file_name_+" is not found\n",3);
		};

		string E_distr_file_name_summary=E_distr_file_name_+"_summary";
		fE_summary.open(E_distr_file_name_summary.data());
		if(!fE_summary)
		{
			throw error("Error - the file "+E_distr_file_name_summary+" is not found\n",3);
		};
	};




	long int k;


	array_positive<double> *diff=new array_positive<double>[number_of_sets_+1];//statistic
	array_positive<double> *diff_errors=new array_positive<double>[number_of_sets_+1];//statistic
	


	array_positive<long int> *distance_along_direction_1=new array_positive<long int>[number_of_realizations_max];//statistic
	array_positive<long int> *distance_along_direction_2=new array_positive<long int>[number_of_realizations_max];//statistic

	array_positive<double> *ALP_weight=new array_positive<double>[number_of_realizations_max];//statistic
	array_positive<long int> *ALP_edge_max=new array_positive<long int>[number_of_realizations_max];//statistic

	//--------------------

	
	double *sum_of_weights=new double[number_of_sets_+1];//statistic;
	double *sum_of_weights_error=new double[number_of_sets_+1];//statistic;

	for(k1=0;k1<=number_of_sets_;k1++)
	{
		sum_of_weights[k1]=0;//statistic
		sum_of_weights_error[k1]=0;//statistic
	};

	long int M_threshold=(long int)ceil(log((*eps_K_)*(1.0-ungapped_lambda_))/ungapped_lambda_);
	M_threshold*=2;


	if(M_threshold>=0)
	{
		M_threshold=-1;
	};

	two_dim_layer_alignment_algorithm_test_.d_M_threshold=M_threshold;
	
	long int M_max=0;

	if(save_states_flag_)
	{
		states_new_->d_average_ALP_pos1=0;
		states_new_->d_average_ALP_pos2=0;
		states_new_->d_average_ALP_pos1_mult_ALP_pos2=0;

		states_new_->d_average_expanding_length1=0;
		states_new_->d_average_expanding_length2=0;
		states_new_->d_average_expanding_length1_mult_expanding_length2=0;
	};


	//variables for saving the state at the moment when M_thr_ is reached

	two_dim_layer_alignment_algorithm<long int> *two_dim_layer_alignment_algorithm_test_M_thr=NULL;
	IS1_general_simulation *IS1_general_simulation_M_thr=NULL;

	bool M_thr_flag=(M_thr_>0)&&futher_expanding_;

	long int seq1_length_M_thr=0;
	long int seq2_length_M_thr=0;

	long int not_all_ALP_are_sampled=0;

	bool stop_loop_flag=false;
	for(k=0;k<number_of_realizations_;k++)
	{

		//allocate memory
		if(save_states_flag_)
		{
			if(k>=number_of_realizations_old)
			{
				states_new_->d_states[k].d_ALP_weights_array=new array_positive<double>;
				states_new_->d_states[k].d_M_array=new array_positive<long int>;
				states_new_->d_states[k].d_seq1_length_array=new array_positive<long int>;
				states_new_->d_states[k].d_seq2_length_array=new array_positive<long int>;
				states_new_->d_states[k].d_distance_along_direction_1_array=new array_positive<long int>;
				states_new_->d_states[k].d_distance_along_direction_2_array=new array_positive<long int>;

				
			};
		};


		long int sn=realization_number_into_set_number(k,number_of_realizations_set);


		IS1_general_simulation_test_.init();

		IS1_general_simulation_test_.d_W1_seq1_current_length=-1;
		IS1_general_simulation_test_.d_W1_seq2_current_length=-1;

		two_dim_layer_alignment_algorithm_test_.init();

		long int current_ALP_number=0;
		long int current_M=0;
		distrE[0][current_ALP_number]->increase_elem_by_x(current_M,1.0);

		distrE[sn][current_ALP_number]->increase_elem_by_x(current_M,1.0);
		
		distance_along_direction_1[k].set_elem(current_ALP_number,0);
		distance_along_direction_2[k].set_elem(current_ALP_number,0);

		ALP_weight[k].set_elem(current_ALP_number,1.0);
		ALP_edge_max[k].set_elem(current_ALP_number,0);

		score_max[0][current_ALP_number]=current_M;
		score_max[sn][current_ALP_number]=current_M;


		long int seq1_length=0;
		long int seq2_length=0;

		if(save_states_flag_)
		{
			if(!states_old_||k>=number_of_realizations_old)
			{
				states_new_->d_states[k].d_ALP_weights_array->set_elem(0,1.0);
				states_new_->d_states[k].d_M_array->set_elem(0,0);
				states_new_->d_states[k].d_seq1_length_array->set_elem(0,0);
				states_new_->d_states[k].d_seq2_length_array->set_elem(0,0);
				states_new_->d_states[k].d_distance_along_direction_1_array->set_elem(0,0);
				states_new_->d_states[k].d_distance_along_direction_2_array->set_elem(0,0);

			};
		};


		//weights calculation
		double weight2=1;
		double tmp_w=1;
		double tmp_w2=1;

		double tmp_w_M_thr=0;
		double tmp_w2_M_thr=0;


		bool M_thr_var_are_saved=false;

		bool restored_flag=false;
		bool flag_tmp3=false;

		
		while((current_ALP_number<target_ALP_)||(M_thr_flag&&!M_thr_var_are_saved))
		{
			bool ALP_stat_flag=(current_ALP_number<target_ALP_);

			
			if(ALP_stat_flag&&M_thr_var_are_saved)
			{
				flag_tmp3=true;
			};


			seq1_length+=two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth1;
			seq2_length+=two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth2;


			if(seq1_length>two_dim_layer_alignment_algorithm_test_.d_max_ind1||
				seq2_length>two_dim_layer_alignment_algorithm_test_.d_max_ind2)
			{
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			};


			bool new_calculation_flag=true;

			if(states_old_)
			{
				if(k<number_of_realizations_old)
				{
					if(current_ALP_number<states_old_->d_states[k].d_ALP_number)
					{
						new_calculation_flag=false;
					};

				};
			};
			

			bool ALP_flag=false;


			if(new_calculation_flag)
			{
				if(states_old_)
				{
					if(k<number_of_realizations_old)
					{
						if(current_ALP_number==states_old_->d_states[k].d_ALP_number&&!restored_flag)
						{//restore the state; 				seq1_length, seq2_length was increased, require correction

							restored_flag=true;

							restore_arrays(
							k,
							states_old_,
							two_dim_layer_alignment_algorithm_test_,
							IS1_general_simulation_test_);


						};
					};
				};

				if(save_states_flag_)
				{
					long int seq1_length_prev=IS1_general_simulation_test_.d_W1_seq1_current_length;
					long int seq2_length_prev=IS1_general_simulation_test_.d_W1_seq2_current_length;

					states_new_->d_total_number_of_ALP_cells+=
						seq1_length*seq2_length-seq1_length_prev*seq2_length_prev;
				};




				IS1_general_simulation_test_.simulate_upto_target_lengths(
				seq1_length,
				seq2_length);


				IS1_general_simulation_test_.calculate_weight_W1_upto_target_lengths(
				seq1_length,//target length of the sequence #1
				seq2_length,//target length of the sequence #2
				seq1_length,//the weights are calculated upto this length for the sequence #1
				seq2_length);//the weights are calculated upto this length for the sequence #2



				two_dim_layer_alignment_algorithm_test_.align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
				//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
				seq1_length,//target length of the sequence #1
				seq2_length);//target length of the sequence #2


				ALP_flag=(two_dim_layer_alignment_algorithm_test_.d_M>current_M);

				if(save_states_flag_)
				{
					if(ALP_flag)
					{
						states_new_->d_total_number_of_ALPs++;
					};
				};
			};


			if(states_old_)
			{
				if(!new_calculation_flag)
				{
					ALP_flag=(seq1_length==states_old_->d_states[k].d_seq1_length_array->get_elem(current_ALP_number+1))&&
						(seq2_length==states_old_->d_states[k].d_seq2_length_array->get_elem(current_ALP_number+1));

					

					if(ALP_flag)
					{
						//restore the state

						two_dim_layer_alignment_algorithm_test_.d_M=states_old_->d_states[k].d_M_array->get_elem(current_ALP_number+1);

						weight2=states_old_->d_states[k].d_ALP_weights_array->get_elem(current_ALP_number+1);

						two_dim_layer_alignment_algorithm_test_.d_distance_along_direction_1=
							states_old_->d_states[k].d_distance_along_direction_1_array->get_elem(current_ALP_number+1);

						two_dim_layer_alignment_algorithm_test_.d_distance_along_direction_2=
							states_old_->d_states[k].d_distance_along_direction_2_array->get_elem(current_ALP_number+1);


					};

				};
			};

			if(ALP_flag)
			{//new ALP found

				current_M=two_dim_layer_alignment_algorithm_test_.d_M;
				current_ALP_number++;



				if(new_calculation_flag)
				{
					IS1_general_simulation_test_.calculate_weight_W1_for_fixed_lengths_using_recursions(
					seq1_length,//target length of the sequence #1
					seq2_length,//target length of the sequence #2
					weight2);//the resulted weight

					if(save_states_flag_)
					{
						states_new_->d_states[k].d_ALP_weights_array->set_elem(current_ALP_number,weight2);
						states_new_->d_states[k].d_distance_along_direction_1_array->set_elem(current_ALP_number,two_dim_layer_alignment_algorithm_test_.d_distance_along_direction_1);
						states_new_->d_states[k].d_distance_along_direction_2_array->set_elem(current_ALP_number,two_dim_layer_alignment_algorithm_test_.d_distance_along_direction_2);
						states_new_->d_states[k].d_M_array->set_elem(current_ALP_number,two_dim_layer_alignment_algorithm_test_.d_M);
						states_new_->d_states[k].d_seq1_length_array->set_elem(current_ALP_number,seq1_length);
						states_new_->d_states[k].d_seq2_length_array->set_elem(current_ALP_number,seq2_length);
					};



				};




				if(ALP_stat_flag)
				{
					score_max[0][current_ALP_number]=FSA_utils::Tmax(score_max[0][current_ALP_number],current_M);
					score_max[sn][current_ALP_number]=FSA_utils::Tmax(score_max[sn][current_ALP_number],current_M);
				};

				if(weight2<=0)
				{
					inside_simulation_flag_=false;
					tmp_w=0;
				}
				else
				{
					tmp_w=1/weight2;
				};

				tmp_w2=tmp_w*tmp_w;

			

				if(ALP_stat_flag)
				{
					distrE[0][current_ALP_number]->increase_elem_by_x(current_M,tmp_w);
					distrE[sn][current_ALP_number]->increase_elem_by_x(current_M,tmp_w);
					
					distrE_errors[0][current_ALP_number]->increase_elem_by_x(current_M,tmp_w2);
					distrE_errors[sn][current_ALP_number]->increase_elem_by_x(current_M,tmp_w2);

					distance_along_direction_1[k].set_elem(current_ALP_number,two_dim_layer_alignment_algorithm_test_.d_distance_along_direction_1);
					distance_along_direction_2[k].set_elem(current_ALP_number,two_dim_layer_alignment_algorithm_test_.d_distance_along_direction_2);

					ALP_weight[k].set_elem(current_ALP_number,tmp_w);
					ALP_edge_max[k].set_elem(current_ALP_number,current_M);


					M_max=FSA_utils::Tmax(current_M,M_max);
				};

				

		
			};

			if(seq2_length>limit2_)
			{
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			};


			

			
			bool save_flag_tmp=((current_ALP_number==target_ALP_)&&!M_thr_flag)||
			(M_thr_flag&&(two_dim_layer_alignment_algorithm_test_.d_M>=M_thr_)&&(current_ALP_number>=target_ALP_));
			
			if(save_flag_tmp)
			{//save the matrices for weights and alingment matrix


				if(states_old_)
				{
					if(k<number_of_realizations_old)
					{
						if(current_ALP_number==states_old_->d_states[k].d_ALP_number&&!restored_flag)
						{//restore the state; 				

							restored_flag=true;

							restore_arrays(
							k,
							states_old_,
							two_dim_layer_alignment_algorithm_test_,
							IS1_general_simulation_test_);


						};
					};
				};

				//save states
				if(save_states_flag_)
				{

					states_new_->d_states[k].d_ALP_number=current_ALP_number;

					states_new_->d_states[k].d_IS1_general_simulation=
					new IS1_general_simulation(//creates a copy of IS1_general_simulation_ with smallest possible memory allocation
					&IS1_general_simulation_test_);//maximum sequence length

					states_new_->d_states[k].d_two_dim_layer_alignment_algorithm=
					new two_dim_layer_alignment_algorithm<long int>(
					seq1_length,//target length of the sequence #1
					seq2_length,//target length of the sequence #2
					&two_dim_layer_alignment_algorithm_test_);

					states_new_->d_average_ALP_pos1+=seq1_length;
					states_new_->d_average_ALP_pos2+=seq2_length;
					states_new_->d_average_ALP_pos1_mult_ALP_pos2+=seq1_length*seq2_length;
				};


			};

			//check wheter we reached the score threshold M_thr_

			if(M_thr_flag&&(two_dim_layer_alignment_algorithm_test_.d_M>=M_thr_)&&!M_thr_var_are_saved)
			{

				M_thr_var_are_saved=true;
				tmp_w_M_thr=tmp_w;
				tmp_w2_M_thr=tmp_w2;

				delete IS1_general_simulation_M_thr;

				IS1_general_simulation_M_thr=
				new IS1_general_simulation(//creates a copy of IS1_general_simulation_ with smallest possible memory allocation
				&IS1_general_simulation_test_);//maximum sequence length

				delete two_dim_layer_alignment_algorithm_test_M_thr;

				two_dim_layer_alignment_algorithm_test_M_thr=
				new two_dim_layer_alignment_algorithm<long int>(
				seq1_length,//target length of the sequence #1
				seq2_length,//target length of the sequence #2
				&two_dim_layer_alignment_algorithm_test_);

				two_dim_layer_alignment_algorithm_test_.d_scores_counts_flag=false;

				seq1_length_M_thr=seq1_length;
				seq2_length_M_thr=seq2_length;

			};



			if(ending_time_>=0)
			{
				if(new_calculation_flag)
				{
					double time1_current;
					FSA_utils::get_current_time(time1_current);
					if(time1_current>=ending_time_)
					{
						stop_loop_flag=true;
						stopped_by_ending_time_flag_=true;
					};
				};
			};

			if(ending_time_to_test_logarithmic_regime_>0)
			{
				double time1_current;
				FSA_utils::get_current_time(time1_current);
				if(time1_current>=ending_time_to_test_logarithmic_regime_)
				{
					throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
				};
			};

		};

		
		if(flag_tmp3)
		{
			not_all_ALP_are_sampled++;
		};

		if(screen_output_k_flag_)
		{
			if((k+1)%100==0)
			{
				cout<<k+1<<endl;
			};
		};



		//expanding

		if(M_thr_flag)
		{
			if(!M_thr_var_are_saved)
			{
				throw error("Unexpected error in test::collect_and_output_E_distr_upto_fixed_ALP_with_saving_of_states\n",1);
			};

			two_dim_layer_alignment_algorithm_test_.d_scores_counts_flag=true;

			seq1_length=seq1_length_M_thr;
			seq2_length=seq2_length_M_thr;

			tmp_w=tmp_w_M_thr;
			tmp_w2=tmp_w2_M_thr;



			test::restore_arrays(
			two_dim_layer_alignment_algorithm_test_M_thr,
			IS1_general_simulation_M_thr,
			two_dim_layer_alignment_algorithm_test_,
			IS1_general_simulation_test_);

			two_dim_layer_alignment_algorithm_test_.d_M=two_dim_layer_alignment_algorithm_test_M_thr->d_M;

		};

		
		long int seq1_length_before_expanding=seq1_length;
		long int seq2_length_before_expanding=seq2_length;
		
		if(futher_expanding_)
		{
			bool expanding_new_flag=true;
			//restore state

			if(states_old_)
			{
				if(k<number_of_realizations_old)
				{
					if(current_ALP_number==states_old_->d_states[k].d_ALP_number&&(!restored_flag||states_old_->d_states[k].d_expanding_finished_flag))
					{
						if(states_old_->d_states[k].d_expanding_finished_flag)
						{
							expanding_new_flag=false;
							two_dim_layer_alignment_algorithm_test_.d_M=states_old_->d_states[k].d_M_after_expanding;
							array_v<double>::copy_x2_into_x1(
								two_dim_layer_alignment_algorithm_test_.d_scores_counts,
								states_old_->d_states[k].d_scores_counts_after_expanding);

							seq1_length=states_old_->d_states[k].d_seq1_length_after_expanding;
							seq2_length=states_old_->d_states[k].d_seq2_length_after_expanding;

						}
						else
						{
							restore_arrays(
							k,
							states_old_,
							two_dim_layer_alignment_algorithm_test_,
							IS1_general_simulation_test_);
						};
					};
				};
			};

			if(expanding_new_flag)
			{
				bool tmp_bool=two_dim_layer_alignment_algorithm_test_.d_E_flag;
				two_dim_layer_alignment_algorithm_test_.d_E_flag=true;

				long int current_state=0;//the current state
				long int new_state;

				bool success_flag=false;

				long int j;
				for(j=1;j<=*number_of_cs_steps_for_expanding_;j++)
				{
					IS1_general_cs_->one_step_of_the_importance_samping(//an initial state when sequences are empty must be defined outside
					current_state,//the current state
					IS1_general_simulation_test_.d_seq1+seq1_length,//letters for the sequence #1; the array must be allocated
					IS1_general_simulation_test_.d_seq2+seq2_length,//letters for the sequence #2; the array must be allocated
					new_state);//a new state

					seq1_length+=two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth1;
					seq2_length+=two_dim_layer_alignment_algorithm_test_.d_edge_maximum_calculation_depth2;


					two_dim_layer_alignment_algorithm_test_.align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
					//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
					seq1_length,//target length of the sequence #1
					seq2_length);//target length of the sequence #2

					long int &M=two_dim_layer_alignment_algorithm_test_.d_M;
					long int &E=two_dim_layer_alignment_algorithm_test_.d_E;


					if(E-M<M_threshold)
					{
						success_flag=true;
						break;
					};


					if(ending_time_to_test_logarithmic_regime_>0)
					{
						double time1_current;
						FSA_utils::get_current_time(time1_current);
						if(time1_current>=ending_time_to_test_logarithmic_regime_)
						{
							throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
						};
					};


				};


				if(!success_flag)
				{
					throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
				};

				two_dim_layer_alignment_algorithm_test_.d_E_flag=tmp_bool;

			};

				

			IS1_general_simulation_test_.d_seq1_current_length=seq1_length;
			IS1_general_simulation_test_.d_seq2_current_length=seq2_length;

			if(save_states_flag_)
			{
				states_new_->d_average_expanding_length1+=seq1_length-seq1_length_before_expanding;
				states_new_->d_average_expanding_length2+=seq2_length-seq2_length_before_expanding;
				states_new_->d_average_expanding_length1_mult_expanding_length2+=
					(seq1_length-seq1_length_before_expanding)*(seq2_length-seq2_length_before_expanding);

				states_new_->d_total_number_of_exp_cells+=
					seq1_length*seq2_length-seq1_length_before_expanding*seq2_length_before_expanding;

				states_new_->d_states[k].d_expanding_finished_flag=true;
				states_new_->d_states[k].d_M_after_expanding=two_dim_layer_alignment_algorithm_test_.d_M;
				states_new_->d_states[k].d_seq1_length_after_expanding=seq1_length;
				states_new_->d_states[k].d_seq2_length_after_expanding=seq2_length;

				array_v<double>::copy_x2_into_x1(
				states_new_->d_states[k].d_scores_counts_after_expanding,
				two_dim_layer_alignment_algorithm_test_.d_scores_counts);



			};

			

		};
		


		//collect differences
		if(two_dim_layer_alignment_algorithm_test_.d_scores_counts_flag&&futher_expanding_)
		{
			//sum of weigths
			sum_of_weights[0]+=tmp_w;
			sum_of_weights[sn]+=tmp_w;
			sum_of_weights_error[0]+=tmp_w2;
			sum_of_weights_error[sn]+=tmp_w2;

			long int &M=two_dim_layer_alignment_algorithm_test_.d_M;
			long int ii;
			
			for(ii=FSA_utils::Tmax(M_threshold,two_dim_layer_alignment_algorithm_test_.d_scores_counts.d_ind0);ii<=FSA_utils::Tmin(M,two_dim_layer_alignment_algorithm_test_.d_scores_counts.d_dim_plus_d_ind0);ii++)
			{
				long int diff_tmp=(M-ii);
				double diff_tmp_weight=two_dim_layer_alignment_algorithm_test_.d_scores_counts.get_elem(ii)*tmp_w;
				double diff_tmp_weight2=diff_tmp_weight*diff_tmp_weight;

				diff[0].increase_elem_by_x(diff_tmp,diff_tmp_weight);
				diff[sn].increase_elem_by_x(diff_tmp,diff_tmp_weight);

				diff_errors[0].increase_elem_by_x(diff_tmp,diff_tmp_weight2);
				diff_errors[sn].increase_elem_by_x(diff_tmp,diff_tmp_weight2);
			};

		};


	
		if(ending_time_>=0)
		{
			if(stop_loop_flag)
			{
				//the loop is stopped due to time limitations
				number_of_realizations_with_ending_time_=k+1;
				if(number_of_realizations_with_ending_time_%number_of_realizations_set==0)
				{
					number_of_sets_=sn;
				}
				else
				{
					number_of_sets_=sn-1;
				};

				number_of_realizations_=number_of_realizations_with_ending_time_;

				break;
			};
		};
	};


	
	
	//cout<<"Number of real when not all ALP are sampled\t"<<not_all_ALP_are_sampled<<"\t"<<number_of_realizations_<<endl;

	delete IS1_general_simulation_M_thr;
	delete two_dim_layer_alignment_algorithm_test_M_thr;

	if(stop_loop_flag)
	{
		//the loop is stopped due to time limitations
		//number of realizations and number of subsets is different
		if(save_states_flag_)
		{
			mult_states_type *states_new_loop=new mult_states_type;
			*states_new_loop=*states_new_;

			states_new_loop->d_number_of_realizations=number_of_realizations_with_ending_time_;

			if(states_old_)
			{
				states_new_loop->d_number_of_realizations=
					FSA_utils::Tmax(number_of_realizations_with_ending_time_,states_old_->d_number_of_realizations);
			};

			states_new_loop->d_states=new state_type[states_new_loop->d_number_of_realizations];
			long int ii;
			for(ii=0;ii<states_new_loop->d_number_of_realizations;ii++)
			{


				states_new_loop->d_states[ii].d_ALP_number=states_new_->d_states[ii].d_ALP_number;
				states_new_loop->d_states[ii].d_ALP_weights_array=states_new_->d_states[ii].d_ALP_weights_array;
				states_new_loop->d_states[ii].d_M_array=states_new_->d_states[ii].d_M_array;
				states_new_loop->d_states[ii].d_seq1_length_array=states_new_->d_states[ii].d_seq1_length_array;
				states_new_loop->d_states[ii].d_seq2_length_array=states_new_->d_states[ii].d_seq2_length_array;
				states_new_loop->d_states[ii].d_distance_along_direction_1_array=states_new_->d_states[ii].d_distance_along_direction_1_array;
				states_new_loop->d_states[ii].d_distance_along_direction_2_array=states_new_->d_states[ii].d_distance_along_direction_2_array;

				states_new_loop->d_states[ii].d_IS1_general_simulation=
					states_new_->d_states[ii].d_IS1_general_simulation;
				states_new_loop->d_states[ii].d_two_dim_layer_alignment_algorithm=
					states_new_->d_states[ii].d_two_dim_layer_alignment_algorithm;

				states_new_loop->d_states[ii].d_expanding_finished_flag=
					states_new_->d_states[ii].d_expanding_finished_flag;
				if(states_new_->d_states[ii].d_expanding_finished_flag)
				{
					states_new_loop->d_states[ii].d_M_after_expanding=
						states_new_->d_states[ii].d_M_after_expanding;
					states_new_loop->d_states[ii].d_seq1_length_after_expanding=
						states_new_->d_states[ii].d_seq1_length_after_expanding;
					states_new_loop->d_states[ii].d_seq2_length_after_expanding=
						states_new_->d_states[ii].d_seq2_length_after_expanding;

					array_v<double>::copy_x2_into_x1(
					states_new_loop->d_states[ii].d_scores_counts_after_expanding,
						states_new_->d_states[ii].d_scores_counts_after_expanding);
				};
			};




			if(states_old_)
			{
				for(ii=number_of_realizations_with_ending_time_;ii<FSA_utils::Tmin(states_old_->d_number_of_realizations,number_of_realizations_before_loop);ii++)
				{


					states_new_loop->d_states[ii].d_IS1_general_simulation=
					new IS1_general_simulation(//creates a copy of IS1_general_simulation_ with smallest possible memory allocation
					states_old_->d_states[ii].d_IS1_general_simulation);//maximum sequence length

					states_new_loop->d_states[ii].d_two_dim_layer_alignment_algorithm=
					new two_dim_layer_alignment_algorithm<long int>(
					//states_old_->d_states[ii].d_seq1_length[states_old_->d_states[ii].d_ALP_number],//target length of the sequence #1
					states_old_->d_states[ii].d_seq1_length_array->get_elem(states_old_->d_states[ii].d_ALP_number),//target length of the sequence #1
					//states_old_->d_states[ii].d_seq2_length[states_old_->d_states[ii].d_ALP_number],//target length of the sequence #2
					states_old_->d_states[ii].d_seq2_length_array->get_elem(states_old_->d_states[ii].d_ALP_number),//target length of the sequence #2
					states_old_->d_states[ii].d_two_dim_layer_alignment_algorithm);


					states_new_loop->d_states[ii].d_ALP_number=states_old_->d_states[ii].d_ALP_number;

					states_new_loop->d_states[ii].d_expanding_finished_flag=
						states_old_->d_states[ii].d_expanding_finished_flag;

					if(states_old_->d_states[ii].d_expanding_finished_flag)
					{
						states_new_loop->d_states[ii].d_M_after_expanding=
							states_old_->d_states[ii].d_M_after_expanding;
						states_new_loop->d_states[ii].d_seq1_length_after_expanding=
							states_old_->d_states[ii].d_seq1_length_after_expanding;
						states_new_loop->d_states[ii].d_seq2_length_after_expanding=
							states_old_->d_states[ii].d_seq2_length_after_expanding;

						array_v<double>::copy_x2_into_x1(
						states_new_loop->d_states[ii].d_scores_counts_after_expanding,
							states_old_->d_states[ii].d_scores_counts_after_expanding);
					};

				};
			};


			for(ii=states_new_loop->d_number_of_realizations;ii<states_new_->d_number_of_realizations;ii++)
			{
				delete states_new_->d_states[ii].d_ALP_weights_array;
				delete states_new_->d_states[ii].d_M_array;
				delete states_new_->d_states[ii].d_seq1_length_array;
				delete states_new_->d_states[ii].d_seq2_length_array;
				delete states_new_->d_states[ii].d_distance_along_direction_1_array;
				delete states_new_->d_states[ii].d_distance_along_direction_2_array;

				delete states_new_->d_states[ii].d_IS1_general_simulation;
				delete states_new_->d_states[ii].d_two_dim_layer_alignment_algorithm;
			};

			delete []states_new_->d_states;
			delete states_new_;
			states_new_=states_new_loop;
		};
	};


	
	double time2;
	if(save_states_flag_)
	{
		FSA_utils::get_current_time(time2);
		double total_time_tmp=time2-time1;
		if(states_old_&&!futher_expanding_)
		{
			total_time_tmp+=states_old_->d_total_calculation_time;
		};
		states_new_->d_total_calculation_time=total_time_tmp;
		states_new_->d_number_of_ALP=target_ALP_;

		states_new_->d_average_ALP_pos1/=(double)(states_new_->d_number_of_realizations);
		states_new_->d_average_ALP_pos2/=(double)(states_new_->d_number_of_realizations);
		states_new_->d_average_ALP_pos1_mult_ALP_pos2/=(double)(states_new_->d_number_of_realizations);

		states_new_->d_average_expanding_length1/=(double)(states_new_->d_number_of_realizations);
		states_new_->d_average_expanding_length2/=(double)(states_new_->d_number_of_realizations);
		states_new_->d_average_expanding_length1_mult_expanding_length2/=(double)(states_new_->d_number_of_realizations);

	};

	//calculation of parameters
	for(k1=0;k1<=number_of_sets_;k1++)
	{
		long int number_of_realizations_tmp=number_of_realizations_set;
		if(k1==0)
		{
			number_of_realizations_tmp=number_of_realizations_;
		};

		for(k=0;k<=target_ALP_;k++)
		{
			long int j;
			for(j=distrE[k1][k]->d_ind0;j<=score_max[k1][k];j++)
			{
				distrE[k1][k]->divide_elem_by_x(j,(double)number_of_realizations_tmp);
				distrE_errors[k1][k]->divide_elem_by_x(j,(double)number_of_realizations_tmp);
				distrE_errors[k1][k]->increase_elem_by_x(j,-distrE[k1][k]->get_elem(j)*distrE[k1][k]->get_elem(j));
				distrE_errors[k1][k]->divide_elem_by_x(j,(double)number_of_realizations_tmp);
			};
		};
	};


	//calculate M_thr_estimation_

	{
		double M_thr_estimation_tmp=0;
		double sum_of_weigths_tmp=0;

		long int j;
		for(j=distrE[0][target_ALP_]->d_ind0;j<=score_max[0][target_ALP_];j++)
		{
			M_thr_estimation_tmp+=distrE[0][target_ALP_]->get_elem(j)*j;
			sum_of_weigths_tmp+=distrE[0][target_ALP_]->get_elem(j);
		};

		if(sum_of_weigths_tmp>=0)
		{
			M_thr_estimation_tmp/=sum_of_weigths_tmp;
			M_thr_estimation_=(long int)FSA_utils::round(M_thr_estimation_tmp);
		}
		else
		{
			throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
		};
	};

	if(compute_allocated_memory_flag_)
	{
		double states_static_size=(sizeof(IS1_general_simulation)+sizeof(array_positive<double>)+5*sizeof(array_positive<long int>));
		allocated_memory_+=states_new_->d_number_of_realizations*states_static_size;

		double double_dim=0;
		double long_dim=0;
		long int ii;
		for(ii=0;ii<states_new_->d_number_of_realizations;ii++)
		{
			double_dim+=states_new_->d_states[ii].d_ALP_weights_array->d_dim;
			long_dim+=states_new_->d_states[ii].d_M_array->d_dim;
			long_dim+=states_new_->d_states[ii].d_seq1_length_array->d_dim;
			long_dim+=states_new_->d_states[ii].d_seq2_length_array->d_dim;
			long_dim+=states_new_->d_states[ii].d_distance_along_direction_1_array->d_dim;
			long_dim+=states_new_->d_states[ii].d_distance_along_direction_2_array->d_dim;


			long int n_states=states_new_->d_states[ii].d_IS1_general_simulation->d_IS1_general_obj->d_number_of_states;
			allocated_memory_+=sizeof(two_dim_layer<double>)*(double)n_states;
			for(k=0;k<n_states;k++)
			{
				two_dim_layer<double> *&tmp1=states_new_->d_states[ii].d_IS1_general_simulation->d_W1[k];

				double_dim+=
				tmp1->d_layers_number1 *
					tmp1->d_max_ind2+
				tmp1->d_layers_number2 *
					tmp1->d_max_ind1;

			};

			long_dim+=states_new_->d_states[ii].d_IS1_general_simulation->d_max_seq1_length+
				states_new_->d_states[ii].d_IS1_general_simulation->d_max_seq2_length;

			long int var_n=states_new_->d_states[ii].d_two_dim_layer_alignment_algorithm->d_number_of_variables;
			allocated_memory_+=sizeof(two_dim_layer<long int>)*(double)var_n;

			for(k=0;k<var_n;k++)
			{
				two_dim_layer<long int> *&tmp1=states_new_->d_states[ii].d_two_dim_layer_alignment_algorithm->d_vars[k];
				
					long_dim+=
					tmp1->d_layers_number1 *
						tmp1->d_max_ind2+
					tmp1->d_layers_number2 *
						tmp1->d_max_ind1;
			};
		};


		allocated_memory_+=double_dim*sizeof(double);
		allocated_memory_+=long_dim*sizeof(long);

		if(states_new_->d_number_of_realizations>0)
		{
			allocated_memory_/=(double)states_new_->d_number_of_realizations;
		};


	};

	if(parameters_are_not_calculated_)
	{
		if(score_max)
		{
			for(k1=0;k1<=number_of_sets_;k1++)
			{
				delete[]score_max[k1];
			};
			delete []score_max;
		};

		if(distrE)
		{
			for(k1=0;k1<=number_of_sets_;k1++)
			{
				for(k=0;k<=target_ALP_;k++)
				{
					delete distrE[k1][k];
				};
				delete []distrE[k1];
			};
			delete []distrE;
			distrE=NULL;
		};

		if(distrE_errors)
		{
			for(k1=0;k1<=number_of_sets_;k1++)
			{
				for(k=0;k<=target_ALP_;k++)
				{
					delete distrE_errors[k1][k];
				};
				delete []distrE_errors[k1];
			};
			delete []distrE_errors;
			distrE_errors=NULL;
		};

		delete[]diff;
		delete[]diff_errors;

		delete[]distance_along_direction_1;
		delete[]distance_along_direction_2;

		delete[]ALP_weight;
		delete[]ALP_edge_max;

		delete[]sum_of_weights;
		delete[]sum_of_weights_error;

		return;
	};


	

	par_.realizations_number=number_of_realizations_;

	double *lambda_out_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(lambda_out_vect);

	double *lambda_out_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(lambda_out_error_vect);


	double *C_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(C_vect);

	double *C_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(C_error_vect);



	{

		long int nalp=target_ALP_;

		for(k1=0;k1<=number_of_sets_;k1++)
		{
			bool check_the_criteria=false;
			long int nalp_thr;
			bool inside_simulation_flag;
			void **alp_distr=(void**)distrE[k1];
			void **alp_distr_errors=(void**)distrE_errors[k1];
			double ungapped_lambda=ungapped_lambda_;
			double lambda;
			double lambda_error;

			fsa_par::calculate_lambda(
			check_the_criteria,
			nalp,
			nalp_thr,
			inside_simulation_flag,
			alp_distr,
			alp_distr_errors,
			ungapped_lambda,
			lambda,
			lambda_error);


			if(inside_simulation_flag)
			{
				lambda_out_vect[k1]=lambda;
				lambda_out_error_vect[k1]=lambda_error;

				if(k1==0)
				{
					struct_for_lambda_calculation tmp_struct;
					tmp_struct.d_alp_distr=alp_distr;
					tmp_struct.d_alp_distr_errors=alp_distr_errors;
					tmp_struct.d_nalp=nalp;
					tmp_struct.d_calculate_alp_number=false;

					void * data_tmp=(void*)(&tmp_struct);
					fsa_par::function_for_lambda_calculation(
					lambda,
					data_tmp);

					double error2=FSA_utils::Tmax(fabs(tmp_struct.d_last_sum-tmp_struct.d_before_last_sum),fabs(tmp_struct.d_last_sum_error));
					double fabs_tmp=fabs(FSA_utils::Tmin(tmp_struct.d_last_sum,tmp_struct.d_before_last_sum));


					if(fabs_tmp<=0)
					{
						throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
					};

					error2/=fabs_tmp;

					par_.lambda_last_ALP_relative_error=error2;

				};
			}
			else
			{
				inside_simulation_flag_=false;
				lambda_out_vect[k1]=-1;
				lambda_out_error_vect[k1]=-1;
			};
		};

		
		lambda_out_=lambda_out_vect[0];
		lambda_out_error_=lambda_out_error_vect[0];

		if(screen_output_flag_)
		{
			//cout<<endl;
			cout<<"Parameter\tvalue\terror\n";
		};


		if(screen_output_flag_)
		{
			cout<<"Lambda\t"<<lambda_out_;
		};

		if(number_of_sets_>1)
		{
			lambda_out_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			lambda_out_vect+1);

			if(average_set_results_flag)
			{
				lambda_out_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				lambda_out_vect+1);


				if(screen_output_flag_)
				{
					cout<<"\t"<<lambda_out_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<lambda_out_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_LambdaSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_LambdaSbs.push_back(lambda_out_vect[k1]);
		};


		for(k1=0;k1<=number_of_sets_;k1++)
		{

			void **alp_distr=(void**)distrE[k1];
			void **alp_distr_errors=(void**)distrE_errors[k1];

			long int starting_point=0;

			bool inside_simulation_flag;

			fsa_par::calculate_C(
			starting_point,
			nalp,
			alp_distr,
			alp_distr_errors,
			lambda_out_vect[k1],
			lambda_out_error_vect[k1],
			C_vect[k1],
			C_error_vect[k1],
			inside_simulation_flag);


			if(!inside_simulation_flag)
			{
				inside_simulation_flag_=false;
				C_vect[k1]=-1;
				C_error_vect[k1]=-1;
			};
		};

		C_=C_vect[0];
		C_error_=C_error_vect[0];

		if(screen_output_flag_)
		{
			cout<<"C\t"<<C_;
		};

		if(number_of_sets_>1)
		{
			C_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			C_vect+1);

			if(average_set_results_flag)
			{
				C_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				C_vect+1);

				if(screen_output_flag_)
				{
					cout<<"\t"<<C_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<C_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_CSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_CSbs.push_back(C_vect[k1]);
		};

	};



	for(k=0;k<=target_ALP_;k++)
	{

		if(output_flag)
		{
			long int s;

			double sum_tmp=0;
			double sum_tmp_error=0;

			fE<<k<<endl;
			fE<<score_max[0][k]+1<<endl;

			for(s=0;s<=score_max[0][k];s++)
			{
				fE<<s<<"\t"<<distrE[0][k]->d_elem[s-distrE[0][k]->d_ind0]<<"\t"<<FSA_utils::sqrt_plus(distrE_errors[0][k]->d_elem[s-distrE_errors[0][k]->d_ind0])<<endl;
				sum_tmp+=distrE[0][k]->d_elem[s-distrE[0][k]->d_ind0]*exp((double)s*ungapped_lambda_);
				sum_tmp_error+=FSA_utils::sqrt_plus(distrE_errors[0][k]->d_elem[s-distrE_errors[0][k]->d_ind0])*exp((double)s*ungapped_lambda_);
			};
			fE<<endl;

			fE_summary<<k<<"\t"<<sum_tmp<<"\t"<<sum_tmp_error<<endl;
		};
	};


	double *K_C_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(K_C_vect);

	double *K_C_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(K_C_error_vect);


	//calculation of the ratio K/C
	if(two_dim_layer_alignment_algorithm_test_.d_scores_counts_flag&&futher_expanding_)
	{
		for(k1=0;k1<=number_of_sets_;k1++)
		{
			long int number_of_realizations_tmp=number_of_realizations_set;
			if(k1==0)
			{
				number_of_realizations_tmp=number_of_realizations_;
			};

			long int ii;
			for(ii=0;ii<=diff[k1].d_dim;ii++)
			{
				diff[k1].divide_elem_by_x(ii,(double)number_of_realizations_tmp);
				diff_errors[k1].divide_elem_by_x(ii,(double)number_of_realizations_tmp);
				diff_errors[k1].increase_elem_by_x(ii,-diff[k1].get_elem(ii)*diff[k1].get_elem(ii));
				diff_errors[k1].divide_elem_by_x(ii,(double)number_of_realizations_tmp);
				double tmp=FSA_utils::sqrt_plus(diff_errors[k1].get_elem(ii));
				diff_errors[k1].set_elem(ii,tmp);
			};

			

			sum_of_weights[k1]/=(double)number_of_realizations_tmp;
			sum_of_weights_error[k1]/=(double)number_of_realizations_tmp;
			sum_of_weights_error[k1]-=sum_of_weights[k1]*sum_of_weights[k1];
			sum_of_weights_error[k1]/=(double)number_of_realizations_tmp;
			sum_of_weights_error[k1]=FSA_utils::sqrt_plus(sum_of_weights_error[k1]);
			

			double den=0;
			double den_error=0;
			for(ii=0;ii<=diff[k1].d_dim;ii++)
			{
				double tmp=exp(-lambda_out_vect[k1]*(double)ii);
				den+=tmp*diff[k1].get_elem(ii);
				den_error+=tmp*tmp*diff_errors[k1].get_elem(ii)*diff_errors[k1].get_elem(ii);

			};


			den_error=FSA_utils::sqrt_plus(den_error);

			if(den!=0)
			{
				K_C_vect[k1]=sum_of_weights[k1]/den;
				K_C_error_vect[k1]=FSA_utils::error_of_the_ratio(sum_of_weights[k1],sum_of_weights_error[k1],den,den_error);
			}
			else
			{
				inside_simulation_flag_=false;
				K_C_vect[k1]=-1;
				K_C_error_vect[k1]=-1;
			};

		};

		K_C_=K_C_vect[0];
		K_C_error_=K_C_error_vect[0];


		if(output_flag)
		{
			fE<<endl;
			fE<<"K/C distributions\n";
			fE<<diff[0].d_dim+1<<endl;
			long int ii;
			for(ii=0;ii<=diff[0].d_dim;ii++)
			{
				//double tmp=exp(-lambda_out_*(double)ii);
				fE<<ii<<"\t"<<diff[0].get_elem(ii)/sum_of_weights[0]<<"\t"<<diff_errors[0].get_elem(ii)/sum_of_weights[0]<<endl;
			};

		};

		if(screen_output_flag_)
		{
			cout<<"K/C\t"<<K_C_;
		};

		if(number_of_sets_>1)
		{
			K_C_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			K_C_vect+1);

			if(average_set_results_flag)
			{
				K_C_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				K_C_vect+1);

				if(screen_output_flag_)
				{
					cout<<"\t"<<K_C_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<K_C_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_K_CSbs.clear();
		par_.m_KSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_K_CSbs.push_back(K_C_vect[k1]);
			par_.m_KSbs.push_back(K_C_vect[k1]*C_vect[k1]);
		};




	};


	double *a_I_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(a_I_vect);
	double *a_I_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(a_I_error_vect);
	double *a_J_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(a_J_vect);
	double *a_J_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(a_J_error_vect);
	double *sigma_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(sigma_vect);
	double *sigma_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(sigma_error_vect);
	double *alpha_I_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(alpha_I_vect);
	double *alpha_I_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(alpha_I_error_vect);
	double *alpha_J_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(alpha_J_vect);
	double *alpha_J_error_vect=new double[number_of_sets_+1];
	FSA_utils::assert_mem(alpha_J_error_vect);


	//FSC
	{
		long int k1;
		for(k1=0;k1<=number_of_sets_;k1++)
		{

			pair<long int,long int> boundaries;
			if(k1==0)
			{
				boundaries.first=0;
				boundaries.second=number_of_realizations_-1;
			}
			else
			{
				boundaries=set_number_boundaries(//boundaries of realization order numbers for the set; numeration of realizations starts from 0
				k1,//set order number; the numeration starts from 1
				number_of_realizations_set);//total number of sets
			};

			
			bool test_FSC_flag=false;
			double *test_FSC_vect=NULL;
			if(test_FSC_flag)
			{
				test_FSC_vect=new double[target_ALP_+1];
			};

			bool inside_simulation_flag;

			fsa_par::calculate_FSC(
			target_ALP_,
			boundaries.first,
			boundaries.second,
			

			lambda_out_vect[k1],

			M_max,

			distance_along_direction_1,
			distance_along_direction_2,

			ALP_weight,
			ALP_edge_max,

			a_I_vect[k1],
			a_I_error_vect[k1],
			a_J_vect[k1],
			a_J_error_vect[k1],
			sigma_vect[k1],
			sigma_error_vect[k1],
			alpha_I_vect[k1],
			alpha_I_error_vect[k1],
			alpha_J_vect[k1],
			alpha_J_error_vect[k1],
			inside_simulation_flag,
			test_FSC_vect);


			if(test_FSC_flag)
			{

				string E_distr_file_name_summary=E_distr_file_name_+"_summary_II";
				fE_summary.open(E_distr_file_name_summary.data());
				if(!fE_summary)
				{
					throw error("Error - the file "+E_distr_file_name_summary+" is not found\n",3);
				};

				fE_summary<<target_ALP_<<endl;
				long int j;
				for(j=0;j<=target_ALP_;j++)
				{
					fE_summary<<j<<"\t"<<test_FSC_vect[j]<<"\t0.0\n";
				};


				fE_summary.close();
				exit(1);
			};

			if(!inside_simulation_flag)
			{
				a_I_vect[k1]=-1;
				a_I_error_vect[k1]=-1;
				a_J_vect[k1]=-1;
				a_J_error_vect[k1]=-1;
				sigma_vect[k1]=-1;
				sigma_error_vect[k1]=-1;
				alpha_I_vect[k1]=-1;
				alpha_I_error_vect[k1]=-1;
				alpha_J_vect[k1]=-1;
				alpha_J_error_vect[k1]=-1;

				inside_simulation_flag_=false;
			};
		};

		a_I_=a_I_vect[0];
		a_I_error_=a_I_error_vect[0];
		a_J_=a_J_vect[0];
		a_J_error_=a_J_error_vect[0];
		sigma_=sigma_vect[0];
		sigma_error_=sigma_error_vect[0];
		alpha_I_=alpha_I_vect[0];
		alpha_I_error_=alpha_I_error_vect[0];
		alpha_J_=alpha_J_vect[0];
		alpha_J_error_=alpha_J_error_vect[0];


		if(screen_output_flag_)
		{
			cout<<"a_I\t"<<a_I_;
		};

		if(number_of_sets_>1)
		{
			a_I_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			a_I_vect+1);

			if(average_set_results_flag)
			{
				a_I_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				a_I_vect+1);


				if(screen_output_flag_)
				{
					cout<<"\t"<<a_I_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<a_I_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_AISbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_AISbs.push_back(a_I_vect[k1]);
		};



		if(screen_output_flag_)
		{
			cout<<"a_J\t"<<a_J_;
		};

		if(number_of_sets_>1)
		{
			a_J_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			a_J_vect+1);

			if(average_set_results_flag)
			{
				a_J_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				a_J_vect+1);


				if(screen_output_flag_)
				{
					cout<<"\t"<<a_J_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<a_J_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_AJSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_AJSbs.push_back(a_J_vect[k1]);
		};


		if(screen_output_flag_)
		{
			cout<<"sigma\t"<<sigma_;
		};

		if(number_of_sets_>1)
		{
			sigma_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			sigma_vect+1);

			if(average_set_results_flag)
			{
				sigma_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				sigma_vect+1);


				if(screen_output_flag_)
				{
					cout<<"\t"<<sigma_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<sigma_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_SigmaSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_SigmaSbs.push_back(sigma_vect[k1]);
		};



		if(screen_output_flag_)
		{
			cout<<"alpha_I\t"<<alpha_I_;
		};

		if(number_of_sets_>1)
		{
			alpha_I_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			alpha_I_vect+1);

			if(average_set_results_flag)
			{
				alpha_I_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				alpha_I_vect+1);

				if(screen_output_flag_)
				{
					cout<<"\t"<<alpha_I_;
				};
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<alpha_I_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_AlphaISbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_AlphaISbs.push_back(alpha_I_vect[k1]);
		};


		if(screen_output_flag_)
		{
			cout<<"alpha_J\t"<<alpha_J_;
		};

		if(number_of_sets_>1)
		{
			alpha_J_error_=FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
			number_of_sets_,
			alpha_J_vect+1);

			if(average_set_results_flag)
			{
				alpha_J_=FSA_utils::average(//average of elements of vect_
				number_of_sets_,
				alpha_J_vect+1);

				cout<<"\t"<<alpha_J_;
			};

			if(screen_output_flag_)
			{
				cout<<"\t"<<alpha_J_error_;
			};
		};

		if(screen_output_flag_)
		{
			cout<<endl;
		};

		par_.m_AlphaJSbs.clear();
		for(k1=1;k1<=number_of_sets_;k1++)
		{
			par_.m_AlphaJSbs.push_back(alpha_J_vect[k1]);
		};

	};



	if(output_flag)
	{
		fE.close();
		fE_summary.close();
	};


	if(score_max)
	{
		for(k1=0;k1<=number_of_sets_;k1++)
		{
			delete[]score_max[k1];
		};
		delete []score_max;
	};

	if(distrE)
	{
		for(k1=0;k1<=number_of_sets_;k1++)
		{
			for(k=0;k<=target_ALP_;k++)
			{
				delete distrE[k1][k];
			};
			delete []distrE[k1];
		};
		delete []distrE;
		distrE=NULL;
	};

	if(distrE_errors)
	{
		for(k1=0;k1<=number_of_sets_;k1++)
		{
			for(k=0;k<=target_ALP_;k++)
			{
				delete distrE_errors[k1][k];
			};
			delete []distrE_errors[k1];
		};
		delete []distrE_errors;
		distrE_errors=NULL;
	};
	


	delete[]distance_along_direction_1;
	delete[]distance_along_direction_2;
	delete[]ALP_weight;
	delete[]ALP_edge_max;

	delete[]diff;
	delete[]diff_errors;

	delete[]sum_of_weights;
	delete[]sum_of_weights_error;



	delete[]lambda_out_vect;
	delete[]lambda_out_error_vect;
	delete[]C_vect;
	delete[]C_error_vect;
	delete[]K_C_vect;
	delete[]K_C_error_vect;
	delete[]a_I_vect;
	delete[]a_I_error_vect;
	delete[]a_J_vect;
	delete[]a_J_error_vect;
	delete[]sigma_vect;
	delete[]sigma_error_vect;
	delete[]alpha_I_vect;
	delete[]alpha_I_error_vect;
	delete[]alpha_J_vect;
	delete[]alpha_J_error_vect;

	//end4

}

//tests
void test::classical_IS(

unsigned rand_,//randomization number
long int open1_,//gap opening penalty for the nucleotide sequence #1
long int open2_,//gap opening penalty for the amino acid sequence #2

long int epen1_,//gap extension penalty for the nucleotide sequence #1
long int epen2_,//gap extension penalty for the amino acid sequence #2


string smatr_file_name_,//scoring matrix file name
string RR1_file_name_,//background frequencies file name for the sequence #1
string RR2_file_name_)//background frequencies file name for the sequence #2

//long int seq1_length_,//length of sequence #1
//long int seq2_length_,//length of sequence #2

//long int seq_number_)//number of tested alignments
{
	FSA_utils::srand2(rand_);

	long int alphabet_letters_number1;//number of letters in the sequence #1
	long int alphabet_letters_number2;//number of letters in the sequence #2
	double *RR1=NULL;//background probability for the sequence #1
	double *RR2=NULL;//background probability for the sequence #2
	long int number_of_states=3;//number of states
	double **transition_probabilities=NULL;//transition probabilities between states; matrix d_number_of_states x d_number_of_states
	pair<long int, long int> *states_description=NULL;//description of the states; the index is a state number
	double ***states_distr=NULL;//distributions of the states; the index is a state number
								//the second and the third indexes correspond to an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second

	double *RR1_sum=NULL;
	long int *RR1_sum_elements=NULL;
	double *RR2_sum=NULL;
	long int *RR2_sum_elements=NULL;

	long int **smatr=NULL;

	double ungapped_lambda;

	d_IS1=NULL;
	d_IS1_general_simulation=NULL;

	long int *seq1=NULL;
	long int *seq2=NULL;

	{
		long int number_of_AA_smatr;
		long int smatr_min;

		FSA_utils::read_smatr(
		smatr_file_name_,
		smatr,
		number_of_AA_smatr,
		smatr_min);

		

		FSA_utils::read_RR(
		RR1_file_name_,
		RR1,
		RR1_sum,
		RR1_sum_elements,
		alphabet_letters_number1);


		FSA_utils::read_RR(
		RR2_file_name_,
		RR2,
		RR2_sum,
		RR2_sum_elements,
		alphabet_letters_number2);

		if(number_of_AA_smatr!=alphabet_letters_number1||number_of_AA_smatr!=alphabet_letters_number2)
		{
			throw error("Error - different numbers of letters in the files "+smatr_file_name_+", "+RR1_file_name_+", "+RR2_file_name_+"\n",1);
		};

		calculate_ungapped_lambda(
		alphabet_letters_number1,
		RR1,
		RR2,
		smatr,
		ungapped_lambda);

	};

	map<string, long int> state_name_into_number;

	state_name_into_number["S"]=0;
	state_name_into_number["D"]=1;
	state_name_into_number["I"]=2;

	map<long int,string> state_number_into_name;

	state_number_into_name[0]="S";
	state_number_into_name[1]="D";
	state_number_into_name[2]="I";


	FSA_utils::get_memory_for_matrix(number_of_states,number_of_states,transition_probabilities);

	long int min_open=FSA_utils::Tmin(open1_,open2_);
	long int min_epen=FSA_utils::Tmin(epen1_,epen2_);

	double mu=exp(-ungapped_lambda*(min_open+min_epen));
	double v=exp(-ungapped_lambda*min_epen);



	transition_probabilities[state_name_into_number["S"]][state_name_into_number["S"]]=
		((1-v)*(1-v))/((1+mu-v)*(1+mu-v));

	transition_probabilities[state_name_into_number["S"]][state_name_into_number["D"]]=
		mu/(1+mu-v);

	transition_probabilities[state_name_into_number["S"]][state_name_into_number["I"]]=
		(mu*(1-v))/((1+mu-v)*(1+mu-v));




	transition_probabilities[state_name_into_number["D"]][state_name_into_number["S"]]=
		((1-v)*(1-v))/(1+mu-v);

	transition_probabilities[state_name_into_number["D"]][state_name_into_number["D"]]=
		v;

	transition_probabilities[state_name_into_number["D"]][state_name_into_number["I"]]=
		(mu*(1-v))/(1+mu-v);




	transition_probabilities[state_name_into_number["I"]][state_name_into_number["S"]]=
		1-v;

	transition_probabilities[state_name_into_number["I"]][state_name_into_number["D"]]=
		0;

	transition_probabilities[state_name_into_number["I"]][state_name_into_number["I"]]=
		v;

	//FSA_utils::print_matrix("mtp.out",3,3,transition_probabilities);

//--------------------
	states_description=new pair<long int, long int>[number_of_states];
	FSA_utils::assert_mem(states_description);

	states_description[state_name_into_number["S"]]=make_pair(1,1);
	states_description[state_name_into_number["D"]]=make_pair(1,0);
	states_description[state_name_into_number["I"]]=make_pair(0,1);

//--------------------

	states_distr=new double **[number_of_states];

	long int s;
	for(s=0;s<number_of_states;s++)
	{
		states_distr[s]=IS1_general::allocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
			alphabet_letters_number1,//number of letters in the sequence #1
			alphabet_letters_number2,//number of letters in the sequence #2
			states_description[s]);//state description

	};

	double tmp_sum=0;
	long int i,j;
	for(i=0;i<alphabet_letters_number1;i++)
	{
		for(j=0;j<alphabet_letters_number2;j++)
		{
			states_distr[state_name_into_number["S"]][i][j]=RR1[i]*RR2[j]*exp(ungapped_lambda*smatr[i][j]);
			tmp_sum+=states_distr[state_name_into_number["S"]][i][j];
		};
	};

	//cout<<"Total sum of joint distribution of (S_A,S_B) is "<<tmp_sum<<endl;

	for(i=0;i<alphabet_letters_number1;i++)
	{
		states_distr[state_name_into_number["D"]][i][0]=RR1[i];
	};

	for(j=0;j<alphabet_letters_number2;j++)
	{
		states_distr[state_name_into_number["I"]][0][j]=RR2[j];
	};
	


	d_IS1=new IS1_general(
	alphabet_letters_number1,//number of letters in the sequence #1
	alphabet_letters_number2,//number of letters in the sequence #2
	RR1,//background probability for the sequence #1
	RR2,//background probability for the sequence #2
	number_of_states,//number of states
	transition_probabilities,//transition probabilities between states; matrix d_number_of_states x d_number_of_states
	states_description,//description of the states; the index is a state number
	states_distr);//distributions of the states; the index is a state number

//-----------------------------------
//crude sampling object

	long int number_of_states_cs=1;
	double **transition_probabilities_cs=NULL;//transition probabilities between states; matrix d_number_of_states x d_number_of_states
	pair<long int, long int> *states_description_cs=NULL;//description of the states; the index is a state number
	double ***states_distr_cs=NULL;//distributions of the states; the index is a state number

	FSA_utils::get_memory_for_matrix(number_of_states_cs,number_of_states_cs,transition_probabilities_cs);
	transition_probabilities_cs[state_name_into_number["S"]][state_name_into_number["S"]]=1;

	states_description_cs=new pair<long int, long int>[number_of_states_cs];
	FSA_utils::assert_mem(states_description_cs);

	states_description_cs[state_name_into_number["S"]]=make_pair(1,1);


	states_distr_cs=new double **[number_of_states_cs];

	for(s=0;s<number_of_states_cs;s++)
	{
		states_distr_cs[s]=IS1_general::allocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
			alphabet_letters_number1,//number of letters in the sequence #1
			alphabet_letters_number2,//number of letters in the sequence #2
			states_description_cs[s]);//state description

	};


	double tmp_sum_cs=0;
	for(i=0;i<alphabet_letters_number1;i++)
	{
		for(j=0;j<alphabet_letters_number2;j++)
		{
			states_distr_cs[state_name_into_number["S"]][i][j]=RR1[i]*RR2[j];
			tmp_sum_cs+=states_distr_cs[state_name_into_number["S"]][i][j];
		};
	};


	IS1_general* IS1_cs=new IS1_general(
	alphabet_letters_number1,//number of letters in the sequence #1
	alphabet_letters_number2,//number of letters in the sequence #2
	RR1,//background probability for the sequence #1
	RR2,//background probability for the sequence #2
	number_of_states_cs,//number of states
	transition_probabilities_cs,//transition probabilities between states; matrix d_number_of_states x d_number_of_states
	states_description_cs,//description of the states; the index is a state number
	states_distr_cs);//distributions of the states; the index is a state number


//-----------------------------------------------------------
//test for the layers object
	{

		long int max_ind1=100;//max of the index #1
		long int max_ind2=200;//max of the index #2



		long int layers_number1=6;//number of layers for the index #1
		long int layers_number2=20;//number of layers for the index #2

		two_dim_layer<long int> two_dim_layer_test(
		max_ind1,//max of the index #1
		max_ind2,//max of the index #2
		layers_number1,//number of layers for the index #1
		layers_number2,//number of layers for the index #2
		0);

	
	};
	//------------------------------------------------------
	//test for the two_dim_layer_alignment_algorithm object
	{
		long int depth1=1;//the maximum difference of the first index in the dynamic equations
		long int depth2=1;//the maximum difference of the second index in the dynamic equations

		long int max_ind1=100;//max of the index #1
		long int max_ind2=200;//max of the index #2

		long int var_num_dim=1000;
		long int *var_num=new long int[var_num_dim];
		long int i;
		for(i=0;i<var_num_dim;i++)
		{
			var_num[i]=-1;
		};

		var_num['S']=0;
		var_num['D']=1;
		var_num['I']=2;

		data_for_classical_alignment data_test;

		data_test.d_open1=open1_;
		data_test.d_open2=open2_;

		data_test.d_epen1=epen1_;
		data_test.d_epen2=epen2_;

		data_test.d_smatr=smatr;



		long int initial_state=state_name_into_number["S"];

		d_IS1_general_simulation=new IS1_general_simulation(
		d_IS1,
		initial_state,//initial state for the IS
		max_ind1*5,//maximum sequence length
		max_ind2*5);//maximum sequence length


		data_test.d_seq1=d_IS1_general_simulation->d_seq1;//sequence #1
		data_test.d_seq2=d_IS1_general_simulation->d_seq2;//sequence #2


		



		two_dim_layer_alignment_algorithm<long int> two_dim_layer_alignment_algorithm_test(
		number_of_states,//total number of variables in dynamic equations
		depth1,//the maximum difference of the first index in the dynamic equations
		depth2,//the maximum difference of the second index in the dynamic equations
		max_ind1,//max of the index #1 (minimum index is 0 by default)
		max_ind2,//max of the index #2 (minimum index is 0 by default)
		0,//null element of T
		var_num);//d_variable_name_to_number[c] converts a letter 'c' (of type char) to order number of variable

		two_dim_layer_alignment_algorithm_test.d_two_dim_layer_alignment_function=test::two_dim_layer_alignment_function_classical_global;
		two_dim_layer_alignment_algorithm_test.d_two_dim_layer_boundary_function=test::two_dim_layer_alignment_function_classical_global;
		two_dim_layer_alignment_algorithm_test.d_par=&data_test;

		long int target_seq1_length=100;
		long int target_seq2_length=200;

		d_IS1_general_simulation->simulate_upto_target_lengths(
			target_seq1_length,
			target_seq2_length);



		two_dim_layer_alignment_algorithm_test.align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
		//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
		target_seq1_length,//target length of the sequence #1
		target_seq2_length);//target length of the sequence #2


		delete[]var_num;
	};



	//test for the two_dim_layer_alignment_algorithm object
	//crude sampling
	{
		long int number_of_realizations=100;

		long int depth1=1;//the maximum difference of the first index in the dynamic equations
		long int depth2=1;//the maximum difference of the second index in the dynamic equations

		long int max_ind1=100;//max of the index #1
		long int max_ind2=100;//max of the index #2

		long int var_num_dim=1000;
		long int *var_num=new long int[var_num_dim];
		long int i;
		for(i=0;i<var_num_dim;i++)
		{
			var_num[i]=-1;
		};

		var_num['S']=0;
		var_num['D']=1;
		var_num['I']=2;

		data_for_classical_alignment data_test;

		data_test.d_open1=open1_;
		data_test.d_open2=open2_;

		data_test.d_epen1=epen1_;
		data_test.d_epen2=epen2_;

		data_test.d_smatr=smatr;



		long int initial_state=state_name_into_number["S"];

		IS1_general_simulation *IS1_general_simulation_cs=new IS1_general_simulation(
		IS1_cs,
		initial_state,//initial state for the IS
		max_ind1,//maximum sequence length
		max_ind2);//maximum sequence length


		data_test.d_seq1=IS1_general_simulation_cs->d_seq1;//sequence #1
		data_test.d_seq2=IS1_general_simulation_cs->d_seq2;//sequence #2


		



		two_dim_layer_alignment_algorithm<long int> two_dim_layer_alignment_algorithm_test_cs(
		number_of_states,//total number of variables in dynamic equations
		depth1,//the maximum difference of the first index in the dynamic equations
		depth2,//the maximum difference of the second index in the dynamic equations
		max_ind1,//max of the index #1 (minimum index is 0 by default)
		max_ind2,//max of the index #2 (minimum index is 0 by default)
		0,//null element of T
		var_num);//d_variable_name_to_number[c] converts a letter 'c' (of type char) to order number of variable

		two_dim_layer_alignment_algorithm_test_cs.d_two_dim_layer_alignment_function=test::two_dim_layer_alignment_function_classical_global;
		two_dim_layer_alignment_algorithm_test_cs.d_two_dim_layer_boundary_function=test::two_dim_layer_alignment_function_classical_global;
		two_dim_layer_alignment_algorithm_test_cs.d_par=&data_test;

		two_dim_layer_alignment_algorithm_test_cs.d_M_flag=true;

		long int target_seq1_length=max_ind1;
		long int target_seq2_length=max_ind2;

		ofstream fM;

		array_v<double> *distrM=NULL;
		long int score_max=-inf;

		if(two_dim_layer_alignment_algorithm_test_cs.d_M_flag)
		{
			//distribution of M
			string M_distr_file_name="distr_M_tmp.out";

			fM.open(M_distr_file_name.data());
			if(!fM)
			{
				throw error("Error - the file "+M_distr_file_name+" is not found\n",3);
			};


			distrM=new array_v<double>(NULL);

			

		};

		//--------------------

		long int k;
		for(k=1;k<=number_of_realizations;k++)
		{

			IS1_general_simulation_cs->init();
			IS1_general_simulation_cs->simulate_upto_target_lengths(
			target_seq1_length,
			target_seq2_length);

			two_dim_layer_alignment_algorithm_test_cs.init();
			two_dim_layer_alignment_algorithm_test_cs.align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
			//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
			target_seq1_length,//target length of the sequence #1
			target_seq2_length);//target length of the sequence #2

			if(two_dim_layer_alignment_algorithm_test_cs.d_M_flag)
			{
				score_max=FSA_utils::Tmax(score_max,two_dim_layer_alignment_algorithm_test_cs.d_M);
				distrM->increase_elem_by_1(two_dim_layer_alignment_algorithm_test_cs.d_M);
			};

			if(k%1000==0)
			{
				cout<<k<<endl;
			};

		};

		if(two_dim_layer_alignment_algorithm_test_cs.d_M_flag)
		{
			long int s;

			double sum_tmp=0;
			for(s=0;s<=score_max;s++)
			{
				sum_tmp+=distrM->d_elem[s-distrM->d_ind0];
			};

			if(sum_tmp<=0)
			{
				throw error("Unexpected error - sum_tmp<=0\n",1);
			};
			

			fM<<score_max+1<<endl;
			for(s=0;s<=score_max;s++)
			{
				distrM->d_elem[s-distrM->d_ind0]/=sum_tmp;
				fM<<s<<"\t"<<distrM->d_elem[s-distrM->d_ind0]<<endl;
			};

			fM.close();
		};



		delete[]var_num;
		delete IS1_general_simulation_cs;
		delete distrM;
	};


	delete[]RR1;
	delete[]RR2;

	delete[]RR1_sum;
	delete[]RR2_sum;

	delete[]RR1_sum_elements;
	delete[]RR2_sum_elements;

	FSA_utils::delete_memory_for_matrix(alphabet_letters_number1,smatr);
	FSA_utils::delete_memory_for_matrix(number_of_states,transition_probabilities);
	FSA_utils::delete_memory_for_matrix(number_of_states_cs,transition_probabilities_cs);

	for(s=0;s<number_of_states;s++)
	{
		IS1_general::deallocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
			states_distr[s],
			alphabet_letters_number1,//number of letters in the sequence #1
			states_description[s]);//state description
	};

	for(s=0;s<number_of_states_cs;s++)
	{
		IS1_general::deallocate_state_distribution(//allocates an array of dimensions d_A_letters^state_description_type.first x d_B_letters^state_description_type.second
			states_distr_cs[s],
			alphabet_letters_number1,//number of letters in the sequence #1
			states_description_cs[s]);//state description
	};


	delete[]states_description;

	delete[]states_description_cs;
	

	delete[]seq1;
	delete[]seq2;

	delete d_IS1;

	delete IS1_cs;

	delete d_IS1_general_simulation;
	
}

double test::lambda_equation(double x_,void* func_number_)
{
	data_for_lambda_equation *data=(data_for_lambda_equation*)func_number_;
	long int d_number_of_AA=data->d_number_of_AA;
	long int** d_smatr=data->d_smatr;
	double *d_RR1=data->d_RR1;
	double *d_RR2=data->d_RR2;

	double res=0;
	long int i,j;

	for(i=0;i<d_number_of_AA;i++)
	{
		for(j=0;j<d_number_of_AA;j++)
		{
			res+=d_RR1[i]*d_RR2[j]*exp(x_*d_smatr[i][j]);
		};
	};

	return res-1.0;
}

double test::lambda_equation_FSA(double x_,void* func_number_)
{
	data_for_lambda_equation_FSA *data=(data_for_lambda_equation_FSA*)func_number_;
	long int number_of_letters1=data->d_number_of_letters1;
	long int number_of_letters2=data->d_number_of_letters2;

	long int** d_smatr=data->d_smatr;
	double *d_RR1=data->d_RR1;
	double *d_RR2=data->d_RR2;

	long int codon_length=data->d_codon_length;//codon length 
	long int *codon_AA=data->d_codon_AA;//<codon code,AA number>

	long int *codon=new long int [codon_length];
	FSA_utils::assert_mem(codon);

	long int i,j;

	long int code_n=1;
	for(i=1;i<=codon_length;i++)
	{
		code_n*=number_of_letters1;
	};

	double res=0;
	

	for(i=0;i<code_n;i++)
	{
		FSA_utils::convert_code_into_codon(
		i,//the input code
		codon_length,//codon length 
		number_of_letters1,//number of letters for the sequence 1
		codon);//must be allocated

		long int AA1=codon_AA[i];

		double RR1_mult=1;
		
		for(j=0;j<codon_length;j++)
		{
			RR1_mult*=d_RR1[codon[j]];
		};

		for(j=0;j<number_of_letters2;j++)
		{
			res+=RR1_mult*d_RR2[j]*exp(x_*d_smatr[AA1][j]);
		};
	};

	delete[]codon;
	return res-1.0;
}


void test::calculate_ungapped_lambda(
long int number_of_AA_,
double *RR1_,
double *RR2_,
long int**smatr_,
double &ungapped_lambda_)
{
	//calculation of the importance sampling theta

	data_for_lambda_equation tmp_ptr;
	tmp_ptr.d_number_of_AA=number_of_AA_;
	tmp_ptr.d_RR1=RR1_;
	tmp_ptr.d_RR2=RR2_;
	tmp_ptr.d_smatr=smatr_;

	//calculate maximum of smatr_ elements
	long int smatr_max=smatr_[0][0];
	long int smatr_max_i=0;
	long int smatr_max_j=0;
	long int smatr_min=smatr_[0][0];

	long int smatr_pos_max=LONG_MIN;
	long int smatr_neg_min=LONG_MAX;

	double eps=0.0000001;
	double threshold=DBL_MIN*10.0;

	double aver_score=0;
	long int i,j;
	for(i=0;i<number_of_AA_;i++)
	{
		for(j=0;j<number_of_AA_;j++)
		{
			if(RR1_[i]*RR2_[j]<=threshold)
			{
				continue;
			};
								


			aver_score+=RR1_[i]*RR2_[j]*smatr_[i][j];

			if(smatr_max<smatr_[i][j])
			{
				smatr_max=smatr_[i][j];
				smatr_max_i=i;
				smatr_max_j=j;
			};
			smatr_min=FSA_utils::Tmin(smatr_min,smatr_[i][j]);
			

			if(smatr_[i][j]>0)
			{
				smatr_pos_max=FSA_utils::Tmax(smatr_pos_max,smatr_[i][j]);
			};

			if(smatr_[i][j]<0)
			{
				smatr_neg_min=FSA_utils::Tmin(smatr_neg_min,smatr_[i][j]);
			};

		};
	};

	if(aver_score>=-threshold)
	{
		throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
	};

	if(smatr_max<=0)
	{
		throw error("Error - at least one element of the scoring matrix must be positive\n",3);
	};

	

	double a=eps;

	while(test::lambda_equation(a,(void*)(&tmp_ptr))>0)
	{
		a/=2.0;

		if(a<threshold*100.0)
		{
			throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
		};
	};

	if(a<threshold*100.0)
	{
		throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
	};

	eps=a/10.0;


	double tmp_pr=RR1_[smatr_max_i]*RR2_[smatr_max_j];
	double b=(log(1+10*eps)-log(tmp_pr))/(double)smatr_max;

	
	long int n_partition=2;
	std::vector<double> res_lambda;
	
	
	alp_reg::find_tetta_general(
	test::lambda_equation,
	(void*)(&tmp_ptr),
	a,
	b,
	n_partition,
	eps,
	res_lambda);

	sort(res_lambda.begin(),res_lambda.end());

	if(res_lambda.size()==0)
	{
		throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
	};

	
	ungapped_lambda_=res_lambda[res_lambda.size()-1];

	cout<<"Ungapped lambda is "<<ungapped_lambda_<<endl;


}

void test::calculate_ungapped_lambda_general(
function_type *func_,
void* func_pointer_,
double &ungapped_lambda_)
{

	double eps=1e-7;
	double threshold=DBL_MIN*10.0;

	double threshold2=DBL_MAX/4;


	double a=eps;
	while(func_(a,func_pointer_)>0)
	{
		a/=2.0;

		if(a<threshold*100.0)
		{
			throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
		};
	};

	if(a<threshold*100.0)
	{
		throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
	};

	eps=a/10.0;


	double b=a;


	while(func_(b,func_pointer_)<=0)
	{
		b*=2.0;

		if(b>threshold2)
		{
			throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
		};
	};


	
	long int n_partition=2;
	std::vector<double> res_lambda;
	
	
	alp_reg::find_tetta_general(
	func_,
	func_pointer_,
	a,
	b,
	n_partition,
	eps,
	res_lambda);

	sort(res_lambda.begin(),res_lambda.end());

	if(res_lambda.size()==0)
	{
		throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
	};

	
	ungapped_lambda_=res_lambda[res_lambda.size()-1];

	//cout<<setprecision(12)<<"Ungapped lambda is "<<ungapped_lambda_<<endl;
	//cout<<"Ungapped lambda is "<<ungapped_lambda_<<endl;


}

//======================================================

void test::two_dim_layer_alignment_function_FSA(
long int i1_,
long int i2_,
two_dim_layer_alignment_algorithm<long int> *obj_,
void*par_)
{
	data_for_FSA_alignment &data=*(data_for_FSA_alignment*)par_;

	long int &alphabet_letters_number1=data.d_alphabet_letters_number1;//number of letters in the sequence #1

	long int &open1=data.d_open1;//gap opening penalty for the nucleotide sequence #1
	long int &open2=data.d_open2;//gap opening penalty for the amino acid sequence #2

	long int &epen1=data.d_epen1;//gap extension penalty for the nucleotide sequence #1
	long int &epen2=data.d_epen2;//gap extension penalty for the amino acid sequence #2

	long int &gamma=data.d_gamma;//frame shift penalty

	long int**&smatr=data.d_smatr;//the scoring matrix
	long int *&seq1=data.d_seq1;//element of sequence #1
	long int *&seq2=data.d_seq2;//element of sequence #2

	long int *&codon_AA=data.d_codon_AA;//<codon code,AA number>

	bool &insertions_after_deletions=data.d_insertions_after_deletions;//if true, then insertions after deletions are allowed

	string &alignment_type=data.d_alignment_type;//possible values are "global" or "local"

	if(obj_->d_number_of_variables!=3)
	{
		throw error("Error - obj_->d_number_of_variables!=3 in two_dim_layer_alignment_function_FSA\n",1);
	};

	long int *d_var_num=obj_->d_var_num;

	two_dim_layer<long int> **vars=obj_->d_vars;

	long int Sn=d_var_num['S'];
	long int Dn=d_var_num['D'];
	long int In=d_var_num['I'];

	long int codon_length=3;//codon length 

//----------------------------------------

	bool local_flag=(alignment_type=="local");

	long int init_for_S=-inf;
	if(local_flag)
	{
		init_for_S=0;
	};


	long int i=i1_;
	long int j=i2_;

	if(j==0)
	{
		if(local_flag)
		{
			vars[Dn]->set_element(i,0,-inf);
		}
		else
		{
			vars[Dn]->set_element(i,0,-inf);
			
		};


		
		vars[In]->set_element(i,0,-inf);

		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//if(i<=2)
		if(i<=0)
		{
			vars[Sn]->set_element(i,0,0);
		}
		else
		{
			vars[Sn]->set_element(i,0,init_for_S);
		};


		return;

	};

	if(i==0)
	{

		vars[Sn]->set_element(0,j,init_for_S);
		vars[Dn]->set_element(0,j,-inf);
		if(local_flag)
		{
			vars[In]->set_element(0,j,-inf);
		}
		else
		{
			vars[In]->set_element(0,j,-inf);
		};

		return;
	};

	if(i==1)
	{
		//check these conditions
		vars[Sn]->set_element(1,j,init_for_S);
		vars[Dn]->set_element(1,j,-inf);
		if(local_flag)
		{
			vars[In]->set_element(1,j,-inf);
		}
		else
		{
			vars[In]->set_element(1,j,-inf);
		};

		return;
	};

	if(i==2)
	{
		vars[Sn]->set_element(2,j,init_for_S);
		vars[Dn]->set_element(2,j,-inf);
		if(local_flag)
		{
			vars[In]->set_element(2,j,-inf);
		}
		else
		{
			vars[In]->set_element(2,j,-inf);
		};

		return;
	};


	{

		{

			if(insertions_after_deletions)
			{
				vars[In]->set_element(i,j,
					FSA_utils::Tmax(
					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					vars[Dn]->get_element(i,j-1)-open2-epen2,
					vars[Sn]->get_element(i,j-1)-open2-epen2,vars[In]->get_element(i,j-1)-epen2));
			}
			else
			{
				vars[In]->set_element(i,j,
					FSA_utils::Tmax(
					vars[Sn]->get_element(i,j-1)-open2-epen2,vars[In]->get_element(i,j-1)-epen2));

			};
			
			if(i==3)
			{

				long int AA1=FSA_utils::convert_codon_into_AA(
				codon_length,//codon length 
				codon_AA,//<codon code,AA number>
				alphabet_letters_number1,//number of letters for the sequence 1
				seq1+i-3);


				long int smart_score=smatr[AA1][seq2[j-1]];
				vars[Sn]->set_element(i,j,FSA_utils::Tmax(

					//for local alignment
					init_for_S,
					
					FSA_utils::Tmax(vars[Sn]->get_element(i-3,j-1),vars[Dn]->get_element(i-3,j-1),vars[In]->get_element(i-3,j-1))+smart_score,
					FSA_utils::Tmax(vars[Sn]->get_element(i-2,j-1),vars[Dn]->get_element(i-2,j-1),vars[In]->get_element(i-2,j-1))+smart_score-gamma
					
					)
				);

				vars[Dn]->set_element(i,j,FSA_utils::Tmax(
					
					FSA_utils::Tmax(vars[Sn]->get_element(i-3,j)-open1-epen1,vars[Dn]->get_element(i-3,j)-epen1),
					FSA_utils::Tmax(vars[Sn]->get_element(i-2,j)-open1-epen1,vars[Dn]->get_element(i-2,j)-epen1)-gamma

					)
					);
				return;
				
			};
			if(i>=4)
			{

				long int AA1=FSA_utils::convert_codon_into_AA(
				codon_length,//codon length 
				codon_AA,//<codon code,AA number>
				alphabet_letters_number1,//number of letters for the sequence 1
				seq1+i-3);

				long int smart_score=smatr[AA1][seq2[j-1]];
				vars[Sn]->set_element(i,j,FSA_utils::Tmax(
					
					//for local alignment
					init_for_S,
					
					FSA_utils::Tmax(vars[Sn]->get_element(i-3,j-1),vars[Dn]->get_element(i-3,j-1),vars[In]->get_element(i-3,j-1))+smart_score,
					
					FSA_utils::Tmax(
					FSA_utils::Tmax(vars[Sn]->get_element(i-2,j-1),vars[Dn]->get_element(i-2,j-1),vars[In]->get_element(i-2,j-1)),
					FSA_utils::Tmax(vars[Sn]->get_element(i-4,j-1),vars[Dn]->get_element(i-4,j-1),vars[In]->get_element(i-4,j-1))
					)+smart_score-gamma
					
					)
					);



				vars[Dn]->set_element(i,j,FSA_utils::Tmax(
					
					FSA_utils::Tmax(vars[Sn]->get_element(i-3,j)-open1-epen1,vars[Dn]->get_element(i-3,j)-epen1),
					FSA_utils::Tmax(vars[Sn]->get_element(i-2,j)-open1-epen1,vars[Dn]->get_element(i-2,j)-epen1)-gamma,
					FSA_utils::Tmax(vars[Sn]->get_element(i-4,j)-open1-epen1,vars[Dn]->get_element(i-4,j)-epen1)-gamma

					)
					);

			};



		};
	};


}

//======================================================


void test::two_dim_layer_alignment_function_classical_global(
long int i1_,
long int i2_,
two_dim_layer_alignment_algorithm<long int> *obj_,
void*par_)
{
	data_for_classical_alignment data=*(data_for_classical_alignment*)par_;

	long int &open1=data.d_open1;//gap opening penalty for the nucleotide sequence #1
	long int &open2=data.d_open2;//gap opening penalty for the amino acid sequence #2

	long int &epen1=data.d_epen1;//gap extension penalty for the nucleotide sequence #1
	long int &epen2=data.d_epen2;//gap extension penalty for the amino acid sequence #2

	long int**&smatr=data.d_smatr;//the scoring matrix
	long int *seq1=data.d_seq1;//element of sequence #1
	long int *seq2=data.d_seq2;//element of sequence #2

	if(obj_->d_number_of_variables!=3)
	{
		throw error("Error - obj_->d_number_of_variables!=3 in two_dim_layer_alignment_function_classical_global\n",1);
	};

	long int *d_var_num=obj_->d_var_num;

	two_dim_layer<long int> **vars=obj_->d_vars;

	long int Sn=d_var_num['S'];
	long int Dn=d_var_num['D'];
	long int In=d_var_num['I'];

	//boundary conditions
	if(i1_==0)
	{
		if(i2_==0)
		{
			vars[Sn]->set_element(0,0,0);
			vars[Dn]->set_element(0,0,-open1);
		}
		else
		{
			vars[Sn]->set_element(0,i2_,-inf);
			vars[Dn]->set_element(0,i2_,-inf);
		};

		vars[In]->set_element(0,i2_,-open2-i2_*epen2);
		return;
	};

	if(i2_==0)
	{
		if(i1_==0)
		{
			vars[Sn]->set_element(0,0,0);
			vars[In]->set_element(0,0,-open2);
		}
		else
		{
			vars[Sn]->set_element(i1_,0,-inf);
			vars[In]->set_element(i1_,0,-inf);
		};

		vars[Dn]->set_element(i1_,0,-open1-i1_*epen1);
		return;
	};

	//here i1_>0, i2_>0
	long int S_1_1=vars[Sn]->get_element(i1_-1,i2_-1);
	long int D_1_1=vars[Dn]->get_element(i1_-1,i2_-1);
	long int I_1_1=vars[In]->get_element(i1_-1,i2_-1);


	long int S_1_0=vars[Sn]->get_element(i1_-1,i2_);
	long int D_1_0=vars[Dn]->get_element(i1_-1,i2_);
	//long int I_1_0=vars[In]->get_element(i1_-1,i2_);

	long int S_0_1=vars[Sn]->get_element(i1_,i2_-1);
	long int D_0_1=vars[Dn]->get_element(i1_,i2_-1);
	long int I_0_1=vars[In]->get_element(i1_,i2_-1);

	long int S_0_0=FSA_utils::Tmax(S_1_1,D_1_1,I_1_1)+smatr[seq1[i1_-1]][seq2[i2_-1]];
	long int D_0_0=FSA_utils::Tmax(S_1_0-open1-epen1,D_1_0-epen1);
	long int I_0_0=FSA_utils::Tmax(S_0_1-open2-epen2,I_0_1-epen2,D_0_1-open2-epen2);

	vars[Sn]->set_element(i1_,i2_,S_0_0);
	vars[Dn]->set_element(i1_,i2_,D_0_0);
	vars[In]->set_element(i1_,i2_,I_0_0);


}

//--------------------------------------
template<typename T> 
two_dim_layer_alignment_algorithm<T>::two_dim_layer_alignment_algorithm(
long int number_of_variables_,//total number of variables in dynamic equations
long int depth1_,//the maximum difference of the first index in the dynamic equations
long int depth2_,//the maximum difference of the second index in the dynamic equations
long int max_ind1_,//max of the index #1 (minimum index is 0 by default)
long int max_ind2_,//max of the index #2 (minimum index is 0 by default)
T null_,//null element of T
long int *var_num_)//d_variable_name_to_number[c] converts a letter 'c' (of type char) to order number of variable
{

	d_max_ind1=max_ind1_;
	d_max_ind2=max_ind2_;

	d_number_of_variables=number_of_variables_;
	d_depth1=depth1_;
	d_depth2=depth2_;
	d_var_num=var_num_;
	d_null=null_;
	d_step1=depth1_*2+1;
	d_step2=depth2_*2+1;

	//d_step1=depth1_;
	//d_step2=depth2_;

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	d_edge_maximum_calculation_depth1=3;
	d_edge_maximum_calculation_depth2=1;

	long int layers_number1=depth1_+d_step1;//number of layers for the index #1
	long int layers_number2=depth2_+d_step2;//number of layers for the index #2

	d_vars=new two_dim_layer<T> *[number_of_variables_];

	long int i;
	for(i=0;i<d_number_of_variables;i++)
	{
		d_vars[i]=new two_dim_layer<T>(
		max_ind1_,//max of the index #1 (minimum index is 0 by default)
		max_ind2_,//max of the index #2 (minimum index is 0 by default)
		layers_number1,//number of layers for the index #1
		layers_number2,//number of layers for the index #2
		null_);
	};

	d_current_align_ind1=-1;
	d_current_align_ind2=-1;

	//statistics 
	d_M_flag=false;
	d_M=-inf;

	d_E_flag=false;
	d_E=-inf;

}


template<typename T> 
two_dim_layer_alignment_algorithm<T>::two_dim_layer_alignment_algorithm(
long int max_ind1_,//max of the index #1
long int max_ind2_,//max of the index #2
two_dim_layer_alignment_algorithm<T> *two_dim_layer_alignment_algorithm_)
{

	d_max_ind1=two_dim_layer_alignment_algorithm_->d_max_ind1;
	d_max_ind2=two_dim_layer_alignment_algorithm_->d_max_ind2;

	d_number_of_variables=two_dim_layer_alignment_algorithm_->d_number_of_variables;
	d_depth1=two_dim_layer_alignment_algorithm_->d_depth1;
	d_depth2=two_dim_layer_alignment_algorithm_->d_depth2;
	d_var_num=two_dim_layer_alignment_algorithm_->d_var_num;
	d_null=two_dim_layer_alignment_algorithm_->d_null;
	d_step1=two_dim_layer_alignment_algorithm_->d_step1;
	d_step2=two_dim_layer_alignment_algorithm_->d_step2;


	d_edge_maximum_calculation_depth1=two_dim_layer_alignment_algorithm_->d_edge_maximum_calculation_depth1;
	d_edge_maximum_calculation_depth2=two_dim_layer_alignment_algorithm_->d_edge_maximum_calculation_depth2;

	d_vars=new two_dim_layer<T> *[d_number_of_variables];

	long int i;
	for(i=0;i<d_number_of_variables;i++)
	{
		d_vars[i]=new two_dim_layer<T>(
		max_ind1_,//max of the index #1 (minimum index is 0 by default)
		max_ind2_,//max of the index #2 (minimum index is 0 by default)
		two_dim_layer_alignment_algorithm_->d_vars[i]);
	};

	d_current_align_ind1=two_dim_layer_alignment_algorithm_->d_current_align_ind1;
	d_current_align_ind2=two_dim_layer_alignment_algorithm_->d_current_align_ind2;

	//statistics 
	d_M_flag=two_dim_layer_alignment_algorithm_->d_M_flag;
	d_M=two_dim_layer_alignment_algorithm_->d_M;

	d_E_flag=two_dim_layer_alignment_algorithm_->d_E_flag;
	d_E=two_dim_layer_alignment_algorithm_->d_E;

}


template<typename T> 
two_dim_layer_alignment_algorithm<T>::~two_dim_layer_alignment_algorithm()
{
	long int i;
	for(i=0;i<d_number_of_variables;i++)
	{
		d_vars[i]->~two_dim_layer<T>();
	};

	delete[]d_vars;
}

template<typename T> 
void two_dim_layer_alignment_algorithm<T>::init()//initialization for a new alignment
{
	d_M=-inf;
	d_E=-inf;
	d_current_align_ind1=-1;
	d_current_align_ind2=-1;

	d_scores_counts.clear();
}

template<typename T> 
void two_dim_layer_alignment_algorithm<T>::align_upto_target_lengths_step(//aligns the sequences upto the target sequence lengths including 
//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
long int target_align1_length_,//target length for the sequence #1
long int target_align2_length_)//target length for the sequence #2
{
	long int i;
	for(i=0;i<d_number_of_variables;i++)
	{
		d_vars[i]->set_max_ind(target_align1_length_,target_align2_length_);
	};

	if(target_align1_length_-d_current_align_ind1>d_step1||target_align2_length_-d_current_align_ind2>d_step2)
	{
		throw error("Unexpected error in two_dim_layer_alignment_algorithm<T>::align_upto_target_lengths: target_align1_length_-d_current_align_ind1>d_step1||target_align2_length_-d_current_align_ind2>d_step2\n",1);
	};

	long int i1,i2;
	for(i1=0;i1<=d_current_align_ind1;i1++)
	{
		for(i2=d_current_align_ind2+1;i2<=target_align2_length_;i2++)
		{
			two_dim_layer_alignment_function_type *tmp_func=d_two_dim_layer_alignment_function;
			if(i1<d_depth1||i2<d_depth2)
			{
				tmp_func=d_two_dim_layer_boundary_function;
			};

			tmp_func(i1,i2,this,d_par);
		};
	};

	for(i2=0;i2<=d_current_align_ind2;i2++)
	{
		for(i1=d_current_align_ind1+1;i1<=target_align1_length_;i1++)
		{
			two_dim_layer_alignment_function_type *tmp_func=d_two_dim_layer_alignment_function;
			if(i1<d_depth1||i2<d_depth2)
			{
				tmp_func=d_two_dim_layer_boundary_function;
			};

			tmp_func(i1,i2,this,d_par);
		};
	};

	for(i1=d_current_align_ind1+1;i1<=target_align1_length_;i1++)
	{
		for(i2=d_current_align_ind2+1;i2<=target_align2_length_;i2++)
		{
			two_dim_layer_alignment_function_type *tmp_func=d_two_dim_layer_alignment_function;
			if(i1<d_depth1||i2<d_depth2)
			{
				tmp_func=d_two_dim_layer_boundary_function;
			};

			tmp_func(i1,i2,this,d_par);
		};
	};

	d_current_align_ind1=FSA_utils::Tmax(target_align1_length_,d_current_align_ind1);
	d_current_align_ind2=FSA_utils::Tmax(target_align2_length_,d_current_align_ind2);

}

template<typename T> 
void two_dim_layer_alignment_algorithm<T>::align_upto_target_lengths(//aligns the sequences upto the target sequence lengths including 
//the sequences #1, #2 must be defined upto target_seq1_length_-1 and target_seq2_length_-1 lengths 
long int target_align1_length_,//target length for the sequence #1
long int target_align2_length_)//target length for the sequence #2
{

	while(d_current_align_ind1<target_align1_length_||d_current_align_ind2<target_align2_length_)
	{
		long int current_align_ind1_old=d_current_align_ind1;
		long int current_align_ind2_old=d_current_align_ind2;

		align_upto_target_lengths_step(
		FSA_utils::Tmin(target_align1_length_,d_current_align_ind1+d_step1),
		FSA_utils::Tmin(target_align2_length_,d_current_align_ind2+d_step2));

		long int M_tmp=d_M;

		//statistics
		if(d_M_flag||d_FSC_flag)
		{
			long int i1,i2,j;
			for(i1=0;i1<=d_current_align_ind1;i1++)
			{
				for(i2=current_align_ind2_old+1;i2<=d_current_align_ind2;i2++)
				{
					for(j=0;j<d_number_of_variables;j++)
					{
						d_M=FSA_utils::Tmax(d_M,d_vars[j]->get_element(i1,i2));
					};
				};
			};

			for(i1=current_align_ind1_old+1;i1<=d_current_align_ind1;i1++)
			{
				for(i2=0;i2<=current_align_ind2_old;i2++)
				{
					for(j=0;j<d_number_of_variables;j++)
					{
						d_M=FSA_utils::Tmax(d_M,d_vars[j]->get_element(i1,i2));
					};
				};
			};

		};

		if(d_scores_counts_flag)
		{
			long int i1,i2,j;
			for(i1=0;i1<=d_current_align_ind1;i1++)
			{
				for(i2=current_align_ind2_old+1;i2<=d_current_align_ind2;i2++)
				{
					long int M_local=d_vars[0]->get_element(i1,i2);
					for(j=1;j<d_number_of_variables;j++)
					{
						M_local=FSA_utils::Tmax(M_local,d_vars[j]->get_element(i1,i2));
					};


					if(M_local>=d_M_threshold)
					{
						d_scores_counts.increase_elem_by_1(M_local);
					};
				};
			};

			for(i1=current_align_ind1_old+1;i1<=d_current_align_ind1;i1++)
			{
				for(i2=0;i2<=current_align_ind2_old;i2++)
				{
					long int M_local=d_vars[0]->get_element(i1,i2);
					for(j=1;j<d_number_of_variables;j++)
					{
						M_local=FSA_utils::Tmax(M_local,d_vars[j]->get_element(i1,i2));
					};

					if(M_local>=d_M_threshold)
					{
						d_scores_counts.increase_elem_by_1(M_local);
					};
				};
			};

		};




		if(d_FSC_flag&&M_tmp<d_M)
		{//calculated only along the edge containing the next ALP

			d_distance_along_direction_1=-1;

			d_distance_along_direction_2=-1;

			long int i1,i2,j;
			for(i1=0;i1<=d_current_align_ind1;i1++)
			{
				for(i2=current_align_ind2_old+1;i2<=d_current_align_ind2;i2++)
				{
					long int M_local=d_vars[0]->get_element(i1,i2);
					for(j=1;j<d_number_of_variables;j++)
					{
						M_local=FSA_utils::Tmax(M_local,d_vars[j]->get_element(i1,i2));
					};

					if(M_local==d_M)
					{
						if(d_distance_along_direction_1==-1)
						{
							d_distance_along_direction_1=i1;
						}
						else
						{
							d_distance_along_direction_1=FSA_utils::Tmin(d_distance_along_direction_1,i1);
						};

						if(d_distance_along_direction_2==-1)
						{
							d_distance_along_direction_2=i2;
						}
						else
						{
							d_distance_along_direction_2=FSA_utils::Tmin(d_distance_along_direction_2,i2);
						};
					};

				};
			};

			for(i1=current_align_ind1_old+1;i1<=d_current_align_ind1;i1++)
			{
				for(i2=0;i2<=d_current_align_ind2;i2++)
				{
					long int M_local=d_vars[0]->get_element(i1,i2);
					for(j=1;j<d_number_of_variables;j++)
					{
						M_local=FSA_utils::Tmax(M_local,d_vars[j]->get_element(i1,i2));
					};

					if(M_local==d_M)
					{
						if(d_distance_along_direction_1==-1)
						{
							d_distance_along_direction_1=i1;
						}
						else
						{
							d_distance_along_direction_1=FSA_utils::Tmin(d_distance_along_direction_1,i1);
						};

						if(d_distance_along_direction_2==-1)
						{
							d_distance_along_direction_2=i2;
						}
						else
						{
							d_distance_along_direction_2=FSA_utils::Tmin(d_distance_along_direction_2,i2);
						};
					};
				};
			};


			if(d_distance_along_direction_1==-1||d_distance_along_direction_2==-1)
			{
				throw error("Unexpected error - d_distance_along_direction_1==-1||d_distance_along_direction_2==-1\n",1);
			};
		};



	};


	//statistics
	if(d_E_flag||d_scores_counts_flag)
	{

		d_E=-inf;
		long int i1,i2,j;
		for(i1=0;i1<=d_current_align_ind1;i1++)
		{
			for(i2=d_current_align_ind2-d_edge_maximum_calculation_depth2+1;i2<=d_current_align_ind2;i2++)
			{
				for(j=0;j<d_number_of_variables;j++)
				{
					d_E=FSA_utils::Tmax(d_E,d_vars[j]->get_element(i1,i2));
				};
			};
		};

		for(i1=d_current_align_ind1-d_edge_maximum_calculation_depth1+1;i1<=d_current_align_ind1;i1++)
		{
			for(i2=0;i2<=d_current_align_ind2-d_edge_maximum_calculation_depth2;i2++)
			{
				for(j=0;j<d_number_of_variables;j++)
				{
					d_E=FSA_utils::Tmax(d_E,d_vars[j]->get_element(i1,i2));
				};
			};
		};




	};


}
//---------------------------------------------------------

