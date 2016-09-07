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

File name: sls_fsa1_utils.cpp

Author: Sergey Sheetlin

Contents: Frameshift alignment algorithms 

******************************************************************************/

#include "sls_fsa1_utils.hpp"

using namespace Sls;
using namespace std;

static long int length_max=1000;

void FSA_utils::read_RR(
string RR_file_name_,
double *&RR_,
double *&RR_sum_,
long int *&RR_sum_elements_,
long int &number_of_AA_RR_,
long int number_of_AA_RR_default_)
{
	read_RR(
	RR_file_name_,
	RR_,
	number_of_AA_RR_,
	number_of_AA_RR_default_);

	calculate_RR_sum(
	RR_,
	number_of_AA_RR_,
	RR_sum_,
	RR_sum_elements_);

}

void FSA_utils::check_RR_sum(
double sum_tmp_,
long int number_of_AA_RR_,
string RR_file_name_)
{


	if(number_of_AA_RR_<=0)
	{
		throw error("Error - number of letters in the probabilities file must be greater than 0\n",3);
	};

	double diff_tmp=fabs(sum_tmp_-1.0);
	if(diff_tmp>0)
	{
		double lg_diff=-(log(diff_tmp)-log((double)number_of_AA_RR_))/log(10.0);
		double lg_eps=-log(DBL_EPSILON)/log(10.0)-1;
		if(lg_diff<lg_eps)
		{

			if(sum_tmp_<=0)
			{
				if(RR_file_name_!="")
				{
					throw error("Error: the sum of the probabilities from the file "+RR_file_name_+" is non-positive\n",3);
				}
				else
				{
					throw error("Error: the sum of the probabilities is non-positive\n",3);
				};

			};

			if(RR_file_name_!="")
			{
				static map<string, bool> flag_RR;

				if(!flag_RR[RR_file_name_])
				{
					cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
					cout<<"Warning: the sum of the probabilities from the file "<<RR_file_name_<<" is not equal to 1\n";
					cout<<"The probabilities will be normalized for the computation\n";
					cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n";

					flag_RR[RR_file_name_]=true;
				};
			}
			else
			{
				//no messages if called from the library functions
			};

		};

	};

}

void FSA_utils::calculate_RR_sum(
double *RR_,
long int number_of_AA_RR_,
double *&RR_sum_,
long int *&RR_sum_elements_)
{

	RR_sum_=NULL;
	RR_sum_elements_=NULL;


	try
	{

		long int i;
		if(number_of_AA_RR_<=0)
		{
			throw error("Error - number of letters in the probabilities file must be greater than 0\n",3);
		};
		
		RR_sum_=new double[number_of_AA_RR_];
		assert_mem(RR_sum_);

		RR_sum_elements_=new long int [number_of_AA_RR_];
		assert_mem(RR_sum_elements_);


		for(i=0;i<number_of_AA_RR_;i++)
		{
			if(RR_[i]<0)
			{
				throw error("Error - the frequencies must be non-negative\n",3);
			};

			if(i!=0)
			{
				RR_sum_[i]=RR_sum_[i-1]+RR_[i];
			}
			else
			{
				RR_sum_[i]=RR_[i];
			};
			RR_sum_elements_[i]=i;
		};

		double sum_tmp=RR_sum_[number_of_AA_RR_-1];

		check_RR_sum(
		sum_tmp,
		number_of_AA_RR_,
		"");

		if(sum_tmp>0)
		{
			long int i;
			for(i=0;i<number_of_AA_RR_;i++)
			{
				RR_[i]/=sum_tmp;
				RR_sum_[i]/=sum_tmp;
			};
		};

	}
	catch (...)
	{ 
		delete[]RR_sum_;RR_sum_=NULL;
		delete[]RR_sum_elements_;RR_sum_elements_=NULL;
		throw;
	};


}

void FSA_utils::read_RR(
string RR_file_name_,
double *&RR_,
long int &number_of_AA_RR_,
long int number_of_AA_RR_default_)
{
	ifstream f;

	RR_=NULL;

	try
	{

		if(RR_file_name_=="")
		{
		
			if(number_of_AA_RR_default_==4)
			{
				//default for RR1
				number_of_AA_RR_=4;

				RR_=new double[number_of_AA_RR_];
				assert_mem(RR_);


				RR_[0]=0.25; RR_[1]=0.25; RR_[2]=0.25; RR_[3]=0.25; 
				return;
			};

			if(number_of_AA_RR_default_==25)
			{
				//default for RR2
				number_of_AA_RR_=25;

				RR_=new double[number_of_AA_RR_];
				assert_mem(RR_);

				RR_[0]=0.07805; RR_[1]=0.05129; RR_[2]=0.04487; RR_[3]=0.05364; RR_[4]=0.01925; RR_[5]=0.04264; RR_[6]=0.06295; RR_[7]=0.07377; RR_[8]=0.02199; RR_[9]=0.05142; RR_[10]=0.09019; RR_[11]=0.05744; RR_[12]=0.02243; RR_[13]=0.03856; RR_[14]=0.05203; RR_[15]=0.0712; RR_[16]=0.05841; RR_[17]=0.0133; RR_[18]=0.03216; RR_[19]=0.06441; RR_[20]=0; RR_[21]=0; RR_[22]=0; RR_[23]=0; RR_[24]=0; 
				return;
			};

			throw error("Unexpected parameters in void FSA_utils::read_RR\n",1);

		};

		long int i;
		f.open(RR_file_name_.data(),ios::in);
		if(!f)
		{
			throw error("Error - file "+RR_file_name_+" is not found\n",3);
		};

		f>>number_of_AA_RR_;

		if(number_of_AA_RR_<=0)
		{
			throw error("Error - number of letters in the probabilities file must be greater than 0\n",3);
		};
		
		RR_=new double[number_of_AA_RR_];
		assert_mem(RR_);


		double sum_tmp=0;
		for(i=0;i<number_of_AA_RR_;i++)
		{
			f>>RR_[i];

			if(RR_[i]<0)
			{
				throw error("Error - the frequencies defined in the file "+RR_file_name_+" must be non-negative\n",3);
			};

			sum_tmp+=RR_[i];

		};

		check_RR_sum(
		sum_tmp,
		number_of_AA_RR_,
		RR_file_name_);

		f.close();
	}
	catch (...)
	{ 
		if(f.is_open())
		{
			f.close();
		};
		delete[]RR_;RR_=NULL;
		throw;
	};

}

void FSA_utils::read_smatr(
string smatr_file_name_,
long int **&smatr_,
long int &number_of_AA_smatr_,
long int &smatr_min_)
{
	ifstream f;

	try
	{

		if(smatr_file_name_=="")
		{
			//default for the scoring matrix (BLOSUM80)
			number_of_AA_smatr_=25;
			FSA_utils::get_memory_for_matrix(number_of_AA_smatr_,number_of_AA_smatr_,smatr_);
			smatr_min_=-6;

			smatr_[0][0]=5; smatr_[0][1]=-2; smatr_[0][2]=-2; smatr_[0][3]=-2; smatr_[0][4]=-1; smatr_[0][5]=-1; smatr_[0][6]=-1; smatr_[0][7]=0; smatr_[0][8]=-2; smatr_[0][9]=-2; smatr_[0][10]=-2; smatr_[0][11]=-1; smatr_[0][12]=-1; smatr_[0][13]=-3; smatr_[0][14]=-1; smatr_[0][15]=1; smatr_[0][16]=0; smatr_[0][17]=-3; smatr_[0][18]=-2; smatr_[0][19]=0; smatr_[0][20]=-2; smatr_[0][21]=-2; smatr_[0][22]=-1; smatr_[0][23]=-1; smatr_[0][24]=-6; 
			smatr_[1][0]=-2; smatr_[1][1]=6; smatr_[1][2]=-1; smatr_[1][3]=-2; smatr_[1][4]=-4; smatr_[1][5]=1; smatr_[1][6]=-1; smatr_[1][7]=-3; smatr_[1][8]=0; smatr_[1][9]=-3; smatr_[1][10]=-3; smatr_[1][11]=2; smatr_[1][12]=-2; smatr_[1][13]=-4; smatr_[1][14]=-2; smatr_[1][15]=-1; smatr_[1][16]=-1; smatr_[1][17]=-4; smatr_[1][18]=-3; smatr_[1][19]=-3; smatr_[1][20]=-1; smatr_[1][21]=-3; smatr_[1][22]=0; smatr_[1][23]=-1; smatr_[1][24]=-6; 
			smatr_[2][0]=-2; smatr_[2][1]=-1; smatr_[2][2]=6; smatr_[2][3]=1; smatr_[2][4]=-3; smatr_[2][5]=0; smatr_[2][6]=-1; smatr_[2][7]=-1; smatr_[2][8]=0; smatr_[2][9]=-4; smatr_[2][10]=-4; smatr_[2][11]=0; smatr_[2][12]=-3; smatr_[2][13]=-4; smatr_[2][14]=-3; smatr_[2][15]=0; smatr_[2][16]=0; smatr_[2][17]=-4; smatr_[2][18]=-3; smatr_[2][19]=-4; smatr_[2][20]=5; smatr_[2][21]=-4; smatr_[2][22]=0; smatr_[2][23]=-1; smatr_[2][24]=-6; 
			smatr_[3][0]=-2; smatr_[3][1]=-2; smatr_[3][2]=1; smatr_[3][3]=6; smatr_[3][4]=-4; smatr_[3][5]=-1; smatr_[3][6]=1; smatr_[3][7]=-2; smatr_[3][8]=-2; smatr_[3][9]=-4; smatr_[3][10]=-5; smatr_[3][11]=-1; smatr_[3][12]=-4; smatr_[3][13]=-4; smatr_[3][14]=-2; smatr_[3][15]=-1; smatr_[3][16]=-1; smatr_[3][17]=-6; smatr_[3][18]=-4; smatr_[3][19]=-4; smatr_[3][20]=5; smatr_[3][21]=-5; smatr_[3][22]=1; smatr_[3][23]=-1; smatr_[3][24]=-6; 
			smatr_[4][0]=-1; smatr_[4][1]=-4; smatr_[4][2]=-3; smatr_[4][3]=-4; smatr_[4][4]=9; smatr_[4][5]=-4; smatr_[4][6]=-5; smatr_[4][7]=-4; smatr_[4][8]=-4; smatr_[4][9]=-2; smatr_[4][10]=-2; smatr_[4][11]=-4; smatr_[4][12]=-2; smatr_[4][13]=-3; smatr_[4][14]=-4; smatr_[4][15]=-2; smatr_[4][16]=-1; smatr_[4][17]=-3; smatr_[4][18]=-3; smatr_[4][19]=-1; smatr_[4][20]=-4; smatr_[4][21]=-2; smatr_[4][22]=-4; smatr_[4][23]=-1; smatr_[4][24]=-6; 
			smatr_[5][0]=-1; smatr_[5][1]=1; smatr_[5][2]=0; smatr_[5][3]=-1; smatr_[5][4]=-4; smatr_[5][5]=6; smatr_[5][6]=2; smatr_[5][7]=-2; smatr_[5][8]=1; smatr_[5][9]=-3; smatr_[5][10]=-3; smatr_[5][11]=1; smatr_[5][12]=0; smatr_[5][13]=-4; smatr_[5][14]=-2; smatr_[5][15]=0; smatr_[5][16]=-1; smatr_[5][17]=-3; smatr_[5][18]=-2; smatr_[5][19]=-3; smatr_[5][20]=0; smatr_[5][21]=-3; smatr_[5][22]=4; smatr_[5][23]=-1; smatr_[5][24]=-6; 
			smatr_[6][0]=-1; smatr_[6][1]=-1; smatr_[6][2]=-1; smatr_[6][3]=1; smatr_[6][4]=-5; smatr_[6][5]=2; smatr_[6][6]=6; smatr_[6][7]=-3; smatr_[6][8]=0; smatr_[6][9]=-4; smatr_[6][10]=-4; smatr_[6][11]=1; smatr_[6][12]=-2; smatr_[6][13]=-4; smatr_[6][14]=-2; smatr_[6][15]=0; smatr_[6][16]=-1; smatr_[6][17]=-4; smatr_[6][18]=-3; smatr_[6][19]=-3; smatr_[6][20]=1; smatr_[6][21]=-4; smatr_[6][22]=5; smatr_[6][23]=-1; smatr_[6][24]=-6; 
			smatr_[7][0]=0; smatr_[7][1]=-3; smatr_[7][2]=-1; smatr_[7][3]=-2; smatr_[7][4]=-4; smatr_[7][5]=-2; smatr_[7][6]=-3; smatr_[7][7]=6; smatr_[7][8]=-3; smatr_[7][9]=-5; smatr_[7][10]=-4; smatr_[7][11]=-2; smatr_[7][12]=-4; smatr_[7][13]=-4; smatr_[7][14]=-3; smatr_[7][15]=-1; smatr_[7][16]=-2; smatr_[7][17]=-4; smatr_[7][18]=-4; smatr_[7][19]=-4; smatr_[7][20]=-1; smatr_[7][21]=-5; smatr_[7][22]=-3; smatr_[7][23]=-1; smatr_[7][24]=-6; 
			smatr_[8][0]=-2; smatr_[8][1]=0; smatr_[8][2]=0; smatr_[8][3]=-2; smatr_[8][4]=-4; smatr_[8][5]=1; smatr_[8][6]=0; smatr_[8][7]=-3; smatr_[8][8]=8; smatr_[8][9]=-4; smatr_[8][10]=-3; smatr_[8][11]=-1; smatr_[8][12]=-2; smatr_[8][13]=-2; smatr_[8][14]=-3; smatr_[8][15]=-1; smatr_[8][16]=-2; smatr_[8][17]=-3; smatr_[8][18]=2; smatr_[8][19]=-4; smatr_[8][20]=-1; smatr_[8][21]=-4; smatr_[8][22]=0; smatr_[8][23]=-1; smatr_[8][24]=-6; 
			smatr_[9][0]=-2; smatr_[9][1]=-3; smatr_[9][2]=-4; smatr_[9][3]=-4; smatr_[9][4]=-2; smatr_[9][5]=-3; smatr_[9][6]=-4; smatr_[9][7]=-5; smatr_[9][8]=-4; smatr_[9][9]=5; smatr_[9][10]=1; smatr_[9][11]=-3; smatr_[9][12]=1; smatr_[9][13]=-1; smatr_[9][14]=-4; smatr_[9][15]=-3; smatr_[9][16]=-1; smatr_[9][17]=-3; smatr_[9][18]=-2; smatr_[9][19]=3; smatr_[9][20]=-4; smatr_[9][21]=3; smatr_[9][22]=-4; smatr_[9][23]=-1; smatr_[9][24]=-6; 
			smatr_[10][0]=-2; smatr_[10][1]=-3; smatr_[10][2]=-4; smatr_[10][3]=-5; smatr_[10][4]=-2; smatr_[10][5]=-3; smatr_[10][6]=-4; smatr_[10][7]=-4; smatr_[10][8]=-3; smatr_[10][9]=1; smatr_[10][10]=4; smatr_[10][11]=-3; smatr_[10][12]=2; smatr_[10][13]=0; smatr_[10][14]=-3; smatr_[10][15]=-3; smatr_[10][16]=-2; smatr_[10][17]=-2; smatr_[10][18]=-2; smatr_[10][19]=1; smatr_[10][20]=-4; smatr_[10][21]=3; smatr_[10][22]=-3; smatr_[10][23]=-1; smatr_[10][24]=-6; 
			smatr_[11][0]=-1; smatr_[11][1]=2; smatr_[11][2]=0; smatr_[11][3]=-1; smatr_[11][4]=-4; smatr_[11][5]=1; smatr_[11][6]=1; smatr_[11][7]=-2; smatr_[11][8]=-1; smatr_[11][9]=-3; smatr_[11][10]=-3; smatr_[11][11]=5; smatr_[11][12]=-2; smatr_[11][13]=-4; smatr_[11][14]=-1; smatr_[11][15]=-1; smatr_[11][16]=-1; smatr_[11][17]=-4; smatr_[11][18]=-3; smatr_[11][19]=-3; smatr_[11][20]=-1; smatr_[11][21]=-3; smatr_[11][22]=1; smatr_[11][23]=-1; smatr_[11][24]=-6; 
			smatr_[12][0]=-1; smatr_[12][1]=-2; smatr_[12][2]=-3; smatr_[12][3]=-4; smatr_[12][4]=-2; smatr_[12][5]=0; smatr_[12][6]=-2; smatr_[12][7]=-4; smatr_[12][8]=-2; smatr_[12][9]=1; smatr_[12][10]=2; smatr_[12][11]=-2; smatr_[12][12]=6; smatr_[12][13]=0; smatr_[12][14]=-3; smatr_[12][15]=-2; smatr_[12][16]=-1; smatr_[12][17]=-2; smatr_[12][18]=-2; smatr_[12][19]=1; smatr_[12][20]=-3; smatr_[12][21]=2; smatr_[12][22]=-1; smatr_[12][23]=-1; smatr_[12][24]=-6; 
			smatr_[13][0]=-3; smatr_[13][1]=-4; smatr_[13][2]=-4; smatr_[13][3]=-4; smatr_[13][4]=-3; smatr_[13][5]=-4; smatr_[13][6]=-4; smatr_[13][7]=-4; smatr_[13][8]=-2; smatr_[13][9]=-1; smatr_[13][10]=0; smatr_[13][11]=-4; smatr_[13][12]=0; smatr_[13][13]=6; smatr_[13][14]=-4; smatr_[13][15]=-3; smatr_[13][16]=-2; smatr_[13][17]=0; smatr_[13][18]=3; smatr_[13][19]=-1; smatr_[13][20]=-4; smatr_[13][21]=0; smatr_[13][22]=-4; smatr_[13][23]=-1; smatr_[13][24]=-6; 
			smatr_[14][0]=-1; smatr_[14][1]=-2; smatr_[14][2]=-3; smatr_[14][3]=-2; smatr_[14][4]=-4; smatr_[14][5]=-2; smatr_[14][6]=-2; smatr_[14][7]=-3; smatr_[14][8]=-3; smatr_[14][9]=-4; smatr_[14][10]=-3; smatr_[14][11]=-1; smatr_[14][12]=-3; smatr_[14][13]=-4; smatr_[14][14]=8; smatr_[14][15]=-1; smatr_[14][16]=-2; smatr_[14][17]=-5; smatr_[14][18]=-4; smatr_[14][19]=-3; smatr_[14][20]=-2; smatr_[14][21]=-4; smatr_[14][22]=-2; smatr_[14][23]=-1; smatr_[14][24]=-6; 
			smatr_[15][0]=1; smatr_[15][1]=-1; smatr_[15][2]=0; smatr_[15][3]=-1; smatr_[15][4]=-2; smatr_[15][5]=0; smatr_[15][6]=0; smatr_[15][7]=-1; smatr_[15][8]=-1; smatr_[15][9]=-3; smatr_[15][10]=-3; smatr_[15][11]=-1; smatr_[15][12]=-2; smatr_[15][13]=-3; smatr_[15][14]=-1; smatr_[15][15]=5; smatr_[15][16]=1; smatr_[15][17]=-4; smatr_[15][18]=-2; smatr_[15][19]=-2; smatr_[15][20]=0; smatr_[15][21]=-3; smatr_[15][22]=0; smatr_[15][23]=-1; smatr_[15][24]=-6; 
			smatr_[16][0]=0; smatr_[16][1]=-1; smatr_[16][2]=0; smatr_[16][3]=-1; smatr_[16][4]=-1; smatr_[16][5]=-1; smatr_[16][6]=-1; smatr_[16][7]=-2; smatr_[16][8]=-2; smatr_[16][9]=-1; smatr_[16][10]=-2; smatr_[16][11]=-1; smatr_[16][12]=-1; smatr_[16][13]=-2; smatr_[16][14]=-2; smatr_[16][15]=1; smatr_[16][16]=5; smatr_[16][17]=-4; smatr_[16][18]=-2; smatr_[16][19]=0; smatr_[16][20]=-1; smatr_[16][21]=-1; smatr_[16][22]=-1; smatr_[16][23]=-1; smatr_[16][24]=-6; 
			smatr_[17][0]=-3; smatr_[17][1]=-4; smatr_[17][2]=-4; smatr_[17][3]=-6; smatr_[17][4]=-3; smatr_[17][5]=-3; smatr_[17][6]=-4; smatr_[17][7]=-4; smatr_[17][8]=-3; smatr_[17][9]=-3; smatr_[17][10]=-2; smatr_[17][11]=-4; smatr_[17][12]=-2; smatr_[17][13]=0; smatr_[17][14]=-5; smatr_[17][15]=-4; smatr_[17][16]=-4; smatr_[17][17]=11; smatr_[17][18]=2; smatr_[17][19]=-3; smatr_[17][20]=-5; smatr_[17][21]=-3; smatr_[17][22]=-3; smatr_[17][23]=-1; smatr_[17][24]=-6; 
			smatr_[18][0]=-2; smatr_[18][1]=-3; smatr_[18][2]=-3; smatr_[18][3]=-4; smatr_[18][4]=-3; smatr_[18][5]=-2; smatr_[18][6]=-3; smatr_[18][7]=-4; smatr_[18][8]=2; smatr_[18][9]=-2; smatr_[18][10]=-2; smatr_[18][11]=-3; smatr_[18][12]=-2; smatr_[18][13]=3; smatr_[18][14]=-4; smatr_[18][15]=-2; smatr_[18][16]=-2; smatr_[18][17]=2; smatr_[18][18]=7; smatr_[18][19]=-2; smatr_[18][20]=-3; smatr_[18][21]=-2; smatr_[18][22]=-3; smatr_[18][23]=-1; smatr_[18][24]=-6; 
			smatr_[19][0]=0; smatr_[19][1]=-3; smatr_[19][2]=-4; smatr_[19][3]=-4; smatr_[19][4]=-1; smatr_[19][5]=-3; smatr_[19][6]=-3; smatr_[19][7]=-4; smatr_[19][8]=-4; smatr_[19][9]=3; smatr_[19][10]=1; smatr_[19][11]=-3; smatr_[19][12]=1; smatr_[19][13]=-1; smatr_[19][14]=-3; smatr_[19][15]=-2; smatr_[19][16]=0; smatr_[19][17]=-3; smatr_[19][18]=-2; smatr_[19][19]=4; smatr_[19][20]=-4; smatr_[19][21]=2; smatr_[19][22]=-3; smatr_[19][23]=-1; smatr_[19][24]=-6; 
			smatr_[20][0]=-2; smatr_[20][1]=-1; smatr_[20][2]=5; smatr_[20][3]=5; smatr_[20][4]=-4; smatr_[20][5]=0; smatr_[20][6]=1; smatr_[20][7]=-1; smatr_[20][8]=-1; smatr_[20][9]=-4; smatr_[20][10]=-4; smatr_[20][11]=-1; smatr_[20][12]=-3; smatr_[20][13]=-4; smatr_[20][14]=-2; smatr_[20][15]=0; smatr_[20][16]=-1; smatr_[20][17]=-5; smatr_[20][18]=-3; smatr_[20][19]=-4; smatr_[20][20]=5; smatr_[20][21]=-4; smatr_[20][22]=0; smatr_[20][23]=-1; smatr_[20][24]=-6; 
			smatr_[21][0]=-2; smatr_[21][1]=-3; smatr_[21][2]=-4; smatr_[21][3]=-5; smatr_[21][4]=-2; smatr_[21][5]=-3; smatr_[21][6]=-4; smatr_[21][7]=-5; smatr_[21][8]=-4; smatr_[21][9]=3; smatr_[21][10]=3; smatr_[21][11]=-3; smatr_[21][12]=2; smatr_[21][13]=0; smatr_[21][14]=-4; smatr_[21][15]=-3; smatr_[21][16]=-1; smatr_[21][17]=-3; smatr_[21][18]=-2; smatr_[21][19]=2; smatr_[21][20]=-4; smatr_[21][21]=3; smatr_[21][22]=-3; smatr_[21][23]=-1; smatr_[21][24]=-6; 
			smatr_[22][0]=-1; smatr_[22][1]=0; smatr_[22][2]=0; smatr_[22][3]=1; smatr_[22][4]=-4; smatr_[22][5]=4; smatr_[22][6]=5; smatr_[22][7]=-3; smatr_[22][8]=0; smatr_[22][9]=-4; smatr_[22][10]=-3; smatr_[22][11]=1; smatr_[22][12]=-1; smatr_[22][13]=-4; smatr_[22][14]=-2; smatr_[22][15]=0; smatr_[22][16]=-1; smatr_[22][17]=-3; smatr_[22][18]=-3; smatr_[22][19]=-3; smatr_[22][20]=0; smatr_[22][21]=-3; smatr_[22][22]=5; smatr_[22][23]=-1; smatr_[22][24]=-6; 
			smatr_[23][0]=-1; smatr_[23][1]=-1; smatr_[23][2]=-1; smatr_[23][3]=-1; smatr_[23][4]=-1; smatr_[23][5]=-1; smatr_[23][6]=-1; smatr_[23][7]=-1; smatr_[23][8]=-1; smatr_[23][9]=-1; smatr_[23][10]=-1; smatr_[23][11]=-1; smatr_[23][12]=-1; smatr_[23][13]=-1; smatr_[23][14]=-1; smatr_[23][15]=-1; smatr_[23][16]=-1; smatr_[23][17]=-1; smatr_[23][18]=-1; smatr_[23][19]=-1; smatr_[23][20]=-1; smatr_[23][21]=-1; smatr_[23][22]=-1; smatr_[23][23]=-1; smatr_[23][24]=-6; 
			smatr_[24][0]=-6; smatr_[24][1]=-6; smatr_[24][2]=-6; smatr_[24][3]=-6; smatr_[24][4]=-6; smatr_[24][5]=-6; smatr_[24][6]=-6; smatr_[24][7]=-6; smatr_[24][8]=-6; smatr_[24][9]=-6; smatr_[24][10]=-6; smatr_[24][11]=-6; smatr_[24][12]=-6; smatr_[24][13]=-6; smatr_[24][14]=-6; smatr_[24][15]=-6; smatr_[24][16]=-6; smatr_[24][17]=-6; smatr_[24][18]=-6; smatr_[24][19]=-6; smatr_[24][20]=-6; smatr_[24][21]=-6; smatr_[24][22]=-6; smatr_[24][23]=-6; smatr_[24][24]=1; 

			return;
		};


		long int i,j;
		f.open(smatr_file_name_.data(),ios::in);
		if(!f)
		{
			throw error("Error - file "+smatr_file_name_+" is not found\n",3);
		};

		f>>number_of_AA_smatr_;

		if(number_of_AA_smatr_<=0)
		{
			throw error("Error - number of letters in the scoring matrix file must be greater than 0\n",3);
		};

		get_memory_for_matrix(number_of_AA_smatr_,number_of_AA_smatr_,smatr_);


		for(i=0;i<number_of_AA_smatr_;i++)
		{
			for(j=0;j<number_of_AA_smatr_;j++)
			{
				if(f.eof())
				{
					throw error("Error - file "+smatr_file_name_+" is not correct (please check dimensions of the scoring matrix)\n",3);
				};

				f>>smatr_[i][j];
			};
		};

		f.close();


		smatr_min(
		smatr_,
		number_of_AA_smatr_,
		smatr_min_);

		
	}
	catch (...)
	{ 
		if(f.is_open())
		{
			f.close();
		};
		throw;
	};

}

void FSA_utils::smatr_min(
long int **smatr_,
long int number_of_AA_smatr_,
long int &smatr_min_)
{
	long int small_number=numeric_limits<long int>::min()/10000;
	long int i,j;

	for(i=0;i<number_of_AA_smatr_;i++)
	{
		for(j=0;j<number_of_AA_smatr_;j++)
		{
			if(smatr_[i][j]<small_number)
			{
				smatr_[i][j]=small_number;
			};

			if(i==0&&j==0)
			{
				smatr_min_=smatr_[0][0];
			}
			else
			{
				if(smatr_min_>smatr_[i][j])
				{
					smatr_min_=smatr_[i][j];
				};
			};
		};
	};

}

void FSA_utils::remove_zero_probabilities(
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
long int *&codon_AA_)//<codon code,AA number>
{

	vector<bool> codon_flag_RR2;
	
	{
		long int number_of_codons_old=FSA_utils::power_long(alphabet_letters_number1_,codon_length_);
		codon_flag_RR2.resize(number_of_codons_old,false);
		long int i;
		for(i=0;i<number_of_codons_old;i++)
		{
			codon_flag_RR2[codon_AA_[i]]=true;
		};
	};


	map<long int,bool> zero_flag_RR1;
	map<long int,bool> zero_flag_RR2;

	map<long int,long int> old1_numbering_to_the_new1;
	map<long int,long int> new1_numbering_to_the_old1;
	long int count1=-1;
	long int i;
	for(i=0;i<alphabet_letters_number1_;i++)
	{
		if(RR1_[i]==0)
		{
			zero_flag_RR1[i]=true;
			old1_numbering_to_the_new1[i]=-1;
		}
		else
		{
			zero_flag_RR1[i]=false;
			count1++;
			old1_numbering_to_the_new1[i]=count1;
			new1_numbering_to_the_old1[count1]=i;
		};

	};

	map<long int,long int> old2_numbering_to_the_new2;
	map<long int,long int> new2_numbering_to_the_old2;
	long int count2=-1;
	for(i=0;i<alphabet_letters_number2_;i++)
	{
		if(RR2_[i]==0&&!codon_flag_RR2[i])
		{
			zero_flag_RR2[i]=true;
			old2_numbering_to_the_new2[i]=-1;
		}
		else
		{
			zero_flag_RR2[i]=false;
			count2++;
			old2_numbering_to_the_new2[i]=count2;
			new2_numbering_to_the_old2[count2]=i;
		};
	};
//-------------------------------------
	long int alphabet_letters_number1_new=count1+1;
	double *RR1_new=new double[alphabet_letters_number1_new];
	double *RR1_sum_new=new double[alphabet_letters_number1_new];
	long int *RR1_sum_elements_new=new long int[alphabet_letters_number1_new];
	
	long int alphabet_letters_number2_new=count2+1;
	double *RR2_new=new double[alphabet_letters_number2_new];
	double *RR2_sum_new=new double[alphabet_letters_number2_new];
	long int *RR2_sum_elements_new=new long int[alphabet_letters_number2_new];


	long int c1;
	for(c1=0;c1<alphabet_letters_number1_new;c1++)
	{
		RR1_new[c1]=RR1_[new1_numbering_to_the_old1[c1]];
		RR1_sum_new[c1]=RR1_sum_[new1_numbering_to_the_old1[c1]];
		RR1_sum_elements_new[c1]=c1;
	};

	long int c2;
	for(c2=0;c2<alphabet_letters_number2_new;c2++)
	{
		RR2_new[c2]=RR2_[new2_numbering_to_the_old2[c2]];
		RR2_sum_new[c2]=RR2_sum_[new2_numbering_to_the_old2[c2]];
		RR2_sum_elements_new[c2]=c2;
	};



	//reallocation
	
	delete[]RR1_;
	RR1_=RR1_new;
	delete[]RR1_sum_;
	RR1_sum_=RR1_sum_new;
	delete[]RR1_sum_elements_;
	RR1_sum_elements_=RR1_sum_elements_new;

	
	delete[]RR2_;
	RR2_=RR2_new;
	delete[]RR2_sum_;
	RR2_sum_=RR2_sum_new;
	delete[]RR2_sum_elements_;
	RR2_sum_elements_=RR2_sum_elements_new;


//----------------------------

	long int **smatr_new=NULL;

	FSA_utils::get_memory_for_matrix(alphabet_letters_number2_new,alphabet_letters_number2_new,smatr_new);

	long int smatr_min_new=smatr_[new2_numbering_to_the_old2[0]][new2_numbering_to_the_old2[0]];
	for(c1=0;c1<alphabet_letters_number2_new;c1++)
	{
		for(c2=0;c2<alphabet_letters_number2_new;c2++)
		{
			smatr_new[c1][c2]=
				smatr_[new2_numbering_to_the_old2[c1]][new2_numbering_to_the_old2[c2]];

			smatr_min_new=FSA_utils::Tmin(smatr_min_new,
				smatr_new[c1][c2]);

		};
	};

	FSA_utils::delete_memory_for_matrix(number_of_AA_smatr_,smatr_);
	smatr_=smatr_new;

	char *alphabet1_new=new char[alphabet_letters_number1_new];
	char *alphabet2_new=new char[alphabet_letters_number2_new];

	for(c1=0;c1<alphabet_letters_number1_new;c1++)
	{
		alphabet1_new[c1]=alphabet1_[new1_numbering_to_the_old1[c1]];
	};

	for(c2=0;c2<alphabet_letters_number2_new;c2++)
	{
		alphabet2_new[c2]=alphabet2_[new2_numbering_to_the_old2[c2]];
	};

	delete[]alphabet1_;
	alphabet1_=alphabet1_new;
	delete[]alphabet2_;
	alphabet2_=alphabet2_new;

	//-----------------------------

	
	{
		long int k;
		long int*alphabet1_to_long_new=new long int [length_max];
		FSA_utils::assert_mem(alphabet1_to_long_new);
		long int*alphabet2_to_long_new=new long int [length_max];
		FSA_utils::assert_mem(alphabet2_to_long_new);

		for(k=0;k<length_max;k++)
		{
			alphabet1_to_long_new[k]=-1;
			alphabet2_to_long_new[k]=-1;
		};

		for(k=0;k<alphabet_letters_number1_new;k++)
		{
			alphabet1_to_long_new[(size_t)alphabet1_[k]]=k;
		};
		for(k=0;k<alphabet_letters_number2_new;k++)
		{
			alphabet2_to_long_new[(size_t)alphabet2_[k]]=k;
		};

		delete[]alphabet1_to_long_;
		alphabet1_to_long_=alphabet1_to_long_new;
		delete[]alphabet2_to_long_;
		alphabet2_to_long_=alphabet2_to_long_new;

	};



	//-----------------------------

	long int number_of_codons_new=FSA_utils::power_long(alphabet_letters_number1_new,codon_length_);
	long int *codon_AA_new = new long int[number_of_codons_new];//<codon code,AA number>

	long int *codon=new long int [codon_length_];
	FSA_utils::assert_mem(codon);

	for(c1=0;c1<number_of_codons_new;c1++)
	{
		FSA_utils::convert_code_into_codon(
		c1,//the input code
		codon_length_,//codon length 
		alphabet_letters_number1_new,//number of letters for the sequence 1
		codon);//must be allocated

		long int k;
		for(k=0;k<codon_length_;k++)
		{
			codon[k]=new1_numbering_to_the_old1[codon[k]];
		};

		long int code_old=FSA_utils::convert_codon_into_code(
		codon_length_,//codon length 
		alphabet_letters_number1_,//number of letters for the sequence 1
		codon);//input codon


		codon_AA_new[c1]=old2_numbering_to_the_new2[codon_AA_[code_old]];
	};

	delete[]codon;

	delete[]codon_AA_;
	codon_AA_=codon_AA_new;


	//-----------------------------
	alphabet_letters_number1_=alphabet_letters_number1_new;
	alphabet_letters_number2_=alphabet_letters_number2_new;
	number_of_AA_smatr_=alphabet_letters_number2_new;
	smatr_min_=smatr_min_new;
	number_of_letters1_=alphabet_letters_number1_new;
	number_of_letters2_=alphabet_letters_number2_new;

}


void FSA_utils::reverse_codons(
long int *codon_AA_,//<codon code,AA number>; original codons
long int alphabet_letters_number1_,//number of letters for the sequence #1
long int codon_length_,//codon length 
long int *&codon_AA_reversed_)//<codon code,AA number>; reversed codons
{

	long int number_of_codons=FSA_utils::power_long(alphabet_letters_number1_,codon_length_);

	codon_AA_reversed_ = new long int[number_of_codons];//<codon code,AA number>
	FSA_utils::assert_mem(codon_AA_reversed_);

	long int *codon=new long int [codon_length_];
	FSA_utils::assert_mem(codon);
	long int *codon2=new long int [codon_length_];
	FSA_utils::assert_mem(codon2);

	long int c1;
	for(c1=0;c1<number_of_codons;c1++)
	{
		FSA_utils::convert_code_into_codon(
		c1,//the input code
		codon_length_,//codon length 
		alphabet_letters_number1_,//number of letters for the sequence 1
		codon);//must be allocated

		long int k;
		for(k=0;k<codon_length_;k++)
		{
			codon2[k]=codon[codon_length_-1-k];
		};

		long int code_reversed=FSA_utils::convert_codon_into_code(
		codon_length_,//codon length 
		alphabet_letters_number1_,//number of letters for the sequence 1
		codon2);//input codon


		codon_AA_reversed_[code_reversed]=codon_AA_[c1];


	};

	delete[]codon;
	delete[]codon2;

}


long int FSA_utils::convert_codon_into_code(
long int codon_length_,//codon length 
long int number_of_letters1_,//number of letters for the sequence 1
long int *codon_)//input codon
{
	long int codon_code=codon_[0];
	long int i;
	for(i=1;i<codon_length_;i++)
	{
		codon_code=codon_code*number_of_letters1_+codon_[i];
	};

	return codon_code;

}


long int FSA_utils::convert_codon_into_AA(
long int codon_length_,//codon length 
long int *codon_AA_,//<codon code,AA number>
long int number_of_letters1_,//number of letters for the sequence 1
long int *codon_)
{
	long int codon_code=codon_[0];
	long int i;
	for(i=1;i<codon_length_;i++)
	{
		codon_code=codon_code*number_of_letters1_+codon_[i];
	};

	return codon_AA_[codon_code];

}

void FSA_utils::convert_code_into_codon(
long int code_,//the input code
long int codon_length_,//codon length 
long int number_of_letters1_,//number of letters for the sequence 1
long int *codon_)//must be allocated
{
	long int i;
	for(i=codon_length_-1;i>=0;i--)
	{
		codon_[i]=code_%number_of_letters1_;
		code_-=codon_[i];
		code_/=number_of_letters1_;
	};

}

void FSA_utils::read_codon_AA_file(
string file_name_,
long int &number_of_letters1_,//number of letters for the sequence 1
long int &number_of_letters2_,//number of letters for the sequence 2
char *&alphabet1_,//alphabet letters for the sequence #1
char *&alphabet2_,//alphabet letters for the sequence #2

long int *&alphabet1_to_long_,//d_alphabet1_to_long[c] returns order number of the letter c (for the sequence #1)
long int *&alphabet2_to_long_,//d_alphabet2_to_long[c] returns order number of the letter c (for the sequence #2)

long int &codon_length_,//codon length 
long int *&codon_AA_,//<codon code,AA number>
bool reverse_codons_flag_)//if true, then the codons are reversed
{
	ifstream f;

	alphabet1_=NULL;
	alphabet2_=NULL;
	alphabet1_to_long_=NULL;
	alphabet2_to_long_=NULL;
	codon_AA_=NULL;
	try
	{


		if(file_name_=="")
		{
			number_of_letters1_=4;
			alphabet1_=new char[4];
			FSA_utils::assert_mem(alphabet1_);
			alphabet1_[0]='A'; alphabet1_[1]='C'; alphabet1_[2]='G'; alphabet1_[3]='T'; 
			number_of_letters2_=25;
			alphabet2_=new char[25];
			FSA_utils::assert_mem(alphabet2_);
			alphabet2_[0]='A'; alphabet2_[1]='R'; alphabet2_[2]='N'; alphabet2_[3]='D'; alphabet2_[4]='C'; alphabet2_[5]='Q'; alphabet2_[6]='E'; alphabet2_[7]='G'; alphabet2_[8]='H'; alphabet2_[9]='I'; alphabet2_[10]='L'; alphabet2_[11]='K'; alphabet2_[12]='M'; alphabet2_[13]='F'; alphabet2_[14]='P'; alphabet2_[15]='S'; alphabet2_[16]='T'; alphabet2_[17]='W'; alphabet2_[18]='Y'; alphabet2_[19]='V'; alphabet2_[20]='B'; alphabet2_[21]='J'; alphabet2_[22]='Z'; alphabet2_[23]='X'; alphabet2_[24]='*'; 
			alphabet1_to_long_=new long int [length_max];
			FSA_utils::assert_mem(alphabet1_to_long_);
			alphabet2_to_long_=new long int [length_max];
			FSA_utils::assert_mem(alphabet2_to_long_);

			long int k;
			for(k=0;k<length_max;k++)
			{
				alphabet1_to_long_[k]=-1;
				alphabet2_to_long_[k]=-1;
			};

			for(k=0;k<number_of_letters1_;k++)
			{
				alphabet1_to_long_[(size_t)alphabet1_[k]]=k;
			};
			for(k=0;k<number_of_letters2_;k++)
			{
				alphabet2_to_long_[(size_t)alphabet2_[k]]=k;
			};


			codon_length_=3;
			codon_AA_=new long int [64];
			FSA_utils::assert_mem(codon_AA_);
			codon_AA_[63]=13; codon_AA_[61]=13; codon_AA_[60]=10; codon_AA_[62]=10; codon_AA_[31]=10; codon_AA_[29]=10; codon_AA_[28]=10; codon_AA_[30]=10; codon_AA_[15]=9; codon_AA_[13]=9; codon_AA_[12]=9; codon_AA_[14]=12; codon_AA_[47]=19; codon_AA_[45]=19; codon_AA_[44]=19; codon_AA_[46]=19; codon_AA_[55]=15; codon_AA_[53]=15; codon_AA_[52]=15; codon_AA_[54]=15; codon_AA_[23]=14; codon_AA_[21]=14; codon_AA_[20]=14; codon_AA_[22]=14; codon_AA_[7]=16; codon_AA_[5]=16; codon_AA_[4]=16; codon_AA_[6]=16; codon_AA_[39]=0; codon_AA_[37]=0; codon_AA_[36]=0; codon_AA_[38]=0; codon_AA_[51]=18; codon_AA_[49]=18; codon_AA_[48]=24; codon_AA_[50]=24; codon_AA_[19]=8; codon_AA_[17]=8; codon_AA_[16]=5; codon_AA_[18]=5; codon_AA_[3]=2; codon_AA_[1]=2; codon_AA_[0]=11; codon_AA_[2]=11; codon_AA_[35]=3; codon_AA_[33]=3; codon_AA_[32]=6; codon_AA_[34]=6; codon_AA_[59]=4; codon_AA_[57]=4; codon_AA_[56]=24; codon_AA_[58]=17; codon_AA_[27]=1; codon_AA_[25]=1; codon_AA_[24]=1; codon_AA_[26]=1; codon_AA_[11]=15; codon_AA_[9]=15; codon_AA_[8]=1; codon_AA_[10]=1; codon_AA_[43]=7; codon_AA_[41]=7; codon_AA_[40]=7; codon_AA_[42]=7; 

			return;
		};

		f.open(file_name_.data(),ios::in);
		if(!f)
		{
			throw error("Error - the file "+file_name_+" is not found\n",3);
		};

		//reading alphabet #1
		f>>number_of_letters1_;

		if(number_of_letters1_<=0)
		{
			throw error("Error - the file "+file_name_+" is not correct: the number of letters must be positive\n",3);
		};

		alphabet1_=new char[number_of_letters1_];
		FSA_utils::assert_mem(alphabet1_);

		string line_tmp="";

		f>>line_tmp;

		if((long int)line_tmp.length()<number_of_letters1_)
		{
			throw error("Error - the file "+file_name_+" is not correct\n",3);
		};

		long int k;
		for(k=0;k<(long int)line_tmp.length();k++)
		{
			if(!isalpha(line_tmp[k]))
			{
				throw error("Error - the file "+file_name_+" is not correct\n",3);
			};
			alphabet1_[k]=line_tmp[k];
		};


		//reading alphabet #2
		f>>number_of_letters2_;

		if(number_of_letters2_<=0)
		{
			throw error("Error - the file "+file_name_+" is not correct\n",3);
		};

		alphabet2_=new char[number_of_letters2_];
		FSA_utils::assert_mem(alphabet2_);

		line_tmp="";

		f>>line_tmp;

		if((long int)line_tmp.length()<number_of_letters2_)
		{
			throw error("Error - the file "+file_name_+" is not correct\n",3);
		};

		for(k=0;k<(long int)line_tmp.length();k++)
		{
			//if(!isalpha(line_tmp[k]))
			if(!isalpha(line_tmp[k])&&line_tmp[k]!='*')
			{
				throw error("Error - the file "+file_name_+" is not correct\n",3);
			};
			alphabet2_[k]=line_tmp[k];
		};


		alphabet1_to_long_=new long int [length_max];
		FSA_utils::assert_mem(alphabet1_to_long_);
		alphabet2_to_long_=new long int [length_max];
		FSA_utils::assert_mem(alphabet2_to_long_);

		for(k=0;k<length_max;k++)
		{
			alphabet1_to_long_[k]=-1;
			alphabet2_to_long_[k]=-1;
		};

		for(k=0;k<number_of_letters1_;k++)
		{
			alphabet1_to_long_[(size_t)alphabet1_[k]]=k;
		};
		for(k=0;k<number_of_letters2_;k++)
		{
			alphabet2_to_long_[(size_t)alphabet2_[k]]=k;
		};

		long int AA_length;
		//f>>codon_length_>>AA_length;
		codon_length_=3;
		AA_length=1;
		if(AA_length!=1||codon_length_<=0)
		{
			throw error("Error - the file "+file_name_+" is not correct\n",3);
		};
		//read the table
		long int number_of_codons=1;
		for(k=0;k<codon_length_;k++)
		{
			number_of_codons*=number_of_letters1_;
		};

		codon_AA_=new long int [number_of_codons];
		FSA_utils::assert_mem(codon_AA_);

		for(k=0;k<number_of_codons;k++)
		{
			string line_tmp_codon,line_tmp_AA;

			f>>line_tmp_codon;
			if(reverse_codons_flag_)
			{
				line_tmp_codon=string(line_tmp_codon.rbegin(),line_tmp_codon.rend()); 
			};
			f>>line_tmp_AA;

			if((codon_length_!=(long int)line_tmp_codon.length())||(AA_length!=(long int)line_tmp_AA.length()))
			{
				throw error("Error - the file "+file_name_+" is not correct\n",3);
			};

			long int codon_code=alphabet1_to_long_[(size_t)line_tmp_codon[0]];
			if(alphabet1_to_long_[(size_t)line_tmp_codon[0]]<0)
			{
				throw error("Error - the file "+file_name_+" is not correct\n",3);
			};
			long int i;
			for(i=1;i<(long int)line_tmp_codon.length();i++)
			{
				if(alphabet1_to_long_[(size_t)line_tmp_codon[i]]<0)
				{
					throw error("Error - the file "+file_name_+" is not correct\n",3);
				};
				codon_code=codon_code*number_of_letters1_+alphabet1_to_long_[(size_t)line_tmp_codon[i]];
			};

			if(alphabet2_to_long_[(size_t)line_tmp_AA[0]]<0)
			{
				throw error("Error - the file "+file_name_+" is not correct\n",3);
			};


			codon_AA_[codon_code]=alphabet2_to_long_[(size_t)line_tmp_AA[0]];
		};

		


		f.close();
	}
	catch (...)
	{ 
		delete[]alphabet1_;alphabet1_=NULL;
		delete[]alphabet2_;alphabet2_=NULL;
		delete[]alphabet1_to_long_;alphabet1_to_long_=NULL;
		delete[]alphabet2_to_long_;alphabet2_to_long_=NULL;
		delete[]codon_AA_;codon_AA_=NULL;
		if(f.is_open())
		{
			f.close();
		};
		throw;
	};


}

string FSA_utils::long_to_string(//convert interer ot string
long int number_)
{
	string res_="";
	string tmp_string;
	if(number_>0)
	{
		tmp_string="";
	}
	else
	{
		if(number_==0)
		{
			tmp_string="";
		}
		else
		{
			tmp_string="-";
		};
	};
	number_=labs(number_);

	for( ; ; )
	{
		long int reminder=number_%10;
		number_=(number_-reminder)/10;
		res_=digit_to_string(reminder)+res_;
		if (number_==0)
		{
			break;
		};
	};

	return tmp_string+res_;
}

char FSA_utils::digit_to_string(//convert interer ot string
long int digit_)
{
	switch(digit_)
	{
	case 0:return '0';
	case 1:return '1';
	case 2:return '2';
	case 3:return '3';
	case 4:return '4';
	case 5:return '5';
	case 6:return '6';
	case 7:return '7';
	case 8:return '8';
	case 9:return '9';
	default:return '?';
	};
}

bool FSA_utils::the_value_is_double(
string str_,
double &val_)
{
	if(str_=="")
	{
		return false;
	};

	bool res=false;
	errno=0;
	char *p;
	val_=strtod(str_.c_str(),&p);
	if(errno!=0)
	{
		res=false;
	}
	else
	{
		res=(*p==0);
	};
	return res;
}

bool FSA_utils::the_value_is_long(
string str_,
long int &val_)
{

	if(str_.length()==0)
	{
		return false;
	};
	if(!(str_[0]=='+'||str_[0]=='-'||isdigit(str_[0])))
	{
		return false;
	};

	long int start_digit=0;

	if(str_[0]=='+'||str_[0]=='-')
	{
		start_digit=1;
	};


	long int i;
	for(i=start_digit;i<(long int)str_.size();i++)
	{
		if(!isdigit(str_[i]))
		{
			return false;
		};
	};

	if(((long int)str_.size()-start_digit)<=0)
	{
		return false;
	};

	if(((long int)str_.size()-start_digit)>1)
	{
		while(str_[start_digit]=='0')
		{
			string::iterator it=str_.begin()+start_digit;


			str_.erase(it);
			if((long int)str_.size()<=start_digit+1)
			{
				break;
			};
		};
	};

	if(((long int)str_.size()-start_digit>10)||((long int)str_.size()-start_digit)<=0)
	{
		return false;
	};


	if((long int)str_.size()-start_digit==10)
	{
		if(!(str_[start_digit]=='1'||str_[start_digit]=='2'))
		{
			return false;
		};

		if(str_[start_digit]=='2')
		{

			long int val2;
			string str2=str_.substr(start_digit+1,9);
			int flag=sscanf(str2.c_str(),"%ld",&val2);
			if(flag!=1)
			{
				return false;
			};

			bool positive=true;
			if(start_digit>0)
			{
				if(str_[0]=='-')
				{
					positive=false;
				};
			};

			if(positive)
			{
				if(val2>147483647)
				{
					return false;
				};
			}
			else
			{
				if(val2>147483648)
				{
					return false;
				};
			};

		};
	};

	int flag=sscanf(str_.c_str(),"%ld",&val_);
	if(flag!=1)
	{
		return false;
	};

	return true;
}

double FSA_utils::sqrt_plus(
double x_)
{
	if(x_>=0)
	{
		return sqrt(x_);
	}
	else
	{
		return 0;
	};
}

double FSA_utils::error_of_the_sum(//v1_+v2_
double v1_error_,
double v2_error_)
{
	if(v1_error_>=1e100||v2_error_>=1e100)
	{
		return 1e100;
	};

	return sqrt(v1_error_*v1_error_+v2_error_*v2_error_);
}

double FSA_utils::error_of_the_product(//v1_*v2_
double v1_,
double v1_error_,
double v2_,
double v2_error_)
{
	if(v1_error_>=1e100||v2_error_>=1e100)
	{
		return 1e100;
	};

	double a1=(v1_+v1_error_)*(v2_+v2_error_);
	double a2=(v1_-v1_error_)*(v2_+v2_error_);
	double a3=(v1_+v1_error_)*(v2_-v2_error_);
	double a4=(v1_-v1_error_)*(v2_-v2_error_);

	double a=v1_*v2_;

	return FSA_utils::Tmax(fabs(a1-a),fabs(a2-a),fabs(a3-a),fabs(a4-a));

}

double FSA_utils::error_of_the_lg(//lg(v1_)
double v1_,
double v1_error_)
{
	if(v1_error_>=1e100||v1_<=0)
	{
		return 1e100;
	};

	return FSA_utils::Tmin(fabs(log(v1_)/log(10.0)),v1_error_/v1_/log(10.0));
}

double FSA_utils::error_of_the_sqrt(//sqrt(v1_)
double v1_,
double v1_error_)
{
	if(v1_error_>=1e100||v1_<0)
	{
		return 1e100;
	};

	double s=sqrt(v1_);
	double s1=sqrt(FSA_utils::Tmax(0.0,v1_-v1_error_));
	double s2=sqrt(FSA_utils::Tmax(0.0,v1_+v1_error_));

	return FSA_utils::Tmax(fabs(s-s1),fabs(s-s2));
}

double FSA_utils::error_of_the_ratio(//v1_/v2_
double v1_,
double v1_error_,
double v2_,
double v2_error_)
{
	if(v1_error_>=1e100||v2_error_>=1e100)
	{
		return 1e100;
	};


	if(v2_==0)
	{
		return 1e100;
	};

	if(v1_==0&&v1_error_==0)
	{
		return 0.0;
	};

	double a=v1_/v2_;


	if(((v2_+v2_error_)*v2_<=0))
	{
		double a3=(v1_+v1_error_)/(v2_-v2_error_);
		double a4=(v1_-v1_error_)/(v2_-v2_error_);
		return FSA_utils::Tmax(fabs(a-a3),fabs(a-a4));
	};

	if(((v2_-v2_error_)*v2_<=0))
	{
		double a1=(v1_+v1_error_)/(v2_+v2_error_);
		double a2=(v1_-v1_error_)/(v2_+v2_error_);
		return FSA_utils::Tmax(fabs(a-a1),fabs(a-a2));
	};


	double a1=(v1_+v1_error_)/(v2_+v2_error_);
	double a2=(v1_-v1_error_)/(v2_+v2_error_);
	double a3=(v1_+v1_error_)/(v2_-v2_error_);
	double a4=(v1_-v1_error_)/(v2_-v2_error_);

	return FSA_utils::Tmax(fabs(a-a1),fabs(a-a2),fabs(a-a3),fabs(a-a4));
}

double FSA_utils::error_of_the_sum_with_coeff(//c1_*v1_+c2_*v2_
double c1_,
double v1_error_,
double c2_,
double v2_error_)
{
	if(v1_error_>=1e100||v2_error_>=1e100)
	{
		return 1e100;
	};

	return sqrt(c1_*c1_*v1_error_*v1_error_+c2_*c2_*v2_error_*v2_error_);
}


void FSA_utils::read_alphabet(
string alphabet_file_name_,
long int &number_of_AA_alphabet_,
char* &alphabet_)
{
	ifstream a(alphabet_file_name_.data());
	if(!a)
	{
		throw error("Error - the file "+alphabet_file_name_+" is not found\n",1);
	};

	a>>number_of_AA_alphabet_;
	if(number_of_AA_alphabet_<=0)
	{
		throw error("Error - the file "+alphabet_file_name_+" is wrong\n",1);
	};
	alphabet_=new char[number_of_AA_alphabet_];
	FSA_utils::assert_mem(alphabet_);

	long int i;
	for(i=0;i<number_of_AA_alphabet_;i++)
	{
		if(a.eof())
		{
			throw error("Error - the file "+alphabet_file_name_+" is wrong\n",1);
		};
		a>>alphabet_[i];
	};

	a.close();

}

void FSA_utils::calculate_composition_frequencies(
double number_of_AA_alphabet_,
char* alphabet_,
string fasta_file_name_,
double *comp_frequencies_)
{

	long int max_c=1000;
	long int *letter_to_int=new long int[max_c];
	FSA_utils::assert_mem(letter_to_int);

	long int i;
	for(i=0;i<max_c;i++)
	{
		letter_to_int[i]=-1;
	};

	for(i=0;i<number_of_AA_alphabet_;i++)
	{
		letter_to_int[(size_t)alphabet_[i]]=i;
		comp_frequencies_[i]=0;
	};

	ifstream f(fasta_file_name_.data());
	if(!f)
	{
		throw error("Error - the file "+fasta_file_name_+" is not found\n",1);
	};




	
	string header;

	if(f.eof())
	{
		throw error("Error - the file "+fasta_file_name_+" is incorrect\n",1);
	};

	getline(f,header);

	if(header.size()==0)
	{
		throw error("Error - the file "+fasta_file_name_+" is incorrect\n",1);
	};


	if(header[0]!='>')
	{
		throw error("Error - the file "+fasta_file_name_+" is incorrect\n",1);
	};
	

	string st;
	if(f.eof())
	{
		throw error("Error - the file "+fasta_file_name_+" is incorrect\n",1);
	};

	getline(f,st);
	if(st.size()==0)
	{
		throw error("Error - the file "+fasta_file_name_+" is incorrect\n",1);
	};
	for( ; ; )
	{
		long int j;
		for(j=0;j<(long int)st.size();j++)
		{
			long int letter_ind=letter_to_int[(size_t)st[j]];
			if(letter_ind<0||letter_ind>=number_of_AA_alphabet_)
			{
				continue;
			};
			comp_frequencies_[letter_ind]++;
		};

		if(f.eof())
		{
			break;
		};

		getline(f,st);
		if(st.size()==0)
		{
			break;
		};
	};

	f.close();	

	double sum=0;
	for(i=0;i<number_of_AA_alphabet_;i++)
	{
		sum+=comp_frequencies_[i];
	};
	cout<<"Total number of allowed letters = "<<sum<<endl;

	if(sum>0)
	{
		for(i=0;i<number_of_AA_alphabet_;i++)
		{
			comp_frequencies_[i]/=sum;
		};
	}
	else
	{
		throw error("Error - the file "+fasta_file_name_+" does not have letters from the allowed alphabet\n",1);
	};

	delete[]letter_to_int;
}

void FSA_utils::reverse_sequence(//reverse the letters of the sequence
long int *seq_,
long int seq_length_)
{
	long int i;
	for(i=0;i<=(long int)floor((double)seq_length_/2.0)-1;i++)
	{
		long int reverse_ind=seq_length_-i-1;
		long int tmp=seq_[i];
		seq_[i]=seq_[reverse_ind];
		seq_[reverse_ind]=tmp;
	};
}

void FSA_utils::extract_AA_frequencies_for_DNA_sequence(
const long int *codon_AA_,//<codon code,AA number>
long int &codon_length_,//codon length 
long int number_of_letters1_,//number of letters for the sequence 1
long int number_of_letters2_,//number of letters for the sequence 2
const double *RR1_,//nucleotide probabilities
double *&RR1_AA_)//the resulted frequencies
{
	long int *codon_tmp=new long int [codon_length_];
	FSA_utils::assert_mem(codon_tmp);

	RR1_AA_=new double[number_of_letters2_];
	FSA_utils::assert_mem(RR1_AA_);

	long int k;
	for(k=0;k<number_of_letters2_;k++)
	{
		RR1_AA_[k]=0.0;
	};

	long int number_of_codons=FSA_utils::power_long(number_of_letters1_,codon_length_);
	for(k=0;k<number_of_codons;k++)
	{
		convert_code_into_codon(
		k,//the input code
		codon_length_,//codon length 
		number_of_letters1_,//number of letters for the sequence 1
		codon_tmp);//must be allocated

		long int AA1=codon_AA_[k];
		if(AA1<0||AA1>=number_of_letters2_)
		{
			throw error("Unexpected errro in FSA_utils::extract_AA_frequencies_for_DNA_sequence\n",1);
		};

		double prob_tmp=1.0;
		long int i;
		for(i=0;i<codon_length_;i++)
		{
			prob_tmp*=RR1_[codon_tmp[i]];
		};

		RR1_AA_[AA1]+=prob_tmp;

	};

	delete[]codon_tmp;
}


long int FSA_utils::power_long(//returns a_^b_
long int a_,
long int b_)
{
	if(b_<0)
	{
		throw error("Error - unexpected parameter b_<0 in the function FSA_utils::power_long\n",1);
	};

	long int res=1;
	long int i;
	for(i=1;i<=b_;i++)
	{
		res*=a_;
	};

	return res;

}

void FSA_utils::convert_distr_into_sum(//convert distr_[0], distr_[1], distr_[2],... into distr_[0], distr_[0]+distr_[1], distr_[0]+distr_[1]+distr_[2],...
long int dim_,
double *distr_)
{
	long int i;
	for(i=1;i<dim_;i++)
	{
		distr_[i]=distr_[i]+distr_[i-1];
	};
}

bool FSA_utils::Gauss(
std::vector<std::vector<double> > A_,//matrix n*(n+1)
std::vector<double> &x_,//solution
double inside_eps_,
std::vector<std::vector<double> > *inv_A_)
{
	

	long int i,j,jj;
	long int matr_size=(long int)A_.size();
	if(matr_size==0)
	{
		//throw error("Error in FSA_utils::Gauss - sizes of matrix are wrong\n",1);
		return false;
	};

	std::vector<std::vector<double> > E;

	if(inv_A_)
	{
		vector<double> zero(matr_size,0);
		(*inv_A_).resize(matr_size, zero);
		
		E.resize(matr_size, zero);
		long int i;
		for(i=0;i<matr_size;i++)
		{
			E[i][i]=1.0;
		};
	};


	for(i=0;i<matr_size;i++)
	{
		if((long int)A_[i].size()!=matr_size+1)
		{
			throw error("Error in FSA_utils::Gauss - sizes of matrix are wrong\n",1);
		};
	};
	x_.clear();
	x_.resize(matr_size);
	//forward trace
	for(j=0;j<matr_size;j++)
	{
		long int absmax=j;
		for(i=j+1;i<matr_size;i++)
		{
			if(fabs(A_[absmax][j])<fabs(A_[i][j]))
			{
				absmax=i;
			};
		};

		if(j!=absmax)
		{
			for(jj=j;jj<matr_size+1;jj++)
			{
				double tmp=A_[absmax][jj];
				A_[absmax][jj]=A_[j][jj];
				A_[j][jj]=tmp;
			};

			if(inv_A_)
			{
				for(jj=0;jj<matr_size;jj++)
				{
					double tmp=E[absmax][jj];
					E[absmax][jj]=E[j][jj];
					E[j][jj]=tmp;
				};
			};
		};

		if(fabs(A_[j][j])<=inside_eps_)
		{
			throw error("Error in FSA_utils::Gauss - matrix is singular\n",1);
		};

		for(i=j+1;i<matr_size;i++)
		{
			double tmp=A_[i][j]/A_[j][j];
			for(jj=j+1;jj<matr_size+1;jj++)
			{
				A_[i][jj]=A_[i][jj]-tmp*A_[j][jj];
			};

			if(inv_A_)
			{
				for(jj=0;jj<matr_size;jj++)
				{
					E[i][jj]=E[i][jj]-tmp*E[j][jj];
				};
			};
		};
	};

	//reverse trace
	x_[matr_size-1]=A_[matr_size-1][matr_size]/A_[matr_size-1][matr_size-1];
	for(i=matr_size-2;i>=0;i--)
	{
		x_[i]=A_[i][matr_size];
		for(j=i+1;j<matr_size;j++)
		{
			x_[i]-=A_[i][j]*x_[j];
		};
		x_[i]/=A_[i][i];
	};

	if(inv_A_)
	{
		long int k;
		for(k=0;k<matr_size;k++)
		{
			(*inv_A_)[matr_size-1][k]=E[matr_size-1][k]/A_[matr_size-1][matr_size-1];
			long int i;
			for(i=matr_size-2;i>=0;i--)
			{
				(*inv_A_)[i][k]=E[i][k];
				long int j;
				for(j=i+1;j<matr_size;j++)
				{
					(*inv_A_)[i][k]-=A_[i][j]*(*inv_A_)[j][k];
				};
				(*inv_A_)[i][k]/=A_[i][i];
			};
		};

	};

	return true;
}

void FSA_utils::multiply_matrices(
const std::vector<std::vector<double> > &A_,
const std::vector<std::vector<double> > &B_,
std::vector<std::vector<double> > &res_)
{
	long int size1=(long int)A_.size();
	if(size1==0)
	{
		throw error("Error in FSA_utils::multiply_matrices\n",1);
	};
	long int size2=(long int)A_[0].size();
	if(size2==0)
	{
		throw error("Error in FSA_utils::multiply_matrices\n",1);
	};

	long int i,j,k;
	for(i=1;i<size1;i++)
	{
		if((long int)A_[i].size()!=size2)
		{
			throw error("Error in FSA_utils::multiply_matrices\n",1);
		};
	};

	if(size2!=(long int)B_.size())
	{
		throw error("Error in FSA_utils::multiply_matrices\n",1);
	};

	long int size3=(long int)B_[0].size();
	if(size3==0)
	{
		throw error("Error in FSA_utils::multiply_matrices\n",1);
	};

	for(i=1;i<size2;i++)
	{
		if((long int)B_[i].size()!=size3)
		{
			throw error("Error in FSA_utils::multiply_matrices\n",1);
		};
	};

	res_.clear();
	res_.resize(size1);
	for(i=0;i<size1;i++)
	{
		res_[i].resize(size3,0);
	};

	for(i=0;i<size1;i++)
	{
		for(j=0;j<size3;j++)
		{
			for(k=0;k<size2;k++)
			{
				res_[i][j]+=A_[i][k]*B_[k][j];
			};
		};
	};
}

void FSA_utils::multiply_matrix_and_vector(
const std::vector<std::vector<double> > &A_,
const std::vector<double> &y_,
std::vector<double> &res_)
{
	long int size1=(long int)A_.size();
	if(size1==0)
	{
		throw error("Error in FSA_utils::multiply_matrix_and_vector\n",1);
	};
	long int size2=(long int)A_[0].size();
	if(size2==0)
	{
		throw error("Error in FSA_utils::multiply_matrix_and_vector\n",1);
	};

	long int i,k;
	for(i=1;i<size1;i++)
	{
		if((long int)A_[i].size()!=size2)
		{
			throw error("Error in FSA_utils::multiply_matrix_and_vector\n",1);
		};
	};

	if(size2!=(long int)y_.size())
	{
		throw error("Error in FSA_utils::multiply_matrix_and_vector\n",1);
	};


	res_.clear();
	res_.resize(size1,0);

	for(i=0;i<size1;i++)
	{
		for(k=0;k<size2;k++)
		{
			res_[i]+=A_[i][k]*y_[k];
		};
	};
}

void FSA_utils::transpose_matrix(
const std::vector<std::vector<double> > &A_,
std::vector<std::vector<double> > &res_)
{
	long int size1=(long int)A_.size();
	if(size1==0)
	{
		res_.clear();
		return;
	};
	long int size2=(long int)A_[0].size();

	long int i,j;
	for(i=1;i<size1;i++)
	{
		if((long int)A_[i].size()!=size2)
		{
			throw error("Error in FSA_utils::transpose_matrix\n",1);
		};
	};

	res_.clear();
	res_.resize(size2);
	for(i=0;i<size2;i++)
	{
		res_[i].resize(size1);
	};

	for(i=0;i<size2;i++)
	{
		for(j=0;j<size1;j++)
		{
			res_[i][j]=A_[j][i];
		};
	};
}

void FSA_utils::print_matrix(
const std::vector<std::vector<double> > A_)
{
	long int i,j;
	for(i=0;i<(long int)A_.size();i++)
	{
		for(j=0;j<(long int)A_[i].size();j++)
		{
			if(j<(long int)A_[i].size()-1)
			{
				cout<<A_[i][j]<<"\t";
			}
			else
			{
				cout<<A_[i][j]<<"\n";
			};
		};
	};
}

void FSA_utils::process_random_factor(
long int &random_factor_,
bool *rand_flag_)
{
	if(rand_flag_)
	{
		*rand_flag_=true;
	};

	if(random_factor_<0)
	{
		random_factor_=(long int)time(NULL);
		#ifndef _MSC_VER //UNIX program
			struct timeval tv;
			struct timezone tz;
			gettimeofday(&tv, &tz);
			random_factor_+=tv.tv_usec*10000000;
		#else
			struct _timeb timebuffer;
			char *timeline;
			_ftime( &timebuffer );
			timeline = ctime( & ( timebuffer.time ) );
			random_factor_+=timebuffer.millitm*10000000;
		#endif

		random_factor_=abs(random_factor_);

		if(rand_flag_)
		{
			*rand_flag_=false;
		};

	};
}

double FSA_utils::standard_deviation(//standard deviation for average of elements of vect_
long int dim_,
double *vect_)
{
	if(dim_<1)
	{
		throw error("Unexpected error in FSA_utils::standard_deviation\n",1);
	};

	if(dim_==1)
	{
		return 0;
	};

	double E=0;
	long int i;
	for(i=0;i<dim_;i++)
	{
		E+=vect_[i];
	};
	E/=(double)dim_;

	double E2=0;
	for(i=0;i<dim_;i++)
	{
		E2+=(vect_[i]-E)*(vect_[i]-E);
	};

	E2=sqrt(E2/(double)((dim_-1)*dim_));

	return E2;

}

double FSA_utils::average(//average of elements of vect_
long int dim_,
double *vect_)
{
	if(dim_<1)
	{
		throw error("Unexpected error in FSA_utils::average\n",1);
	};

	double E=0;
	long int i;
	for(i=0;i<dim_;i++)
	{
		E+=vect_[i];
	};
	E/=(double)dim_;

	return E;

}

void FSA_utils::read_string(
ifstream &f_,
string &st_,
bool &end_of_file_flag_)
{
	if(f_.eof())
	{
		end_of_file_flag_=true;
		return;
	};

	st_="";
	end_of_file_flag_=false;
	while(!end_of_file_flag_)
	{
		getline(f_,st_);

		if(st_.size()>0)
		{
			bool flag=true;
			long int i;
			for(i=0;i<(long int)st_.size();i++)
			{
				if(!(st_[i]==' '||(long int)st_[i]==13))
				{
					flag=false;
					break;
				};
			};
			if(flag)
			{
				continue;
			};
		};

		if(st_.size()>0)
		{
			if((long int)st_[st_.size()-1]==13)
			{
				st_.resize(st_.size()-1);
			};
			break;
		};
		if(f_.eof())
		{
			end_of_file_flag_=true;
		};
	};
}

void FSA_utils::read_sequences_for_alingment(

string input_file_name_,

long int &number_of_letters1_,//number of letters for the sequence 1
long int &number_of_letters2_,//number of letters for the sequence 2

char *&alphabet1_,//alphabet letters for the sequence #1
char *&alphabet2_,//alphabet letters for the sequence #2

long int& number_of_sequences_,

string *&headers_,
long int *&lengths1_,//lengths of the sequences #1
long int *&lengths2_,//lengths of the sequences #2
long int **&sequences1_,//the first index numerates sequences; the second - sequence letters
long int **&sequences2_)
{

	long int max_c=1000;
	long int *letter_to_int1=new long int[max_c];
	long int *letter_to_int2=new long int[max_c];

	long int i;
	for(i=0;i<max_c;i++)
	{
		letter_to_int1[i]=-1;
		letter_to_int2[i]=-1;
	};

	for(i=0;i<number_of_letters1_;i++)
	{
		letter_to_int1[(size_t)alphabet1_[i]]=i;
	};
	for(i=0;i<number_of_letters2_;i++)
	{
		letter_to_int2[(size_t)alphabet2_[i]]=i;
	};

	long int k;
	ifstream f(input_file_name_.data());
	if(!f)
	{
		throw error("Error - the input file "+input_file_name_+" is not found\n",1);
	};

	

	long int count=0;

	vector<long int *> sequences1_vect;
	vector<long int *> sequences2_vect;
	vector<long int > lengths1_vect;
	vector<long int > lengths2_vect;
	vector<string> headers_vect;

	for( ; ; )
	{

		bool end_of_file_flag;
		string header;

		read_string(
		f,
		header,
		end_of_file_flag);

		
		if(end_of_file_flag)
		{
			if(header!="")
			{
				throw error("Error - the file "+input_file_name_+" is incorrect\n",1);
			};

			if(count==0)
			{
				throw error("Error - the file "+input_file_name_+" is incorrect\n",1);
			};
			break;
		};


		if(header[0]!='>')
		{
			throw error("Error - the file "+input_file_name_+" is incorrect\n",1);
		};
		

		string st1="";
		string st2="";
		if(f.eof())
		{
			throw error("Error - the file "+input_file_name_+" is incorrect\n",1);
		};

		read_string(
		f,
		st1,
		end_of_file_flag);

		if(st1.size()==0||end_of_file_flag)
		{
			throw error("Error - the file "+input_file_name_+" is incorrect\n",1);
		};

		read_string(
		f,
		st2,
		end_of_file_flag);

		if(st2.size()==0||end_of_file_flag)
		{
			throw error("Error - the file "+input_file_name_+" is incorrect\n",1);
		};
		count++;

		//checking the sequences
		long int *sequences1_tmp=new long int [st1.size()];
		FSA_utils::assert_mem(sequences1_tmp);
		long int *sequences2_tmp=new long int [st2.size()];
		FSA_utils::assert_mem(sequences2_tmp);

		headers_vect.push_back(header);
		sequences1_vect.push_back(sequences1_tmp);
		sequences2_vect.push_back(sequences2_tmp);
		lengths1_vect.push_back((long int)st1.size());
		lengths2_vect.push_back((long int)st2.size());

		for(k=0;k<(long int)st1.size();k++)
		{

			long int long_tmp=letter_to_int1[(size_t)st1[k]];
			if(long_tmp>=0)
			{
				sequences1_tmp[k]=long_tmp;
			}
			else
			{
				throw error("Error - the file "+input_file_name_+" is incorrect: a non-alphabet letter with the code "+FSA_utils::long_to_string((long int)st1[k])+"\n",1);
			};
		};
		for(k=0;k<(long int)st2.size();k++)
		{

			long int long_tmp=letter_to_int2[(size_t)st2[k]];
			if(long_tmp>=0)
			{
				sequences2_tmp[k]=long_tmp;
			}
			else
			{
				throw error("Error - the file "+input_file_name_+" is incorrect: a non-alphabet letter with the code "+FSA_utils::long_to_string((long int)st2[k])+"\n",1);
			};
		};

	};

	f.close();

	delete[]letter_to_int1;
	delete[]letter_to_int2;

	//allocate memory for the arrays
	number_of_sequences_=count;
	sequences1_=new long int*[number_of_sequences_];
	FSA_utils::assert_mem(sequences1_);
	sequences2_=new long int*[number_of_sequences_];
	FSA_utils::assert_mem(sequences2_);
	
	lengths1_=new long int[number_of_sequences_];
	FSA_utils::assert_mem(lengths1_);
	lengths2_=new long int[number_of_sequences_];
	FSA_utils::assert_mem(lengths2_);
	headers_=new string[number_of_sequences_];
	FSA_utils::assert_mem(headers_);

	for(k=0;k<number_of_sequences_;k++)
	{
		headers_[k]=headers_vect[k];
		lengths1_[k]=lengths1_vect[k];
		lengths2_[k]=lengths2_vect[k];
		sequences1_[k]=sequences1_vect[k];
		sequences2_[k]=sequences2_vect[k];
	};

}

