#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include<bits/stdc++.h> 

#include "LEGACY.h"
using namespace std;
using namespace NTL;
int main()
{
	LEGACY test;

	long long m = 7;
 	long long modular = test.find_prime(105,8);	
	long long RA_out[m];
	long long RA_in[m];	

	vector<ZZ> input(m);	
	vector<ZZ> output(m);
	ZZ input_tmp[m];	
	
	int gen = test.find_gen(m);
	int gen_inv = test.find_inv(gen, m);	
	int in_idx[m] ={0} ;
	ZZ in_idx_tmp;
	int out_idx[m] ={0};
	ZZ out_tmp[m] ;	
	ZZ in_tmp[m] ;		


//cout << "input = " << endl;
for(int i = 0; i < m ; i++)
{
	//input_tmp[i] = 35*i + 21;
}
input_tmp[0] = 13581 ;
input_tmp[1] = 19940 ;
input_tmp[2] = 99    ;
input_tmp[3] = 7335  ;
input_tmp[4] = 7608  ;
input_tmp[5] = 4154  ;
input_tmp[6] = 8663  ;



for(int i = 0; i < m-1 ; i++)
{
	//PowerMod(in_idx[i+1], (ZZ)gen, i, (ZZ)m);
	in_idx[i+1] = test.prou_power(gen,i,m);	
	//cout << input[i+1] << endl;
}


for(int i = 0; i < m ; i++)
{
	//in_idx_tmp = in_idx[i];
	input[i] =  input_tmp[in_idx[i]] ;  //reindex
	//input[i] =  input_tmp[i] ;	    //no reindex
}


for(int i = 0; i < m ; i++)
{
	//PowerMod(input[i+1], (ZZ)gen, i, (ZZ)m);
	cout << input[i] << endl;
	//input[i] = in_tmp[i];
	
}

//test.FFT_1024_radix2_config(output, input, m, prou,(ZZ)modular);

	
for(int i = 0; i < m-1 ; i++)
{
	//PowerMod(out_idx[i+1], (ZZ)gen_inv, i, (ZZ)m);
	out_idx[i+1] = test.prou_power(gen_inv,i,m);
	//out_idx[i+1] %= m;
	//cout << input[i+1] << endl;
}

//cout << "out_idx = " << endl;
for(int i = 0; i < m ; i++)
{
	out_tmp[out_idx[i]] = output[i];
}


cout << "output = " << endl;
for(int i = 0; i < m; i++)
{
	cout << out_tmp[i] << endl;// correct 
	//cout << output[i] << endl;   // no reindex
}	
	
	
	
	return 0;
}