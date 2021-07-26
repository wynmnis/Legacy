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
/* 	long long amount;
	long long factor[20] = {0};
    amount = test.Factorize(105);
	cout << endl; */


	long long m = 1024;
 	long long modular = test.find_prime(m,6);	
/*	long long RA_out[m];
	long long RA_in[m];	
	long long DFT_out[m];
	long long DFT_in[m]; */	
	
/* 	for(int i = 0; i < m; i++){
		RA_in[i] = i;
		DFT_in[i] = i;
	}
	
	long long prou = test.find_prou(m, modular);
	test.Rader_DFT(RA_out, RA_in, m, prou, modular);
	test.DFT(DFT_out, DFT_in, m, prou, modular);
	for(int i = 0; i < m; i++){
		cout << RA_out[i] << endl;
	}
	cout << endl;
	
	for(int i = 0; i < m; i++){
		cout << DFT_out[i] << endl;
	} */
	vector<ZZ> output(m);
	vector<ZZ> input(m);	
	
	
	//test.Config_PFA_Rader_FFT(output, input, (ZZ)m, (ZZ)modular);
	
	
	test.FFT_1024_radix2(output, input, m, (ZZ)modular);
	int DFT_m = 16;
	
	long long DFT_in[DFT_m] ;
	long long DFT_out[DFT_m] ;
	long long prou = 3/*test.find_prou(4, 17)*/;	
	//cout << "prou=" << prou << endl;
	for (int i = 0; i < DFT_m; i++){
		DFT_in[i] = i+1; 
	}
	
	
	test.DFT(DFT_out, DFT_in, DFT_m, prou, 17);
	
	
	cout << endl;
	for(int i = 0; i < DFT_m; i++){
		//cout << DFT_out[i] << endl;
	}

	
	return 0;
}