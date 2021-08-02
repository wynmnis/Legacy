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

	long long m = 21845;
 	long long modular = test.find_prime(m,8);	
	vector<ZZ> input(m);
	vector<ZZ> output(m);	
	ZZ prou = test.find_prou(m, (ZZ)modular);
	cout << "prou = " << prou << endl;
	
	
	for(int i = 0; i < m; i++){
		input[i] = i ;
	}	
	
	test.Config_PFA_Rader_FFT(output, input, (ZZ)m, prou, (ZZ)modular);

  std::ofstream output_ans("./output.txt");
  std::ofstream golden("./golden.txt");
	
		
	

	
	
	
	
	//cout << "modular = " << modular << endl;
	long long DFT_data_out[m] ; 
	long long DFT_in[m];
	
	for(int i = 0; i < m; i++){
		DFT_in[i] = i ;
		//cout << DFT_data_out[i] << endl ;
	}		
	
	test.DFT(DFT_data_out, DFT_in, m, test.ZZ2int(prou), modular);
	//cout << "golden = " << endl;
	for(int i = 0; i < m; i++){
		//cout << DFT_data_out[i] << endl ;
	}		


	for(int i = 0; i < m; i++){
		//cout << output[i] << endl;;
		output_ans << output[i] << endl;
		golden << DFT_data_out[i] << endl;
	}
	
	ZZ error[m];
	for (int i = 0; i< m ; i++){	
		//cout << DFT_data_out[i] << endl;
		error[i] = output[i] - (ZZ)DFT_data_out[i];
		//std::cout << error[i] << " ";	
	}
	cout <<endl ;	
	int k = 0;
	for (int i = 0; i< m ; i++){	
		if(output[i] != DFT_data_out[i])	{
			cout << "fail" <<endl;
			break;
		}
		else {
			k++;
		}
		if(k == m){
			cout << "done" <<endl;
		}
			
	}
	
	
	
	return 0;
}