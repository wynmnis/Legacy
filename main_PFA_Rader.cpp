#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <iostream>
#include<bits/stdc++.h> 

#include "LEGACY.h"
using namespace std;

int main()
{
    LEGACY test_105_PFA;
	//notice that the primitive root of unity exist
	//when m | modular - 1 (modular -1 can be divided by m)	
	/*
	long long m = 255;	
	long long m1 = 3;	
	long long m2 = 85;	
	long long s_m1 = 5;	
	long long s_m2 = 17;		
	long long modular_n = test_105_PFA.find_prime(m,4);	
	*/
/*
	long long m = 105;	
	long long m1 = 3;	
	long long m2 = 35;	
	long long s_m1 = 5;	
	long long s_m2 = 7;		
	long long modular_n = test_105_PFA.find_prime(m,4);		
*/
	
	long long m = 21845;	
	long long m1 = 5;	
	long long m2 = 4369;	
	long long s_m1 = 17;	
	long long s_m2 = 257;		
	long long modular_n = test_105_PFA.find_prime(m,8);	
		
	long long prou_n;	
	long long DFT_data_out[m];
	long long DFT_data[m];	
	long long PFA_data_out[m];	
	long long IDFT_data_out[m];
	long long error[m];
	long long data_in[m];
	long long data_in2[m]={5,4,3,2,1};
	long long data_tmp[m];
	long long m1_inv,m2_inv,s_m1_inv,s_m2_inv;
	prou_n = test_105_PFA.find_prou( m, modular_n);
	m1_inv = test_105_PFA.find_inv(m1, m2);
	m2_inv = test_105_PFA.find_inv(m2, m1);
	s_m1_inv = test_105_PFA.find_inv(s_m1, s_m2);
	s_m2_inv = test_105_PFA.find_inv(s_m2, s_m1);	
	
/*	
	cout << "modulus = " <<modular_n << endl;
	cout << " m1_inv = " << m1_inv << endl;
	cout << " m2_inv = " << m2_inv << endl;
	cout << " s_m1_inv = " << s_m1_inv << endl;
	cout << " s_m2_inv = " << s_m2_inv << endl;	
*/	
/*
	cout << " 1285_inv(17) = " << test_105_PFA.find_inv(1285, 17) << endl;	
	cout << " 85_inv(257) = " << test_105_PFA.find_inv(85, 257) << endl;
*/	
	//------------------------------------------
	for (int i = 0; i < m ; i++){
		data_in[i] = i*2 ;
	}			



	long long index_tmp[m];
	long long PFA_data_tmp[m];		
	//-------------------------------------------
	test_105_PFA.PFA3_v4(PFA_data_out, data_in, m1,m2,s_m1,s_m2,m1_inv,m2_inv,s_m1_inv,s_m2_inv, prou_n, modular_n);           	
	/*std::cout << "PFA_data_out =  ";	

	
	
	for (int i = 0; i< m ; i++){	
		std::cout << PFA_data_out[i] << endl;	
	}
	*/
	
	std::cout << "\n  ";
	test_105_PFA.DFT(DFT_data_out, data_in, m, prou_n, modular_n);
/*
	std::cout << "DFT_data_out =  ";	
	for (int i = 0; i< m ; i++){	
		std::cout << DFT_data_out[i] <<" ";	
	}
	cout <<endl ;
	
	std::cout << "FFT_data_out =  ";	
	for (int i = 0; i< m ; i++){	
		std::cout << PFA_data_out[i] <<" ";	
	}
	cout <<endl ;
*/

	
	//std::cout << "FFT_data_out =  ";	
	for (int i = 0; i< m ; i++){	
		error[i] = PFA_data_out[i] - DFT_data_out[i];
		//std::cout << error[i] << " ";	
	}
	cout <<endl ;	
	int k = 0;
	for (int i = 0; i< m ; i++){	
		if(PFA_data_out[i] != DFT_data_out[i])	{
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
		
	
	
	
	
	
	std::cout << "\n  ";		
}