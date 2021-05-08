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
    LEGACY test_15_PFA;
	//notice that the primitive root of unity exist
	//when m | modular - 1 (modular -1 can be divided by m)	
	long long m = 15;	
	long long m1 = 3;	
	long long m2 = 5;	
	long long modular_n = 31;
	long long prou_n;	
	long long DFT_data_out[m];
	long long DFT_data[m];	
	long long PFA_data_out[m];	
	long long IDFT_data_out[m];	
	long long data_in[m]/*={1,2,3,4,5}*/;
	long long data_tmp[m];
	
	
	prou_n = test_15_PFA.find_prou(m,modular_n);
	
	test_15_PFA.PFA2(PFA_data_out,data_in,m1,m2,2,2,prou_n,modular_n);

	
	std::cout << "PFA_data_out =  ";	
	for (int i = 0; i< m ; i++){	
		std::cout << PFA_data_out[i] << " ";	
	}
	std::cout << "\n  ";
		
	test_15_PFA.DFT(DFT_data_out, data_in, m, prou_n, modular_n);           	
	std::cout << "DFT_data_out =  ";	
	for (int i = 0; i< m ; i++){	
		std::cout << DFT_data_out[i] << " ";	
	}
	std::cout << "\n  ";
		
}