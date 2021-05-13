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
	LEGACY test;
	//notice that the primitive root of unity exist
	//when m | modular - 1 (modular -1 can be divided by m)	
	long long n = 17;	
	long long modular_n ;
	long long prou_n;	
	long long FFT_data_out[n];
	long long DFT_data_out[n];
	long long data_in[n]/*={0,1,2,3,4,5,6}*/;
	//cout << test.find_prime(11,8) <<endl;	
	modular_n = test.find_prime(n,8);
	
	for (int i = 0; i < n ; i++){
		data_in[i] = i;
	}	
	
	
	prou_n = test.find_prou( n, modular_n) ;
	test.Rader(FFT_data_out , data_in , n , prou_n , modular_n);
	test.DFT(DFT_data_out, data_in, n, prou_n, modular_n);
	/*
	std::cout << "\nDFT_out =  ";	
	for (int i = 0; i< n ; i++){	
		std::cout << DFT_data_out[i] << " ";	
	}
    printf("\nDone \n");	
	*/
	int k = 0;
	for (int i = 0; i< n ; i++){	
		if(FFT_data_out[i] != DFT_data_out[i])	{
			cout << "fail" <<endl;
			break;
		}
		else {
			k++;
		}
		if(k == n){
			cout << "done" <<endl;
		}
			
	}	
	

}