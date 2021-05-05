#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <iostream>
#include "LEGACY.h"
using namespace std;

int main()
{
    LEGACY test_5_RA;
	//notice that the primitive root of unity exist
	//when m | modular - 1 (modular -1 can be divided by m)	
	long long n = 5;	
	long long modular_n = 41;
	long long prou_n;	
	long long FFT_data_out[n];
	long long DFT_data_out[n];
	long long data_in[n]/*={1,2,3,4}*/;
	
	prou_n = test_5_RA.find_prou( n, modular_n) ;
	
	for (int i = 0; i < n ; i++){
		data_in[i] = i ;
	}
	std::cout << "data_in =  ";	
	for (int i = 0; i< 5 ; i++){	
		std::cout << data_in[i] << " ";	
	}
//---------Test Rader 5-point with 4-FFT--------------------------------------------------------	
	
    long long data_in_RA[4] = {data_in[1],data_in[2],data_in[4],data_in[3]};
	long long data_out_RA[4];
	long long prou_n_RA;
	
	prou_n_RA = test_5_RA.find_prou( 4, modular_n) ;	
	test_5_RA.FFT(data_out_RA , data_in_RA , 4 , prou_n_RA , modular_n);
	std::cout << "\nFFT_out_RA =  ";	
	for (int i = 0; i< 4 ; i++){	
		std::cout << data_out_RA[i] << " ";	
	}
	//---------pre-compute tw factor-----------------
	long long prou_set_RA[4];
	
	prou_set_RA[0] = test_5_RA.prou_power(prou_n, 1, modular_n);
	prou_set_RA[1] = test_5_RA.prou_power(prou_n, 3, modular_n);
	prou_set_RA[2] = test_5_RA.prou_power(prou_n, 4, modular_n);	
	prou_set_RA[3] = test_5_RA.prou_power(prou_n, 2, modular_n);			
	
	std::cout <<"\nprou_set_RA = ";		
	for (int i = 0; i< 4 ; i++){	
		std::cout << prou_set_RA[i] << " ";	
	}
	
	long long prou_FFT_out[4];
	test_5_RA.FFT(prou_FFT_out , prou_set_RA , 4 , prou_n_RA , modular_n);
	std::cout << "\nprou_FFT_out =  ";	
	for (int i = 0; i< 4 ; i++){	
		std::cout << prou_FFT_out[i] << " ";	
	}	
	
    //----------------elementwise-mul--------------	
	long long ele_mul_RA[4];
	for (int i = 0; i< 4 ; i++){	
		ele_mul_RA[i] = (data_out_RA[i]*prou_FFT_out[i]) % modular_n ;
	}
	std::cout <<"\nele_mul_RA = ";		
	for (int i = 0; i< 4 ; i++){	
		std::cout << ele_mul_RA[i] << " ";	
	}	
	
	//----------------IFFT--------------------------
	long long IFFT_out_RA[4];
	test_5_RA.IFFT( IFFT_out_RA , ele_mul_RA , 4 , prou_n_RA , modular_n);
	std::cout << "\nIFFT_out_RA =  ";	
	for (int i = 0; i< 4 ; i++){	
		std::cout << IFFT_out_RA[i] << " ";	
	}		
	//--------------add d0--------------------------
	long long RA_out_tmp[5];
	RA_out_tmp[0] = data_in[0]+data_in[1]+data_in[2]+data_in[3]+data_in[4] ;
	RA_out_tmp[1] = IFFT_out_RA[0] + data_in[0];
	RA_out_tmp[2] = IFFT_out_RA[1] + data_in[0];
	RA_out_tmp[3] = IFFT_out_RA[2] + data_in[0];
	RA_out_tmp[4] = IFFT_out_RA[3] + data_in[0];
	std::cout << "\nRA_out_tmp =  ";
	for (int i = 0; i< n ; i++){	
		std::cout << RA_out_tmp[i] << " ";	
	}
	//---------------re-index-----------------------
	long long RA_out[5];	
	RA_out[0] = RA_out_tmp[0] ;
	RA_out[1] = RA_out_tmp[1] ;
	RA_out[2] = RA_out_tmp[4] ;
	RA_out[3] = RA_out_tmp[2] ;	
	RA_out[4] = RA_out_tmp[3] ;
	std::cout << "\nRA_out =  ";	
	for (int i = 0; i< n ; i++){	
		std::cout << RA_out[i] << " ";	
	}	
	
//-----------Golden DFT----------------------------------------------------
	test_5_RA.DFT(DFT_data_out, data_in, n, prou_n, modular_n);
	std::cout << "\nDFT_out =  ";	
	for (int i = 0; i< n ; i++){	
		std::cout << DFT_data_out[i] << " ";	
	}
    printf("\nDone \n");
    return 0;
}