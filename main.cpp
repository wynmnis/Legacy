#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <iostream>
#include<bits/stdc++.h> 
#include <chrono>
#include <thread>

#include "LEGACY.h"
using namespace std;

int main()
{
		
	ZZ n ;
	n = 2821; // 7*13*31
	LEGACY test_2821(n);
////////////////////////////////////////////////////////////	
	long long N1 = 7;
	long long N2 = 13;
	long long N3 = 31;
	long long N1N2_inv = test_2821.find_inv(N1*N2,N3);
	long long N1N3_inv = test_2821.find_inv(N1*N3,N2);
	long long N2N3_inv = test_2821.find_inv(N3*N2,N1);	
	long long N = N1*N2*N3;
	long long powerof2 = ceil(log(2*N-3)/log(2.0));
	long long blue_m_prime = pow(2, powerof2);
	cout << blue_m_prime << endl;
	int bank_num = N3;
	int addr_num = N1*N2;
	long long memory[bank_num][addr_num];
	long long data_in[N];
	long long data_out[N];	
	long long LCM = test_2821.LCM(test_2821.LCM(N1-1,N2-1),N3-1);
	long long modular_n = test_2821.find_prime(N*LCM, 1);
	cout << " modular_n = " << modular_n << endl;	
	long long prou_N = test_2821.find_prou(N, modular_n);
	//cout << " prou = " << prou << endl;	
	long long prou_N1 = PowerMod(prou_N, N/N1, modular_n);
	long long prou_N2 = PowerMod(prou_N, N/N2, modular_n);
	long long prou_N3 = PowerMod(prou_N, N/N3, modular_n);	
	//long long prou_blue_m_prime = test_2821.find_prou(blue_m_prime, modular_n);
	long long prou_blue_m_prime = 80346;
	cout << " prou_blue_m_prime = " << prou_blue_m_prime << endl;	
	long long rand_num;	
	
	srand(5);
	
	for (int i = 0; i < N ; i++){
		rand_num = rand();
		data_in[i] = rand_num % modular_n ;
		//cout << data_in[i] << endl;
	}	
	
	
	
	
	
//-----preprocessing-----//
long long prou_N1_neg_1 = test_2821.find_prou( N1-1 , modular_n);
long long prou_N2_neg_1 = test_2821.find_prou( N2-1 , modular_n);
long long prou_N3_neg_1 = test_2821.find_prou( N3-1 , modular_n);
long long Rader_N1_In_index[N1];
long long Rader_N2_In_index[N2];
long long Rader_N3_In_index[N3];
long long Rader_N1_Out_index[N1];
long long Rader_N2_Out_index[N2];
long long Rader_N3_Out_index[N3];
long long N1_gen = test_2821.find_gen(N1);
long long N2_gen = test_2821.find_gen(N2);
long long N3_gen = test_2821.find_gen(N3);

long long tmp = 1;


Rader_N1_In_index[0] = 0;
Rader_N2_In_index[0] = 0;
Rader_N3_In_index[0] = 0;
Rader_N1_Out_index[0] = 0;
Rader_N2_Out_index[0] = 0;
Rader_N3_Out_index[0] = 0;
Rader_N1_Out_index[1] = 1;
Rader_N2_Out_index[1] = 1;
Rader_N3_Out_index[1] = 1;
long long Tf_In_N1_neg_1[N1-1];
long long Tf_In_N2_neg_1[N2-1];
long long Tf_In_N3_neg_1[N3-1];
long long Tf_Out_N1_neg_1[N1-1];
long long Tf_Out_N2_neg_1[N2-1];
long long Tf_Out_N3_neg_1[N3-1];

	for (int i = 0; i < N1-1; i++){
		tmp = PowerMod(N1_gen, i, N1);
		Rader_N1_In_index[i+1] = tmp ;
	}
	tmp = 1;
	for (int i = 0; i < N2-1; i++){
		tmp = PowerMod(N2_gen, i, N2);
		Rader_N2_In_index[i+1] = tmp ;
	}
	tmp = 1;
	for (int i = 0; i < N3-1; i++){
		tmp = PowerMod(N3_gen, i, N3);
		Rader_N3_In_index[i+1] = tmp ;
	}	

	for (int i = 2; i < N1; i++){
		Rader_N1_Out_index[i] = Rader_N1_In_index[N1-i+1];
	}	
	for (int i = 2; i < N2; i++){
		Rader_N2_Out_index[i] = Rader_N2_In_index[N2-i+1];
	}	
	for (int i = 2; i < N3; i++){
		Rader_N3_Out_index[i] = Rader_N3_In_index[N3-i+1];
	}

	for (int i = 1; i < N1; i++){
		Tf_In_N1_neg_1[i-1] = PowerMod(prou_N1, Rader_N1_Out_index[i], modular_n);
	}
	for (int i = 1; i < N2; i++){
		Tf_In_N2_neg_1[i-1] = PowerMod(prou_N2, Rader_N2_Out_index[i], modular_n);
	}
	for (int i = 1; i < N3; i++){
		Tf_In_N3_neg_1[i-1] = PowerMod(prou_N3, Rader_N3_Out_index[i], modular_n);
	}	
	
	test_2821.DFT(Tf_Out_N1_neg_1, Tf_In_N1_neg_1, N1-1, prou_N1_neg_1, modular_n) ;	
	test_2821.DFT(Tf_Out_N2_neg_1, Tf_In_N2_neg_1, N2-1, prou_N2_neg_1, modular_n) ;
	test_2821.DFT(Tf_Out_N3_neg_1, Tf_In_N3_neg_1, N3-1, prou_N3_neg_1, modular_n) ;
	
/*
long long RA_out[N1];
long long DFT_out[N1];
long long RA_in[N1];

	for (int i = 0; i < N1; i++){
		RA_in[i] = i;
		cout << RA_in[i] << endl;		
	}


	test_2821.Rader_DFT_0118(RA_out,RA_in,N1,Tf_Out_N1_neg_1,Rader_N1_In_index,Rader_N1_Out_index,prou_N1_neg_1,modular_n);
	
	cout << "RA_out" << endl;	
	for (int i = 0; i < N1; i++){
		cout << RA_out[i] << endl;
	}	
	
	test_2821.DFT(DFT_out, RA_in, N1, prou_N1, modular_n) ;	
	
	cout << "DFT_out" << endl;	
	for (int i = 0; i < N1; i++){
		cout << DFT_out[i] << endl;
	}		
*/	
	
/////////////input reindex/////////////////////////
// n = n1(N2*N3) + n2(N1N3) + n3(N1N2)
// n = n1*403 + n2*217 + n3*91
// n ---> (n1,n2,n3)

auto start = std::chrono::high_resolution_clock::now();

int bank,addr;
long long index;
int n1_,n2_,n3_;
for(int n3 = 0; n3 < N3; n3++){
	for(int n2 = 0; n2 < N2; n2++){
		for(int n1 = 0; n1 < N1; n1++){
			//index = (n1*(N2*N3) + n2*(N1*N3) + n3*(N1*N2))%N;
			if(n1 == 0)
				n1_ = 0;
			else	
				n1_ = PowerMod(N1_gen, n1-1, N1);
			if(n2 == 0)
				n2_ = 0;
			else	
				n2_ = PowerMod(N2_gen, n2-1, N2);
			if(n3 == 0)
				n3_ = 0;
			else	
				n3_ = PowerMod(N3_gen, n3-1, N3);	
				
			index = (n1_*(N2*N3) + n2_*(N1*N3) + n3_*(N1*N2))%N;
			bank = (n1 + n2 + n3) % bank_num;
			addr = N2*n1 + n2;
			memory[bank][addr] = data_in[index];
		}
	}
}



long long N1_tmp_in[N1];
long long N1_tmp_out[N1];
// N1-point FFT (0 ~ (N1 - 1) , c , c ) 
for(int n3 = 0; n3 < N3; n3++){
	for(int n2 = 0; n2 < N2; n2++){
		for(int n1 = 0; n1 < N1; n1++){
			bank = (n1 + n2 + n3) % bank_num;
			addr = N2*n1 + n2;			
			N1_tmp_in[n1] = memory[bank][addr];
		}
		// NTT 
		//test_2821.DFT(N1_tmp_out, N1_tmp_in, N1, prou_N1, modular_n);
		//test_2821.Rader(N1_tmp_out, N1_tmp_in, N1, prou_N1, modular_n);
		//test_2821.Rader_DFT(N1_tmp_out, N1_tmp_in, N1, prou_N1, modular_n);
		test_2821.Rader_DFT_0118(N1_tmp_out,N1_tmp_in,N1,Tf_Out_N1_neg_1,Rader_N1_In_index,Rader_N1_Out_index,prou_N1_neg_1,modular_n);
		for(int n1 = 0; n1 < N1; n1++){
			bank = (n1 + n2 + n3) % bank_num;
			addr = N2*n1 + n2;			
			memory[bank][addr] = N1_tmp_out[n1];
		}		
	}
}


long long N2_tmp_in[N2];
long long N2_tmp_out[N2];
// N2-point FFT (c , 0 ~ (N2 - 1) , c ) 
for(int n3 = 0; n3 < N3; n3++){
	for(int n1 = 0; n1 < N1; n1++){
		for(int n2 = 0; n2 < N2; n2++){
			bank = (n1 + n2 + n3) % bank_num;
			addr = N2*n1 + n2;			
			N2_tmp_in[n2] = memory[bank][addr];
		}
		// NTT 
		//test_2821.DFT(N2_tmp_out, N2_tmp_in, N2, prou_N2, modular_n);
		test_2821.Rader_DFT_0118(N2_tmp_out,N2_tmp_in,N2,Tf_Out_N2_neg_1,Rader_N2_In_index,Rader_N2_Out_index,prou_N2_neg_1,modular_n);
		for(int n2 = 0; n2 < N2; n2++){
			bank = (n1 + n2 + n3) % bank_num;
			addr = N2*n1 + n2;			
			memory[bank][addr] = N2_tmp_out[n2];
		}		
	}
}

long long N3_tmp_in[N3];
long long N3_tmp_out[N3];
// N2-point FFT (c , 0 ~ (N2 - 1) , c ) 
for(int n2 = 0; n2 < N2; n2++){
	for(int n1 = 0; n1 < N1; n1++){
		for(int n3= 0; n3 < N3; n3++){
			bank = (n1 + n2 + n3) % bank_num;
			addr = N2*n1 + n2;			
			N3_tmp_in[n3] = memory[bank][addr];
		}
		// NTT 
		//test_2821.DFT(N3_tmp_out, N3_tmp_in, N3, prou_N3, modular_n);
		test_2821.Rader_DFT_0118(N3_tmp_out,N3_tmp_in,N3,Tf_Out_N3_neg_1,Rader_N3_In_index,Rader_N3_Out_index,prou_N3_neg_1,modular_n);
		for(int n3 = 0; n3 < N3; n3++){
			bank = (n1 + n2 + n3) % bank_num;
			addr = N2*n1 + n2;			
			memory[bank][addr] = N3_tmp_out[n3];
		}		
	}
}



// output re-index
long long N1_gen_inv = test_2821.find_inv(N1_gen,N1);
long long N2_gen_inv = test_2821.find_inv(N2_gen,N2);
long long N3_gen_inv = test_2821.find_inv(N3_gen,N3);
int k1_,k2_,k3_;
for(int k3 = 0; k3 < N3; k3++){
	for(int k2 = 0; k2 < N2; k2++){
		for(int k1 = 0; k1 < N1; k1++){
			
			if(k1 == 0)
				k1_ = 0;
			else	
				k1_ = PowerMod(N1_gen_inv, k1-1, N1);
			if(k2 == 0)
				k2_ = 0;
			else	
				k2_ = PowerMod(N2_gen_inv, k2-1, N2);
			if(k3 == 0)
				k3_ = 0;
			else	
				k3_ = PowerMod(N3_gen_inv, k3-1, N3);
				
			//cout << "k1_" << k1_ << endl;	
			//index = (k1*(N2*N3)*N2N3_inv + k2*(N1*N3)*N1N3_inv + k3*(N1*N2)*N1N2_inv) % N;	
			index = (k1_*(N2*N3)*N2N3_inv + k2_*(N1*N3)*N1N3_inv + k3_*(N1*N2)*N1N2_inv) % N;	
			bank = (k1 + k2 + k3) % bank_num;
			addr = N2*k1 + k2;
			data_out[index] = memory[bank][addr];
		}
	}
}

auto end = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::milli> float_ms = end - start;
std::cout << "PFA_NTT time is " << float_ms.count() << " milliseconds" << std::endl;



long long golden_out[N];
//Golden 

auto start2 = std::chrono::high_resolution_clock::now();
test_2821.DFT(golden_out, data_in, N, prou_N, modular_n);
auto end2 = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::milli> float_ms2 = end2 - start2;
std::cout << "NTT time is " << float_ms2.count() << " milliseconds" << std::endl;

long long blue_data_in[blue_m_prime];
	for (int i = 0; i < blue_m_prime ; i++){
		rand_num = rand();
		blue_data_in[i] = rand_num % modular_n ;
		//cout << data_in[i] << endl;
	}	

// bluestein
long long blue_out[blue_m_prime];
auto start3 = std::chrono::high_resolution_clock::now();
//test_2821.FFT(blue_out, blue_data_in, blue_m_prime, prou_blue_m_prime, modular_n);
auto end3 = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::milli> float_ms3 = end3 - start3;
std::cout << "Blue time is " << float_ms3.count() << " milliseconds" << std::endl;




for(int i = 0; i < N ; i++){
	if(golden_out[i] != data_out[i]){
		cout << "num " << i << " error !" << endl;
		break;
	}
	if(i == N1){		
		cout << golden_out[i] << endl;
		cout << data_out[i] << endl;
	}
	else if(i == N2){
		cout << golden_out[i] << endl;
		cout << data_out[i] << endl;		
	} 
	else if(i == N3){
		cout << golden_out[i] << endl;
		cout << data_out[i] << endl;		
	} 
	
	if(i == N-1){
		cout << golden_out[i] << endl;
		cout << data_out[i] << endl;
		cout << "success !" << endl;		
	} 
}
/*
cout << "data_out" << endl;  
for(int i = 0; i < N ; i++){
	cout << data_out[i] << endl;
}
cout << "golden_out" << endl;  
for(int i = 0; i < N ; i++){
	cout << golden_out[i] << endl;
}*/	

}