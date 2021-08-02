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
	//notice that the primitive root of unity exist
	//when m | modular - 1 (modular -1 can be divided by m)	
	int test_time = 1;
//---------setting-----------	
	long long m = 21845;	
	long long m1 = 5;	
	long long m2 = 4369;	
	long long s_m1 = 17;	
	long long s_m2 = 257;		
	long long modular_n = test.find_prime(m,8);	
	//long long modular_n = 68303552103055361 ;
	//long long modular_n = 65699948296929281 ;
	
	ZZ modular_n_ZZ;
	modular_n_ZZ = modular_n;

	//cout << "modular_n_ZZ = " << modular_n_ZZ << endl;
		
	long long prou_n;
	ZZ prou_n_ZZ;
	long long DFT_data_out[m];
	long long DFT_data_out2[m];	
	long long C_K[m];	
	long long DFT_data[m];	
	long long PFA_data_out[m];	
	long long IDFT_data_out[m];
	long long error[m];
	long long data_in[m]={0};
	ZZ data_in_ZZ[m];
	long long debug_data[m];
	long long data_in2[m];
	long long data_tmp[m];
	long long m1_inv,m2_inv,s_m1_inv,s_m2_inv;
	long long phi_m ;
	prou_n_ZZ = test.find_prou( m, modular_n_ZZ);
	//prou_n = test.find_prou( m, modular_n);	
	//prou_n_ZZ = prou_n;
	cout << " prou_n_ZZ = " << prou_n_ZZ << endl;
	m1_inv = test.find_inv(m1, m2);
	m2_inv = test.find_inv(m2, m1);
	s_m1_inv = test.find_inv(s_m1, s_m2);
	s_m2_inv = test.find_inv(s_m2, s_m1);
	
	phi_m = test.Euler(m);
	srand( time(NULL));

	//cout << " 1 " << endl;
//--------Mem-----------------
	vector<vector<ZZ>> Dual_port_mem_4w_4096(4096);
	vector<ZZ> Dual_port_mem_1w_4096(4096);
	vector<vector<ZZ>> Dual_port_mem_4w_273(257+16);
	vector<ZZ> Dual_port_mem_1w_273(257+16);	
	vector<ZZ> four_word_mem(4);
	
	for(int i = 0; i < 4096; i++){	
		Dual_port_mem_4w_4096[i].resize(4);
	}
	for(int i = 0; i < 273; i++){	
		Dual_port_mem_4w_273[i].resize(4);
	}	
	
//------------------------ rader parameter m1 set up-------------------------

	long long tmp;
	long long m1_prime, gen1;
	long long rader_index_in1[m1];
	long long rader_index_out1[m1];	
	gen1 = test.find_gen(m1);
	if( (m1==3) || (m1==5) || (m1==17) || (m1==257)){
		m1_prime = m1-1;
	}
	else{
		m1_prime = test.find_m_prime(m1);
	}
	
	long long tw_FFT_index1[m1_prime];
	long long tw_FFT_in1[m1_prime];	
	ZZ tw_FFT_in1_ZZ[m1_prime];		
	long long m1_prou;
	ZZ m1_prou_ZZ;


	//cout << " 2 " << endl;
	
	//m1_prou = test.prou_power(prou_n, m2, modular_n);
	//PowerMod(m1_prou_ZZ, (ZZ)prou_n, (ZZ)m2, modular_n_ZZ);
	m1_prou_ZZ = test.find_prou( m1, (ZZ)modular_n);
	cout << " m1_prou_ZZ = " << m1_prou_ZZ << endl;	
	rader_index_in1[0] = 0;
	rader_index_out1[0] = 0;
	rader_index_in1[1] = 1;
	rader_index_out1[1] = 1;
	tmp = 1;
	for (int i = 2; i < m1 ; i++){
		tmp *= gen1;
		tmp %= m1;
		rader_index_in1[i] = tmp;
	}	
	for (int i = 2; i < m1 ; i++){
		rader_index_out1[i] = rader_index_in1[m1-i+1];
	}


	//cout << "m1_prou = " << m1_prou << endl;
	//cout << "m1_prou_ZZ = " << m1_prou_ZZ << endl;	


	//cout << " 3 " << endl;

		//cout << "tw_FFT_index[i]" <<endl;	
	for (int i = 0; i < m1_prime ; i++){
		if(i < (m1 - 1) )//0-5
			tw_FFT_index1[i] = rader_index_out1[i+1];
		else if(i > (m1_prime - m1 + 1))
			tw_FFT_index1[i] = rader_index_out1[i + m1 - m1_prime];			
		else		
			tw_FFT_index1[i] = 0; //00000
		//cout << tw_FFT_index1[i] << endl;	
	}
	
	for (int i = 0; i < m1_prime ; i++){
		if(tw_FFT_index1[i] == 0){
			//tw_FFT_in1[i] = 0;
			tw_FFT_in1_ZZ[i] = 0;
		}
		else{
			//tw_FFT_in1[i] = test.prou_power(m1_prou,tw_FFT_index1[i],modular_n);
			PowerMod(tw_FFT_in1_ZZ[i], m1_prou_ZZ, (ZZ)tw_FFT_index1[i], (ZZ)modular_n);
		}
	}	
	//cout << "prou = " << prou << endl;
	//cout << "m1_prou = " << m1_prou << endl;	
	//cout << "tw_FFT_in1 = " << endl;
	for (int i = 0; i < m1_prime ; i++){
			//cout << tw_FFT_in1[i] <<endl;
			//cout << tw_FFT_in1_ZZ[i] <<endl;			
	}
	
	long long tw_FFT_out1[m1_prime];
	ZZ tw_FFT_out1_ZZ[m1_prime];	
	long long m1_prime_prou;
	ZZ m1_prime_prou_ZZ;	
	//cout << " 4 " << endl;	
	//m1_prime_prou = test.find_prou(m1_prime, modular_n);
	m1_prime_prou_ZZ = test.find_prou( m1_prime, modular_n_ZZ);
	//m1_prime_prou_ZZ = m1_prime_prou;
	cout << " m1_prime_prou_ZZ = " << m1_prime_prou_ZZ << endl;		
	//cout << "m1_prime_prou = " << m1_prime_prou << endl;
	//cout << "m1_prime_prou_ZZ = " << m1_prime_prou_ZZ << endl;		
	
	//test.DFT(tw_FFT_out1, tw_FFT_in1, m1_prime, m1_prime_prou, modular_n) ;	
	test.DFT(tw_FFT_out1_ZZ, tw_FFT_in1_ZZ, m1_prime, m1_prime_prou_ZZ, modular_n_ZZ) ;	
	//cout << "tw_FFT_out1 = " << endl;
	for (int i = 0; i < m1_prime ; i++){
			//cout << tw_FFT_out1[i] <<endl;
			//cout << tw_FFT_out1_ZZ[i] <<endl;			
	}	
	

	//cout << " 4 " << endl;	
//-----------------------
	int m2_prime, gen2;
	long long rader_index_in2[s_m1];
	long long rader_index_out2[s_m1];	
	gen2 = test.find_gen(s_m1);	
	if( (s_m1==3) || (s_m1==5) || (s_m1==17) || (s_m1==257)){
		m2_prime = s_m1-1;
	}
	else{
		m2_prime = test.find_m_prime(s_m1);
	}
	long long tw_FFT_index2[m2_prime];	
	long long tw_FFT_in2[m2_prime];	
	long long m2_prou;	
	ZZ tw_FFT_in2_ZZ[m2_prime];	
	ZZ m2_prou_ZZ;	
	
	//m2_prou = test.prou_power(prou_n, m1*s_m2, modular_n);
	//PowerMod(m2_prou_ZZ, (ZZ)prou_n, (ZZ)m1*s_m2, modular_n_ZZ);	
	m2_prou_ZZ = test.find_prou( s_m1, (ZZ)modular_n);	
	cout << " m2_prou_ZZ = " << m2_prou_ZZ << endl;		
	rader_index_in2[0] = 0;
	rader_index_out2[0] = 0;
	rader_index_in2[1] = 1;
	rader_index_out2[1] = 1;
	tmp = 1;
	for (int i = 2; i < s_m1 ; i++){
		tmp *= gen2;
		tmp %= s_m1;
		rader_index_in2[i] = tmp;	
	}	
	for (int i = 2; i < s_m1 ; i++){
		rader_index_out2[i] = rader_index_in2[s_m1-i+1];
	}
	for (int i = 0; i < m2_prime ; i++){
		if(i < (s_m1 - 1) )//0-5
			tw_FFT_index2[i] = rader_index_out2[i+1];
		else if(i > (m2_prime - s_m1 + 1))
			tw_FFT_index2[i] = rader_index_out2[i + s_m1 - m2_prime];			
		else		
			tw_FFT_index2[i] = 0; //00000
		//cout << tw_FFT_index2[i] << endl;	
	}
	
	for (int i = 0; i < m2_prime ; i++){
		if(tw_FFT_index2[i] == 0){
			//tw_FFT_in2[i] = 0;
			tw_FFT_in2_ZZ[i] = 0;
		}			
		else{
			//tw_FFT_in2[i] = test.prou_power(m2_prou,tw_FFT_index2[i],modular_n);
			PowerMod(tw_FFT_in2_ZZ[i], m2_prou_ZZ, (ZZ)tw_FFT_index2[i], (ZZ)modular_n);			
		}
	}
	long long tw_FFT_out2[m2_prime];
	ZZ tw_FFT_out2_ZZ[m2_prime];	
	long long m2_prime_prou;
	ZZ m2_prime_prou_ZZ;	
	//cout << " 4 " << endl;	
	//m2_prime_prou = test.find_prou(m2_prime, modular_n);
	m2_prime_prou_ZZ = test.find_prou( m2_prime, modular_n_ZZ);
	cout << " m2_prime_prou_ZZ = " << m2_prime_prou_ZZ << endl;		
	//m2_prime_prou_ZZ = m2_prime_prou;
	
	//test.FFT(tw_FFT_out2, tw_FFT_in2, m2_prime, m2_prime_prou, modular_n) ;	
	test.DFT(tw_FFT_out2_ZZ, tw_FFT_in2_ZZ, m2_prime, m2_prime_prou_ZZ, modular_n_ZZ) ;	
	
	//cout << "tw_FFT_out2 = " << endl;
	for (int i = 0; i < m2_prime ; i++){
			//cout << tw_FFT_out2[i] <<endl;
			//cout << tw_FFT_out2_ZZ[i] <<endl;			
	}	
	//cout << endl;	
//----------	
	int m3_prime, gen3;
	long long rader_index_in3[s_m2];
	long long rader_index_out3[s_m2];
	gen3 = test.find_gen(s_m2);
	if( (s_m2==3) || (s_m2==5) || (s_m2==17) || (s_m2==257)){
		m3_prime = s_m2-1;
	}
	else{
		m3_prime = test.find_m_prime(s_m2);
	}
	long long tw_FFT_index3[m3_prime];
	long long tw_FFT_in3[m3_prime];
	ZZ tw_FFT_in3_ZZ[m3_prime];		
	long long m3_prou;	
	ZZ m3_prou_ZZ;
	
	//m3_prou = test.prou_power(prou_n, m1*s_m1, modular_n);
	//PowerMod(m3_prou_ZZ, (ZZ)prou_n, (ZZ)m1*s_m1, modular_n_ZZ);
	m3_prou_ZZ = test.find_prou( s_m2, (ZZ)modular_n);
	cout << " m3_prou_ZZ = " << m3_prou_ZZ << endl;		
	rader_index_in3[0] = 0;
	rader_index_out3[0] = 0;
	rader_index_in3[1] = 1;
	rader_index_out3[1] = 1;
	tmp = 1;
	for (int i = 2; i < s_m2 ; i++){
		tmp *= gen3;
		tmp %= s_m2;
		rader_index_in3[i] = tmp;	
	}
	for (int i = 2; i < s_m2 ; i++){
		rader_index_out3[i] = rader_index_in3[s_m2-i+1];
	}
	
	for (int i = 0; i < s_m2 ; i++){
		//cout << rader_index_out3[i] << endl;
	}	
	//cout << endl;
	
	
	for (int i = 0; i < m3_prime ; i++){
		if(i < (s_m2 - 1) )//0-5
			tw_FFT_index3[i] = rader_index_out3[i+1];
		else if(i > (m3_prime - s_m2 + 1))
			tw_FFT_index3[i] = rader_index_out3[i + s_m2 - m3_prime];			
		else		
			tw_FFT_index3[i] = 0; //00000
		//cout << tw_FFT_index3[i] << endl;	
	}	
	
	for (int i = 0; i < m3_prime ; i++){
		if(tw_FFT_index3[i] == 0){
			//tw_FFT_in3[i] = 0;
			tw_FFT_in3_ZZ[i] = 0;
		}
		else{
			//tw_FFT_in3[i] = test.prou_power(m3_prou,tw_FFT_index3[i],modular_n);
			PowerMod(tw_FFT_in3_ZZ[i], m3_prou_ZZ, (ZZ)tw_FFT_index3[i], (ZZ)modular_n);
		}
	}	
	
	
	
	long long tw_FFT_out3[m3_prime];
	ZZ tw_FFT_out3_ZZ[m3_prime];	
	long long m3_prime_prou;
	ZZ m3_prime_prou_ZZ;
	
	//m3_prime_prou = test.find_prou(m3_prime, modular_n);
	m3_prime_prou_ZZ = test.find_prou( m3_prime, modular_n_ZZ);
	//m3_prime_prou_ZZ = m3_prime_prou;
	cout << " m3_prime_prou_ZZ = " << m3_prime_prou_ZZ << endl;			
	//test.FFT(tw_FFT_out3, tw_FFT_in3, m3_prime, m3_prime_prou, modular_n) ;
	test.DFT(tw_FFT_out3_ZZ, tw_FFT_in3_ZZ, m3_prime, m3_prime_prou_ZZ, modular_n_ZZ) ;	
	
	for (int i = 0; i < m3_prime ; i++){
			//cout << tw_FFT_out3[i] <<endl;
			//cout << tw_FFT_out3_ZZ[i] <<endl;			
	}	
//-----------------------------------------------------------------------------
//cout << " 1 " << endl;

for(int time = 0; time < test_time; time++){
//---------I/O Re-Index-------------------
//------------------------------------------	
	for(int i = 0; i < phi_m; i++){	
		data_in[i] = 1/*(rand() % modular_n)*/ ;
		data_in_ZZ[i] =  data_in[i] ;
	}
	for(int i = 0; i < m; i++){	
		debug_data[i] = i ;
	}	
	
    long long n1,n2,n3,k1,k2,k3;
    long long index_in;
    long long index_out;
    long long index_m;
	long long Re_Index_Input_Data[m];	
	long long Re_Index_Output_Data[m];
	for(n3 = 0; n3 < s_m2; n3++)//257     17
	{ 
		for(n2 = 0; n2 < s_m1; n2++)//17   5
		{
			for(n1 = 0; n1 < m1; n1++)//5  3
			{
				index_m = n1 + n2 * m1 + n3 * m1*s_m1; // 1 2 3 ... m
				k1 = rader_index_in1[n1] ;//5
				k2 = rader_index_in2[n2] ;//17
				k3 = rader_index_in3[n3] ;//257
				index_in = (k1 * s_m1 * s_m2 + k2 * m1 * s_m2 + k3 * m1 * s_m1) % (m1 * s_m1 *s_m2);// 1
				
				//index_m = n1 + n2 * m1; // 1 2 3 ... m
				//k3 = (index_m / (m1*m2)) % m3 ;
				//k2 = (index_m / m1) % m2 ;
				//k1 = index_m % m1;
				//index_in = (k1 * m2 + k2 * m3*m1 + k3 * m2*m1) % (m1 * m2);// 1
				Re_Index_Input_Data[index_m] = data_in[index_in];	
				//cout<< DFT_data[index_m] << endl;
			}    
		}
	}


//------------------- put into memory ----------------------------
	for(int i = 0; i < 273 ; i++){
		if(i<257){
			Dual_port_mem_1w_273[i] = Re_Index_Input_Data[85*i];
			//Dual_port_mem_1w_273[i] = debug_data[85*i];		
		}
		else {
			Dual_port_mem_1w_273[i] = Re_Index_Input_Data[5*(i-256)];
			//Dual_port_mem_1w_273[i] = debug_data[5*(i-256)];
		}
		//cout << Dual_port_mem_1w_273[i] << endl;
	}

	for(int i = 0; i < 273 ; i++){
		for(int j = 0; j < 4 ; j++){
			if(i<257){
				Dual_port_mem_4w_273[i][j] = Re_Index_Input_Data[85*i+1+j];
				//Dual_port_mem_4w_273[i][j] = debug_data[85*i+1+j];		
			}
			else {
				Dual_port_mem_4w_273[i][j] = Re_Index_Input_Data[5*(i-256)+1+j];
				//Dual_port_mem_4w_273[i][j] = debug_data[5*(i-256)+1+j];
			}
			//cout << Dual_port_mem_4w_273[i][j] << " ";
		}	
		//cout << endl;
	}


	for(int k1 = 0; k1 < 256 ; k1++){
		for(int k2 = 0; k2 < 16 ; k2++){
			Dual_port_mem_1w_4096[16*k1+k2] = Re_Index_Input_Data[90+85*k1+5*k2];	
			//Dual_port_mem_1w_4096[16*k1+k2] = debug_data[90+85*k1+5*k2];	
			//cout << Dual_port_mem_1w_4096[16*k1+k2]<< endl;
		}	
	}

	for(int k1 = 0; k1 < 256 ; k1++){
		for(int k2 = 0; k2 < 16 ; k2++){
			for(int j = 0; j < 4 ; j++){
				Dual_port_mem_4w_4096[16*k1+k2][j] = Re_Index_Input_Data[90+85*k1+5*k2+1+j];
				//Dual_port_mem_4w_4096[16*k1+k2][j] = debug_data[90+85*k1+5*k2+1+j];	
				//cout << Dual_port_mem_4w_4096[16*k1+k2][j]<<" ";				
			}	
			//cout << endl;
		}	
	}


//cout << " 2 " << endl;

//------------------stage1 5-FFT---------------------------		
	//long long m1_prime_prou_inv = test.find_inv(m1_prime_prou,modular_n);
	ZZ m1_prime_prou_inv_ZZ = test.find_inv(m1_prime_prou_ZZ, (ZZ)modular_n);
	//m1_prime_prou_inv_ZZ = test.find_inv(m1_prime_prou_ZZ, (ZZ)modular_n);
	
	//cout << "m1_prime_prou_inv = " << m1_prime_prou_inv << endl;
	//cout << "m1_prime_prou_inv_ZZ = " << m1_prime_prou_inv_ZZ << endl;	
	
	//long long inv_4 = test.find_inv(4, modular_n);
	ZZ inv_4_ZZ = test.find_inv((ZZ)4, (ZZ)modular_n);	
	
	ZZ tw_1_4, tw_1_4_inv;
	tw_1_4 = m1_prime_prou_ZZ;
	tw_1_4_inv = m1_prime_prou_inv_ZZ;

/* cout << "5-point input = " << endl;
cout << Dual_port_mem_1w_273[1] << endl;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_273[1][i] << " ";
}
cout << endl ;  */

//cout << " mid FFT" << endl ; 

 	for(int i = 0; i < 273 ; i++){
		//first 4-FFT
		test.Radix_4_BU(Dual_port_mem_4w_273[i], Dual_port_mem_4w_273[i], tw_1_4, 4, (ZZ)modular_n);
/*  if(i ==1){
	for(int j = 0; j < 4; j++)
		cout << tw_FFT_out1[j] << " ";
} */
		for(int j = 0; j < 4 ; j++){ // point-wise mul tw
			//cout << Dual_port_mem_4w_273[0][j] << " ";
			MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], tw_FFT_out1_ZZ[j], (ZZ)modular_n );	
			//cout << Dual_port_mem_4w_273[0][j] << " ";
			//MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], inv_4, (ZZ)modular_n );	

		}
/*  if(i ==1){
	for(int j = 0; j < 4; j++)
		cout << Dual_port_mem_4w_273[i][j] << " ";
} */
		
		for(int j = 0; j < 4 ; j++){ // point-wise mul tw
			//cout << Dual_port_mem_4w_273[0][j] << " ";
			//MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], tw_FFT_out1[j], (ZZ)modular_n );	
			//cout << Dual_port_mem_4w_273[0][j] << " ";
			MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], inv_4_ZZ, (ZZ)modular_n );	

		}
	 				
		test.Radix_4_BU(Dual_port_mem_4w_273[i], Dual_port_mem_4w_273[i], tw_1_4_inv, 4, (ZZ)modular_n);
		
/*  if(i ==1){
	for(int j = 0; j < 4; j++)
		cout << Dual_port_mem_4w_273[i][j] << " ";
}	 */ 

		for(int j = 0; j < 4 ; j++){ // add d[0]
			AddMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], Dual_port_mem_1w_273[i], (ZZ)modular_n );	
		}
	} 
	
/* cout << endl ; 

cout << "5-point output = " << endl;
//cout << Dual_port_mem_1w_273[1] << endl;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_273[1][i] << endl;
}
cout << endl ;  */
	
 	for(int i = 257; i <= 260 ; i++){	// Last 16 relocation
		test.Relocation_4(Dual_port_mem_4w_273[i],Dual_port_mem_4w_273[i+4],Dual_port_mem_4w_273[i+8],Dual_port_mem_4w_273[i+12]);
	}
	
	for(int i = 0; i < 4096 ; i++){
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_4, 4, (ZZ)modular_n);
		for(int j = 0; j < 4 ; j++){ // point-wise mul tw
			//cout << Dual_port_mem_4w_273[0][j] << " ";
			MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], tw_FFT_out1_ZZ[j], (ZZ)modular_n );	
			//cout << Dual_port_mem_4w_273[0][j] << " ";
			MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], inv_4_ZZ, (ZZ)modular_n );				
		}
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_4_inv, 4, (ZZ)modular_n);
		for(int j = 0; j < 4 ; j++){ // add d[0]
			AddMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], Dual_port_mem_1w_4096[i], (ZZ)modular_n );	
		}			
	}	

 	for(int i = 0; i < 256 ; i++){	// relocation 256*16
	 	for(int j = 0; j <4  ; j++){
			test.Relocation_4(Dual_port_mem_4w_4096[16*i+j],Dual_port_mem_4w_4096[16*i+j+4],Dual_port_mem_4w_4096[16*i+j+8],Dual_port_mem_4w_4096[16*i+j+12]);
		}	
	}

//cout << " 3 " << endl;
//-----------------------------------------------End first stage----------------------------------------------//
//---------------/////////////////////-----------Second stage---------------///////////////////////-----------//
//17-RAFFT ---> 1 + 16-FFT ---> 2 stages radix-4 

//FFT    first stage radix4 ---> first relocation ---> second stage radix4 ---> point-wise mul 
//IFFT   first stage radix4 ---> first relocation ---> second stage radix4 ---> add d[0]
	//long long m2_prime_prou_inv = test.find_inv(m2_prime_prou,modular_n);
	ZZ m2_prime_prou_inv_ZZ = test.find_inv(m2_prime_prou_ZZ,(ZZ)modular_n);
	//long long inv_16 = test.find_inv(16, modular_n);
	ZZ inv_16_ZZ = test.find_inv((ZZ)16, (ZZ)modular_n);	
	
	ZZ tw_1_16, tw_1_16_inv;
	tw_1_16 = m2_prime_prou_ZZ;
	tw_1_16_inv = m2_prime_prou_inv_ZZ;
	vector<vector<ZZ>> tw_16p(4), tw_16p_inv(4);
	for(int i = 0; i < 4; i++){
		tw_16p[i].resize(4);
		tw_16p_inv[i].resize(4);		
	}	
	for(int i = 0; i < 4; i++){
		PowerMod(tw_16p[i][0], tw_1_16, 0  , (ZZ)modular_n);
		PowerMod(tw_16p[i][1], tw_1_16, i, (ZZ)modular_n);
		PowerMod(tw_16p[i][2], tw_1_16, 2*i  , (ZZ)modular_n);
		PowerMod(tw_16p[i][3], tw_1_16, 3*i, (ZZ)modular_n);		
	}	
	for(int i = 0; i < 4; i++){
		PowerMod(tw_16p_inv[i][0], tw_1_16_inv, 0  , (ZZ)modular_n);
		PowerMod(tw_16p_inv[i][1], tw_1_16_inv, i, (ZZ)modular_n);
		PowerMod(tw_16p_inv[i][2], tw_1_16_inv, 2*i  , (ZZ)modular_n);
		PowerMod(tw_16p_inv[i][3], tw_1_16_inv, 3*i, (ZZ)modular_n);		
	}		

/*   cout << "input = " << endl;

cout << Dual_port_mem_4w_273[0][1] << endl;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_273[261][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_273[262][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_273[263][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_273[264][i] << " ";
}
cout << endl ;  */ 

/* cout << "input = " << endl;

cout << Dual_port_mem_4w_273[256][3] << endl;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[4092][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[4093][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[4094][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[4095][i] << " ";
}
cout << endl ;  */

{
 	for(int i = 257; i < 273 ; i++){	// first stage radix4
		test.Radix_4_BU(Dual_port_mem_4w_273[i], Dual_port_mem_4w_273[i], tw_1_16, 16, (ZZ)modular_n);
	}

 	for(int i = 0; i < 4096 ; i++){	// first stage radix4
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_16, 16, (ZZ)modular_n);
	}
	
 	for(int i = 257; i < 273 ; i++){	// mul 16FFT middle tw 
		if((i - 257) % 4 == 0){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], tw_16p[0][j], (ZZ)modular_n);
			}
		}
		else if((i - 257) % 4 == 1){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], tw_16p[1][j], (ZZ)modular_n);
			}
		}
		else if((i - 257) % 4 == 2){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], tw_16p[2][j], (ZZ)modular_n);
			}
		}
		else if((i - 257) % 4 == 3){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], tw_16p[3][j], (ZZ)modular_n);
			}
		}		
	}	

 	for(int i = 0; i < 4096 ; i++){	// mul 16FFT middle tw 
		if((i) % 4 == 0){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], tw_16p[0][j], (ZZ)modular_n);
			}
		}
		else if((i) % 4 == 1){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], tw_16p[1][j], (ZZ)modular_n);
			}
		}
		else if((i) % 4 == 2){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], tw_16p[2][j], (ZZ)modular_n);
			}
		}
		else if((i) % 4 == 3){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], tw_16p[3][j], (ZZ)modular_n);
			}
		}		
	}

 	for(int i = 257; i < 273 ; i=i+4){	// first relocation
		test.Relocation_4(Dual_port_mem_4w_273[i],Dual_port_mem_4w_273[i+1],Dual_port_mem_4w_273[i+2],Dual_port_mem_4w_273[i+3]);
	}

 	for(int i = 0; i < 4096 ; i=i+4){	// first relocation
		test.Relocation_4(Dual_port_mem_4w_4096[i],Dual_port_mem_4w_4096[i+1],Dual_port_mem_4w_4096[i+2],Dual_port_mem_4w_4096[i+3]);
	}
	
	
 	for(int i = 257; i < 273 ; i++){	// second stage radix4
		test.Radix_4_BU(Dual_port_mem_4w_273[i], Dual_port_mem_4w_273[i], tw_1_16, 16, (ZZ)modular_n);
	}

 	for(int i = 0; i < 4096 ; i++){	// second stage radix4
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_16, 16, (ZZ)modular_n);
	}
	
	
/*   cout << "after FFT = " << endl;

for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[0][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[1][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[2][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[3][i] << " ";
}
cout << endl ;  */ 	
	
		
	//  point-wise mul 
 	for(int i = 257; i < 273 ; i=i+4){ // 4 group	and each 16 word
		for(int j = 0; j < 4; j++){  // 4 group and each 4 word
			for(int k = 0; k < 4; k++){  // 4 group and each 1 word		
				MulMod(Dual_port_mem_4w_273[i+j][k], Dual_port_mem_4w_273[i+j][k], tw_FFT_out2_ZZ[j+4*k], (ZZ)modular_n );
				//cout <<"["<<i+j<<"]"<<"["<<k<<"]"<<"-->"<<"["<<j+4*k<<"]"<< endl;
				MulMod(Dual_port_mem_4w_273[i+j][k], Dual_port_mem_4w_273[i+j][k], inv_16_ZZ, (ZZ)modular_n );
			}
		}
	}

	//  point-wise mul 
 	for(int i = 0; i < 1024 ; i++){ // 1024 group 
		for(int j = 0; j < 4; j++){  // each group 4 addr
			for(int k = 0; k < 4; k++){  // each addr 4 word	 0 4 8 12, 1 5 9 13, 2 6 10 14, 3 7 11 15
				MulMod(Dual_port_mem_4w_4096[4*i+j][k], Dual_port_mem_4w_4096[4*i+j][k], tw_FFT_out2_ZZ[j+4*k], (ZZ)modular_n );
				//cout <<"["<<i+j<<"]"<<"["<<k<<"]"<<"-->"<<"["<<j+4*k<<"]"<< endl;
				MulMod(Dual_port_mem_4w_4096[4*i+j][k], Dual_port_mem_4w_4096[4*i+j][k], inv_16_ZZ, (ZZ)modular_n );
			}
		}
	}
	
//IFFT
 	for(int i = 257; i < 273 ; i++){	// first stage radix4
		test.Radix_4_BU(Dual_port_mem_4w_273[i], Dual_port_mem_4w_273[i], tw_1_16_inv, 16, (ZZ)modular_n);
	}

 	for(int i = 0; i < 4096 ; i++){	// first stage radix4
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_16_inv, 16, (ZZ)modular_n);
	}

 	for(int i = 257; i < 273 ; i++){	// mul 16FFT middle tw 
		if((i - 257) % 4 == 0){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], tw_16p_inv[0][j], (ZZ)modular_n);
			}
		}
		else if((i - 257) % 4 == 1){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], tw_16p_inv[1][j], (ZZ)modular_n);
			}
		}
		else if((i - 257) % 4 == 2){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], tw_16p_inv[2][j], (ZZ)modular_n);
			}
		}
		else if((i - 257) % 4 == 3){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], tw_16p_inv[3][j], (ZZ)modular_n);
			}
		}		
	}

 	for(int i = 0; i < 4096 ; i++){	// mul 16FFT middle tw 
		if((i) % 4 == 0){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], tw_16p_inv[0][j], (ZZ)modular_n);
			}
		}
		else if((i) % 4 == 1){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], tw_16p_inv[1][j], (ZZ)modular_n);
			}
		}
		else if((i) % 4 == 2){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], tw_16p_inv[2][j], (ZZ)modular_n);
			}
		}
		else if((i) % 4 == 3){
			for(int j = 0; j < 4; j++){
				MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], tw_16p_inv[3][j], (ZZ)modular_n);
			}
		}		
	}

 	for(int i = 257; i < 273 ; i=i+4){	// first relocation
		test.Relocation_4(Dual_port_mem_4w_273[i],Dual_port_mem_4w_273[i+1],Dual_port_mem_4w_273[i+2],Dual_port_mem_4w_273[i+3]);
	}

 	for(int i = 0; i < 4096 ; i=i+4){	// first relocation
		test.Relocation_4(Dual_port_mem_4w_4096[i],Dual_port_mem_4w_4096[i+1],Dual_port_mem_4w_4096[i+2],Dual_port_mem_4w_4096[i+3]);
	}

 	for(int i = 257; i < 273 ; i++){	// second stage radix4
		test.Radix_4_BU(Dual_port_mem_4w_273[i], Dual_port_mem_4w_273[i], tw_1_16_inv, 16, (ZZ)modular_n);
	}

 	for(int i = 0; i < 4096 ; i++){	// second stage radix4
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_16_inv, 16, (ZZ)modular_n);
	}

 	for(int i = 257; i < 273 ; i++){
		for(int j = 0; j < 4 ; j++){	// add d0
			AddMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[0][(i-257)/4], (ZZ)modular_n );
		}
	}
	
 	for(int i = 0; i < 256 ; i++){   //4096 --> 256 group
		for(int j = 0; j < 16 ; j++){	// each group 16 addr
			for(int k = 0; k < 4 ; k++){	// each addr 4 word
				AddMod(Dual_port_mem_4w_4096[16*i+j][k], Dual_port_mem_4w_4096[16*i+j][k], Dual_port_mem_4w_273[i+1][j/4], (ZZ)modular_n );
			}
		}
	}	
}

/*  cout << "output = " << endl;

for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_273[261][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_273[262][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_273[263][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_273[264][i] << " ";
}
cout << endl ;  */

/* cout << "output = " << endl;

//cout << Dual_port_mem_4w_273[256][3] << endl;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[4092][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[4093][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[4094][i] << " ";
}
cout << endl ;
for(int i = 0; i < 4; i++){
	cout << Dual_port_mem_4w_4096[4095][i] << " ";
}
cout << endl ;  */


//---next stage relocation 17->257
// p11
 	for(int i = 0; i < 1024 ; i++){	// relocation 256*16
		test.Relocation_4(Dual_port_mem_4w_4096[i],Dual_port_mem_4w_4096[i+1024],Dual_port_mem_4w_4096[i+2048],Dual_port_mem_4w_4096[i+3072]);	
	}
	
	
//--------------------------------------------------------------------------------------------//
//test data
/* long long test_256_FFT_input[256]; 
long long test_256_FFT_output[256]; 
long long test_257_FFT_input[257]; 
long long test_257_FFT_output[257]; 
long long test_256_FFT_prou = test.find_prou(256, modular_n);
long long test_257_FFT_prou = test.find_prou(257, modular_n);
//Dual_port_mem_4w_273[257][0] = 0;
//cout << Dual_port_mem_4w_273[257][0] << endl;
for (int i = 0; i < 64 ; i++){
	for (int j = 0; j < 4 ; j++){
		//Dual_port_mem_4w_4096[16*i][j] = rader_index_in3[i + 64*j + 1];
		//Dual_port_mem_4w_4096[16*i][j] =  i + 64*j  ;
		//Dual_port_mem_4w_4096[16*i][j] =  i ;
		//cout << Dual_port_mem_4w_4096[16*i][j] << " ";
	}
	//cout << endl;
}

for (int i = 0; i < 257 ; i++){
	//for (int j = 0; j < 4 ; j++){
		test_257_FFT_input[i] =  i ; // 0 - 256
		//cout << test_257_FFT_input[i] << endl;
	//}
} 

//test.DFT(test_257_FFT_output, test_257_FFT_input, 257, test_257_FFT_prou , modular_n);
for (int i = 0; i < 257 ; i++){
	//for (int j = 0; j < 4 ; j++){
		//test_257_FFT_input[i] =  i ; // 0 - 256
		//cout << test_257_FFT_output[i] << endl;
	//}
}  */
	//cout << endl;
//cout << " 4 " << endl;
//--------------------------------------------------------------------------------------------//
//-----------------------------third stage 256 FFT -------------------------------------------// 
	//long long m3_prime_prou_inv = test.find_inv(m3_prime_prou,modular_n);
	ZZ m3_prime_prou_inv_ZZ = test.find_inv(m3_prime_prou_ZZ,(ZZ)modular_n);	
	
	//long long inv_256 = test.find_inv(256, modular_n);
	ZZ inv_256_ZZ = test.find_inv((ZZ)256, (ZZ)modular_n);	
	
	ZZ tw_1_256, tw_1_256_inv;	
	tw_1_256 = m3_prime_prou_ZZ;
	tw_1_256_inv = m3_prime_prou_inv_ZZ;
	
	vector<vector<ZZ>> tw_256p_s1(64), tw_256p_inv_s1(64);
	vector<vector<ZZ>> tw_256p_s2(16), tw_256p_inv_s2(16);
	vector<vector<ZZ>> tw_256p_s3(4), tw_256p_inv_s3(4);	
	
	for(int i = 0; i < 64; i++){
		tw_256p_s1[i].resize(4);
		tw_256p_inv_s1[i].resize(4);		
	}
	for(int i = 0; i < 16; i++){
		tw_256p_s2[i].resize(4);
		tw_256p_inv_s2[i].resize(4);		
	}
	for(int i = 0; i < 4; i++){
		tw_256p_s3[i].resize(4);
		tw_256p_inv_s3[i].resize(4);		
	}	
	for(int i = 0; i < 64; i++){
		PowerMod(tw_256p_s1[i][0], tw_1_256, 0  , (ZZ)modular_n);
		PowerMod(tw_256p_s1[i][1], tw_1_256, i, (ZZ)modular_n);
		PowerMod(tw_256p_s1[i][2], tw_1_256, 2*i  , (ZZ)modular_n);
		PowerMod(tw_256p_s1[i][3], tw_1_256, 3*i, (ZZ)modular_n);
		PowerMod(tw_256p_inv_s1[i][0], tw_1_256_inv, 0  , (ZZ)modular_n);
		PowerMod(tw_256p_inv_s1[i][1], tw_1_256_inv, i, (ZZ)modular_n);
		PowerMod(tw_256p_inv_s1[i][2], tw_1_256_inv, 2*i  , (ZZ)modular_n);
		PowerMod(tw_256p_inv_s1[i][3], tw_1_256_inv, 3*i, (ZZ)modular_n);
	}
	for(int i = 0; i < 16; i++){
		tw_256p_s2[i] = tw_256p_s1[4*i];
		tw_256p_inv_s2[i] = tw_256p_inv_s1[4*i];
	}
	for(int i = 0; i < 4; i++){
		tw_256p_s3[i] = tw_256p_s1[16*i];
		tw_256p_inv_s3[i] = tw_256p_inv_s1[16*i];
	}

// FFT
{
//cout << "stage1" << endl;
 	for(int i = 0; i < 4096 ; i++){	// first stage radix4 s1
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_256, 256, (ZZ)modular_n);
	}

	// mul s1_tw factor
	long long mem_addr;
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			for(int t = 0; t < 64 ; t++){  //each FFT 64 addr
				mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i + 16 * t; //64 to 64
				for(int k = 0; k < 4 ; k++){  //each addr 4 words
					MulMod(Dual_port_mem_4w_4096[mem_addr][k], Dual_port_mem_4w_4096[mem_addr][k], tw_256p_s1[t][k], (ZZ)modular_n);				
				}			
			}		
		}
	}

	//s1 relacation
	//excel 256-relocation		
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i ;
			for(int t1 = 0; t1 < 16 ; t1++){  //each FFT 64 addr 
				test.Relocation_4(Dual_port_mem_4w_4096[mem_addr + 16 * t1],Dual_port_mem_4w_4096[mem_addr + 16 * (t1+16)],Dual_port_mem_4w_4096[mem_addr + 16 * (t1+32)],Dual_port_mem_4w_4096[mem_addr + 16 * (t1+48)]);											
			}		
		}
	}
//cout << "stage2" << endl;
 	for(int i = 0; i < 4096 ; i++){	// first stage radix4 s2
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_256, 256, (ZZ)modular_n);
	}

	// mul s2_tw factor
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			for(int t1 = 0; t1 < 4 ; t1++){  //each FFT 64 addr
				for(int t2 = 0; t2 < 16 ; t2++){  //each FFT 64 addr
					mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i + 16 * (16 * t1 + t2); //64 to 16 * 4
					for(int k = 0; k < 4 ; k++){  //each addr 4 words
						MulMod(Dual_port_mem_4w_4096[mem_addr][k], Dual_port_mem_4w_4096[mem_addr][k], tw_256p_s2[t2][k], (ZZ)modular_n);				
					}			
				}				
			}		
		}
	}

	//s2 relacation
	//excel 256-relocation		
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i ; 		
			for(int t1 = 0; t1 < 4 ; t1++){  //each FFT 64 addr
				for(int t2 = 0; t2 < 4 ; t2++){  //each FFT 64 addr
					test.Relocation_4(Dual_port_mem_4w_4096[mem_addr + 16 * (16 * t1 + t2 + 0)],Dual_port_mem_4w_4096[mem_addr + 16 * (16 * t1 + t2 + 4)],Dual_port_mem_4w_4096[mem_addr + 16 * (16 * t1 + t2 + 8)],Dual_port_mem_4w_4096[mem_addr + 16 * (16 * t1 + t2 + 12)]);										
				}						
			}		
		}
	}
//cout << "stage3" << endl;
 	for(int i = 0; i < 4096 ; i++){	// first stage radix4 s3
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_256, 256, (ZZ)modular_n);
	}


	// mul s3_tw factor
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			for(int t1 = 0; t1 < 16 ; t1++){  //each FFT 64 addr
				for(int t2 = 0; t2 < 4 ; t2++){  //each FFT 64 addr
					mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i + 16 * (4 * t1 + t2); //64 to 4 * 16
					for(int k = 0; k < 4 ; k++){  //each addr 4 words
						MulMod(Dual_port_mem_4w_4096[mem_addr][k], Dual_port_mem_4w_4096[mem_addr][k], tw_256p_s3[t2][k], (ZZ)modular_n);				
					}			
				}				
			}		
		}
	}


	//s3 relacation
	//excel 256-relocation		
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i ; 		
			for(int t1 = 0; t1 < 16 ; t1++){  //each FFT 64 addr
				test.Relocation_4(Dual_port_mem_4w_4096[mem_addr + 16 * (4 * t1 + 0)],Dual_port_mem_4w_4096[mem_addr + 16 * (4 * t1 + 1)],Dual_port_mem_4w_4096[mem_addr + 16 * (4 * t1 + 2)],Dual_port_mem_4w_4096[mem_addr + 16 * (4 * t1 + 3)]);																
			}		
		}
	}
//cout << "stage4" << endl;
 	for(int i = 0; i < 4096 ; i++){	// first stage radix4 s4
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_256, 256, (ZZ)modular_n);
	}


}
//end FFT
//cout << " 5 " << endl;

//--------------------------------------------------------------------------------------------//
//test data
 for (int i = 0; i < 64 ; i++){
	for (int j = 0; j < 4 ; j++){
		//cout << Dual_port_mem_4w_4096[16*i][j] << " ";
	}
	//cout << endl;
}

//test.DFT(test_256_FFT_output, (long long *)test_256_FFT_input, 256, test_256_FFT_prou , modular_n);
for (int i = 0; i < 256 ; i++){
	//cout << test_256_FFT_output[i] << endl;
}
 
//--------------------------------------------------------------------------------------------//
	long long mem_addr;
//-------------------------------- 256-point-wise mul-----------------------------------------//
 	for(int i = 0; i < 4 ; i++){ // 4 group 
		for(int j = 0; j < 16; j++){  // each group 16 256-FFT
			mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i;
			for(int t1 = 0; t1 < 4; t1++){  // each addr 4 word	 0 64 128 192 16 80 144 208 
				for(int t2 = 0; t2 < 4; t2++){  // each addr 4 word	 0 64 128 192 16 80 144 208 
					for(int t3 = 0; t3 < 4; t3++){  // each addr 4 word	 0 64 128 192 16 80 144 208 
						for(int k = 0; k < 4; k++){  // each addr 4 word	 0 64 128 192 16 80 144 208 
							MulMod(Dual_port_mem_4w_4096[mem_addr+16*(16*t1+4*t2+t3)][k], Dual_port_mem_4w_4096[mem_addr+16*(16*t1+4*t2+t3)][k], tw_FFT_out3_ZZ[(t1+4*t2+16*t3)+64*k], (ZZ)modular_n );
							//cout <<"["<<16*t1+4*t2+t3<<"]"<<"["<<k<<"]"<<"-->"<<"["<<(t1+4*t2+16*t3)+64*k<<"]"<< endl;
						}
					}
				}
			}
		}
	}

//--------------------------------------------------------------------------------------------//

// IFFT
{
//cout << "stage1" << endl;
 	for(int i = 0; i < 4096 ; i++){	// first stage radix4 s1
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_256_inv, 256, (ZZ)modular_n);
	}

	// mul s1_tw factor

	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			for(int t = 0; t < 64 ; t++){  //each FFT 64 addr
				mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i + 16 * test.Bit_convert(t); //64 to 64
				for(int k = 0; k < 4 ; k++){  //each addr 4 words
					MulMod(Dual_port_mem_4w_4096[mem_addr][k], Dual_port_mem_4w_4096[mem_addr][k], tw_256p_inv_s1[t][k], (ZZ)modular_n);				
				}			
			}		
		}
	}

	//s1 relacation
	//excel 256-relocation		
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i ;
			for(int t1 = 0; t1 < 16 ; t1++){  //each FFT 64 addr 
				test.Relocation_4(Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(t1)],Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(t1+16)],Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(t1+32)],Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(t1+48)]);											
			}		
		}
	}
//cout << "stage2" << endl;
 	for(int i = 0; i < 4096 ; i++){	// first stage radix4 s2
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_256_inv, 256, (ZZ)modular_n);
	}

	// mul s2_tw factor
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			for(int t1 = 0; t1 < 4 ; t1++){  //each FFT 64 addr
				for(int t2 = 0; t2 < 16 ; t2++){  //each FFT 64 addr
					mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i + 16 * test.Bit_convert(16 * t1 + t2); //64 to 16 * 4
					for(int k = 0; k < 4 ; k++){  //each addr 4 words
						MulMod(Dual_port_mem_4w_4096[mem_addr][k], Dual_port_mem_4w_4096[mem_addr][k], tw_256p_inv_s2[t2][k], (ZZ)modular_n);				
					}			
				}				
			}		
		}
	}

	//s2 relacation
	//excel 256-relocation		
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i ; 		
			for(int t1 = 0; t1 < 4 ; t1++){  //each FFT 64 addr
				for(int t2 = 0; t2 < 4 ; t2++){  //each FFT 64 addr
					test.Relocation_4(Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(16 * t1 + t2 + 0)],Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(16 * t1 + t2 + 4)],Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(16 * t1 + t2 + 8)],Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(16 * t1 + t2 + 12)]);										
				}						
			}		
		}
	}
//cout << "stage3" << endl;
 	for(int i = 0; i < 4096 ; i++){	// first stage radix4 s3
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_256_inv, 256, (ZZ)modular_n);
	}


	// mul s3_tw factor
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			for(int t1 = 0; t1 < 16 ; t1++){  //each FFT 64 addr
				for(int t2 = 0; t2 < 4 ; t2++){  //each FFT 64 addr
					mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i + 16 * test.Bit_convert(4 * t1 + t2); //64 to 4 * 16
					for(int k = 0; k < 4 ; k++){  //each addr 4 words
						MulMod(Dual_port_mem_4w_4096[mem_addr][k], Dual_port_mem_4w_4096[mem_addr][k], tw_256p_inv_s3[t2][k], (ZZ)modular_n);				
					}			
				}				
			}		
		}
	}


	//s3 relacation
	//excel 256-relocation		
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i ; 		
			for(int t1 = 0; t1 < 16 ; t1++){  //each FFT 64 addr
				test.Relocation_4(Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(4 * t1 + 0)],Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(4 * t1 + 1)],Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(4 * t1 + 2)],Dual_port_mem_4w_4096[mem_addr + 16 * test.Bit_convert(4 * t1 + 3)]);																
			}		
		}
	}
//cout << "stage4" << endl;
 	for(int i = 0; i < 4096 ; i++){	// first stage radix4 s4
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_1_256_inv, 256, (ZZ)modular_n);
		for(int k = 0; k < 4; k++){
			MulMod(Dual_port_mem_4w_4096[i][k], Dual_port_mem_4w_4096[i][k], inv_256_ZZ, (ZZ)modular_n );
		}
	}
}
//end IFFT
//cout << "end IFFT" << endl;
//--------------------------add d[0]---------------------------------------------//
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i ; 		
			for(int t1 = 0; t1 < 64 ; t1++){  //each FFT 64 addr
				for(int k = 0; k < 4 ; k++){  //each 4 words
					AddMod(Dual_port_mem_4w_4096[mem_addr + 16*(t1) ][k], Dual_port_mem_4w_4096[mem_addr + 16*(t1)][k], Dual_port_mem_4w_273[257+4*i+(j%4)][(j/4)], (ZZ)modular_n );														
				}																														
			}		
		}
	}

//-------------------------------------------------------------------------------//
//cout << "end add d[0]" << endl;




//--------------------output array----------------------------------//
  std::ofstream o_o_16384("./output_16384.txt");
  std::ofstream o_o("./output_21845.txt");
  std::ofstream g_o("./golden.txt");  



ZZ output_array_tmp[16384];
ZZ output_array_tmp2[m];
ZZ output_array[m];
ZZ golden_output[m];
ZZ golden_output2[m];
	for(int i = 0; i < 4 ; i++){  // 4 group
		for(int j = 0; j < 16 ; j++){  //each group 16 256-FFT
			mem_addr = 1024 * (j / 4) + (j % 4) + 4 * i ; 		
			for(int t1 = 0; t1 < 4 ; t1++){  //each FFT 64 addr
				for(int k = 0; k < 64 ; k++){  //each 4 words
					output_array_tmp[4096 * i + 256 * j + 64 * t1 + k] = Dual_port_mem_4w_4096[mem_addr + 16*k][t1];												
				}																														
			}		
		}
	}

for(int i = 1; i < 5; i++){
	for(int j = 1; j < 17 ; j++){
		for(int k = 1; k < 257 ; k++){
			output_array_tmp2[4369*i + 257*j + k] = output_array_tmp[4096*(i-1) + 256*(j-1) + (k-1)];
		}		
	}
}




	for(n1 = 0; n1 < m1; n1++) //5    3  
	{ 
		for(n2 = 0; n2 < s_m1; n2++)//17    5
		{
			for(n3 = 0; n3 < s_m2; n3++)//257     7
			{
				index_m = n3 + n2 * s_m2 + n1 * s_m2*s_m1; // 1 2 3 ... m
				k1 = rader_index_out1[n1] ;//5
				k2 = rader_index_out2[n2] ;//17
				k3 = rader_index_out3[n3] ;//257
				index_out = (k1 * m2_inv * m2 + (k2 * s_m2 * s_m2_inv + k3 * s_m1 * s_m1_inv)*m1*m1_inv ) % (m1 * m2);
                output_array[index_out] = output_array_tmp2[index_m];	
				
				//index_m = n1 + n2 * m1; // 1 2 3 ... m
				//k3 = (index_m / (m1*m2)) % m3 ;
				//k2 = (index_m / m1) % m2 ;
				//k1 = index_m % m1;
				//index_in = (k1 * m2 + k2 * m3*m1 + k3 * m2*m1) % (m1 * m2);// 1
				//DFT_data[index_m] = data_in[index_in];	
				//test << k1 << endl;
			}    
		}
	}





long long test_21845_FFT_output[m];
ZZ test_21845_FFT_output_ZZ[m];
test.DFT(test_21845_FFT_output, data_in, 21845, prou_n , modular_n);
test.DFT(test_21845_FFT_output_ZZ, data_in_ZZ, 21845, prou_n_ZZ , modular_n_ZZ);


for( int i = 0; i < 15; i++){
	//cout << " test_21845_FFT_output = " << test_21845_FFT_output[i] << endl;
	//cout << " test_21845_FFT_output_ZZ = " << test_21845_FFT_output_ZZ[i] << endl;	
}


long long zmstar[phi_m];
test.find_zmstar(zmstar, m);
for(int i = 0; i < phi_m; i++){
	golden_output[zmstar[i]] = test_21845_FFT_output_ZZ[zmstar[i]];
}

	for(int i = 0; i < m; i++){
		o_o << output_array[i] << endl;
		g_o << golden_output[i] << endl;
		//cout << output_array[i]<< endl;
		//o_o << golden_output2[i] << endl;
		
	}


	for(int i = 0; i < 16384; i++){
		o_o_16384 << output_array_tmp[i] << endl;
		//g_o << test_21845_FFT_output[i] << endl;
		//cout << output_array[i]<< endl;
		
	}


	int k = 0;
	for (int i = 0; i< m ; i++){	
		if(output_array[i] != golden_output[i])	{
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


	//test.PFA3_v4(PFA_data_out, data_in, m1,m2,s_m1,s_m2,m1_inv,m2_inv,s_m1_inv,s_m2_inv, prou_n, modular_n); 











  for (int i = 0; i < 64 ; i++){
	for (int j = 0; j < 4 ; j++){
		//cout << Dual_port_mem_4w_4096[16*i][j] << " ";
	}
	//cout << endl;
} 

ZZ tmp_256[256];
ZZ ans_256[256];
 for (int i = 0; i < 4 ; i++){
	for (int j = 0; j < 64 ; j++){
		//cout << Dual_port_mem_4w_4096[16*j][i] << endl;
		//tmp_256[64 * i + j] =  Dual_port_mem_4w_4096[16*j][i] ;
		//cout << tmp_256[64 * i + j] << endl;
	}
	//cout << endl;
}

for (int j = 0; j < 256 ; j++){
	//cout << Dual_port_mem_4w_4096[16*j][i] << endl;
	//ans_256[rader_index_out3[j+1]-1] =  tmp_256[j] ;
}


for (int j = 0; j < 256 ; j++){
	//cout << Dual_port_mem_4w_4096[16*j][i] << endl;
	//cout << ans_256[j] - test_257_FFT_output[j+1] << endl;
}



//---------------/////////////////////-----------End Second stage---------------///////////////////////-----------//
//cout << "golden = " << endl;
//-----golden rader-----
//long long rader_in[17] = {14577791,10629372,6680953,11561170,16649673,13251367,10360485,12714571,4365733,12793856,10103857,643192,16340040,15243222,11238477,11237192,11235907};
//long long rader_out[17];
//long long rader_prou = test.find_prou(17, modular_n);
// long long rader_in[5] = {85 ,4454 ,8823 ,13192,0};
// long long rader_out[5];
// long long rader_prou = test.find_prou(5, modular_n);
//test.Rader(rader_out, rader_in, 17, rader_prou , modular_n);
//cout << " prou_n = " <<  rader_prou << endl;	
//test.DFT(rader_out, rader_in, 17, rader_prou , modular_n);
 	//for(int i = 0; i < 17 ; i++){
		//cout << rader_out[i] << " ";
	//}
	

//long long in_16[16] = {10629372 ,11561170 ,12793856 ,10103857,15243222 ,13251367 ,11237192 ,643192   ,11235907 ,11238477 ,4365733  ,12714571 ,16649673,16340040,6680953,10360485};
//test.DFT(rader_out, in_16, 16, m2_prime_prou , modular_n);
 	//for(int i = 0; i < 16 ; i++){
		//cout << rader_out[i] << " ";
	//}



/* 2-D vector example 
	vector<vector<ZZ>> Dual_port_mem(25);
	vector<ZZ> four_word_mem(4);

	for(int i = 0; i < 25; i++){	
		Dual_port_mem[i].resize(4);
	}
	
	for(int i = 0; i < 100; i++){
		four_word_mem[i%4] = i;
		if(i%4 == 3){
			Dual_port_mem[i/4] = four_word_mem;
		}
	}
	
	for(int i = 0; i < 25; i++){
		for(int j = 0; j < 4; j++){
			cout << Dual_port_mem[i][j] << " "; 
		}
		cout << endl;
	}
*/
/*
	vector<ZZ> input;
	input.resize(4);
	input[0] = 1;
	input[1] = 3;
	input[2] = 5;
	input[3] = 7;	
	vector<ZZ> Output;	
	Output.resize(4);
	vector<ZZ> tw;	
	tw.resize(4);	
    PowerMod(tw[0], (ZZ)prou_n, 0, (ZZ)modular_n);	
    PowerMod(tw[1], (ZZ)prou_n, 1, (ZZ)modular_n);
    PowerMod(tw[2], (ZZ)prou_n, 2, (ZZ)modular_n);
    PowerMod(tw[3], (ZZ)prou_n, 3, (ZZ)modular_n);	
	test.Radix_4_BU(Output, input, tw, (ZZ)modular_n);
	
	
	cout << endl;		
	for (int i = 0; i < m ; i++){	
		//cout << Output[i] <<endl;
	}
	
	cout <<endl;	
*/	












}

}