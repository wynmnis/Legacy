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
using namespace NTL;
int main()
{
	LEGACY test;
	//notice that the primitive root of unity exist
	//when m | modular - 1 (modular -1 can be divided by m)	

//---------setting-----------	
	long long m = 21845;	
	long long m1 = 5;	
	long long m2 = 4369;	
	long long s_m1 = 17;	
	long long s_m2 = 257;		
	long long modular_n = test.find_prime(m,8);	
		
	long long prou_n;	
	long long DFT_data_out[m];
	long long DFT_data_out2[m];	
	long long C_K[m];	
	long long DFT_data[m];	
	long long PFA_data_out[m];	
	long long IDFT_data_out[m];
	long long error[m];
	long long data_in[m]={0};
	long long debug_data[m];
	long long data_in2[m];
	long long data_tmp[m];
	long long m1_inv,m2_inv,s_m1_inv,s_m2_inv;
	long long phi_m ;
	prou_n = test.find_prou( m, modular_n);
	m1_inv = test.find_inv(m1, m2);
	m2_inv = test.find_inv(m2, m1);
	s_m1_inv = test.find_inv(s_m1, s_m2);
	s_m2_inv = test.find_inv(s_m2, s_m1);
	
	prou_n = test.find_prou( m, modular_n);
	phi_m = test.Euler(m);
	//------------------------------------------	
	for(int i = 0; i < phi_m; i++){	
		data_in[i] = i ;
	}
	for(int i = 0; i < m; i++){	
		debug_data[i] = i ;
	}	
	
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
	long long m1_prou;	
	
	m1_prou = test.prou_power(prou_n, m2, modular_n);
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
		if(tw_FFT_index1[i] == 0)
			tw_FFT_in1[i] = 0;
		else
			tw_FFT_in1[i] = test.prou_power(m1_prou,tw_FFT_index1[i],modular_n);
	}
	//cout << "prou = " << prou << endl;
	//cout << "m1_prou = " << m1_prou << endl;	
	//cout << "tw_FFT_in1 = " << endl;
	/*for (int i = 0; i < m1_prime ; i++){
			//cout << tw_FFT_in1[i] <<endl;
	}*/
	
	long long tw_FFT_out1[m1_prime];
	long long m1_prime_prou;

	m1_prime_prou = test.find_prou(m1_prime, modular_n);
	test.FFT(tw_FFT_out1, tw_FFT_in1, m1_prime, m1_prime_prou, modular_n) ;	
	//cout << "tw_FFT_out1 = " << endl;
	/*for (int i = 0; i < m1_prime ; i++){
			//cout << tw_FFT_out1[i] <<endl;
	}*/	
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
	
	m2_prou = test.prou_power(prou_n, m1*s_m2, modular_n);	
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
		if(tw_FFT_index2[i] == 0)
			tw_FFT_in2[i] = 0;
		else
			tw_FFT_in2[i] = test.prou_power(m2_prou,tw_FFT_index2[i],modular_n);
	}
	long long tw_FFT_out2[m2_prime];
	long long m2_prime_prou;
	m2_prime_prou = test.find_prou(m2_prime, modular_n);
	test.FFT(tw_FFT_out2, tw_FFT_in2, m2_prime, m2_prime_prou, modular_n) ;	
	
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
	long long m3_prou;	
	
	m3_prou = test.prou_power(prou_n, m1*s_m1, modular_n);	
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
		if(tw_FFT_index3[i] == 0)
			tw_FFT_in3[i] = 0;
		else
			tw_FFT_in3[i] = test.prou_power(m3_prou,tw_FFT_index3[i],modular_n);
	}	
	long long tw_FFT_out3[m3_prime];
	long long m3_prime_prou;
	m3_prime_prou = test.find_prou(m3_prime, modular_n);
	test.FFT(tw_FFT_out3, tw_FFT_in3, m3_prime, m3_prime_prou, modular_n) ;
//-----------------------------------------------------------------------------


//---------I/O Re-Index-------------------
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
                Re_Index_Output_Data[index_out] = data_tmp[index_m];	
				
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

//------------------stage1 5-FFT---------------------------
	vector<ZZ> tw_4;	
	tw_4.resize(4);	
    PowerMod(tw_4[0], (ZZ)m1_prime_prou, 0, (ZZ)modular_n);	
    PowerMod(tw_4[1], (ZZ)m1_prime_prou, 1, (ZZ)modular_n);
    PowerMod(tw_4[2], (ZZ)m1_prime_prou, 2, (ZZ)modular_n);
    PowerMod(tw_4[3], (ZZ)m1_prime_prou, 3, (ZZ)modular_n);	
	vector<ZZ> tw_4_inv;	
	long long m1_prime_prou_inv = test.find_inv(m1_prime_prou,modular_n);
	tw_4_inv.resize(4);	
    PowerMod(tw_4_inv[0], (ZZ)m1_prime_prou_inv, 0, (ZZ)modular_n);	
    PowerMod(tw_4_inv[1], (ZZ)m1_prime_prou_inv, 1, (ZZ)modular_n);
    PowerMod(tw_4_inv[2], (ZZ)m1_prime_prou_inv, 2, (ZZ)modular_n);
    PowerMod(tw_4_inv[3], (ZZ)m1_prime_prou_inv, 3, (ZZ)modular_n);	
	long long inv_4 = test.find_inv(4, modular_n);
	//cout << " modulus = " <<  modular_n << endl;
	//cout << " prou_n = " <<  m1_prime_prou << endl;	

	
	cout << endl ; 
	
 	for(int i = 0; i < 273 ; i++){
		//first 4-FFT
		test.Radix_4_BU(Dual_port_mem_4w_273[i], Dual_port_mem_4w_273[i], tw_4, (ZZ)modular_n);
		for(int j = 0; j < 4 ; j++){ // point-wise mul tw
			//cout << Dual_port_mem_4w_273[0][j] << " ";
			MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], tw_FFT_out1[j], (ZZ)modular_n );	
			//cout << Dual_port_mem_4w_273[0][j] << " ";
			MulMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], inv_4, (ZZ)modular_n );				
		}
		test.Radix_4_BU(Dual_port_mem_4w_273[i], Dual_port_mem_4w_273[i], tw_4_inv, (ZZ)modular_n);
		for(int j = 0; j < 4 ; j++){ // add d[0]
			AddMod(Dual_port_mem_4w_273[i][j], Dual_port_mem_4w_273[i][j], Dual_port_mem_1w_273[i], (ZZ)modular_n );	
		}
	} 
	
 	for(int i = 257; i <= 260 ; i++){	// Last 16 relocation
		test.Relocation_4(Dual_port_mem_4w_273[i],Dual_port_mem_4w_273[i+4],Dual_port_mem_4w_273[i+8],Dual_port_mem_4w_273[i+12]);
	}
	
	for(int i = 0; i < 4096 ; i++){
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_4, (ZZ)modular_n);
		for(int j = 0; j < 4 ; j++){ // point-wise mul tw
			//cout << Dual_port_mem_4w_273[0][j] << " ";
			MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], tw_FFT_out1[j], (ZZ)modular_n );	
			//cout << Dual_port_mem_4w_273[0][j] << " ";
			MulMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], inv_4, (ZZ)modular_n );				
		}
		test.Radix_4_BU(Dual_port_mem_4w_4096[i], Dual_port_mem_4w_4096[i], tw_4_inv, (ZZ)modular_n);
		for(int j = 0; j < 4 ; j++){ // add d[0]
			AddMod(Dual_port_mem_4w_4096[i][j], Dual_port_mem_4w_4096[i][j], Dual_port_mem_1w_4096[i], (ZZ)modular_n );	
		}			
	}	

 	for(int i = 0; i < 256 ; i++){	// relocation 256*16
	 	for(int j = 0; j <4  ; j++){
			test.Relocation_4(Dual_port_mem_4w_4096[16*i+j],Dual_port_mem_4w_4096[16*i+j+4],Dual_port_mem_4w_4096[16*i+j+8],Dual_port_mem_4w_4096[16*i+j+12]);
		}	
	}
//-----------------------------------------------end first stage----------------------------------------------//



//-----golden rader-----
long long rader_in[5] = {0,1,2,4,3};
long long rader_out[5];
long long rader_prou = test.find_prou(5, modular_n);
//test.Rader(rader_out, rader_in, 5, m1_prou , modular_n);
//cout << " prou_n = " <<  rader_prou << endl;	
test.DFT(rader_out, rader_in, 5, rader_prou , modular_n);
 	for(int i = 0; i < 5 ; i++){
		//cout << rader_out[i] << " ";
	}



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