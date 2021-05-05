#ifndef _LEGACY_H_
#define _LEGACY_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace NTL;
using namespace std;
class LEGACY{
public:	
    void init();	
	long long find_inv(long long data_in, long long modular);
	long long find_prou(long long m, long long modular);
	long long prou_power(long long data_in, long long power, long long modular);
	long long DFT(long long *DFT_data, long long *data_in, long long m, long long prou, long long modular);
	//long long IDFT(long long *DFT_data, long long *data_in, long long m, long long prou, long long modular);	
	long long FFT(long long *DFT_data, long long *data_in, long long n, long long prou, long long modular);
	long long IFFT(long long *DFT_data, long long *data_in, long long n, long long prou, long long modular);
    long long PFA2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long inv1, long long inv2, long long prou, long long modular);	
	long long PFA3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) ;
	long long cyc_DFT(long long *DFT_data, long long *data_in, long long m, long long prou, long long modular);
	long long cyc_PFA2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long inv1, long long inv2, long long prou, long long modular);
	long long cyc_PFA3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular);
	int Gray(int index,int group);
	int RR(int BC, int shift_bit, int Bit_WIDTH);
	int unary_xor(int data_in, int Bit_WIDTH);
	void int2vec(int integer, int Bit_WIDTH, vector<int> &bit_array);
	int vec2int(vector<int> &bit_array, int Bit_WIDTH);
	
	
};

#endif
