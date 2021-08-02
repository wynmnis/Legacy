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
	long ZZ2int(ZZ n)	;
	bool coprime(long long a, long long b);
	long long Euler(long long data_in);
	long long find_m_prime(long long m);
	ZZ find_m_prime(ZZ m)	;
	long long Prime(long long i);	
	ZZ Prime(ZZ i)	;
	long long Factorize(long long *factor ,long long num);
	long long Factorize(ZZ *factor, ZZ num);
	long long find_gen(long long n);
	long long find_gen(ZZ n);	
	long long find_inv(long long data_in, long long modular);
	ZZ find_inv(ZZ data_in, ZZ modular);
	ZZ exgcd(ZZ a, ZZ b, ZZ &x, ZZ &y);
	ZZ find_inv_exgcd(ZZ a, ZZ m) ;		
	bool isPowerBy2(long long n);
	bool isPowerBy2(ZZ n)	;
	bool isPrime(long long n);
	bool isPrime(ZZ n);	
	long long find_prime(long long m, long long powerof2);
	long long find_prou(long long m, long long modular);
	ZZ find_n_rou(ZZ base, long long m, ZZ modular)	;
	bool check_prou(ZZ n_rou, long long m, ZZ modular)	;
	ZZ find_prou(long long m, ZZ modular);
	void find_zmstar(long long *zmstar, long long m);	
	long long prou_power(long long data_in, long long power, long long modular);
	long long DFT(long long *DFT_data, long long *data_in, long long m, long long prou, long long modular);
	void DFT(ZZ *DFT_data, ZZ *data_in, long long m, ZZ prou, ZZ modular);
	void Rader(long long *RA_out, long long *data_in, long long n, long long prou, long long modular);
	void Rader_DFT(long long *RA_out, long long *data_in, long long n, long long prou, long long modular);
	void Rader_v3(long long *RA_out, long long *data_in, long long *tw_FFT_out, long long *index_in, long long *index_out, long long m, long long m_prime, long long m_prime_prou, long long modular)	;
	void Rader_v4(long long *RA_out, long long *data_in, long long *tw_FFT_out, long long *index_in, long long *index_out, long long m, long long m_prime, long long m_prime_prou, long long modular)	;
	long long IDFT(long long *IDFT_data, long long *data_in, long long n, long long prou, long long modular);	
	long long FFT(long long *DFT_data, long long *data_in, long long n, long long prou, long long modular);
	//void FFT(ZZ *DFT_data, ZZ *data_in, ZZ n, ZZ prou, ZZ modular);
	long long IFFT(long long *DFT_data, long long *data_in, long long n, long long prou, long long modular);
    long long PFA2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long inv1, long long inv2, long long prou, long long modular);	
	long long PFA3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) ;
	long long PFA2_v2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long inv1, long long inv2, long long prou, long long modular);	
	long long PFA2_v3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long m1_prime, long long m2_prime, long long m1_prime_prou, long long m2_prime_prou, long long inv1, long long inv2, long long *tw_FFT_out1, long long *tw_FFT_out2, long long *index_out1, long long *index_out2, long long modular,long long m1_prou,long long m2_prou,long long *index_in1,long long *index_in2);
	long long PFA2_v4(long long *DFT_data, long long *data_in, long long m1, long long m2, long long m1_prime, long long m2_prime, long long m1_prime_prou, long long m2_prime_prou, long long inv1, long long inv2, long long *tw_FFT_out1, long long *tw_FFT_out2, long long *index_out1, long long *index_out2, long long modular,long long m1_prou,long long m2_prou,long long *index_in1,long long *index_in2);
	long long PFA3_v2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) ;
	//long long PFA3_v3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long m3, long long inv1, long long inv2, long long inv3, long long prou, long long modular) ;
	long long PFA3_v3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular);	
	long long PFA3_v4(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular);	
	long long cyc_DFT(long long *DFT_data, long long *data_in, long long m, long long prou, long long modular);
	long long cyc_PFA2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long inv1, long long inv2, long long prou, long long modular);
	long long cyc_PFA3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular);
	int Gray(int index,int group);
	int RR(int BC, int shift_bit, int Bit_WIDTH);
	int bit_reverse(int num, int Bit_WIDTH)	;
	int Bit_convert(int addr);	
	int unary_xor(int data_in, int Bit_WIDTH);
	void int2vec(int integer, int Bit_WIDTH, vector<int> &bit_array);
	int vec2int(vector<int> &bit_array, int Bit_WIDTH);
	void Radix_4_BU(vector<ZZ> &output, vector<ZZ> &input, ZZ tw_1_N, long N, ZZ modular);	
	void Relocation_4(vector<ZZ> &v0, vector<ZZ> &v1, vector<ZZ> &v2, vector<ZZ> &v3);
	void Relocation_2(vector<ZZ> &v0, vector<ZZ> &v1);
	void Radix_2_BU(vector<ZZ> &output, vector<ZZ> &input, ZZ modular);
	void Config_PFA_Rader_FFT(vector<ZZ> &output, vector<ZZ> &input, ZZ m, ZZ prou_m, ZZ modular);
	ZZ expand_point_RA(ZZ m);
	void FFT_1024_radix2(vector<ZZ> &output, vector<ZZ> &input, int point, ZZ modular);	
	void FFT_1024_radix2_config(vector<ZZ> &output, vector<ZZ> &input, int point, ZZ m, ZZ prou_m, ZZ modular);
	void Rader_precompute_data(vector<ZZ> &precompute_data, int n, ZZ m, ZZ prou_m, ZZ modular);
};

#endif
