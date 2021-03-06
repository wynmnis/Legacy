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
	ZZ m;
	ZZ modular_m;
	ZZ prou_m;
	int first_cnt;
	ZZ LCM_ ;
	ZZ prou_LCM;
	ZZ prou_295680, prou_1155;
	vector<ZZ> second_cnt;
	vector<vector<ZZ>> second_decompose;
	vector<ZZ> first_decompose;
	vector<vector<ZZ>> input_index_first_decompose;
	vector<vector<ZZ>> output_index_first_decompose;
	vector<vector<vector<ZZ>>> input_index_second_decompose;
	vector<vector<vector<ZZ>>> output_index_second_decompose;
	vector<vector<ZZ>> first_decompose_precompute_data;
	vector<vector<vector<ZZ>>> second_decompose_precompute_data;
	ZZ prou_3p, prou_5p, prou_7p, prou_11p;
	ZZ prou_6, prou_10;	
	ZZ prou_2, prou_4, prou_8, prou_16, prou_32, prou_64, prou_128, prou_256;
	ZZ prou_9, prou_27;
	ZZ prou_2_inv, prou_4_inv, prou_8_inv, prou_16_inv, prou_32_inv, prou_64_inv, prou_128_inv, prou_256_inv, prou_9_inv, prou_27_inv;	
	ZZ inv_3, inv_5, inv_7, inv_11;
	ZZ inv_2, inv_4, inv_6, inv_10;
	ZZ inv_8, inv_16, inv_32, inv_64, inv_128, inv_256, inv_9, inv_27;
	vector<ZZ> RA_3P_precompute;
	vector<ZZ> RA_5P_precompute;
	vector<ZZ> RA_7P_precompute;
	vector<ZZ> RA_11P_precompute;

	vector<ZZ> RA_3P_inv_precompute;
	vector<ZZ> RA_5P_inv_precompute;
	vector<ZZ> RA_7P_inv_precompute;
	vector<ZZ> RA_11P_inv_precompute;
	
    LEGACY(ZZ m_);	 // constructor
	long ZZ2int(ZZ n)	;
	bool coprime(long long a, long long b);
	long long Euler(long long data_in);
	long long find_m_prime(long long m);
	ZZ find_m_prime(ZZ m)	;
	long long Prime(long long i);	
	ZZ Prime(ZZ i)	;
	long long Factorize(long long *factor ,long long num);
	long long Factorize_no_power(vector<ZZ> &factor, ZZ num);
	long long Factorize(ZZ *factor, ZZ num);	
	long long Factorize(vector<ZZ> &factor, ZZ num);	
	long long Factorize_fine(vector<ZZ> &first_decompose, vector<vector<ZZ>> &second_decompose, vector<ZZ> &second_cnt ,ZZ num)	;
	void Factorize_2(vector<ZZ> &factor, ZZ num);
	long long find_gen(long long n);
	//long long find_gen(ZZ n);	
	ZZ find_gen(ZZ n);		
	long long find_inv(long long data_in, long long modular);
	ZZ find_inv(ZZ data_in, ZZ modular);
	ZZ exgcd(ZZ a, ZZ b, ZZ &x, ZZ &y);
	long long exgcd(long long a, long long b, long long &x, long long &y);	
	ZZ LCM(ZZ a, ZZ b);
	long long LCM(long long a, long long b);	
	ZZ LCM(vector<ZZ> &a);
	long long LCM(vector<long long> &a);	
	ZZ find_inv_exgcd(ZZ a, ZZ m) ;		
	bool isPowerBy2(long long n);
	bool isPowerBy2(ZZ n);
	bool isPowerBy3(ZZ n);
	bool isPowerBy5(ZZ n);
	bool isPowerBy7(ZZ n);	
	bool isPrime(long long n);
	bool isPrime(ZZ n);	
	long long find_prime(long long m, long long powerof2);
    long long find_prime_conditional(long long m, long long lowerbound)	;
	
	ZZ find_prime(ZZ m, long long powerof2);	
	long long find_prou(long long m, long long modular);
	ZZ find_n_rou(ZZ base, long long m, ZZ modular)	;
	bool check_prou(ZZ n_rou, long long m, ZZ modular)	;
	bool check_prou2(vector<ZZ> &factor, long long cnt, ZZ n_rou, long long m, ZZ modular)	;	
	ZZ find_prou(long long m, ZZ modular);
	ZZ find_prou2(long long m, ZZ modular);	
	void C(int n,int r, int a[], int product_out[], int m);
	int C(int n, int r);
	int Combin2(int A[], int m, int n, int r);
	int allprint(int iAllLen, int n, vector<ZZ> &a, vector<ZZ> &out)	;
	long long find_AllProduct(vector<ZZ> &product, vector<ZZ> &factor, long long cnt);	
	void find_zmstar(long long *zmstar, long long m);	
	long long prou_power(long long data_in, long long power, long long modular);
	long long DFT(long long *DFT_data, long long *data_in, long long m, long long prou, long long modular);
	void DFT(ZZ *DFT_data, ZZ *data_in, long long m, ZZ prou, ZZ modular);
	void DFT(vector<ZZ> &DFT_data, vector<ZZ> &data_in, long long m, ZZ prou, ZZ modular);	
	void Rader(long long *RA_out, long long *data_in, long long n, long long prou, long long modular);
	void Rader_DFT(long long *RA_out, long long *data_in, long long n, long long prou, long long modular);
	void Rader_DFT_0118(long long *RA_out, long long *data_in, long long n, long long *precompute_data, long long *In_Reindex, long long *Out_Reindex, long long prou_n_neg_1, long long modular);	
	void Rader_DFT_optimize(long long *RA_out, long long *data_in, long long n, long long prou, long long modular);	
	void Rader_v3(long long *RA_out, long long *data_in, long long *tw_FFT_out, long long *index_in, long long *index_out, long long m, long long m_prime, long long m_prime_prou, long long modular)	;
	void Rader_v4(long long *RA_out, long long *data_in, long long *tw_FFT_out, long long *index_in, long long *index_out, long long m, long long m_prime, long long m_prime_prou, long long modular)	;
	void Rader_v4_optimize(long long *RA_out, long long *data_in, long long *tw_FFT_out, long long *index_in, long long *index_out, long long m, long long m_prime, long long m_prime_prou, long long modular)	;	
	long long IDFT(long long *IDFT_data, long long *data_in, long long n, long long prou, long long modular);	
	void IDFT(ZZ *IDFT_data, ZZ *data_in, long long n, ZZ prou, ZZ modular);	
	void IDFT(vector<ZZ> &IDFT_data, vector<ZZ> &data_in, long long n, ZZ prou, ZZ modular);	
	long long FFT(long long *DFT_data, long long *data_in, long long n, long long prou, long long modular);
	//void FFT(ZZ *DFT_data, ZZ *data_in, ZZ n, ZZ prou, ZZ modular);
	long long IFFT(long long *DFT_data, long long *data_in, long long n, long long prou, long long modular);
    long long PFA2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long inv1, long long inv2, long long prou, long long modular);	
	long long PFA2_optimize(long long *DFT_data, long long *data_in, long long m1, long long m2, long long inv1, long long inv2, long long prou, long long modular);		
	long long PFA3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) ;
	//long long PFA3_optimize(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) ;	
	long long PFA2_v2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long inv1, long long inv2, long long prou, long long modular);	
	long long PFA2_v3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long m1_prime, long long m2_prime, long long m1_prime_prou, long long m2_prime_prou, long long inv1, long long inv2, long long *tw_FFT_out1, long long *tw_FFT_out2, long long *index_out1, long long *index_out2, long long modular,long long m1_prou,long long m2_prou,long long *index_in1,long long *index_in2);
	long long PFA2_v4(long long *DFT_data, long long *data_in, long long m1, long long m2, long long m1_prime, long long m2_prime, long long m1_prime_prou, long long m2_prime_prou, long long inv1, long long inv2, long long *tw_FFT_out1, long long *tw_FFT_out2, long long *index_out1, long long *index_out2, long long modular,long long m1_prou,long long m2_prou,long long *index_in1,long long *index_in2);
	long long PFA2_v4_optimize(long long *DFT_data, long long *data_in, long long m1, long long m2, long long m1_prime, long long m2_prime, long long m1_prime_prou, long long m2_prime_prou, long long inv1, long long inv2, long long *tw_FFT_out1, long long *tw_FFT_out2, long long *index_out1, long long *index_out2, long long modular,long long m1_prou,long long m2_prou,long long *index_in1,long long *index_in2);	
	long long PFA3_v2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) ;
	//long long PFA3_v3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long m3, long long inv1, long long inv2, long long inv3, long long prou, long long modular) ;
	long long PFA3_v3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular);	
	long long PFA3_v4(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular);
	long long PFA3_v4_optimize(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular);	
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
	void RA_2P_FFT(vector<ZZ> &out, vector<ZZ> &in, bool inverse);
	void RA_3P_FFT(vector<ZZ> &out, vector<ZZ> &in, bool inverse);	
	void pointwise_mul(vector<ZZ> &out, vector<ZZ> &in1, vector<ZZ> &in2);
	void RA_5P_FFT(vector<ZZ> &out, vector<ZZ> &in, bool inverse);
	void RA_powerof2_FFT(vector<ZZ> &out, vector<ZZ> &in, ZZ point, bool inverse);
	void RA_powerof3_FFT(vector<ZZ> &out, vector<ZZ> &in, ZZ point, bool inverse);	
	void PFA_6P_FFT(vector<ZZ> &out, vector<ZZ> &in, bool inverse);
	void RA_7P_FFT(vector<ZZ> &out, vector<ZZ> &in, bool inverse);
	void PFA_10P_FFT(vector<ZZ> &out, vector<ZZ> &in, bool inverse);	
	void RA_11P_FFT(vector<ZZ> &out, vector<ZZ> &in, bool inverse);
	void PFA_FFT(vector<ZZ> &out, vector<ZZ> &in, ZZ m1, ZZ m2, ZZ m1_inv, ZZ m2_inv, ZZ prou, bool inverse);
	void Winograd_small_FFT(int N, vector<ZZ> &FFT_out, vector<ZZ> &FFT_in, ZZ m, ZZ prou_m, ZZ modular, bool inverse);
};

#endif
