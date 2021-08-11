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
/* 	long long amount;
	long long factor[20] = {0};
    amount = test.Factorize(105);
	cout << endl; */


	long long m = 7;
	long long modular = test.find_prime(m*6,6);	
	long long RA_out[m];
	long long RA_in[m];	
	long long DFT_out[m];
	long long DFT_in[m];	
	long long inv_6 = test.find_inv(6, modular);
	long long inv_2 = test.find_inv(2, modular);	
	long long prou_6 = test.find_prou(6,modular);	
	
	cout << "prou_6 = " << endl;
	for(int i = 0; i < m; i++){
		RA_in[i] = i;
		DFT_in[i] = i;
	}
	
	ZZ prou_7_1 = test.find_prou(m, (ZZ)modular);
	cout << "prou_7_1 = " << prou_7_1 << endl;
	//ZZ prou_6 = test.find_prou(6, (ZZ)modular);	
	ZZ prou_3 = test.find_prou(3, (ZZ)modular);
	ZZ prou_2 = test.find_prou(2, (ZZ)modular);	
	ZZ precompute_3p[2];
	ZZ precompute_3p_inv[2];	
	ZZ precompute_7p[6];
	ZZ prou_7[6], prou_7_inv[6];
	ZZ prou_tw_out[6];
	ZZ prou_3_2;
	PowerMod(prou_3_2, (ZZ)prou_3, 2,  (ZZ)modular);
	precompute_3p[0] = AddMod(prou_3, prou_3_2, (ZZ)modular);
	precompute_3p[1] = SubMod(prou_3, prou_3_2, (ZZ)modular);	
	precompute_3p_inv[0] = AddMod(prou_3_2, prou_3, (ZZ)modular);
	precompute_3p_inv[1] = SubMod(prou_3_2, prou_3, (ZZ)modular);	

	


	PowerMod(prou_7[0], (ZZ)prou_7_1, 1,  (ZZ)modular);
	PowerMod(prou_7[1], (ZZ)prou_7_1, 5,  (ZZ)modular);
	PowerMod(prou_7[2], (ZZ)prou_7_1, 4,  (ZZ)modular);
	PowerMod(prou_7[3], (ZZ)prou_7_1, 6,  (ZZ)modular);
	PowerMod(prou_7[4], (ZZ)prou_7_1, 2,  (ZZ)modular);
	PowerMod(prou_7[5], (ZZ)prou_7_1, 3,  (ZZ)modular);
	
	test.DFT(prou_tw_out, prou_7, 6, (ZZ)prou_6 , (ZZ)modular);	
	

	PowerMod(prou_7_inv[0], (ZZ)prou_7_1, -1,  (ZZ)modular);
	PowerMod(prou_7_inv[1], (ZZ)prou_7_1, -3,  (ZZ)modular);
	PowerMod(prou_7_inv[2], (ZZ)prou_7_1, -2,  (ZZ)modular);
	PowerMod(prou_7_inv[3], (ZZ)prou_7_1, -6,  (ZZ)modular);
	PowerMod(prou_7_inv[4], (ZZ)prou_7_1, -4,  (ZZ)modular);
	PowerMod(prou_7_inv[5], (ZZ)prou_7_1, -5,  (ZZ)modular);
	
	//Rader input reindex
	ZZ testin[7] = {(ZZ)0,(ZZ)1,(ZZ)6,(ZZ)2,(ZZ)5,(ZZ)4,(ZZ)3};
	vector<ZZ> BU_in1(2);
	vector<ZZ> BU_out1(2);	
	vector<ZZ> BU_in2(2);
	vector<ZZ> BU_out2(2);	
	vector<ZZ> BU_in3(2);
	vector<ZZ> BU_out3(2);	
	//PFA input reindex
	BU_in1[0] = testin[1];
	BU_in1[1] = testin[2];	
	BU_in2[0] = testin[3];
	BU_in2[1] = testin[4];
	BU_in3[0] = testin[5];
	BU_in3[1] = testin[6];	
	test.Radix_2_BU(BU_out1, BU_in1,(ZZ)modular);
	test.Radix_2_BU(BU_out2, BU_in2,(ZZ)modular);	
	test.Radix_2_BU(BU_out3, BU_in3,(ZZ)modular);	
	
	
	vector<ZZ> BU_in4(2);
	vector<ZZ> BU_out4(2);	
	vector<ZZ> BU_in5(2);
	vector<ZZ> BU_out5(2);	
	BU_in4[0] = BU_out2[0];
	BU_in4[1] = BU_out3[0];	
	BU_in5[0] = BU_out2[1];
	BU_in5[1] = BU_out3[1];			
	
	test.Radix_2_BU(BU_out4, BU_in4, (ZZ)modular);
	test.Radix_2_BU(BU_out5, BU_in5, (ZZ)modular);	

	MulMod(BU_out4[0], BU_out4[0], precompute_3p[0], (ZZ)modular);
	MulMod(BU_out4[1], BU_out4[1], precompute_3p[1], (ZZ)modular);	
	MulMod(BU_out5[0], BU_out5[0], precompute_3p[0], (ZZ)modular);
	MulMod(BU_out5[1], BU_out5[1], precompute_3p[1], (ZZ)modular);


	vector<ZZ> BU_out6(2);	
	vector<ZZ> BU_out7(2);	
	
	test.Radix_2_BU(BU_out6, BU_out4,(ZZ)modular);
	test.Radix_2_BU(BU_out7, BU_out5,(ZZ)modular);	

	MulMod(BU_out6[0], BU_out6[0], inv_2, (ZZ)modular);
	MulMod(BU_out6[1], BU_out6[1], inv_2, (ZZ)modular);	
	MulMod(BU_out7[0], BU_out7[0], inv_2, (ZZ)modular);
	MulMod(BU_out7[1], BU_out7[1], inv_2, (ZZ)modular);


	
	ZZ tmp[6];
	// PFA output reindex
	AddMod(tmp[0], (BU_in4[0]+BU_in4[1]), BU_out1[0], (ZZ)modular);
	AddMod(tmp[4], BU_out6[0], BU_out1[0], (ZZ)modular);
	AddMod(tmp[2], BU_out6[1], BU_out1[0], (ZZ)modular);
	AddMod(tmp[3], (BU_in5[0]+BU_in5[1]), BU_out1[1], (ZZ)modular);
	AddMod(tmp[1], BU_out7[0], BU_out1[1], (ZZ)modular);
	AddMod(tmp[5], BU_out7[1], BU_out1[1], (ZZ)modular);	






	for(int i = 0; i < 6; i++){
		cout << " tmp = " << tmp[i] << endl;	
	}
	
	ZZ golden_in6[6]={(ZZ)1,(ZZ)3,(ZZ)2,(ZZ)6,(ZZ)4,(ZZ)5};
	ZZ golden_out6[6];

	test.DFT(golden_out6, golden_in6, 6, (ZZ)prou_6 , (ZZ)modular);	

	for(int i = 0; i < 6; i++){
		cout << " golden_out6 = " << golden_out6[i] << endl;	
	}

	//long long golden_in6_2[6]={1,3,2,6,4,5};
	ZZ golden_out6_2[6];
	//long long prou_6 = test.find_prou(6,modular);
	
	MulMod(golden_out6[0], golden_out6[0], prou_tw_out[0], (ZZ)modular);
	MulMod(golden_out6[1], golden_out6[1], prou_tw_out[1], (ZZ)modular);	
	MulMod(golden_out6[2], golden_out6[2], prou_tw_out[2], (ZZ)modular);
	MulMod(golden_out6[3], golden_out6[3], prou_tw_out[3], (ZZ)modular);	
	MulMod(golden_out6[4], golden_out6[4], prou_tw_out[4], (ZZ)modular);
	MulMod(golden_out6[5], golden_out6[5], prou_tw_out[5], (ZZ)modular);		
	
	
		
	
	test.IDFT(golden_out6_2 , golden_out6, 6, (ZZ)prou_6 , (ZZ)modular);	

	for(int i = 0; i < 6; i++){
		cout << " golden_out6_2 = " << golden_out6_2[i] << endl;	
	}




//--------------------------------------------------------------------------------------------//


	MulMod(tmp[0], tmp[0], prou_tw_out[0], (ZZ)modular);
	MulMod(tmp[1], tmp[1], prou_tw_out[1], (ZZ)modular);	
	MulMod(tmp[2], tmp[2], prou_tw_out[2], (ZZ)modular);
	MulMod(tmp[3], tmp[3], prou_tw_out[3], (ZZ)modular);	
	MulMod(tmp[4], tmp[4], prou_tw_out[4], (ZZ)modular);
	MulMod(tmp[5], tmp[5], prou_tw_out[5], (ZZ)modular);	
	//PFA IFFT input index
	BU_in1[0] =tmp[0];
	BU_in1[1] =tmp[3];	
	BU_in2[0] =tmp[2];
	BU_in2[1] =tmp[5];
	BU_in3[0] =tmp[4];
	BU_in3[1] =tmp[1];	
	
	for(int i = 0; i < 6; i++){
		//cout << " golden_out6 = " << golden_out6[i] << endl;	
	}	
//---------------------------------//

	test.Radix_2_BU(BU_out1, BU_in1,(ZZ)modular);
	test.Radix_2_BU(BU_out2, BU_in2,(ZZ)modular);	
	test.Radix_2_BU(BU_out3, BU_in3,(ZZ)modular);		
	
	//cout << test.find_inv(53, 97) << endl;
	
	BU_in4[0] = BU_out2[0];
	BU_in4[1] = BU_out3[0];	
	BU_in5[0] = BU_out2[1];
	BU_in5[1] = BU_out3[1];			
	
	test.Radix_2_BU(BU_out4, BU_in4, (ZZ)modular);
	test.Radix_2_BU(BU_out5, BU_in5, (ZZ)modular);	

	MulMod(BU_out4[0], BU_out4[0], precompute_3p_inv[0], (ZZ)modular);
	MulMod(BU_out4[1], BU_out4[1], precompute_3p_inv[1], (ZZ)modular);	
	MulMod(BU_out5[0], BU_out5[0], precompute_3p_inv[0], (ZZ)modular);
	MulMod(BU_out5[1], BU_out5[1], precompute_3p_inv[1], (ZZ)modular);	
	
	test.Radix_2_BU(BU_out6, BU_out4, (ZZ)modular);
	test.Radix_2_BU(BU_out7, BU_out5, (ZZ)modular);		

	MulMod(BU_out6[0], BU_out6[0], inv_2, (ZZ)modular);
	MulMod(BU_out6[1], BU_out6[1], inv_2, (ZZ)modular);	
	MulMod(BU_out7[0], BU_out7[0], inv_2, (ZZ)modular);
	MulMod(BU_out7[1], BU_out7[1], inv_2, (ZZ)modular);

	//PFA IDFT output index
	AddMod(tmp[0], (BU_in4[0]+BU_in4[1]), BU_out1[0], (ZZ)modular);
	AddMod(tmp[1], BU_out6[0], BU_out1[0], (ZZ)modular);
	AddMod(tmp[2], BU_out6[1], BU_out1[0], (ZZ)modular);
	AddMod(tmp[3], (BU_in5[0]+BU_in5[1]), BU_out1[1], (ZZ)modular);
	AddMod(tmp[4], BU_out7[0], BU_out1[1], (ZZ)modular);
	AddMod(tmp[5], BU_out7[1], BU_out1[1], (ZZ)modular);	

//IFFT inverse 1/N	           
	MulMod(tmp[0], tmp[0], inv_6, (ZZ)modular);
	MulMod(tmp[1], tmp[1], inv_6, (ZZ)modular);	
	MulMod(tmp[2], tmp[2], inv_6, (ZZ)modular);
	MulMod(tmp[3], tmp[3], inv_6, (ZZ)modular);	
	MulMod(tmp[4], tmp[4], inv_6, (ZZ)modular);
	MulMod(tmp[5], tmp[5], inv_6, (ZZ)modular);		
	
	
	vector<ZZ> output(7);		
	//rader output reindex or mix
	AddMod(output[0], (testin[1]+testin[2]+testin[3]+testin[4]+testin[5]+testin[6]), testin[0], (ZZ)modular);
	AddMod(output[1], tmp[0], testin[0], (ZZ)modular);
	AddMod(output[2], tmp[1], testin[0], (ZZ)modular);
	AddMod(output[4], tmp[2], testin[0], (ZZ)modular);
	AddMod(output[6], tmp[3], testin[0], (ZZ)modular);
	AddMod(output[5], tmp[4], testin[0], (ZZ)modular);	
	AddMod(output[3], tmp[5], testin[0], (ZZ)modular);	

	long long golden_in[7]={0,1,2,3,4,5,6};
	long long golden_out[7];
	
	//cout << " prou_7_1 = " << prou_7_1 << endl;
	
	for(int i = 0; i < 7; i++){
		cout << " output = " << output[i] << endl;	
	}
	
	cout <<endl;
	
	test.DFT(golden_out, golden_in, 7, 1454 , modular);
	
	for(int i = 0; i < 7; i++){
		cout << " golden_out = " << golden_out[i] << endl;	
	}

	
	return 0;
}