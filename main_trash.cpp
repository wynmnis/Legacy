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
	
	for(int i = 0; i < m; i++){
		RA_in[i] = i;
		DFT_in[i] = i;
	}
	
	ZZ prou_7 = test.find_prou(m, (ZZ)modular);
	ZZ prou_6 = test.find_prou(6, (ZZ)modular);	
	ZZ prou_3 = test.find_prou(3, (ZZ)modular);
	ZZ prou_2 = test.find_prou(2, (ZZ)modular);	
	ZZ precompute_3p[2];
	ZZ precompute_7p[6];
	ZZ prou_7[6];
	ZZ prou_3_2;
	PowerMod(prou_3_2, (ZZ)prou_3, 2,  (ZZ)modular);
	precompute_3p[0] = AddMod(prou_3, prou_3_2, (ZZ)modular);
	precompute_3p[1] = SubMod(prou_3, prou_3_2, (ZZ)modular);	
	

	PowerMod(prou_7[0], (ZZ)prou_7, 1,  (ZZ)modular);
	PowerMod(prou_7[1], (ZZ)prou_7, 3,  (ZZ)modular);
	PowerMod(prou_7[2], (ZZ)prou_7, 2,  (ZZ)modular);
	PowerMod(prou_7[3], (ZZ)prou_7, 6,  (ZZ)modular);
	PowerMod(prou_7[4], (ZZ)prou_7, 4,  (ZZ)modular);
	PowerMod(prou_7[5], (ZZ)prou_7, 5,  (ZZ)modular);

	
	ZZ testin[6] = {1,3,2,6,4,5};
	vector<ZZ> BU_in1(2);
	vector<ZZ> BU_out1(2);	
	vector<ZZ> BU_in1(2);
	vector<ZZ> BU_out1(2);	
	vector<ZZ> BU_in1(2);
	vector<ZZ> BU_out1(2);		
	BU_in1[0] = 1;
	BU_in1[1] = 3;	
	BU_in2[0] = 2;
	BU_in2[1] = 6;
	BU_in3[0] = 4;
	BU_in3[1] = 5;	
	test.Radix_2_BU(BU_out1, BU_in1, modular);
	test.Radix_2_BU(BU_out2, BU_in2, modular);	
	test.Radix_2_BU(BU_out3, BU_in3, modular);	
	
	
	vector<ZZ> BU_in4(2);
	vector<ZZ> BU_out4(2);	
	vector<ZZ> BU_in5(2);
	vector<ZZ> BU_out5(2);	
	BU_in4[0] = BU_out2[0];
	BU_in4[1] = BU_out3[0];	
	BU_in5[0] = BU_out2[1];
	BU_in5[1] = BU_out3[1];			
	
	test.Radix_2_BU(BU_out4, BU_in4, modular);
	test.Radix_2_BU(BU_out5, BU_in5, modular);	

	MulMod(BU_out4[0], BU_out4[0], precompute_3p[0], (ZZ)modular);
	MulMod(BU_out4[1], BU_out4[1], precompute_3p[1], (ZZ)modular);	
	MulMod(BU_out5[0], BU_out5[0], precompute_3p[0], (ZZ)modular);
	MulMod(BU_out5[1], BU_out5[1], precompute_3p[1], (ZZ)modular);


	vector<ZZ> BU_out6(2);	
	vector<ZZ> BU_out7(2);	
	
	test.Radix_2_BU(BU_out6, BU_out4, modular);
	test.Radix_2_BU(BU_out7, BU_out5, modular);	


	
	ZZ tmp[6];
	AddMod(tmp[0], (BU_in4[0]+BU_in4[1]), BU_out1[0], (ZZ)modular);
	AddMod(tmp[1], BU_out6[0], BU_out1[0], (ZZ)modular);
	AddMod(tmp[2], BU_out6[1], BU_out1[0], (ZZ)modular);
	AddMod(tmp[3], (BU_in5[0]+BU_in5[1]), BU_out1[1], (ZZ)modular);
	AddMod(tmp[4], BU_out7[0], BU_out1[1], (ZZ)modular);
	AddMod(tmp[5], BU_out7[1], BU_out1[1], (ZZ)modular);	

	MulMod(tmp[0], tmp[0], prou_7[0], (ZZ)modular);
	MulMod(tmp[1], tmp[1], prou_7[1], (ZZ)modular);	
	MulMod(tmp[2], tmp[2], prou_7[2], (ZZ)modular);
	MulMod(tmp[3], tmp[3], prou_7[3], (ZZ)modular);	
	MulMod(tmp[4], tmp[4], prou_7[4], (ZZ)modular);
	MulMod(tmp[5], tmp[5], prou_7[5], (ZZ)modular);		
	
	//cout << test.find_inv(53, 97) << endl;
	
	
	
	
	
	
	
	

	return 0;
}