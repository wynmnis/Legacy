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
	ZZ m ;
	m = 1925;
    LEGACY precompute_gen(m);
	ZZ prou_LCM;
	//long long LCM = 4369*1155*256;
	//prou_LCM = (ZZ)2168490747132;
	
	//LCM = 4369*1155*256
	//notice that the primitive root of unity exist
	//when m | modular - 1 (modular -1 can be divided by m)	
	vector<ZZ> a(6) ;
	vector<ZZ> out(100) ;
	a[0] = 2;
	a[1] = 4;
	a[2] = 6;
	a[3] = 8;
	a[4] = 10;
	a[5] = 12;	
	int cnt ;
	//precompute_gen.C(6,3,a,0);
	//precompute_gen.Combin2(a,5,3,3);
	//cout << precompute_gen.C(7,4) << endl;
	//cnt = precompute_gen.allprint(6, 3, a, out);
	//cnt = precompute_gen.find_AllProduct(out, a, 6);
	//cout << "cnt = " << cnt <<endl;
	
	for(int i = 0; i < cnt; i++){
		//cout << out[i] << endl;
	}
	
	//cnt = precompute_gen.Factorize_no_power(out, m*16);
	
	
	for(int i = 0; i < cnt; i++){
		//cout << out[i] << endl;
	}
	//cout << "check " << precompute_gen.check_prou((ZZ)1433201, 105, precompute_gen.modular_m)<< endl;
	
	
	int p = precompute_gen.ZZ2int(m) ;
	vector<ZZ>   RA_in4(p);
	vector<ZZ>  RA_out4(p);	
	vector<ZZ> DFT_out4(p);		
	 for(int i = 0; i < p; i++){
		RA_in4[i] = i+1 ;
	 }	
	ZZ prou_p, inv1, inv2;
	vector<ZZ> factor_p(10);
	prou_p = PowerMod(precompute_gen.prou_LCM, precompute_gen.LCM_/p, precompute_gen.modular_m);
	//prou_p = PowerMod(prou_LCM, LCM/p, precompute_gen.modular_m);
	precompute_gen.Factorize_2(factor_p, (ZZ)p);
	cout << factor_p[0] << " , " << factor_p[1] << endl;
	inv1 = precompute_gen.find_inv_exgcd(factor_p[0], factor_p[1]);
	inv2 = precompute_gen.find_inv_exgcd(factor_p[1], factor_p[0]);
	//precompute_gen.RA_powerof2_FFT(RA_out4, RA_in4, (ZZ)8, 0);
	//precompute_gen.PFA_10P_FFT(RA_out4,RA_in4, 0);
	//precompute_gen.RA_11P_FFT(RA_out4,RA_in4, 0);
	//precompute_gen.PFA_FFT(RA_out4,RA_in4, (ZZ)2, (ZZ)15, (ZZ)8, (ZZ)1,precompute_gen.find_prou(p, precompute_gen.modular_m), 0);
	precompute_gen.PFA_FFT(RA_out4,RA_in4, factor_p[0], factor_p[1], inv1, inv2, prou_p , 0);
	for(int i = 0; i < p; i++){
		//cout << RA_out4[i] << endl;
	}	
	
	//cout << " prou = " << precompute_gen.find_prou(p, precompute_gen.modular_m) << endl;
	precompute_gen.DFT(DFT_out4, RA_in4, p, prou_p, precompute_gen.modular_m);
	for(int i = 0; i < p; i++){
		//cout << DFT_out4[i] << endl;
	}		

	for(int i = 0; i < p; i++){
		if(RA_out4[i] !=  DFT_out4[i]){
			cout << " DFT fail " << endl;
			break;
		}
		else if(i == p - 1)
			cout << " DFT done " << endl;
	}		

	//precompute_gen.RA_powerof2_FFT(DFT_out4, RA_out4, (ZZ)8, 1);
	//precompute_gen.PFA_10P_FFT(DFT_out4,RA_out4, 1);	
	//precompute_gen.RA_11P_FFT(DFT_out4,RA_out4, 1);
	precompute_gen.PFA_FFT(DFT_out4,RA_out4, factor_p[0], factor_p[1], inv1, inv2, prou_p , 1);	
	for(int i = 0; i < p; i++){
		//cout << DFT_out4[i] << endl;
	}	

	for(int i = 0; i < p; i++){
		if(DFT_out4[i] !=  RA_in4[i]){
			cout << " IDFT fail " << endl;
			break;
		}
		else if(i == p - 1)
			cout << " IDFT done " << endl;
	}	

	

}