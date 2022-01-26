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
	m = 2821;
    LEGACY winograd(m);
	long long point = 7;
	
	if(point == 7){
		//  7-WFFT
		long long a[point] ;
		long long DFT_out[point] ;
		long long Winograd_out[point] ;
		long long out(100) ;
		long long modular = 421;
		a[0] = 2;
		a[1] = 14;
		a[2] = 9;
		a[3] = 8;
		a[4] = 13;
		a[5] = 50;
		a[6] = 130;		
		long long prou7;
		long long inv_2, inv_3;
		inv_2 = winograd.find_inv(2, modular);
		inv_3 = winograd.find_inv(3, modular);		
		prou7 = winograd.find_prou(7, modular);
		winograd.DFT(DFT_out, a, 7, prou7, modular);
		long long w1,w2,w3,w4,w5,w6;

		w1 = PowerMod(prou7,1,modular);
		w2 = PowerMod(prou7,2,modular);
		w3 = PowerMod(prou7,3,modular);
		w4 = PowerMod(prou7,4,modular);
		w5 = PowerMod(prou7,5,modular);
		w6 = PowerMod(prou7,6,modular);
		
		long long cosu, cos2u, cos3u, isinu, isin2u, isin3u;
		cosu =  MulMod(AddMod(w1,w6,  modular), inv_2,  modular);
		isinu = MulMod(SubMod(w1,w6, modular), inv_2, modular);
		cos2u = MulMod(AddMod(w2,w5, modular), inv_2, modular);
		isin2u =MulMod(SubMod(w2,w5,modular), inv_2,modular);	
		cos3u = MulMod(AddMod(w3,w4, modular), inv_2, modular);
		isin3u =MulMod(SubMod(w3,w4,modular), inv_2,modular);		

		long long s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31,s32,s33,s34,s35,s36;
		long long m0,m1,m2,m3,m4,m5,m6,m7;
		long long c0,c1,c2,c3,c4,c5,c6,c7;		

		c0 = SubMod(MulMod(AddMod(AddMod(cosu , cos2u, modular), cos3u, modular), inv_3, modular), 1, modular);
		c1 = MulMod(SubMod(SubMod(MulMod(cosu,2,modular) , cos2u, modular), cos3u, modular), inv_3, modular);	
		c3 = MulMod(AddMod(SubMod(cosu , MulMod(cos2u,2,modular), modular), cos3u, modular), inv_3, modular);	
		c2 = MulMod(SubMod(AddMod(cosu , cos2u, modular), MulMod(cos3u,2,modular), modular), inv_3, modular);	
		c4 = MulMod(SubMod(AddMod(isinu , isin2u, modular), isin3u, modular), inv_3, modular);
		c5 = MulMod(AddMod(SubMod(MulMod(isinu,2,modular) , isin2u, modular), isin3u, modular), inv_3, modular);
		c7 = MulMod(SubMod(SubMod(isinu , MulMod(isin2u,2,modular), modular), isin3u, modular), inv_3, modular);
		c6 = MulMod(AddMod(AddMod(isinu , isin2u, modular), MulMod(isin3u,2,modular), modular), inv_3, modular);		
		
		cout << "golden " << endl;
		for(int i = 0; i < point; i++){
			cout << DFT_out[i] << endl;
		}		
		
		s1 = AddMod(a[1] , a[6], modular);
		s2 = SubMod(a[1] , a[6], modular);
		s3 = AddMod(a[4] , a[3], modular);
		s4 = SubMod(a[4] , a[3], modular);
		s5 = AddMod(a[2] , a[5], modular);
		s6 = SubMod(a[2] , a[5], modular);
		s7 =  AddMod(s1 , s3, modular);
		s8 =  SubMod(s1 , s3, modular);
		s9 =  SubMod(s5 , s1, modular);
		s10 = SubMod(s3 , s5, modular);
		s11 = AddMod(s2 , s4, modular);
		s12 = SubMod(s2 , s4, modular);
		s13 = SubMod(s6 , s2, modular);
		s14 = SubMod(s4 , s6, modular);
		s15 = AddMod(s7 , s5, modular);
		s16 = AddMod(s11 ,s6, modular);
		s17 = AddMod(a[0] , s15, modular);
		m0 = MulMod(s15, c0, modular);
		m1 = MulMod(s8,  c1, modular);	
		m2 = MulMod(s9,  c2, modular);
		m3 = MulMod(s10, c3, modular);
		m4 = MulMod(s16, c4, modular);
		m5 = MulMod(s12, c5, modular);
		m6 = MulMod(s13, c6, modular);
		m7 = MulMod(s14, c7, modular);
		s18 = AddMod(m0, s17, modular);	
		s19 = AddMod(s18, m1, modular);
		s20 = SubMod(s18, m1, modular);
		s21 = SubMod(s18, m3, modular);
		s22 = AddMod(m5, m4, modular);
		s23 = SubMod(m4, m5, modular);	
		s24 = SubMod(m4, m7, modular);
		
		s25 = AddMod(s19, m3, modular);
		s26 = SubMod(s20, m2, modular);
		s27 = AddMod(s21, m2, modular);
		s28 = AddMod(s22, m7, modular);
		s29 = SubMod(s23, m6, modular);
		
		s30 = AddMod(m6,  s24, modular);
		s31 = AddMod(s25, s28, modular);
		s32 = AddMod(s26, s29, modular);
		s33 = AddMod(s27, s30, modular);
		s34 = SubMod(s25, s28, modular);
		s35 = SubMod(s26, s29, modular);
		s36 = SubMod(s27, s30, modular);
		Winograd_out[0] = s17;
		Winograd_out[1] = s31;
		Winograd_out[2] = s32;
		Winograd_out[3] = s36;
		Winograd_out[4] = s33;
		Winograd_out[5] = s35;
		Winograd_out[6] = s34;	

		cout << "ans " << endl;
		for(int i = 0; i < point; i++){
			cout << Winograd_out[i] << endl;
		}		
	}
	else if(point == 5){
		//  5-WFFT
		long long a[point] ;
		long long DFT_out[point] ;
		long long Winograd_out[point] ;
		long long out(100) ;
		long long modular = 421;
		a[0] = 2;
		a[1] = 14;
		a[2] = 9;
		a[3] = 8;
		a[4] = 13;
		long long prou5;
		long long inv_2;
		inv_2 = winograd.find_inv(2, modular);
		prou5 = winograd.find_prou(5, modular);
		//cout << "prou5 = " << prou5 << endl;
		winograd.DFT(DFT_out, a, 5, prou5, modular);
		long long w0,w1,w2,w3,w4;
		
		w0 = PowerMod(prou5,0,modular);
		w1 = PowerMod(prou5,1,modular);
		w2 = PowerMod(prou5,2,modular);
		w3 = PowerMod(prou5,3,modular);
		w4 = PowerMod(prou5,4,modular);
		
		long long cosu, cos2u, isinu, isin2u;
		cosu = MulMod(AddMod(w1,w4,  modular), inv_2,  modular);
		isinu = MulMod(SubMod(w1,w4, modular), inv_2, modular);
		cos2u = MulMod(AddMod(w2,w3, modular), inv_2, modular);
		isin2u = MulMod(SubMod(w2,w3,modular), inv_2,modular);	
		
		cout << "golden " << endl;
		for(int i = 0; i < point; i++){
			cout << DFT_out[i] << endl;
		}
		
		// winograd 5-point FFT 
		long long s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17;
		long long m0,m1,m2,m3,m4,m5;
		long long c1,c2,c3,c4,c5;

		c1 = MulMod(AddMod(cosu,cos2u,modular), inv_2, modular) - 1;
		c2 = MulMod(SubMod(cosu,cos2u,modular), inv_2, modular);
		c3 = AddMod(isinu,isin2u,modular);
		c5 = isin2u;
		c4 = SubMod(isinu,isin2u,modular);
		
		s1 = AddMod(a[1], a[4], modular);
		s2 = SubMod(a[1], a[4], modular);
		s3 = AddMod(a[3], a[2], modular);
		s4 = SubMod(a[3], a[2], modular);
		s5 = AddMod(s1, s3,modular);
		s6 = SubMod(s1, s3,modular);
		s7 = AddMod(s2, s4,modular);
		s8 = AddMod(s5, a[0], modular);	

		m0 = s8;
		m1 = MulMod(s5, c1, modular);
		m2 = MulMod(s6, c2, modular);
		m3 = MulMod(s2, c3, modular);
		m4 = MulMod(s4, c4, modular);
		m5 = MulMod(s7, c5, modular);	
		
		s9  = AddMod(m0, m1, modular);
		s10 = SubMod(m3, m5, modular);	 // pic is add !!!! wrong
		s11 = AddMod(m4, m5, modular);	
		s12 = AddMod(s9, m2, modular);
		s13 = SubMod(s9, m2, modular);
		s14 = AddMod(s12, s10, modular);
		s15 = SubMod(s12, s10, modular);	
		s16 = AddMod(s13, s11, modular);
		s17 = SubMod(s13, s11, modular);	
		
		
		
		Winograd_out[0] = m0;
		Winograd_out[1] = s14;
		Winograd_out[2] = s16;
		Winograd_out[3] = s17;
		Winograd_out[4] = s15;
		
		cout << "ans " << endl;
		for(int i = 0; i < point; i++){
			cout << Winograd_out[i] << endl;
		}
	}
	else if(point == 4){
		//  4-WFFT
		long long a[point] ;
		long long DFT_out[point] ;
		long long Winograd_out[point] ;
		long long out(100) ;
		long long modular = 421;
		a[0] = 2;
		a[1] = 14;
		a[2] = 9;
		a[3] = 8;
		long long prou4;
		//long long inv_2;
		//inv_2 = winograd.find_inv(2, modular);
		prou4 = winograd.find_prou(4, modular);
		//cout << "prou5 = " << prou5 << endl;
		winograd.DFT(DFT_out, a, 4, prou4, modular);
		//long long w0,w1,w2,w3,w4;
		//
		//w0 = PowerMod(prou5,0,modular);
		//w1 = PowerMod(prou5,1,modular);
		//w2 = PowerMod(prou5,2,modular);
		//w3 = PowerMod(prou5,3,modular);
		//w4 = PowerMod(prou5,4,modular);
		
		//long long cosu, cos2u, isinu, isin2u;
		//cosu = MulMod(AddMod(w1,w4,  modular), inv_2,  modular);
		//isinu = MulMod(SubMod(w1,w4, modular), inv_2, modular);
		//cos2u = MulMod(AddMod(w2,w3, modular), inv_2, modular);
		//isin2u = MulMod(SubMod(w2,w3,modular), inv_2,modular);	
		
		cout << "golden " << endl;
		for(int i = 0; i < point; i++){
			cout << DFT_out[i] << endl;
		}
		
		// winograd 5-point FFT 
		long long s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17;
		long long m0,m1,m2,m3,m4,m5;
		long long c1,c2,c3,c4,c5;

		c1 = prou4;

		
		s1 = AddMod(a[0], a[2], modular);
		s2 = SubMod(a[0], a[2], modular);
		s3 = AddMod(a[1], a[3], modular);
		s4 = SubMod(a[1], a[3], modular);

		m1 = MulMod(s4, c1, modular);

		
		s5 = AddMod(s1, s3, modular);
		s6 = SubMod(s1, s3, modular);	
		s7 = AddMod(s2, m1, modular);
		s8 = SubMod(s2, m1, modular);	
		
		
		
		Winograd_out[0] = s5;
		Winograd_out[2] = s6;
		Winograd_out[1] = s7;
		Winograd_out[3] = s8;

		
		cout << "ans " << endl;
		for(int i = 0; i < point; i++){
			cout << Winograd_out[i] << endl;
		}	
	}
	else if(point == 3){
		//  4-WFFT
		long long a[point] ;
		long long DFT_out[point] ;
		long long Winograd_out[point] ;
		long long out(100) ;
		long long modular = 421;
		a[0] = 2;
		a[1] = 14;
		a[2] = 9;
		long long prou3;
		long long inv_2;
		inv_2 = winograd.find_inv(2, modular);
		prou3 = winograd.find_prou(3, modular);
		//cout << "prou5 = " << prou5 << endl;
		winograd.DFT(DFT_out, a, 3, prou3, modular);
		long long w0,w1,w2;
		//
		w0 = PowerMod(prou3,0,modular);
		w1 = PowerMod(prou3,1,modular);
		w2 = PowerMod(prou3,2,modular);

		
		long long cosu, cos2u, isinu, isin2u;
		cosu = MulMod(AddMod(w1,w2,  modular), inv_2,  modular);
		isinu = MulMod(SubMod(w1,w2, modular), inv_2, modular);
		//cos2u = MulMod(AddMod(w2,w3, modular), inv_2, modular);
		//isin2u = MulMod(SubMod(w2,w3,modular), inv_2,modular);	
		
		cout << "golden " << endl;
		for(int i = 0; i < point; i++){
			cout << DFT_out[i] << endl;
		}
		
		// winograd 5-point FFT 
		long long s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17;
		long long m0,m1,m2,m3,m4,m5;
		long long c1,c2,c3,c4,c5;

		c1 = SubMod(cosu, 1, modular);
		c2 = isinu;

		
		s1 = AddMod(a[1], a[2], modular);
		s2 = SubMod(a[1], a[2], modular);
		s3 = AddMod(a[0], s1, modular);

		m1 = MulMod(s1, c1, modular);
		m2 = MulMod(s2, c2, modular);
		
		s4 = AddMod(m1, s3, modular);
		s5 = AddMod(s4, m2, modular);	
		s6 = SubMod(s4, m2, modular);
	
		
		
		
		Winograd_out[0] = s3;
		Winograd_out[1] = s5;
		Winograd_out[2] = s6;


		
		cout << "ans " << endl;
		for(int i = 0; i < point; i++){
			cout << Winograd_out[i] << endl;
		}	
	}	
}