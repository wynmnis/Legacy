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
    LEGACY precompute_gen(m);




cout << "*---------------------------------------------------*" << endl;
// reconfig 2,3,4,5,7-WFFT
	ZZ modulus = (ZZ)421;
	ZZ prou420, prou7, prou5, prou4, prou3, prou2;
	ZZ prou7_1, prou7_2, prou7_3, prou7_4, prou7_5, prou7_6;
	bool inverse = 0;
	
	prou420 = precompute_gen.find_prou(420, modulus);
	prou7 = PowerMod(prou420, (ZZ)(420/7), modulus);
	prou5 = PowerMod(prou420, (ZZ)(420/5), modulus);
	prou4 = PowerMod(prou420, (ZZ)(420/4), modulus);
	prou3 = PowerMod(prou420, (ZZ)(420/3), modulus);
	prou2 = PowerMod(prou420, (ZZ)(420/2), modulus);
	
	ZZ inv_3, inv_2, inv_7;
	inv_3 = precompute_gen.find_inv((ZZ)3, modulus);
	inv_2 = precompute_gen.find_inv((ZZ)2, modulus);
	inv_7 = precompute_gen.find_inv((ZZ)7, modulus);
	
	if(inverse)
		prou7_1 = precompute_gen.find_inv(prou7, modulus);
	else 
		prou7_1 = prou7;
	
	prou7_2 = PowerMod(prou7_1, 2, modulus);
	prou7_3 = PowerMod(prou7_1, 3, modulus);
	prou7_4 = PowerMod(prou7_1, 4, modulus);
	prou7_5 = PowerMod(prou7_1, 5, modulus);
	prou7_6 = PowerMod(prou7_1, 6, modulus);	

	ZZ cosu, cos2u, cos3u, isinu, isin2u, isin3u;
	cosu = MulMod(AddMod(prou7_1, prou7_6, modulus), inv_2, modulus);
	isinu = MulMod(SubMod(prou7_1, prou7_6, modulus), inv_2, modulus);
	cos2u = MulMod(AddMod(prou7_2, prou7_5, modulus), inv_2, modulus);
	isin2u = MulMod(SubMod(prou7_2, prou7_5, modulus), inv_2, modulus);
	cos3u = MulMod(AddMod(prou7_3, prou7_4, modulus), inv_2, modulus);
	isin3u = MulMod(SubMod(prou7_3, prou7_4, modulus), inv_2, modulus);	
	
	

	ZZ C70, C71, C72, C73, C74, C75, C76, C77;// 論文有錯 注意Index
	C70 = SubMod(MulMod(AddMod(AddMod(cosu , cos2u, modulus), cos3u, modulus), inv_3, modulus), (ZZ)1, modulus);
	C71 = MulMod(SubMod(SubMod(MulMod(cosu,2,modulus) , cos2u, modulus), cos3u, modulus), inv_3, modulus);	
	C72 = MulMod(AddMod(SubMod(cosu , MulMod(cos2u,2,modulus), modulus), cos3u, modulus), inv_3, modulus);	
	C73 = MulMod(SubMod(AddMod(cosu , cos2u, modulus), MulMod(cos3u,2,modulus), modulus), inv_3, modulus);	
	C74 = MulMod(SubMod(AddMod(isinu , isin2u, modulus), isin3u, modulus), inv_3, modulus);
	C75 = MulMod(AddMod(SubMod(MulMod(isinu,2,modulus) , isin2u, modulus), isin3u, modulus), inv_3, modulus);
	C76 = MulMod(SubMod(SubMod(isinu , MulMod(isin2u,2,modulus), modulus), isin3u, modulus), inv_3, modulus);
	C77 = MulMod(AddMod(AddMod(isinu , isin2u, modulus), MulMod(isin3u,2,modulus), modulus), inv_3, modulus);
	
	//cout << "cosu = " << cosu << endl;
	//cout << "isinu = " << isinu << endl;	
	//cout << "cos2u = " << cos2u << endl;
	//cout << "isin2u = " << isin2u << endl;
	//cout << "cos3u = " << cos3u << endl;
	//cout << "isin3u = " << isin3u << endl;	
	//
	//cout << "prou420 = " << prou420 << endl;
	//cout << "prou7 = " << prou7 << endl;
	//cout << "prou5 = " << prou5 << endl;
	//cout << "prou4 = " << prou4 << endl;
	//cout << "prou3 = " << prou3 << endl;
	//cout << "prou2 = " << prou2 << endl;
	//cout << "inv_3 = " << inv_3 << endl;	
	
	//cout << "C70 = " << C70 << endl;
	//cout << "C71 = " << C71 << endl;
	//cout << "C72 = " << C72 << endl;
	//cout << "C73 = " << C73 << endl;
	//cout << "C74 = " << C74 << endl;
	//cout << "C75 = " << C75 << endl;
	//cout << "C76 = " << C76 << endl;
	//cout << "C77 = " << C77 << endl;	
	ZZ prou5_1, prou5_2, prou5_3, prou5_4;
	prou5_1 = prou5;
	prou5_2 = PowerMod(prou5_1, 2, modulus);
	prou5_3 = PowerMod(prou5_1, 3, modulus);
	prou5_4 = PowerMod(prou5_1, 4, modulus);	
	cosu = MulMod(AddMod(prou5_1, prou5_4, modulus), inv_2, modulus);
	isinu = MulMod(SubMod(prou5_1, prou5_4, modulus), inv_2, modulus);
	cos2u = MulMod(AddMod(prou5_2, prou5_3, modulus), inv_2, modulus);
	isin2u = MulMod(SubMod(prou5_2, prou5_3, modulus), inv_2, modulus);	
	
	ZZ C50,C51,C52,C53,C54;
	C50 = SubMod(MulMod(AddMod(cosu, cos2u, modulus), inv_2, modulus), (ZZ)1, modulus);
	C51 = MulMod(SubMod(cosu, cos2u, modulus), inv_2, modulus);	
	C52 = AddMod(isinu, isin2u, modulus);	
	C53 = isin2u;	
	C54 = SubMod(isinu, isin2u, modulus);
	
	ZZ C41;
	C41 = prou4;
	

	ZZ prou3_1, prou3_2;
	prou3_1 = prou3;
	prou3_2 = PowerMod(prou3_1, 2, modulus);
	cosu = MulMod(AddMod(prou3_1, prou3_2, modulus), inv_2, modulus);
	isinu = MulMod(SubMod(prou3_1, prou3_2, modulus), inv_2, modulus);		
	ZZ C30, C31;
	C30 = SubMod(cosu, (ZZ)1, modulus);
	C31 = isinu;
	

	int point = 5;
	vector<ZZ> x(7);
	vector<ZZ> X(7);
	vector<ZZ> x_5(5), x_4(4), x_3(3), x_2(2);
	ZZ s0,s1,s2,s3;
	ZZ A0,A1,A2,A3,A4,A5,A6,A7;


	if(point == 7){
		x[0] = 1;
		x[1] = 3;
		x[2] = 25;
		x[3] = 5;
		x[4] = 6;
		x[5] = 15;
		x[6] = 91;
		
		s0 = 1;
		s1 = 1;
		s2 = 1;
		s3 = 1;	
		
		A0 = C70;
		A1 = C71;
		A2 = C73;
		A3 = C72;
		A4 = C74;
		A5 = C75;
		A6 = C77;	
		A7 = C76;		
	}
	else if(point == 5){
		x[0] = 1;
		x[1] = 3;
		x[2] = 0;
		x[3] = 5;
		x[4] = 6;
		x[5] = 0;
		x[6] = 91;	
		
		x_5[0] = x[0]; 
		x_5[1] = x[1];
		x_5[2] = x[4];
		x_5[3] = x[3];
		x_5[4] = x[6];
		
		s0 = 0;
		s1 = 1;
		s2 = 1;
		s3 = 1;		
		
		A0 = C50;
		A1 = C53;
		A2 = 0;
		A3 = 0;
		A4 = C52;
		A5 = C51;
		A6 = 0;	
		A7 = C54;			
	}
	else if(point == 4){
		x[0] = 0;
		x[1] = 3;
		x[2] = 0;
		x[3] = 5;
		x[4] = 6;
		x[5] = 0;
		x[6] = 91;	
		
		x_4[0] = x[1]; 
		x_4[3] = x[3];
		x_4[1] = x[4];
		x_4[2] = x[6];	
		
		s0 = 0;
		s1 = 0;
		s2 = 0;
		s3 = 0;	
		
		A0 = 0;
		A1 = 0;
		A2 = 0;
		A3 = 0;
		A4 = 1;
		A5 = C41;
		A6 = 0;
		A7 = 1;			
	}
	else if(point == 3){
		x[0] = 1;
		x[1] = 3;
		x[2] = 0;
		x[3] = 0;
		x[4] = 0;
		x[5] = 0;
		x[6] = 91;	
		
		x_3[0] = x[0]; 
		x_3[1] = x[1];
		x_3[2] = x[6];
		
		s0 = 0;
		s1 = 1;
		s2 = 0;
		s3 = 1;	
		
		A0 = C30;
		A1 = 0;
		A2 = 0;
		A3 = 0;
		A4 = C31;
		A5 = 0;
		A6 = 0;
		A7 = 0;			
	}
	else if(point == 2){
		x[0] = 0;
		x[1] = 3;
		x[2] = 0;
		x[3] = 0;
		x[4] = 0;
		x[5] = 0;
		x[6] = 91;	
		
		x_2[0] = x[1]; 
		x_2[1] = x[6];
		
		s0 = 0;
		s1 = 0;
		s2 = 0;
		s3 = 0;	
		
		A0 = 0;
		A1 = 0;
		A2 = 0;
		A3 = 0;
		A4 = 0;
		A5 = 0;
		A6 = 1;
		A7 = 0;			
	}		
	
	ZZ a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36;
	ZZ m1,m2,m3,m4,m5,m6,m7;
	ZZ mul0,mul1,mul2,mul3,mul4,mul5,mul6,mul7,mul8;


	
	a1 = AddMod(x[1] , x[6], modulus);
	a2 = SubMod(x[1] , x[6], modulus);
	a3 = AddMod(x[4] , x[3], modulus);
	a4 = SubMod(x[4] , x[3], modulus);
	a5 = AddMod(x[2] , x[5], modulus);
	a6 = SubMod(x[2] , x[5], modulus);
	a7 = AddMod(a1 , a3, modulus);
	a8 = SubMod(a1 , a3, modulus);
	a9 = SubMod(a5 , a1, modulus);
	a10 = SubMod(a3 , a5, modulus);
	if(s1==1) m1 = a2;
	else m1 = 0;
	if(s0==1) m2 = a4;
	else m2 = 0;	
	a11 = AddMod(a4 , m1, modulus);///
	//a11 = AddMod(a2 , a4, modulus);
	a12 = SubMod(a2 , m2, modulus);
	//a12 = SubMod(a2 , a4, modulus);
	a13 = SubMod(a6 , a2, modulus);
	a14 = SubMod(a4 , a6, modulus);
	a15 = AddMod(a7 , a5, modulus);
	a16 = AddMod(a11 , a6, modulus);
	a17 = AddMod(x[0] , a15, modulus);
	if(s1==1) m3 = a14;
	else m3 = a8;
	mul0 = MulMod(a15, A0, modulus);
	mul1 = MulMod(a8, A1, modulus);	
	mul2 = MulMod(a9, A2, modulus);
	mul3 = MulMod(a10, A3, modulus);
	mul4 = MulMod(a16, A4, modulus);
	mul5 = MulMod(a12, A5, modulus);
	mul6 = MulMod(a13, A6, modulus);
	mul7 = MulMod(m3, A7, modulus);
	//mul7 = MulMod(a14, A7, modulus);
	if(s3==1) m4 = a17;
	else m4 = 0;	
	a18 = AddMod(mul0, m4, modulus);
	//a18 = AddMod(mul0, a17, modulus);
	if(s2==1) m5 = mul5;
	else m5 = 0;		
	if(s1==1) m6 = mul7;
	else m6 = mul6;		
	a19 = AddMod(a18, mul1, modulus);
	a20 = SubMod(a18, mul1, modulus);
	a21 = SubMod(a18, mul3, modulus);
	a22 = AddMod(m5, mul4, modulus);
	//a22 = AddMod(mul4, mul5, modulus);
	a23 = SubMod(mul4, mul5, modulus);	
	a24 = SubMod(mul4, mul7, modulus);
	if(s0==1) m7 = mul6;
	else m7 = 0;
	a25 = AddMod(a19, mul3, modulus);
	a26 = SubMod(a20, mul2, modulus);
	a27 = AddMod(a21, mul2, modulus);
	a28 = AddMod(a22, m6, modulus);
	//a28 = AddMod(a22, mul7, modulus);
	a29 = SubMod(a23, m7, modulus);
	//a29 = SubMod(a23, mul6, modulus);
	a30 = AddMod(mul6, a24, modulus);
	a31 = AddMod(a25, a28, modulus);
	a32 = AddMod(a26, a29, modulus);
	a33 = AddMod(a27, a30, modulus);
	a34 = SubMod(a25, a28, modulus);
	a35 = SubMod(a26, a29, modulus);
	a36 = SubMod(a27, a30, modulus);
	X[0] = a17;
	X[1] = a31;
	X[2] = a32;
	X[3] = a36;
	X[4] = a33;
	X[5] = a35;
	X[6] = a34;
	
	cout << "a1 = " << a1 << endl;	
	cout << "a2 = " << a2 << endl;	
	cout << "a3 = " << a3 << endl;	
	cout << "a4 = " << a4 << endl;	
	cout << "a7 = " << a7 << endl;	
	cout << "a8 = " << a8 << endl;		
	cout << "a13 = " << a13 << endl;		
	cout << "a18 = " << a18 << endl;
	cout << "mul6 = " << mul6 << endl;
	cout << "a27 = " << a27 << endl;
	cout << "a30 = " << a30 << endl;
	cout << "a34 = " << a34 << endl;	
	cout << "a36 = " << a36 << endl;	
	

	
	cout << "ans = " << endl; 
	if(point == 7){
		for(int i = 0; i < 7; i++){
			if(inverse)		
				cout << MulMod(X[i], inv_7, modulus) << endl;
			else
				cout << X[i] << endl;			
		}
	}
	else if (point == 5){
		for(int i = 0; i < 7; i++){
			if(inverse)		
				cout << MulMod(X[i], inv_7, modulus) << endl;
			else
				cout << X[i] << endl;			
		}		
	}
	
	
	vector<ZZ> DFT_out7(7);
	vector<ZZ> DFT_out5(5);
	vector<ZZ> DFT_out4(4);
	vector<ZZ> DFT_out3(3);	
	vector<ZZ> DFT_out2(2);
	
	cout << "golden = " << endl; 
	if(point == 7){
		if(inverse)
			precompute_gen.IDFT(DFT_out7, x, 7, prou7, modulus);
		else
			precompute_gen.DFT(DFT_out7, x, 7, prou7, modulus);			
		for(int i = 0; i < point; i++){
			cout << DFT_out7[i] << endl;
		}	
	}
	else if (point == 5){
		precompute_gen.DFT(DFT_out5, x_5, 5, prou5, modulus);	
		for(int i = 0; i < point; i++){
			cout << DFT_out5[i] << endl;
		}			
	}
	else if (point == 4){
		precompute_gen.DFT(DFT_out4, x_4, 4, prou4, modulus);	
		for(int i = 0; i < point; i++){
			cout << DFT_out4[i] << endl;
		}			
	}	
	else if (point == 3){
		precompute_gen.DFT(DFT_out3, x_3, 3, prou3, modulus);	
		for(int i = 0; i < point; i++){
			cout << DFT_out3[i] << endl;
		}			
	}	
	else if (point == 2){
		precompute_gen.DFT(DFT_out2, x_2, 2, prou2, modulus);	
		for(int i = 0; i < point; i++){
			cout << DFT_out2[i] << endl;
		}			
	}
	
	
/*	
    //  5-WFFT
	vector<ZZ> a(5) ;
	vector<ZZ> DFT_out5(5) ;
	vector<ZZ> Winograd_out5(5) ;
	vector<ZZ> out(100) ;
	ZZ mod_5 = (ZZ)31;
	a[0] = 2;
	a[1] = 14;
	a[2] = 9;
	a[3] = 8;
	a[4] = 13;
	ZZ prou5;
	ZZ inv_2;
	inv_2 = 16;
	prou5 = precompute_gen.find_prou(5, mod_5);
	cout << "prou5 = " << prou5 << endl;
	precompute_gen.DFT(DFT_out5, a, 5, prou5, mod_5);
	ZZ w0,w1,w2,w3,w4;
	
	w0 = PowerMod(prou5,0,mod_5);
	w1 = PowerMod(prou5,1,mod_5);
	w2 = PowerMod(prou5,2,mod_5);
	w3 = PowerMod(prou5,3,mod_5);
	w4 = PowerMod(prou5,4,mod_5);
	
	ZZ cosu, cos2u, isinu, isin2u;
	cosu = MulMod(AddMod(w1,w4,mod_5), inv_2, mod_5);
	isinu = MulMod(SubMod(w1,w4,mod_5), inv_2, mod_5);
	cos2u = MulMod(AddMod(w2,w3,mod_5), inv_2, mod_5);
	isin2u = MulMod(SubMod(w2,w3,mod_5), inv_2, mod_5);	
	
	
	for(int i = 0; i < 5; i++){
		cout << DFT_out5[i] << endl;
	}
	
	// winograd 5-point FFT 
	ZZ s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17;
	ZZ m0,m1,m2,m3,m4,m5;
	ZZ c1,c2,c3,c4,c5;

	c1 = MulMod(AddMod(cosu,cos2u,mod_5), inv_2, mod_5) - 1;
	c2 = MulMod(SubMod(cosu,cos2u,mod_5), inv_2, mod_5);
	c3 = AddMod(isinu,isin2u,mod_5);
	c5 = isin2u;
	c4 = SubMod(isinu,isin2u,mod_5);
	
	s1 = AddMod(a[1], a[4], mod_5);
	s2 = SubMod(a[1], a[4], mod_5);
	s3 = AddMod(a[3], a[2], mod_5);
	s4 = SubMod(a[3], a[2], mod_5);
	s5 = AddMod(s1, s3, mod_5);
	s6 = SubMod(s1, s3, mod_5);
	s7 = AddMod(s2, s4, mod_5);
	s8 = AddMod(s5, a[0], mod_5);	

	m0 = s8;
	m1 = MulMod(s5, c1, mod_5);
	m2 = MulMod(s6, c2, mod_5);
	m3 = MulMod(s2, c3, mod_5);
	m4 = MulMod(s4, c4, mod_5);
	m5 = MulMod(s7, c5, mod_5);	
	
	s9  = AddMod(m0, m1, mod_5);
    s10 = SubMod(m3, m5, mod_5);	 // pic is add !!!! wrong
	s11 = AddMod(m4, m5, mod_5);	
	s12 = AddMod(s9, m2, mod_5);
	s13 = SubMod(s9, m2, mod_5);
	s14 = AddMod(s12, s10, mod_5);
	s15 = SubMod(s12, s10, mod_5);	
	s16 = AddMod(s13, s11, mod_5);
	s17 = SubMod(s13, s11, mod_5);	
	
	
	
	Winograd_out5[0] = m0;
	Winograd_out5[1] = s14;
	Winograd_out5[2] = s16;
	Winograd_out5[3] = s17;
	Winograd_out5[4] = s15;
	
	cout << "ans " << endl;
	for(int i = 0; i < 5; i++){
		cout << Winograd_out5[i] << endl;
	}
*/	


/*  3-WFFT
	vector<ZZ> a(3) ;
	vector<ZZ> DFT_out3(3) ;
	vector<ZZ> Winograd_out3(3) ;
	ZZ s1,s2,s3,s4,s5,s6;
	ZZ m0,m1,m2,m3;	
	ZZ c1,c2;
	ZZ mod_3 = (ZZ)97;
	a[0] = 2;
	a[1] = 14;
	a[2] = 26;
	c1 = 47; //cosu - 1
	c2 = 84; //


	s1 = AddMod(a[1], a[2], mod_3);
	s2 = SubMod(a[1], a[2], mod_3);
	s3 = AddMod(a[0], s1, mod_3);
	m0 = s3;
	m1 = MulMod(s1, c1, mod_3);
	m2 = MulMod(s2, c2, mod_3);
	s4 = AddMod(m0, m1, mod_3);
	s5 = AddMod(s4, m2, mod_3);
	s6 = SubMod(s4, m2, mod_3);
	Winograd_out3[0] = m0;
	Winograd_out3[1] = s5;
	Winograd_out3[2] = s6;

	cout << "golden " << endl;
	precompute_gen.DFT(DFT_out3, a, 3, (ZZ)35, (ZZ)97);
	for(int i = 0; i < 3; i++){
		cout << DFT_out3[i] << endl;
	}	

	cout << "ans " << endl;
	for(int i = 0; i < 3; i++){
		cout << Winograd_out3[i] << endl;
	}
*/

	/*
	a1 = (x[1] + x[6]) % modulus; 
	a2 = (x[1] - x[6]) % modulus; 
	a3 = (x[4] + x[3]) % modulus; 	
	a4 = (x[4] - x[3]) % modulus; 
	a5 = (x[2] + x[5]) % modulus; 
	a6 = (x[2] - x[5]) % modulus; 
	a7 = (a1 + a3) % modulus; 
	a8 = (a7 + a5) % modulus; 
	a9 = (a8 + x[0]) % modulus; 
	a10 = (a1 - a3) % modulus; 
	a11 = (a3 - a5) % modulus; 
	a12 = (a5 - a1) % modulus; 
	a13 = (a2 + a4) % modulus; 
	a14 = (a13 + a6) % modulus; 
	a15 = (a2 - a4) % modulus; 
	a16 = (a4 - a6) % modulus; 
	a17 = (a6 - a2) % modulus; 
	mul0 = a9;
	mul1 = (A0 * a8)%modulus;
	mul2 = (A1 * a10)%modulus;
	mul3 = (A2 * a11)%modulus;
	mul4 = (A3 * a12)%modulus;
	mul5 = (A4 * a14)%modulus;
	mul6 = (A5 * a15)%modulus;	
	mul7 = (A6 * a16)%modulus;	
	mul8 = (A7 * a17)%modulus;	
	a18 = (mul0 + mul1) % modulus;	
	a19 = (a18 + mul2) % modulus;	
	a20 = (a19 + mul3) % modulus;		
	a21 = (a18 - mul2) % modulus;	
	a22 = (a21 - mul4) % modulus;
	a23 = (a18 - mul3) % modulus;
	a24 = (a23 + mul4) % modulus;
	a25 = (mul5 + mul6) % modulus;
	a26 = (a25 + mul7) % modulus;
	a27 = (mul5 - mul6) % modulus;
	a28 = (a27 - mul8) % modulus;
	a29 = (mul5 - mul7) % modulus;
	a30 = (a29 + mul8) % modulus;
	a31 = (a20 + a26) % modulus;
	a32 = (a20 - a26) % modulus;
	a33 = (a22 + a28) % modulus;
	a34 = (a22 - a28) % modulus;
	a35 = (a24 + a30) % modulus;
	a36 = (a24 - a30) % modulus;	
	X[0] = mul0;
	X[1] = a31;
	X[2] = a33;
	X[3] = a36;
	X[4] = a35;
	X[5] = a34;
	X[6] = a32;
	*/


}