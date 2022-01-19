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
	m = 5;
    LEGACY precompute_gen(m);
	ZZ LCM;
	ZZ prime;
	ZZ prou_LCM;
	long long a[8] = {21845, 30277, 32767, 52429, 53261, 61103, 78881, 92837};
	long long phi_a[8] = {16384, 26112, 27000, 46656, 46080, 49392, 70000, 84672};	
	long long blue_m[8] = {65536, 65536, 65536, 131072, 131072, 131072, 262144, 262144};
	long long rader_m[8] = {256, 512, 512, 256, 512, 128, 256, 256};
	
	long long optimize_m[42] = {2,4,6,10,12,16,18,22,25,28,30,36,40,42,49,60,70,72,88,96,100,108,112,126,132,136,150,240,250,256,330,336,396,432,600,630,672,682,756,880,1800,8190};
	
	vector<long long> optimize_vec(42);
	for(int i = 0; i < 42; i++){
		optimize_vec[i] = optimize_m[i];
	}
	
	long long LCM_optimize;
	
	LCM_optimize = precompute_gen.LCM(optimize_vec);
	
	cout << "LCM = " << LCM_optimize << "(" << ceil(log(LCM_optimize)/log(2.0)) << ")" << endl;
	
	long long p ;
/*	
	for(int i = 0; i < 8; i++){
		p = precompute_gen.find_prime(a[i], 0);
		
		cout << a[i] << ", "<< "prime = " << p << "bit = " << log(p)/log(2.0) << endl;
	}
*/	
/*
	    long long num = 65536;
		p = precompute_gen.find_prime(num, 0);
		cout << num << ", "<< "prime = " << p << "bit = " << log(p)/log(2.0) << endl;
*/
	long long Q;
/*	
	for(int i = 0; i < 8; i++){
		Q = a[i]*a[i]*phi_a[i];
		//cout << "Q = " << Q << "(" << ceil(log(Q)/log(2.0)) << ")" << endl;
		p = precompute_gen.find_prime_conditional(rader_m[i], Q);
		cout << "p = " << p << "(" << ceil(log(p)/log(2.0)) << ")" << endl;
	}
*/
	for(int i = 0; i < 8; i++){
		Q = a[i]*a[i]*phi_a[i];
		//cout << "Q = " << Q << "(" << ceil(log(Q)/log(2.0)) << ")" << endl;
		p = precompute_gen.find_prime_conditional(LCM_optimize, Q);
		cout << "p = " << p << "(" << ceil(log(p)/log(2.0)) << ")" << endl;
	}	
	
	
	
	
	
	
}