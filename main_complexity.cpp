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



	cout << prime << endl;

	vector<ZZ> m_set(8), blue_set(8), Rader_set(8);
	
	m_set[0] = 21845;
	m_set[1] = 30277;
	m_set[2] = 32767;
	m_set[3] = 52429;
	m_set[4] = 53261;
	m_set[5] = 61103;
	m_set[6] = 78881;
	m_set[7] = 92837;	

	blue_set[0] = 65536;
	blue_set[1] = 65536;
	blue_set[2] = 65536;
	blue_set[3] = 131072;
	blue_set[4] = 131072;
	blue_set[5] = 131072;
	blue_set[6] = 262144;
	blue_set[7] = 262144;		

	Rader_set[0] = 256;
	Rader_set[1] = 512;
	Rader_set[2] = 512;
	Rader_set[3] = 256;
	Rader_set[4] = 512;
	Rader_set[5] = 128;
	Rader_set[6] = 256;
	Rader_set[7] = 256;		
	/*
	for(int i = 0; i < 8; i++){
		cout << "m = " << m_set[i] << endl;		
			LCM = precompute_gen.LCM( m_set[i], blue_set[i]);
			prime = precompute_gen.find_prime(LCM, 0);			
		cout << "Blue prime = " << prime << ", bit = " << log(prime)/log(2) << endl;
			LCM = precompute_gen.LCM( m_set[i], Rader_set[i]);
			prime = precompute_gen.find_prime(LCM, 0);					
		cout << "Rader prime = " << prime << ", bit = " << log(prime)/log(2) << endl;	
	}
	*/
	
			vector<ZZ> LCM_set(13);
			ZZ LCM_1, LCM_m;
			LCM_set[0] = 256;
			LCM_set[1] = 12;
			LCM_set[2] = 136;
			LCM_set[3] = 150;
			LCM_set[4] = 36;
			LCM_set[5] = 108;	
			LCM_set[6] = 240;
			LCM_set[7] = 28;
			LCM_set[8] = 42;	
			LCM_set[9] = 49;
			LCM_set[10] = 70;
			LCM_set[11] = 100;	
			LCM_set[12] = 126;
			
			LCM_1 = precompute_gen.LCM(LCM_set);
			cout << "LCM_1 = " << LCM_1 << ", bit = " << log(LCM_1)/log(2) << endl;				
			LCM_m = precompute_gen.LCM(m_set);
			cout << "LCM_m = " << LCM_m << ", bit = " << log(LCM_m)/log(2) << endl;	
			LCM = precompute_gen.LCM( LCM_m, LCM_1);
			prime = precompute_gen.find_prime(LCM, 0);	
			cout << "prime = " << prime << ", bit = " << log(prime)/log(2) << endl;	
			
			
			
			
			
}