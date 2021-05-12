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
	LEGACY test;
	//notice that the primitive root of unity exist
	//when m | modular - 1 (modular -1 can be divided by m)	
	long long n = 11;	
	long long modular_n = 23;
	long long prou_n;	
	long long FFT_data_out[n];
	long long DFT_data_out[n];
	long long data_in[n]/*={0,1,2,3,4,5,6}*/;
	
	for (int i = 0; i < n ; i++){
		data_in[i] = i;
	}	
	
	
	prou_n = test.find_prou( n, modular_n) ;
	test.Rader(DFT_data_out , data_in , n , prou_n , modular_n);
	
}