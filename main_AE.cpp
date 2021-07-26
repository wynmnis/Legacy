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
	//test AE 
	LEGACY AE;
	//varilable definition
	int N = 1024 ;
	int r = 2 ;
	int p,g,s;
	p = log(N)/log(r) ; //4
	g = N/(r*r) ;		//16
	s = log2(r) ;		//8
	int    BC_WIDTH;    //6
	BC_WIDTH = (int)ceil(log2(N/r));	
	//main function
	int BC, BN, MA ;
	int RR_out = 0;
	int xor_out = 0;
	int bit_width = log2(N);
	vector<int> bit_array;	
	int Data ;
	
	for (int t = 0; t < p; t++)
	{
		cout << "stage "<< t << endl ;
		for(int i = 0; i < g; i++)
		{
			for(int j = 0; j < r; j++)
			{
				//-----------------------------
				BC = j*g + AE.Gray(i,g) ;
				//cout << "BC = "<< BC << endl ;
				RR_out = AE.RR(BC, s*t, BC_WIDTH);
				xor_out = AE.unary_xor(RR_out , BC_WIDTH);
				BN = xor_out;
				MA = RR_out >> 1;	
				cout << "(BC, BN, MA) = ";				
				cout << "(" << BC << " , "<< BN <<" , "<< MA << ")";	
				//-------------------------------
				AE.int2vec(BC, bit_width, bit_array);				
				std::rotate(bit_array.begin(), bit_array.begin()+ s*t , bit_array.end());
				Data = AE.vec2int(bit_array, bit_width);
				cout << "Data_index = ";				
				cout << "(" ;	
				for(int k = 0; k < r ; k++ ){
					cout << Data + k*(1<<(bit_width-s-s*t)) <<"  ";	
				}
				cout << ") \n" ;	
			}
		}
	}

	
	//cout << AE.Gray(16,4096) << endl;
	
}