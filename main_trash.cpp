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
	//notice that the primitive root of unity exist
	//when m | modular - 1 (modular -1 can be divided by m)	
/*
	vector<vector<ZZ>> second_decompose(5, vector<ZZ>(5));
	vector<ZZ> second_cnt(5);
	vector<ZZ> first_decompose(5);	
	int first_cnt = precompute_gen.Factorize_fine(first_decompose ,second_decompose, second_cnt ,(ZZ)m); 
*/	
	//cout << "gen = " << precompute_gen.find_gen((ZZ)5);
	vector<ZZ> RA_in(3);
	vector<ZZ> RA_out(3);	
	RA_in[0] = 5;
	RA_in[1] = 7;
	RA_in[2] = 8;	
	
	//precompute_gen.RA_3P_FFT(RA_out,RA_in, 0);
	
	
	for(int i = 0; i < 3; i++){
		//cout << RA_out[i] << endl;
	}
	
	vector<ZZ> DFT_out(3);	
	//precompute_gen.DFT(DFT_out, RA_in, 3, precompute_gen.find_prou(3, precompute_gen.modular_m), precompute_gen.modular_m);
	
	for(int i = 0; i < 3; i++){
		//cout << DFT_out[i] << endl;
	}	
	
	vector<ZZ> RA_inv_out(3);	
	//precompute_gen.RA_3P_FFT(RA_inv_out,RA_out, 1);
	
	for(int i = 0; i < 3; i++){
		//cout << RA_inv_out[i] << endl;
	}
	
	vector<ZZ> RA_in4(8);
	vector<ZZ> RA_out4(8);	
	vector<ZZ> DFT_out4(8);		
	for(int i = 0; i < 8; i++){
		RA_in4[i] = i+1 ;
	}	
	
	precompute_gen.RA_powerof2_FFT(RA_out4, RA_in4, (ZZ)8);
	for(int i = 0; i < 8; i++){
		cout << RA_out4[i] << endl;
	}	
	
	
	precompute_gen.DFT(DFT_out4, RA_in4, 8, precompute_gen.find_prou(8, precompute_gen.modular_m), precompute_gen.modular_m);
	for(int i = 0; i < 8; i++){
		cout << DFT_out4[i] << endl;
	}		

	for(int i = 0; i < precompute_gen.first_cnt; i++){
		//cout << precompute_gen.first_decompose[i] << endl;
		for(int j = 0; j < precompute_gen.second_cnt[i]; j++){
			//cout << " " << precompute_gen.second_decompose[i][j];
			//cout << endl;
		}
		//cout << endl;
	}
		//cout << endl;

	/*
	for(int i = 0; i < first_cnt; i++){
		//cout << first_decompose[i] << endl;
		for(int j = 0; j < second_cnt[i]; j++){
			//cout <<" "<<second_decompose[i][j];
			//cout << endl;
		}
		//cout << endl;
	}
	*/
	//----first stage RA input reindex----//
	//[prime_first][idx]
/*	
	vector<ZZ> gen_first_decompose(first_cnt);
	vector<ZZ> gen_inv_first_decompose(first_cnt);	
	for(int i = 0; i < first_cnt; i++){
		gen_first_decompose[i] = precompute_gen.find_gen(first_decompose[i]);
		gen_inv_first_decompose[i] = precompute_gen.find_inv(gen_first_decompose[i], first_decompose[i] );
		//cout << gen_first_decompose[i] << endl;
	}
	
	vector<vector<ZZ>> input_index_first_decompose(first_cnt);
	vector<vector<ZZ>> output_index_first_decompose(first_cnt);	
	for(int i = 0; i < first_cnt; i++){
		input_index_first_decompose[i].resize(precompute_gen.ZZ2int(first_decompose[i]));
		input_index_first_decompose[i][0] = 0; 
		output_index_first_decompose[i].resize(precompute_gen.ZZ2int(first_decompose[i]));
		output_index_first_decompose[i][0] = 0; 		
		for(int j = 0; j < first_decompose[i] - 1 ; j++){
			input_index_first_decompose[i][j+1] = PowerMod(gen_first_decompose[i], j, first_decompose[i]);
			output_index_first_decompose[i][j+1] = PowerMod(gen_inv_first_decompose[i], j, first_decompose[i]);
		}
	}
*/	
/*
	for(int i = 0; i < precompute_gen.first_cnt; i++){ 
		for(int j = 0; j < precompute_gen.first_decompose[i] ; j++){
			cout << precompute_gen.input_index_first_decompose[i][j] << endl;
			cout << precompute_gen.output_index_first_decompose[i][j] << endl;			
		}
	}	
*/
	//----second stage RA input reindex----//
	//[prime_first][idx]	
/*	
cout << "------------------" << endl;

	vector<vector<ZZ>> gen_second_decompose(first_cnt);
	vector<vector<ZZ>> gen_inv_second_decompose(first_cnt);	
	for(int i = 0; i < first_cnt; i++){
		gen_second_decompose[i].resize(precompute_gen.ZZ2int(second_cnt[i]));
		gen_inv_second_decompose[i].resize(precompute_gen.ZZ2int(second_cnt[i]));
		for(int j = 0; j < second_cnt[i]; j++){
			gen_second_decompose[i][j] = precompute_gen.find_gen(second_decompose[i][j]);	
			gen_inv_second_decompose[i][j] = precompute_gen.find_inv_exgcd(gen_second_decompose[i][j], second_decompose[i][j] );
			//cout << gen_inv_second_decompose[i][j] << endl;
		}
	}	
	
	
	
	vector<vector<vector<ZZ>>> input_index_second_decompose(first_cnt);
	vector<vector<vector<ZZ>>> output_index_second_decompose(first_cnt);	
	
	for(int i = 0; i < first_cnt; i++){
		input_index_second_decompose[i].resize(precompute_gen.ZZ2int(second_cnt[i]));
		output_index_second_decompose[i].resize(precompute_gen.ZZ2int(second_cnt[i]));			
		for(int j = 0; j < second_cnt[i]; j++){
			input_index_second_decompose[i][j].resize(precompute_gen.ZZ2int(second_decompose[i][j]));
			output_index_second_decompose[i][j].resize(precompute_gen.ZZ2int(second_decompose[i][j]));	
			input_index_second_decompose[i][j][0] = 0;
			output_index_second_decompose[i][j][0] = 0;
			for(int k = 0; k < second_decompose[i][j] - 1; k++){	
				if(!precompute_gen.isPowerBy2(second_decompose[i][j])){
					input_index_second_decompose[i][j][k+1] = PowerMod(gen_second_decompose[i][j], k, second_decompose[i][j]);			
					output_index_second_decompose[i][j][k+1] = PowerMod(gen_inv_second_decompose[i][j], k, second_decompose[i][j]);

				}
				else{
					input_index_second_decompose[i][j][k+1] = k+1;
					output_index_second_decompose[i][j][k+1] = k+1;
				}
			}
		}		
	}
*/	
	/*
	for(int i = 0; i < first_cnt; i++){	
		for(int j = 0; j < second_cnt[i]; j++){
			for(int k = 0; k < second_decompose[i][j]; k++){	
				cout << input_index_second_decompose[i][j][k] << endl;
			}
		}		
	}
	*/
	

}