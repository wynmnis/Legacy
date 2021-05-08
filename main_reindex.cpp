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
	vector<int> vec_105 ;
	vec_105.resize(105);
	for(int i = 0; i < 105 ; i++){
		vec_105[i] = i ;
	}
	
	int count = 0;	
	int input_tmp;
	//-------105 PFA input re-index --- 1
	vector<int> vec_3_35_input_105 ;
	vec_3_35_input_105.resize(105);
	for(int n2 = 0; n2 < 35 ; n2++){
		for(int n1 = 0; n1 < 3 ; n1++){
			input_tmp = (n2*3 + n1*35)%105;
			vec_3_35_input_105[count] = input_tmp;
			//cout << "tmp[" << count << "]=" << vec_3_35_input_105[count]<<endl;
			++count ;
		}		
	}
	// 3p -35p middle index --- 2
	vector<int> vec_3_35_middle_105 ;
	vec_3_35_middle_105.resize(105);
	count = 0;	
	for(int n2 = 0; n2 < 3 ; n2++){
		for(int n1 = 0; n1 < 35 ; n1++){
			input_tmp = (n2 + n1*3)%105;
			vec_3_35_middle_105[count] = input_tmp;
			//cout /*<< "tmp[" << count << "]="*/ << vec_3_35_middle_105[count]<<endl;
			++count ;
		}		
	}
	
	// 35 PFA input re-index --- 3
	vector<int> vec_5_7_input_35 ;
	vec_5_7_input_35.resize(35);
	count = 0;	
	for(int n2 = 0; n2 < 7 ; n2++){
		for(int n1 = 0; n1 < 5 ; n1++){
			input_tmp = (n2*5 + n1*7)%35;
			vec_5_7_input_35[count] = input_tmp;
			//cout << "tmp[" << count << "]=" << vec_5_7_input_35[count]<<endl;
			++count ;
		}		
	}	
	vector<int> vec_5_7_input_105 ;
	vec_5_7_input_105.resize(105);
	count = 0;	
	for(int i = 0; i < 3 ; i++){	
		for(int j = 0; j < 35; j++){
			input_tmp = (vec_3_35_middle_105[vec_5_7_input_35[j]+35*i])%105;
			vec_5_7_input_105[count] = input_tmp;
			//cout /*<< "tmp[" << count << "]=" */<< vec_5_7_input_105[count]<<endl;
			++count ;			
		}
	}
	
	// 5p -7p middle index --- 4
	vector<int> vec_5_7_middle_35 ;
	vec_5_7_middle_35.resize(35);
	count = 0;	
	for(int n2 = 0; n2 < 5 ; n2++){
		for(int n1 = 0; n1 < 7 ; n1++){
			input_tmp = (n2*7 + n1*5)%35;
			vec_5_7_middle_35[count] = input_tmp;
			//cout << "tmp[" << count << "]=" << vec_5_7_middle_35[count]<<endl;
			++count ;
		}		
	}	
	
	vector<int> vec_5_7_middle_105 ;
	vec_5_7_middle_105.resize(105);
	count = 0;	
	for(int i = 0; i < 3 ; i++){	
		for(int j = 0; j < 35; j++){
			input_tmp = (vec_5_7_middle_35[j]*3 + i*35)%105;
			vec_5_7_middle_105[count] = input_tmp;
			//cout << vec_5_7_middle_105[count]<<endl;
			++count ;			
		}
	}	


	// 5p -7p output index --- 5
	int output_tmp;
	vector<int> vec_5_7_output_35 ;
	vec_5_7_output_35.resize(35);
	count = 0;	
	for(int k1 = 0; k1 < 5 ; k1++){
		for(int k2 = 0; k2 < 7 ; k2++){
			output_tmp = (k2*15 + k1*21)%35;
			vec_5_7_output_35[count] = output_tmp;
			//cout << "tmp[" << count << "]=" << vec_5_7_output_35[count]<<endl;
			++count ;
		}		
	}
	
	vector<int> vec_5_7_output_105 ;
	vec_5_7_output_105.resize(105);
	count = 0;	
	for(int i = 0; i < 3 ; i++){	
		for(int j = 0; j < 35; j++){
			input_tmp = (vec_5_7_middle_105[vec_5_7_output_35[j]+35*i])%105;
			vec_5_7_output_105[count] = input_tmp;
			//cout << vec_5_7_output_105[count]<<endl;
			++count ;			
		}
	}	

		
}