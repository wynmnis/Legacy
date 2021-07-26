#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "LEGACY.h"
#include <assert.h>
#include <NTL/ZZ.h>
#include<bits/stdc++.h> 
//#define system_m 105            //modular = 211,4
//#define system_phim 48
//#define system_m1 15
//#define system_m2 7
//#define system_small_m1 3
//#define system_small_m2 5
//#define system_modular 211

//#define system_m 1023           //modular = 4093,16
//#define system_phim 600
//#define system_m1 33
//#define system_m2 31
//#define system_small_m1 3
//#define system_small_m2 11
//#define system_modular 4093

//#define system_m 1705           //modular = 436481,256
//#define system_phim 1200
//#define system_m1 55
//#define system_m2 31
//#define system_small_m1 5
//#define system_small_m2 11
//#define system_modular 436481

#define system_m 21845           //modular = 4718521,291
#define system_phim 16384
#define system_m1 85
#define system_m2 257
#define system_small_m1 5
#define system_small_m2 17
#define system_modular 4718521   //use RA-FFT modular : 134215681

//#define system_m 24295           //modular = 4718521,291
//#define system_phim 18816
//#define system_m1 215
//#define system_m2 113
//#define system_small_m1 5
//#define system_small_m2 43
//#define system_modular 145771

//#define system_m 27305          //modular = 4718521,291
//#define system_phim 21168
//#define system_m1 215
//#define system_m2 127
//#define system_small_m1 5
//#define system_small_m2 43
//#define system_modular 327661

//#define system_m 28679          //modular = 4718521,291
//#define system_phim 23040
//#define system_m1 119
//#define system_m2 241
//#define system_small_m1 7
//#define system_small_m2 17
//#define system_modular 229433

//#define system_m 32767          //modular = 1310681,2
//#define system_phim 27000
//#define system_m1 217
//#define system_m2 151
//#define system_small_m1 7
//#define system_small_m2 31
//#define system_modular 50330113 //1068 857736 57277

using namespace std;
using namespace NTL;

long long LEGACY::Euler(long long x){
    if (x < 2) return 0;
    int ret = x;
    int sq = sqrt(x);
    for (int p=2; p<=sq; p++){
        if (x % p == 0){
            while (x % p == 0) x /= p;
            ret -= ret / p;
        }
        if (x == 1) break;
    }
    if (x > 1) ret -= ret / x;
    return ret;
}

// find the power of 2 integer which is m' >= 2m-3
long long LEGACY::find_m_prime(long long m){
	long long n = ceil(log2(2*m-3))  ;

	long long ans = 0; ;
	ans = pow(2,n);
	return ans ;
}

ZZ LEGACY::find_m_prime(ZZ m){
	long long n = ceil(NumBits(2*m-3))  ;
	ZZ ans ;
	power(ans, 2 , n);
	return ans ;	
}

bool LEGACY::coprime(long long a, long long b){

	if(a==1||b==1)    
		return true;
	while(1)
    {       
		int t = a%b;
		if(t == 0) 
        {
            break;
        }
		else
        {
			a = b;
			b = t;
		}
	}
	if(b>1)	
		return false;
	else 
		return true;	 
}

bool LEGACY::isPowerBy2(long long n)
{
    return n > 0 && (n & n - 1) == 0;
}

bool LEGACY::isPowerBy2(ZZ n)
{
    return n > 0 && (n & n - 1) == 0;
}

bool LEGACY::isPrime(long long n)
{
    if(n==1)
        return 0;
    long long i=2;
    for(; i*i<=n; i++)
    {
        if(n%i==0)
        {
            return 0;
        }
    }
    return 1;
}

bool LEGACY::isPrime(ZZ n)
{
    if(n==1)
        return 0;
    long long i=2;
    for(; i*i<=n; i++)
    {
        if(n%i==0)
        {
            return 0;
        }
    }
    return 1;
}

long long LEGACY::Prime(long long i)    /*回傳質數*/ 
{
    if(isPrime(i))
	    return i;
	else
	    return 1;
}

ZZ LEGACY::Prime(ZZ i)    /*回傳質數*/ 
{
    if(isPrime(i))
	    return i;
	else
	    return (ZZ)1;
}

long long LEGACY::Factorize(long long *factor, long long num)   /*列印標準分解式*/
{
	long long a=num,i=2,k;
	long long cnt = 0;
	long long factor_tmp[20];
	
	if(!isPowerBy2(num)) {
		//cout << " initail = " << num << endl;
		while(a!=1)
		{
			if(a%Prime(i)==0 && Prime(i)!=1)
			{
			  k=Prime(i);
			  a/=k;
			  //printf("  %d",k);
			  //printf("\n"); 		  
			  factor_tmp[cnt] = k;
			  factor[cnt] = k;			  
			  cnt++;    
			  //if(a!=1)
				  //printf("\n");    
			}
			else 
			  i++;
		}
		/*
		for (int j = 0; j < cnt; j++){
			if(isPowerBy2(factor_tmp[j]) == 0){
				if(isPrime(factor_tmp[j]))
					Factorize(factor_tmp[j] - 1);
				else 
					Factorize(factor_tmp[j]);
			}
		}*/
	}
	
	

		return cnt;
}

long long LEGACY::Factorize(ZZ *factor, ZZ num)   /*列印標準分解式*/
{
	ZZ a=num,i=(ZZ)2,k;
	long long cnt = 0;
	ZZ factor_tmp[20];
	
	if(!isPowerBy2(num)) {
		//cout << " initail = " << num << endl;
		while(a!=1)
		{
			if(a%Prime(i)==0 && Prime(i)!=1)
			{
			  k=Prime(i);
			  a/=k;
			  //printf("  %d",k);
			  //printf("\n"); 		  
			  factor_tmp[cnt] = k;
			  factor[cnt] = k;			  
			  cnt++;    
			  //if(a!=1)
				  //printf("\n");    
			}
			else 
			  i++;
		}
	}
		return cnt;
}


//   m   | modular - 1
//   2^n | modular - 1

long long LEGACY::find_prime(long long m, long long powerof2){
	bool flag = 0;
	bool powerof2_flag = 0;
	bool prime_flag = 0;	
	long long i = 0;
	long long tmp ;
	long long init = m * pow(2,powerof2);
	while(flag == 0){
		i++ ;
		tmp = init * i;
		powerof2_flag = 0;

		prime_flag = 0;
		prime_flag = isPrime(tmp+1);
		if(prime_flag == 1){		
			flag = 1;
		}		
		
	}
	return tmp+1 ;
	
}

long long LEGACY::find_inv(long long data_in, long long modular)
{
    long long inv;
    
    for(inv=1; inv<modular; inv++)
    {
        if(((data_in * inv) % modular) == 1)
        {
            break;
        }
    }
    
    if(inv == modular)
    {
        return 0;
    }
    
    return inv;
}

ZZ LEGACY::find_inv(ZZ data_in, ZZ modular)
{
    ZZ inv;    
	PowerMod(inv, data_in, (modular-2), modular); //by fermat little theorem
    
    return inv;
}

long long LEGACY::find_prou(long long m, long long modular)
{    
    //output = primitive root of unity
    long long i,j;
    long long prou;
    long long prou_temp;
  
    for(j=2;j<modular;j++)
    {
        if((modular % j) == 0)
        {
            printf("modular is no prime\n");
            return 0;
        }
    }

    for(prou=2;prou<modular;prou++)
    {
        prou_temp = 1;
        for(i=1;i<m;i++)
        {
            prou_temp *= prou;
            prou_temp %= modular;
            if(prou_temp == 1)
            {
                break;
            }
        }
        
        prou_temp *= prou;
        prou_temp %= modular;
        if(prou_temp == 1)
        {
            break;
        }
    }
    
    if(prou == modular)
    {
        return 0;
    }
    else
    {
        return prou;
    }
}


ZZ LEGACY::find_n_rou(ZZ base, long long m, ZZ modular) // a^(p-1) = 1 (mod p)  ---> base^(modular-1) = 1 (mod modular)
{
	assert(( modular % m ) == 1);
	ZZ i;
	ZZ n_rou;
	i = modular/m ;   // base^(modular - 1) = base^( n * i ) = (base^i)^n = 1 (mod modular)
	PowerMod(n_rou, base, i, modular);
	//cout << " n_rou = " << n_rou << endl;
	return n_rou;
}

bool LEGACY::check_prou(ZZ n_rou, long long m, ZZ modular){ //check if n_rou^1, n_rou^2,...,n_rou^(m-1) is not equal 1;
	bool is_prou = true;
	ZZ tmp;
	for(int i = 1; i < m; i++){
		PowerMod(tmp, n_rou, i, modular);
		if(tmp == 1){
			is_prou = false;
			break;
		}
	}
	return is_prou;
}

ZZ LEGACY::find_prou(long long m, ZZ modular)
{   
	bool is_prou = false;
	ZZ i = (ZZ)2 ;
	ZZ n_rou;
	ZZ prou;
	while(is_prou == false)
	{
		n_rou = find_n_rou(i, m, modular);
		is_prou = check_prou(n_rou, m, modular);
		i = i + 1;
	}
	prou = n_rou;
	return prou;
}

void LEGACY::find_zmstar(long long *zmstar, long long m)
{
	int j = 0;
	for(int i = 0; i < m; i++){
		if(coprime(i,m)){
			zmstar[j] = i;
			j++;
		}
	}
}	


long long LEGACY::prou_power(long long data_in, long long power, long long modular)
{
    //output = data_in^power % modular
    long long ans;
    
    ans = 1;
    
    if(power >= 2)
    {
        if((power % 2) == 1)
        {
            ans = prou_power(data_in, (power - 1)/2, modular);
            ans *= ans;
            ans %= modular;
            ans *= data_in;
            ans %= modular;
        }
        else
        {
            ans = prou_power(data_in, power/2, modular);
            ans *= ans;
            ans %= modular;
        }
    }
    else if(power == 1)
    {
        ans = data_in;
    }
    
    return ans;
}

long long LEGACY::DFT(long long *DFT_data, long long *data_in, long long m, long long prou, long long modular) //primitive root of unity in m-point DFT
{
	long long DFT_data_tmp[m];
    long long i, j, prou_tmp;
	
	long long check_m;// check if m | modular - 1
	check_m = (modular-1) % m ;
	assert(check_m == 0) ;
	
	//cout << "DFT_in = " << endl;	
    for(i=0;i<m;i++)
    {
        DFT_data_tmp[i] = 0;
		//cout << data_in[i] << endl ;		
    }
	//cout << endl; 

 
    for(i=0;i<m;i++)
    {
        prou_tmp = prou_power(prou, i, modular);
        for(j=m-1;j>0;j--)
        {
            DFT_data_tmp[i] += data_in[j];
            DFT_data_tmp[i] *= prou_tmp;
            DFT_data_tmp[i] %= modular;
        }
        DFT_data_tmp[i] += data_in[0];
        DFT_data_tmp[i] %= modular;
    } 	
    
	//cout << "DFT_out = " << endl;
    for(i=0;i<m;i++)
    {
        DFT_data[i] = DFT_data_tmp[i];
		//cout << DFT_data[i] << endl ;
    }
	//cout << endl;
    
	return 0;
}

void LEGACY::DFT(ZZ *DFT_data, ZZ *data_in, long long m, ZZ prou, ZZ modular) //primitive root of unity in m-point DFT
{
	ZZ DFT_data_tmp[m];
    ZZ prou_tmp;
	
	ZZ check_m;// check if m | modular - 1
	check_m = (modular-1) % m ;
	assert(check_m == 0) ;
	
	//cout << "DFT_in = " << endl;	
    for(int i = 0; i < m; i++)
    {
        DFT_data_tmp[i] = 0;
		//cout << data_in[i] << endl ;		
    }
	//cout << endl; 

 
    for(int i = 0; i < m; i++)
    {
        //prou_tmp = prou_power(prou, i, modular);
		PowerMod(prou_tmp, prou, i, modular);
        for(int j = m - 1; j > 0; j--)
        {
            // DFT_data_tmp[i] += data_in[j];
            // DFT_data_tmp[i] *= prou_tmp;
            // DFT_data_tmp[i] %= modular;
			AddMod(DFT_data_tmp[i], DFT_data_tmp[i], data_in[j], modular);
			MulMod(DFT_data_tmp[i], DFT_data_tmp[i], prou_tmp, modular);	
        }
        // DFT_data_tmp[i] += data_in[0];
        // DFT_data_tmp[i] %= modular;
		AddMod(DFT_data_tmp[i], DFT_data_tmp[i], data_in[0], modular);		
    } 	
    
	//cout << "DFT_out = " << endl;
    for(int i = 0; i < m ; i++)
    {
        DFT_data[i] = DFT_data_tmp[i];
		//cout << DFT_data[i] << endl ;
    }
	//cout << endl;
    
}

long long LEGACY::IDFT(long long *IDFT_data, long long *data_in, long long n, long long prou, long long modular)
{
	long long prou_inv;
	long long n_inv;
	long long IDFT_tmp[n];
	
	prou_inv = find_inv(prou, modular);
	n_inv = find_inv(n, modular);	
	DFT(IDFT_tmp, data_in, n, prou_inv, modular);
	
	for (int i = 0; i< n ; i++){	
		IDFT_data[i] = ( n_inv * IDFT_tmp[i] ) % modular ;
	}	
	
	return 0;	
	
}

long long LEGACY::FFT(long long *DFT_data, long long *data_in, long long n, long long prou, long long modular) //primitive root of unity in n-point FFT
{
	long long DFT_data_tmp_1[n];
	long long DFT_data_tmp_2[n];
    long long two_to_i, ind_j;
    long long i, j, k;
	
	long long check_n;// check if m | modular - 1
	check_n = (modular-1) % n ;
	assert(check_n == 0) ;    
	
    for(j=0;j<n;j++)
    {
        DFT_data_tmp_1[j] = data_in[j];
    }
        
    for(i=0;i<log2(n);i++)
    {
        two_to_i=1<<i;
        for(k=0;k<two_to_i;k++)
        {
            for(j=0;j<((n/two_to_i)/2);j++)
            {
                ind_j = j + k * (n/two_to_i);
                //BU2 up output
                DFT_data_tmp_2[ind_j] = DFT_data_tmp_1[ind_j] + DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)];
                DFT_data_tmp_2[ind_j] %= modular;
                //BU2 down output
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] = DFT_data_tmp_1[ind_j] - DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)];
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] *= prou_power(prou, j * two_to_i, modular);
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] %= modular;
            }
        }
        for(j=0;j<n;j++)
        {
            DFT_data_tmp_1[j] = DFT_data_tmp_2[j];
            DFT_data_tmp_1[j] %= modular;
        }
    } 	
    
    //output index
    for(i=0;i<n;i++)
    {
        ind_j = 0;
        for(k=0;k<log2(n);k++)
        {
           if(((i >> k) & (long long)1) == (long long)1)
           {
               ind_j |= (1 << (int)(log2(n) - k - 1));
           }
        }

        DFT_data[ind_j] = DFT_data_tmp_1[i]; //deal with negative
        if(DFT_data[ind_j] < 0)
        {
        	DFT_data[ind_j] += modular;
        }
    }
    
	return 0;
}
/*
void LEGACY::FFT(ZZ *DFT_data, ZZ *data_in, ZZ n, ZZ prou, ZZ modular) //primitive root of unity in n-point FFT
{
	ZZ DFT_data_tmp_1[n];
	ZZ DFT_data_tmp_2[n];
    ZZ two_to_i, ind_j;
    ZZ i, j, k;
	
	ZZ check_n;// check if m | modular - 1
	check_n = (modular-1) % n ;
	assert(check_n == 0) ;    
	
    for(j=0;j<n;j++)
    {
        DFT_data_tmp_1[j] = data_in[j];
    }
        
    for(i=0;i<log2(n);i++)
    {
        two_to_i=1<<i;
        for(k=0;k<two_to_i;k++)
        {
            for(j=0;j<((n/two_to_i)/2);j++)
            {
                ind_j = j + k * (n/two_to_i);
                //BU2 up output
                DFT_data_tmp_2[ind_j] = DFT_data_tmp_1[ind_j] + DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)];
                DFT_data_tmp_2[ind_j] %= modular;
                //BU2 down output
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] = DFT_data_tmp_1[ind_j] - DFT_data_tmp_1[ind_j + ((n/two_to_i)/2)];
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] *= prou_power(prou, j * two_to_i, modular);
                DFT_data_tmp_2[ind_j + ((n/two_to_i)/2)] %= modular;
            }
        }
        for(j=0;j<n;j++)
        {
            DFT_data_tmp_1[j] = DFT_data_tmp_2[j];
            DFT_data_tmp_1[j] %= modular;
        }
    } 	
    
    //output index
    for(i=0;i<n;i++)
    {
        ind_j = 0;
        for(k=0;k<log2(n);k++)
        {
           if(((i >> k) & (ZZ)1) == (ZZ)1)
           {
               ind_j |= (1 << (int)(log2(n) - k - 1));
           }
        }

        DFT_data[ind_j] = DFT_data_tmp_1[i]; //deal with negative
        if(DFT_data[ind_j] < 0)
        {
        	DFT_data[ind_j] += modular;
        }
    }
    
	return 0;
}
*/
long long LEGACY::IFFT(long long *IDFT_data, long long *data_in, long long n, long long prou, long long modular) //primitive root of unity in n-point FFT
{
	long long prou_inv;
	long long n_inv;
	long long IFFT_tmp[n];
	
	prou_inv = find_inv(prou, modular);
	FFT(IFFT_tmp, data_in, n, prou_inv, modular);
	n_inv = find_inv(n, modular);
	
	for (int i = 0; i< n ; i++){	
		IDFT_data[i] = ( n_inv * IFFT_tmp[i] ) % modular ;
	}	
	
	return 0;
}
long long LEGACY::PFA2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long inv1, long long inv2, long long prou, long long modular) //primitive root of unity in m(= m1 * m2)-point PFA
{
    long long n1,n2,k1,k2;
    long long index_in;
    long long index_out;
    long long index_m;
    
    long long data_tmp[m1*m2];
    
    //cout << "index_3" <<endl;
    //PFA with m = m1 * m2 -> m2 * m1-DFT + m1 * m2-DFT
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1; // 0 1 2 3 .... m
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (n1 * m2 + n2 * m1) % (m1 * m2);//input re-index
           
            DFT_data[index_m] = data_in[index_in];
			//cout << index_in << endl;
        }    
    }
    
    for(n2=0;n2<m2;n2++)
    {
        //m2 times, m1-point DFT
        DFT(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, prou_power(prou, m2, modular), modular);           
    }    
    //cout << "index_4" <<endl;    
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (k1 + k2 * m1) % (m1 * m2); //re-index between m1 and m2 
            
            DFT_data[index_m] = data_tmp[index_in];
			//cout << index_in << endl;
        }    
    }
   
    {
        for(n1=0;n1<m1;n1++)
        {
            //m2-point DFT
            DFT(data_tmp + m2 * n1, DFT_data + m2 * n1, m2, prou_power(prou, m1, modular), modular);           
        }
    }
         
    //cout << "index_5" <<endl;		 
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_out = (k1 * inv2 * m2 + k2 * inv1 * m1) % (m1 * m2); //output re-index
            
            DFT_data[index_out] = data_tmp[index_m];
			//cout << index_out << endl;
        }    
    }
}


long long LEGACY::PFA2_v2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long inv1, long long inv2, long long prou, long long modular) //primitive root of unity in m(= m1 * m2)-point PFA
{
    long long n1,n2,k1,k2;
    long long index_in;
    long long index_out;
    long long index_m;
    
    long long data_tmp[m1*m2];
    
    
    //PFA with m = m1 * m2 -> m2 * m1-DFT + m1 * m2-DFT
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
			
            index_m = n1 + n2 * m1; // 0 1 2 3 .... m
		    /*
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (n1 * m2 + n2 * m1) % (m1 * m2);//input re-index
            */
            DFT_data[index_m] = data_in[index_m];
        }    
    }
    
    for(n2=0;n2<m2;n2++)
    {
        //m2 times, m1-point DFT
        DFT(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, prou_power(prou, m2, modular), modular);  
		//Rader(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, prou_power(prou, m2, modular), modular); 		
    }    
    
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (k1 + k2 * m1) % (m1 * m2); //re-index between m1 and m2 
            
            DFT_data[index_m] = data_tmp[index_in];
        }    
    }
   
    {
        for(n1=0;n1<m1;n1++)
        {
            //m2-point DFT
            DFT(data_tmp + m2 * n1, DFT_data + m2 * n1, m2, prou_power(prou, m1, modular), modular);    
			//Rader(data_tmp + m2 * n1, DFT_data + m2 * n1, m2, prou_power(prou, m1, modular), modular);    	
        }
    }
            
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_out = (k1 * inv2 * m2 + k2 * inv1 * m1) % (m1 * m2); //output re-index
            
            DFT_data[index_m] = data_tmp[index_m];
			//cout << index_out <<" ";
        }    
    }
}

long long LEGACY::PFA2_v3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long m1_prime, long long m2_prime, long long m1_prime_prou, long long m2_prime_prou, long long inv1, long long inv2, long long *tw_FFT_out1, long long *tw_FFT_out2, long long *index_out1, long long *index_out2, long long modular, long long m1_prou,long long m2_prou,long long *index_in1,long long *index_in2) //primitive root of unity in m(= m1 * m2)-point PFA
{
    long long n1,n2,k1,k2;
    long long index_in;
    long long index_out;
    long long index_m;
    
    long long data_tmp[m1*m2];
    
    
    //PFA with m = m1 * m2 -> m2 * m1-DFT + m1 * m2-DFT
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
			
            index_m = n1 + n2 * m1; // 0 1 2 3 .... m
		    
            //k2 = index_m % m2;
            //k1 = index_m / m2;
            //index_in = (n1 * m2 + n2 * m1) % (m1 * m2);//input re-index
            
            DFT_data[index_m] = data_in[index_m];
        }    
    }
 
/*
    cout << "m1 = " << m1 << endl;
    cout << "m1_prime = " << m1_prime << endl;
    cout << "m1_prime_prou = " << m1_prime_prou << endl;
*/
 
    for(n2=0;n2<m2;n2++)
    {
        //m2 times, m1-point DFT
        //DFT(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, m1_prou, modular);                  
	    //Rader(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, m1_prou, modular); 
		Rader_v3(data_tmp + m1 * n2, DFT_data + m1 * n2, tw_FFT_out1 ,index_in1, index_out1 , m1, m1_prime, m1_prime_prou, modular);           
    }    
    
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (k1 + k2 * m1) % (m1 * m2); //re-index between m1 and m2 
            
            DFT_data[index_m] = data_tmp[index_in];
        }  
    }

        for(n1=0;n1<m1;n1++)
        {
            //m2-point DFT
            //Rader(data_tmp + m2 * n1, DFT_data + m2 * n1, m2, m2_prou, modular);
			//DFT(data_tmp + m2 * n1, DFT_data + m2 * n1, m2, m2_prou, modular);  
			Rader_v3(data_tmp + m2 * n1, DFT_data + m2 * n1, tw_FFT_out2 ,index_in2, index_out2 , m2, m2_prime, m2_prime_prou, modular); 			
        }
    
            
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_out = (k1 * inv2 * m2 + k2 * inv1 * m1) % (m1 * m2); //output re-index
            
            DFT_data[index_m] = data_tmp[index_m];
        }    
    }
}


long long LEGACY::PFA2_v4(long long *DFT_data, long long *data_in, long long m1, long long m2, long long m1_prime, long long m2_prime, long long m1_prime_prou, long long m2_prime_prou, long long inv1, long long inv2, long long *tw_FFT_out1, long long *tw_FFT_out2, long long *index_out1, long long *index_out2, long long modular, long long m1_prou,long long m2_prou,long long *index_in1,long long *index_in2) //primitive root of unity in m(= m1 * m2)-point PFA
{
    long long n1,n2,k1,k2;
    long long index_in;
    long long index_out;
    long long index_m;
    
    long long data_tmp[m1*m2];
    
    
    //PFA with m = m1 * m2 -> m2 * m1-DFT + m1 * m2-DFT
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
			
            index_m = n1 + n2 * m1; // 0 1 2 3 .... m
		    
            //k2 = index_m % m2;
            //k1 = index_m / m2;
            //index_in = (n1 * m2 + n2 * m1) % (m1 * m2);//input re-index
            
            DFT_data[index_m] = data_in[index_m];
        }    
    }
 
/*
    cout << "m1 = " << m1 << endl;
    cout << "m1_prime = " << m1_prime << endl;
    cout << "m1_prime_prou = " << m1_prime_prou << endl;
*/
 
    for(n2=0;n2<m2;n2++)
    {
        //m2 times, m1-point DFT
        //DFT(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, m1_prou, modular);                  
	    //Rader(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, m1_prou, modular); 
		Rader_v4(data_tmp + m1 * n2, DFT_data + m1 * n2, tw_FFT_out1 ,index_in1, index_out1 , m1, m1_prime, m1_prime_prou, modular);           
    }    
    
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (k1 + k2 * m1) % (m1 * m2); //re-index between m1 and m2 
            
            DFT_data[index_m] = data_tmp[index_in];
        }  
    }

        for(n1=0;n1<m1;n1++)
        {
            //m2-point DFT
            //Rader(data_tmp + m2 * n1, DFT_data + m2 * n1, m2, m2_prou, modular);
			//DFT(data_tmp + m2 * n1, DFT_data + m2 * n1, m2, m2_prou, modular);  
			Rader_v4(data_tmp + m2 * n1, DFT_data + m2 * n1, tw_FFT_out2 ,index_in2, index_out2 , m2, m2_prime, m2_prime_prou, modular); 			
        }
    
            
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_out = (k1 * inv2 * m2 + k2 * inv1 * m1) % (m1 * m2); //output re-index
            
            DFT_data[index_m] = data_tmp[index_m];
        }    
    }
}


long long LEGACY::PFA3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) //primitive root of unity in m(= m1(= s_m1 * s_m2) * m2)-point PFA
{
    long long n1,n2,k1,k2;
    long long index_in;
    long long index_out;
    long long index_m;
    
    long long data_tmp[m1*m2];
    
    //cout << "index_1 " << endl;
	
    //PFA with m = m1 * m2 -> m2 * m1-DFT + m1 * m2-DFT
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1; // 1 2 3 ... m
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (n1 * m2 + n2 * m1) % (m1 * m2);// 1
            
            DFT_data[index_m] = data_in[index_in];	
			//cout<< index_in << endl;
        }    
    }


    for(n2=0;n2<m2;n2++)
    {
        //m2 times, m1-point DFT
        DFT(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, prou_power(prou, m2, modular), modular);           
    }    

	//cout << "index_2" <<endl; 
	
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (k1 + k2 * m1) % (m1 * m2);
            
            DFT_data[index_m] = data_tmp[index_in];
			//cout<< index_in << endl;
        }    
    }
 	//cout << "----------------------" <<endl;  
    {
		for(n2=0;n2<m1;n2++)
		{
			//m1(= s_m1 * s_m2)-point PFA
			PFA2(data_tmp + m2 * n2, DFT_data + m2 * n2, s_m1, s_m2, s_inv1, s_inv2, prou_power(prou, m1, modular), modular);         
		}  
    }
        
	//cout << "index_6" <<endl;		
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_out = (k1 * inv2 * m2 + k2 * inv1 * m1) % (m1 * m2);
            
            DFT_data[index_out] = data_tmp[index_m];
			//cout<< index_out << endl;
        }    
    }
	 	//cout << "----------------------" <<endl; 
}

//combine the first and second decomposition input index at the beginnig 

long long LEGACY::PFA3_v2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) //primitive root of unity in m(= m1(= s_m1 * s_m2) * m2)-point PFA
{
    long long n1,n2,k1,k2,k3;
    long long index_in;
    long long index_out;
    long long index_m;
    
    long long data_tmp[m1*m2];
    
    
    //PFA with m = m1 * m2 -> m2 * m1-DFT + m1 * m2-DFT
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1; // 1 2 3 ... m
			k3 = (index_m / (m1*s_m1)) % s_m2 ;
            k2 = (index_m / m1) % s_m1 ;
		    k1 = index_m % m1;
            index_in = (k1 * m2 + k2 * s_m2*m1 + k3 * s_m1*m1) % (m1 * m2);// 1
            
            DFT_data[index_m] = data_in[index_in];	
			//cout<< index_in << endl;
        }    
    }
	//cout << "----------------------" <<endl;

    for(n2=0;n2<m2;n2++)
    {
        //m2 times, m1-point DFT
        DFT(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, prou_power(prou, m2, modular), modular);           
    }    
	
	//cout << "rader ok" << endl;
 
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (k1 + k2 * m1) % (m1 * m2);
            
            DFT_data[index_m] = data_tmp[index_in];
			//cout<< index_in << endl;
        }    
    }
 	//cout << "----------------------" <<endl;  
    {
		for(n2=0;n2<m1;n2++)
		{
			//m1(= s_m1 * s_m2)-point PFA
			PFA2_v2(data_tmp + m2 * n2, DFT_data + m2 * n2, s_m1, s_m2, s_inv1, s_inv2, prou_power(prou, m1, modular), modular);         
		}  
    }
	
    //merge final two reindex
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
			//k2 = index_m % m2;
            //k1 = index_m / m2;
            //index_out = (k1 * inv2 * m2 + k2 * inv1 * m1) % (m1 * m2);
			k1 = (index_m / m2) % m1 ;
            k2 = (index_m / s_m2) % s_m1 ;
		    k3 = index_m % s_m2;
            index_out = (k1 * inv2 * m2 + (k2 * s_m2 * s_inv2 + k3 * s_m1 * s_inv1)*m1*inv1 ) % (m1 * m2);
            DFT_data[index_out] = data_tmp[index_m];
			//cout<< index_out << endl;
        }    
    }
	 	//cout << "----------------------" <<endl; 
}

long long LEGACY::PFA3_v3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) //primitive root of unity in m(= m1(= s_m1 * s_m2) * m2)-point PFA
{
    long long n1,n2,k1,k2,k3;
    long long index_in;
    long long index_out;
    long long index_m;
    
    long long data_tmp[m1*m2];
	long long tmp;
//------------------------ rader parameter m1 set up-------------------------
	long long m1_prime, gen1;
	long long rader_index_in1[m1];
	long long rader_index_out1[m1];	
	gen1 = find_gen(m1);
	if( (m1==3) || (m1==5) || (m1==17) || (m1==257)){
		m1_prime = m1-1;
	}
	else{
		m1_prime = find_m_prime(m1);
	}
	
	long long tw_FFT_index1[m1_prime];
	long long tw_FFT_in1[m1_prime];	
	long long m1_prou;	
	
	m1_prou = prou_power(prou, m2, modular);
	rader_index_in1[0] = 0;
	rader_index_out1[0] = 0;
	rader_index_in1[1] = 1;
	rader_index_out1[1] = 1;
	tmp = 1;
	for (int i = 2; i < m1 ; i++){
		tmp *= gen1;
		tmp %= m1;
		rader_index_in1[i] = tmp;
	}	
	for (int i = 2; i < m1 ; i++){
		rader_index_out1[i] = rader_index_in1[m1-i+1];
	}
	
	for (int i = 0; i < m1_prime ; i++){
		if(i < (m1 - 1) )//0-5
			tw_FFT_index1[i] = rader_index_out1[i+1];
		else if(i > (m1_prime - m1 + 1))
			tw_FFT_index1[i] = rader_index_out1[i + m1 - m1_prime];			
		else		
			tw_FFT_index1[i] = 0; //00000
		//cout << tw_FFT_index1[i] << endl;	
	}
	
	for (int i = 0; i < m1_prime ; i++){
		if(tw_FFT_index1[i] == 0)
			tw_FFT_in1[i] = 0;
		else
			tw_FFT_in1[i] = prou_power(m1_prou,tw_FFT_index1[i],modular);
	}
	//cout << "prou = " << prou << endl;
	//cout << "m1_prou = " << m1_prou << endl;	
	//cout << "tw_FFT_in1 = " << endl;
	/*for (int i = 0; i < m1_prime ; i++){
			//cout << tw_FFT_in1[i] <<endl;
	}*/
	
	long long tw_FFT_out1[m1_prime];
	long long m1_prime_prou;

	m1_prime_prou = find_prou(m1_prime, modular);
	FFT(tw_FFT_out1, tw_FFT_in1, m1_prime, m1_prime_prou, modular) ;	
	//cout << "tw_FFT_out1 = " << endl;
	/*for (int i = 0; i < m1_prime ; i++){
			//cout << tw_FFT_out1[i] <<endl;
	}*/	
//-----------------------
	int m2_prime, gen2;
	long long rader_index_in2[s_m1];
	long long rader_index_out2[s_m1];	
	gen2 = find_gen(s_m1);	
	if( (s_m1==3) || (s_m1==5) || (s_m1==17) || (s_m1==257)){
		m2_prime = s_m1-1;
	}
	else{
		m2_prime = find_m_prime(s_m1);
	}
	long long tw_FFT_index2[m2_prime];	
	long long tw_FFT_in2[m2_prime];	
	long long m2_prou;	
	
	m2_prou = prou_power(prou, m1*s_m2, modular);	
	rader_index_in2[0] = 0;
	rader_index_out2[0] = 0;
	rader_index_in2[1] = 1;
	rader_index_out2[1] = 1;
	tmp = 1;
	for (int i = 2; i < s_m1 ; i++){
		tmp *= gen2;
		tmp %= s_m1;
		rader_index_in2[i] = tmp;	
	}	
	for (int i = 2; i < s_m1 ; i++){
		rader_index_out2[i] = rader_index_in2[s_m1-i+1];
	}
	for (int i = 0; i < m2_prime ; i++){
		if(i < (s_m1 - 1) )//0-5
			tw_FFT_index2[i] = rader_index_out2[i+1];
		else if(i > (m2_prime - s_m1 + 1))
			tw_FFT_index2[i] = rader_index_out2[i + s_m1 - m2_prime];			
		else		
			tw_FFT_index2[i] = 0; //00000
		//cout << tw_FFT_index2[i] << endl;	
	}
	
	for (int i = 0; i < m2_prime ; i++){
		if(tw_FFT_index2[i] == 0)
			tw_FFT_in2[i] = 0;
		else
			tw_FFT_in2[i] = prou_power(m2_prou,tw_FFT_index2[i],modular);
	}
	long long tw_FFT_out2[m2_prime];
	long long m2_prime_prou;
	m2_prime_prou = find_prou(m2_prime, modular);
	FFT(tw_FFT_out2, tw_FFT_in2, m2_prime, m2_prime_prou, modular) ;	
	
//----------	
	int m3_prime, gen3;
	long long rader_index_in3[s_m2];
	long long rader_index_out3[s_m2];
	gen3 = find_gen(s_m2);
	if( (s_m2==3) || (s_m2==5) || (s_m2==17) || (s_m2==257)){
		m3_prime = s_m2-1;
	}
	else{
		m3_prime = find_m_prime(s_m2);
	}
	long long tw_FFT_index3[m3_prime];
	long long tw_FFT_in3[m3_prime];
	long long m3_prou;	
	
	m3_prou = prou_power(prou, m1*s_m1, modular);	
	rader_index_in3[0] = 0;
	rader_index_out3[0] = 0;
	rader_index_in3[1] = 1;
	rader_index_out3[1] = 1;
	tmp = 1;
	for (int i = 2; i < s_m2 ; i++){
		tmp *= gen3;
		tmp %= s_m2;
		rader_index_in3[i] = tmp;	
	}
	for (int i = 2; i < s_m2 ; i++){
		rader_index_out3[i] = rader_index_in3[s_m2-i+1];
	}
	for (int i = 0; i < m3_prime ; i++){
		if(i < (s_m2 - 1) )//0-5
			tw_FFT_index3[i] = rader_index_out3[i+1];
		else if(i > (m3_prime - s_m2 + 1))
			tw_FFT_index3[i] = rader_index_out3[i + s_m2 - m3_prime];			
		else		
			tw_FFT_index3[i] = 0; //00000
		//cout << tw_FFT_index3[i] << endl;	
	}	
	
	for (int i = 0; i < m3_prime ; i++){
		if(tw_FFT_index3[i] == 0)
			tw_FFT_in3[i] = 0;
		else
			tw_FFT_in3[i] = prou_power(m3_prou,tw_FFT_index3[i],modular);
	}	
	long long tw_FFT_out3[m3_prime];
	long long m3_prime_prou;
	m3_prime_prou = find_prou(m3_prime, modular);
	FFT(tw_FFT_out3, tw_FFT_in3, m3_prime, m3_prime_prou, modular) ;
//-----------------------------------------------------------------------------

    
    
    //PFA with m = m1 * m2 -> m2 * m1-DFT + m1 * m2-DFT
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1; // 1 2 3 ... m
			k3 = (index_m / (m1*s_m1)) % s_m2 ;
            k2 = (index_m / m1) % s_m1 ;
		    k1 = index_m % m1;
            index_in = (k1 * m2 + k2 * s_m2*m1 + k3 * s_m1*m1) % (m1 * m2);// 1
            
            DFT_data[index_m] = data_in[index_in];	
			//cout<< index_in << endl;
        }    
    }
	//cout << "----------------------" <<endl;
		//cout << "prou = " << prou << endl;
    for(n2=0;n2<m2;n2++)
    {
        //m2 times, m1-point DFT
        //Rader(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, prou_power(prou, m2, modular), modular); 
		//DFT(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, prou_power(prou, m2, modular), modular); 
		Rader_v3(data_tmp + m1 * n2, DFT_data + m1 * n2, tw_FFT_out1, rader_index_in1, rader_index_out1, m1, m1_prime, m1_prime_prou, modular);
		
    }    
	
	//cout << "rader ok" << endl;
 
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (k1 + k2 * m1) % (m1 * m2);
            
            DFT_data[index_m] = data_tmp[index_in];
			//cout<< index_in << endl;
        }    
    }
 	//cout << "----------------------" <<endl;  
    {
		for(n2=0;n2<m1;n2++)
		{
			//m1(= s_m1 * s_m2)-point PFA
			//PFA2_v2(data_tmp + m2 * n2, DFT_data + m2 * n2, s_m1, s_m2, s_inv1, s_inv2, prou_power(prou, m1, modular), modular); 
			PFA2_v3(data_tmp + m2 * n2, DFT_data + m2 * n2, s_m1, s_m2, m2_prime, m3_prime, m2_prime_prou, m3_prime_prou, s_inv1, s_inv2, tw_FFT_out2, tw_FFT_out3, rader_index_out2, rader_index_out3, modular, m2_prou, m3_prou,rader_index_in2,rader_index_in3); 			
		}  
    }
	
    //merge final two reindex
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
			//k2 = index_m % m2;
            //k1 = index_m / m2;
            //index_out = (k1 * inv2 * m2 + k2 * inv1 * m1) % (m1 * m2);
			k1 = (index_m / m2) % m1 ;
            k2 = (index_m / s_m2) % s_m1 ;
		    k3 = index_m % s_m2;
            index_out = (k1 * inv2 * m2 + (k2 * s_m2 * s_inv2 + k3 * s_m1 * s_inv1)*m1*inv1 ) % (m1 * m2);
            DFT_data[index_out] = data_tmp[index_m];
			//cout<< index_out << endl;
        }    
    }
	 	//cout << "----------------------" <<endl; 
}


long long LEGACY::PFA3_v4(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) //primitive root of unity in m(= m1(= s_m1 * s_m2) * m2)-point PFA
{
    long long n1,n2,n3,k1,k2,k3;
    long long index_in;
    long long index_out;
    long long index_m;
    
    long long data_tmp[m1*m2];
	long long tmp;
	std::ofstream test("./test.txt");
//------------------------ rader parameter m1 set up-------------------------
	long long m1_prime, gen1;
	long long rader_index_in1[m1];
	long long rader_index_out1[m1];	
	gen1 = find_gen(m1);
	if( (m1==3) || (m1==5) || (m1==17) || (m1==257)){
		m1_prime = m1-1;
	}
	else{
		m1_prime = find_m_prime(m1);
	}
	
	long long tw_FFT_index1[m1_prime];
	long long tw_FFT_in1[m1_prime];	
	long long m1_prou;	
	
	m1_prou = prou_power(prou, m2, modular);
	rader_index_in1[0] = 0;
	rader_index_out1[0] = 0;
	rader_index_in1[1] = 1;
	rader_index_out1[1] = 1;
	tmp = 1;
	for (int i = 2; i < m1 ; i++){
		tmp *= gen1;
		tmp %= m1;
		rader_index_in1[i] = tmp;
	}	
	for (int i = 2; i < m1 ; i++){
		rader_index_out1[i] = rader_index_in1[m1-i+1];
	}
	
	for (int i = 0; i < m1_prime ; i++){
		if(i < (m1 - 1) )//0-5
			tw_FFT_index1[i] = rader_index_out1[i+1];
		else if(i > (m1_prime - m1 + 1))
			tw_FFT_index1[i] = rader_index_out1[i + m1 - m1_prime];			
		else		
			tw_FFT_index1[i] = 0; //00000
		//cout << tw_FFT_index1[i] << endl;	
	}
	
	for (int i = 0; i < m1_prime ; i++){
		if(tw_FFT_index1[i] == 0)
			tw_FFT_in1[i] = 0;
		else
			tw_FFT_in1[i] = prou_power(m1_prou,tw_FFT_index1[i],modular);
	}
	//cout << "prou = " << prou << endl;
	//cout << "m1_prou = " << m1_prou << endl;	
	//cout << "tw_FFT_in1 = " << endl;
	/*for (int i = 0; i < m1_prime ; i++){
			//cout << tw_FFT_in1[i] <<endl;
	}*/
	
	long long tw_FFT_out1[m1_prime];
	long long m1_prime_prou;

	m1_prime_prou = find_prou(m1_prime, modular);
	FFT(tw_FFT_out1, tw_FFT_in1, m1_prime, m1_prime_prou, modular) ;	
	//cout << "tw_FFT_out1 = " << endl;
	/*for (int i = 0; i < m1_prime ; i++){
			//cout << tw_FFT_out1[i] <<endl;
	}*/	
//-----------------------
	int m2_prime, gen2;
	long long rader_index_in2[s_m1];
	long long rader_index_out2[s_m1];	
	gen2 = find_gen(s_m1);	
	if( (s_m1==3) || (s_m1==5) || (s_m1==17) || (s_m1==257)){
		m2_prime = s_m1-1;
	}
	else{
		m2_prime = find_m_prime(s_m1);
	}
	long long tw_FFT_index2[m2_prime];	
	long long tw_FFT_in2[m2_prime];	
	long long m2_prou;	
	
	m2_prou = prou_power(prou, m1*s_m2, modular);	
	rader_index_in2[0] = 0;
	rader_index_out2[0] = 0;
	rader_index_in2[1] = 1;
	rader_index_out2[1] = 1;
	tmp = 1;
	for (int i = 2; i < s_m1 ; i++){
		tmp *= gen2;
		tmp %= s_m1;
		rader_index_in2[i] = tmp;	
	}	
	for (int i = 2; i < s_m1 ; i++){
		rader_index_out2[i] = rader_index_in2[s_m1-i+1];
	}
	for (int i = 0; i < m2_prime ; i++){
		if(i < (s_m1 - 1) )//0-5
			tw_FFT_index2[i] = rader_index_out2[i+1];
		else if(i > (m2_prime - s_m1 + 1))
			tw_FFT_index2[i] = rader_index_out2[i + s_m1 - m2_prime];			
		else		
			tw_FFT_index2[i] = 0; //00000
		//cout << tw_FFT_index2[i] << endl;	
	}
	
	for (int i = 0; i < m2_prime ; i++){
		if(tw_FFT_index2[i] == 0)
			tw_FFT_in2[i] = 0;
		else
			tw_FFT_in2[i] = prou_power(m2_prou,tw_FFT_index2[i],modular);
	}
	long long tw_FFT_out2[m2_prime];
	long long m2_prime_prou;
	m2_prime_prou = find_prou(m2_prime, modular);
	FFT(tw_FFT_out2, tw_FFT_in2, m2_prime, m2_prime_prou, modular) ;	
	
//----------	
	int m3_prime, gen3;
	long long rader_index_in3[s_m2];
	long long rader_index_out3[s_m2];
	gen3 = find_gen(s_m2);
	if( (s_m2==3) || (s_m2==5) || (s_m2==17) || (s_m2==257)){
		m3_prime = s_m2-1;
	}
	else{
		m3_prime = find_m_prime(s_m2);
	}
	long long tw_FFT_index3[m3_prime];
	long long tw_FFT_in3[m3_prime];
	long long m3_prou;	
	
	m3_prou = prou_power(prou, m1*s_m1, modular);	
	rader_index_in3[0] = 0;
	rader_index_out3[0] = 0;
	rader_index_in3[1] = 1;
	rader_index_out3[1] = 1;
	tmp = 1;
	for (int i = 2; i < s_m2 ; i++){
		tmp *= gen3;
		tmp %= s_m2;
		rader_index_in3[i] = tmp;	
	}
	for (int i = 2; i < s_m2 ; i++){
		rader_index_out3[i] = rader_index_in3[s_m2-i+1];
	}
	for (int i = 0; i < m3_prime ; i++){
		if(i < (s_m2 - 1) )//0-5
			tw_FFT_index3[i] = rader_index_out3[i+1];
		else if(i > (m3_prime - s_m2 + 1))
			tw_FFT_index3[i] = rader_index_out3[i + s_m2 - m3_prime];			
		else		
			tw_FFT_index3[i] = 0; //00000
		//cout << tw_FFT_index3[i] << endl;	
	}	
	
	for (int i = 0; i < m3_prime ; i++){
		if(tw_FFT_index3[i] == 0)
			tw_FFT_in3[i] = 0;
		else
			tw_FFT_in3[i] = prou_power(m3_prou,tw_FFT_index3[i],modular);
	}	
	long long tw_FFT_out3[m3_prime];
	long long m3_prime_prou;
	m3_prime_prou = find_prou(m3_prime, modular);
	FFT(tw_FFT_out3, tw_FFT_in3, m3_prime, m3_prime_prou, modular) ;
//-----------------------------------------------------------------------------

    
    
    //PFA with m = m1 * m2 -> m2 * m1-DFT + m1 * m2-DFT
	for(n3 = 0; n3 < s_m2; n3++)//257     17
	{ 
		for(n2 = 0; n2 < s_m1; n2++)//17   5
		{
			for(n1 = 0; n1 < m1; n1++)//5  3
			{
				index_m = n1 + n2 * m1 + n3 * m1*s_m1; // 1 2 3 ... m
				k1 = rader_index_in1[n1] ;//5
				k2 = rader_index_in2[n2] ;//17
				k3 = rader_index_in3[n3] ;//257
				index_in = (k1 * s_m1 * s_m2 + k2 * m1 * s_m2 + k3 * m1 * s_m1) % (m1 * s_m1 *s_m2);// 1
				
				//index_m = n1 + n2 * m1; // 1 2 3 ... m
				//k3 = (index_m / (m1*m2)) % m3 ;
				//k2 = (index_m / m1) % m2 ;
				//k1 = index_m % m1;
				//index_in = (k1 * m2 + k2 * m3*m1 + k3 * m2*m1) % (m1 * m2);// 1
				DFT_data[index_m] = data_in[index_in];	
				//cout<< index_in << endl;
			}    
		}
	}
	//cout << "----------------------" <<endl;
		//cout << "prou = " << prou << endl;
    for(n2=0;n2<m2;n2++)
    {
        //m2 times, m1-point DFT
        //Rader(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, prou_power(prou, m2, modular), modular); 
		//DFT(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, prou_power(prou, m2, modular), modular); 
		Rader_v4(data_tmp + m1 * n2, DFT_data + m1 * n2, tw_FFT_out1, rader_index_in1, rader_index_out1, m1, m1_prime, m1_prime_prou, modular);
		
    }    
	
	//cout << "rader ok" << endl;
 
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (k1 + k2 * m1) % (m1 * m2);
            
            DFT_data[index_m] = data_tmp[index_in];
			//cout<< index_in << endl;
        }    
    }
 	//cout << "----------------------" <<endl;  
    {
		for(n2=0;n2<m1;n2++)
		{
			//m1(= s_m1 * s_m2)-point PFA
			//PFA2_v2(data_tmp + m2 * n2, DFT_data + m2 * n2, s_m1, s_m2, s_inv1, s_inv2, prou_power(prou, m1, modular), modular); 
			PFA2_v4(data_tmp + m2 * n2, DFT_data + m2 * n2, s_m1, s_m2, m2_prime, m3_prime, m2_prime_prou, m3_prime_prou, s_inv1, s_inv2, tw_FFT_out2, tw_FFT_out3, rader_index_out2, rader_index_out3, modular, m2_prou, m3_prou,rader_index_in2,rader_index_in3); 			
		}  
    }	
	
    //merge final two input reindex
/*	
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
			//k2 = index_m % m2;
            //k1 = index_m / m2;
            //index_out = (k1 * inv2 * m2 + k2 * inv1 * m1) % (m1 * m2);
			k1 = (index_m / m2) % m1 ;
            k2 = (index_m / s_m2) % s_m1 ;
		    k3 = index_m % s_m2;
            index_out = (k1 * inv2 * m2 + (k2 * s_m2 * s_inv2 + k3 * s_m1 * s_inv1)*m1*inv1 ) % (m1 * m2);
            DFT_data[index_out] = data_tmp[index_m];
        }    
    }
*/


	for(n1 = 0; n1 < m1; n1++) //5    3  
	{ 
		for(n2 = 0; n2 < s_m1; n2++)//17    5
		{
			for(n3 = 0; n3 < s_m2; n3++)//257     7
			{
				index_m = n3 + n2 * s_m2 + n1 * s_m2*s_m1; // 1 2 3 ... m
				k1 = rader_index_out1[n1] ;//5
				k2 = rader_index_out2[n2] ;//17
				k3 = rader_index_out3[n3] ;//257
				index_out = (k1 * inv2 * m2 + (k2 * s_m2 * s_inv2 + k3 * s_m1 * s_inv1)*m1*inv1 ) % (m1 * m2);
                DFT_data[index_out] = data_tmp[index_m];	
				
				//index_m = n1 + n2 * m1; // 1 2 3 ... m
				//k3 = (index_m / (m1*m2)) % m3 ;
				//k2 = (index_m / m1) % m2 ;
				//k1 = index_m % m1;
				//index_in = (k1 * m2 + k2 * m3*m1 + k3 * m2*m1) % (m1 * m2);// 1
				//DFT_data[index_m] = data_in[index_in];	
				//test << data_tmp[index_m] << endl;
			}    
		}
	}

	for(int i = 0; i < m1*m2; i++){
		test << DFT_data[i] << endl;
	}




}



long long LEGACY::find_gen(long long n)
{
	bool flag = 0 ;
	long long ans;
	long long tmp = 1;
	for(int i = 2; i < n; i++)
	{
		for(int j = 1; j < n; j++)
		{
			tmp = tmp * i ;
			tmp %= n;
			if(tmp == 1){
				if(j == n - 1)
					flag = 1;
				else
					break ;
			}
		}
		if(flag == 1){
			ans = i;
			break ; 
		}
	}
	return ans ;
}

// input n must be a prime
void LEGACY::Rader(long long *RA_out, long long *data_in, long long n, long long prou, long long modular)
{
	long long check_m;// check if m | modular - 1
	check_m = (modular-1) % n ;
	assert(check_m == 0) ;
	assert(isPrime(n) == 1);
	assert(isPrime(modular) == 1);	
	
	int m_prime ;
	if( (n==3) || (n==5) || (n==17) || (n==257))
		m_prime = n-1;
	else
		m_prime = find_m_prime(n);

	long long gen;
	long long data_in_reindex[n];
	long long tmp = 1;
	long long index_in[n-1];	
	long long index_out[n-1];	
	gen = find_gen(n);

//----------------input re-index-----------------
	index_in[0] = 1;
	index_out[0] = 1;
	data_in_reindex[0] = data_in[0];
	data_in_reindex[1] = data_in[1];	
	for (int i = 0; i < n-2; i++){
		tmp *= gen;
		tmp %= n;
		index_in[i+1] = tmp ;
		data_in_reindex[i+2] = data_in[tmp];
	}
	
	for (int i = 1; i < n-1 ; i++){
		index_out[i] = index_in[n-1-i];
	}

/*
	for (int i = 0; i < n-1 ; i++){
		//cout << index_in[i] <<endl;
		//cout << index_out[i] <<endl;
	}	
*/	
	//cout <<"rader fft in = " << endl ;
	
	long long FFT_in[m_prime];
	for (int i = 0; i < m_prime ; i++){
		if(i < n - 1)
			FFT_in[i] = data_in_reindex[i+1];
		else 
			FFT_in[i] = 0;
		
		//cout << FFT_in[i] <<" ";
	}
	//cout << endl ;	
//-----------------------------------------------
		//cout <<endl;
//----------------tw input re-index-----------------
	long long tw_FFT_index[m_prime];
	long long tw_FFT_in[m_prime];
	for (int i = 0; i < m_prime ; i++){
		if(i < n-1)
			tw_FFT_index[i] = index_out[i];
		else if(i > (m_prime - n + 1))
			tw_FFT_index[i] = index_out[i + n - m_prime -1];
		else
			tw_FFT_index[i] = 0;
	}
		//cout << tw_FFT_index[i] <<endl;
	for (int i = 0; i < m_prime ; i++){
		//cout << tw_FFT_index[i] <<endl;
	}

		//cout <<endl;
		
	for (int i = 0; i < m_prime ; i++){
		if(tw_FFT_index[i] == 0)
			tw_FFT_in[i] = 0;
		else
			tw_FFT_in[i] = prou_power(prou,tw_FFT_index[i],modular);
	}
	//cout << "prou = " << prou << endl;
	//cout << "tw_FFT_in = " << endl ;
	for (int i = 0; i < m_prime ; i++){
		//cout << tw_FFT_in[i] <<endl;
	}
	
//-------------pointwise mul--------------------------
	long long FFT_out[m_prime];
	long long tw_FFT_out[m_prime];
	long long m_prime_prou;
	

	m_prime_prou = find_prou(m_prime, modular);
	//cout << "m_prime_prou = " << m_prime_prou << endl;	
	//cout << "modular = " << modular << endl;	
	FFT(FFT_out, FFT_in, m_prime, m_prime_prou, modular) ;
	FFT(tw_FFT_out, tw_FFT_in, m_prime, m_prime_prou, modular) ;
	
	//cout << "tw_FFT_out = " << endl ;
	for (int i = 0; i < m_prime ; i++){
		//cout << tw_FFT_out[i] <<endl;
		//cout << FFT_out[i] <<" ";
	}
		//cout <<endl;	
	

	long long ele_mul[m_prime];
	for (int i = 0; i< m_prime ; i++){	
		ele_mul[i] = (FFT_out[i]*tw_FFT_out[i]) % modular ;
	}
	for (int i = 0; i < m_prime ; i++){
		//cout << ele_mul[i] <<endl;
	}
//---------------IFFT------------------------
	long long IFFT_out[m_prime];
	IFFT(IFFT_out , ele_mul , m_prime , m_prime_prou , modular);
	/*
	cout << "IFFT_out" <<endl;	
	for (int i = 0; i< n-1 ; i++){	
		cout << IFFT_out[i] << " ";
	}	
	cout  <<endl;*/
	
//-------------add do------------------------
	//long long RA_out[n];
	RA_out[0] = 0;
	for (int i = 0; i< n ; i++){	
		RA_out[0] =  (RA_out[0] + data_in[i]) %  modular;
	}	
	for (int i = 0; i< n-1 ; i++){	
		RA_out[index_out[i]] = (IFFT_out[i] + data_in[0]) %  modular ;
	}	
	
	//cout << "ok" <<endl;
	/*
	cout << "RA_out" <<endl;
	for (int i = 0; i< n ; i++){	
		cout << RA_out[i] << " ";
	}	
	cout  <<endl;
*/

}


void LEGACY::Rader_DFT(long long *RA_out, long long *data_in, long long n, long long prou, long long modular)
{
	long long check_m;// check if m | modular - 1
	check_m = (modular-1) % n ;
	assert(check_m == 0) ;
	assert(isPrime(n) == 1);
	assert(isPrime(modular) == 1);	
	
	int m_prime = n-1;

	long long gen;
	long long data_in_reindex[n];
	long long tmp = 1;
	long long index_in[n-1];	
	long long index_out[n-1];	
	gen = find_gen(n);

//----------------input re-index-----------------
	index_in[0] = 1;
	index_out[0] = 1;
	data_in_reindex[0] = data_in[0];
	data_in_reindex[1] = data_in[1];	
	for (int i = 0; i < n-2; i++){
		tmp *= gen;
		tmp %= n;
		index_in[i+1] = tmp ;
		data_in_reindex[i+2] = data_in[tmp];
	}
	
	for (int i = 1; i < n-1 ; i++){
		index_out[i] = index_in[n-1-i];
	}

	long long FFT_in[m_prime];
	for (int i = 0; i < m_prime ; i++){
		if(i < n - 1)
			FFT_in[i] = data_in_reindex[i+1];
		else 
			FFT_in[i] = 0;
		
		//cout << FFT_in[i] <<" ";
	}
	//cout << endl ;	
//-----------------------------------------------
		//cout <<endl;
//----------------tw input re-index-----------------
	long long tw_FFT_index[m_prime];
	long long tw_FFT_in[m_prime];
	for (int i = 0; i < m_prime ; i++){
		if(i < n-1)
			tw_FFT_index[i] = index_out[i];
		else if(i > (m_prime - n + 1))
			tw_FFT_index[i] = index_out[i + n - m_prime -1];
		else
			tw_FFT_index[i] = 0;
	}
		//cout << tw_FFT_index[i] <<endl;
	for (int i = 0; i < m_prime ; i++){
		//cout << tw_FFT_index[i] <<endl;
	}

		//cout <<endl;
		
	for (int i = 0; i < m_prime ; i++){
		if(tw_FFT_index[i] == 0)
			tw_FFT_in[i] = 0;
		else
			tw_FFT_in[i] = prou_power(prou,tw_FFT_index[i],modular);
	}
	//cout << "prou = " << prou << endl;
	//cout << "tw_FFT_in = " << endl ;
	for (int i = 0; i < m_prime ; i++){
		//cout << tw_FFT_in[i] <<endl;
	}
	
//-------------pointwise mul--------------------------
	long long FFT_out[m_prime];
	long long tw_FFT_out[m_prime];
	long long m_prime_prou;
	

	m_prime_prou = find_prou(m_prime, modular);
	//cout << "m_prime_prou = " << m_prime_prou << endl;	
	//cout << "modular = " << modular << endl;	
	DFT(FFT_out, FFT_in, m_prime, m_prime_prou, modular) ;
	DFT(tw_FFT_out, tw_FFT_in, m_prime, m_prime_prou, modular) ;
	
	//cout << "tw_FFT_out = " << endl ;
	for (int i = 0; i < m_prime ; i++){
		//cout << tw_FFT_out[i] <<endl;
		//cout << FFT_out[i] <<" ";
	}
		//cout <<endl;	
	

	long long ele_mul[m_prime];
	for (int i = 0; i< m_prime ; i++){	
		ele_mul[i] = (FFT_out[i]*tw_FFT_out[i]) % modular ;
	}
	for (int i = 0; i < m_prime ; i++){
		//cout << ele_mul[i] <<endl;
	}
//---------------IFFT------------------------
	long long IFFT_out[m_prime];
	IDFT(IFFT_out , ele_mul , m_prime , m_prime_prou , modular);

	
//-------------add do------------------------
	//long long RA_out[n];
	RA_out[0] = 0;
	for (int i = 0; i< n ; i++){	
		RA_out[0] =  (RA_out[0] + data_in[i]) %  modular;
	}	
	for (int i = 0; i< n-1 ; i++){	
		RA_out[index_out[i]] = (IFFT_out[i] + data_in[0]) %  modular ;
	}	

}

void LEGACY::Rader_v3(long long *RA_out, long long *data_in, long long *tw_FFT_out, long long *index_in,long long *index_out, long long m, long long m_prime, long long m_prime_prou, long long modular)
{
	/*
	cout << "RA_in = " << data_in << endl;
	for (int i = 0; i < m ; i++){
		cout << data_in[i] << endl ;
	}	
	cout << endl;
	*/	
		
	long long FFT_in[m_prime];
	for (int i = 0; i < m_prime ; i++){
		if(i < m - 1)
			FFT_in[i] = data_in[index_in[i+1]];
		else 
			FFT_in[i] = 0;
		//cout << FFT_in[i] << endl;

	}	
	//cout << endl;
//-------------pointwise mul--------------------------
	long long FFT_out[m_prime];
	FFT(FFT_out, FFT_in, m_prime, m_prime_prou, modular) ;
	
	long long ele_mul[m_prime];
	for (int i = 0; i< m_prime ; i++){	
		ele_mul[i] = (FFT_out[i]*tw_FFT_out[i]) % modular ;
	}

//---------------IFFT------------------------
	long long IFFT_out[m_prime];
	IFFT(IFFT_out , ele_mul , m_prime , m_prime_prou , modular);	
//-------------add do------------------------
	//long long RA_out[n];
	RA_out[0] = 0;
	for (int i = 0; i< m ; i++){	
		RA_out[0] = (RA_out[0] + data_in[i]) % modular ;
	}	
	for (int i = 0; i< m-1 ; i++){	
		RA_out[index_out[i+1]] = (IFFT_out[i] + data_in[0]) % modular;
		//cout << index_out[i+1] << endl;
	}	
	
	//cout << "ok" <<endl;
	/*
	cout << "RA_out = " <<endl;
	for (int i = 0; i< m ; i++){	
		cout << RA_out[i] << " ";
	}	
	cout  <<endl;
	*/
	
}

void LEGACY::Rader_v4(long long *RA_out, long long *data_in, long long *tw_FFT_out, long long *index_in,long long *index_out, long long m, long long m_prime, long long m_prime_prou, long long modular)
{
	/*
	cout << "RA_in = " << data_in << endl;
	for (int i = 0; i < m ; i++){
		cout << data_in[i] << endl ;
	}	
	cout << endl;
	*/	
		
	long long FFT_in[m_prime];
	for (int i = 0; i < m_prime ; i++){
		if(i < m - 1)
			FFT_in[i] = data_in[i+1];
		else 
			FFT_in[i] = 0;
		//cout << FFT_in[i] << endl;

	}	
	//cout << endl;
//-------------pointwise mul--------------------------
	long long FFT_out[m_prime];
	FFT(FFT_out, FFT_in, m_prime, m_prime_prou, modular) ;
	
	long long ele_mul[m_prime];
	for (int i = 0; i< m_prime ; i++){	
		ele_mul[i] = (FFT_out[i]*tw_FFT_out[i]) % modular ;
	}

//---------------IFFT------------------------
	long long IFFT_out[m_prime];
	IFFT(IFFT_out , ele_mul , m_prime , m_prime_prou , modular);	
//-------------add do------------------------
	//long long RA_out[n];
	RA_out[0] = 0;
	for (int i = 0; i< m ; i++){	
		RA_out[0] = (RA_out[0] + data_in[i]) % modular ;
	}	
	for (int i = 0; i< m-1 ; i++){	
		RA_out[i+1] = (IFFT_out[i] + data_in[0]) % modular;
		//cout << index_out[i+1] << endl;
	}	
	
	//cout << "ok" <<endl;
	/*
	cout << "RA_out = " <<endl;
	for (int i = 0; i< m ; i++){	
		cout << RA_out[i] << " ";
	}	
	cout  <<endl;
	*/
	
}


long long LEGACY::cyc_DFT(long long *DFT_data, long long *data_in, long long m, long long prou, long long modular) //primitive root of unity in m-point DFT : m in -> m-1 out
{
	long long DFT_data_tmp[m];
    long long i, j, prou_tmp;
	
    for(i=0;i<m;i++)
    {
        DFT_data_tmp[i] = 0;
    }
    
    //DFT i = 0 to m - 1, cyc_DFT i = 1 to m - 1
    for(i=1;i<m;i++) //cyclotomic gcd(i,m)=1, when m is prime, i = 1~m-1
    {
        prou_tmp = prou_power(prou, i, modular);
        for(j=m-1;j>0;j--)
        {
            DFT_data_tmp[i] += data_in[j];
            DFT_data_tmp[i] *= prou_tmp;
            DFT_data_tmp[i] %= modular;
        }
        DFT_data_tmp[i] += data_in[0];
        DFT_data_tmp[i] %= modular;
    } 	
    
    for(i=0;i<m;i++)
    {
        DFT_data[i] = DFT_data_tmp[i];
    }
    
	return 0;
}

long long LEGACY::cyc_PFA2(long long *DFT_data, long long *data_in, long long m1, long long m2, long long inv1, long long inv2, long long prou, long long modular) //primitive root of unity in m(= m1 * m2)-point PFA
{
    long long n1,n2,k1,k2;
    long long index_in;
    long long index_out;
    long long index_m;
    
    long long data_tmp[system_m1];
    
    
    //PFA with m = m1 * m2 -> m2 * m1-DFT + m1 * m2-DFT
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (n1 * m2 + n2 * m1) % (m1 * m2);
            
            DFT_data[index_m] = data_in[index_in];
        }    
    }
    
    for(n2=0;n2<m2;n2++)
    {
        //m1-point cyc_DFT
        cyc_DFT(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, prou_power(prou, m2, modular), modular);           
    }    
    
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (k1 + k2 * m1) % (m1 * m2);
            
            DFT_data[index_m] = data_tmp[index_in];
        }    
    }

    {
        for(n1=1;n1<m1;n1++) //cyclotomic gcd(n1,m1)=1, when m1 is prime, i = 1~m1-1
        {
            //m2-point cyc_DFT
            cyc_DFT(data_tmp + m2 * n1, DFT_data + m2 * n1, m2, prou_power(prou, m1, modular), modular);           
        }
    }
            
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_out = (k1 * inv2 * m2 + k2 * inv1 * m1) % (m1 * m2);
            
            DFT_data[index_out] = data_tmp[index_m];
        }    
    }
}

long long LEGACY::cyc_PFA3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) //primitive root of unity in m(= m1(= s_m1 * s_m2) * m2)-point PFA
{
    long long n1,n2,k1,k2;
    long long index_in;
    long long index_out;
    long long index_m;
    
    long long data_tmp[system_m];
    
    
    //PFA with m = m1 * m2 -> m2 * m1-DFT + m1 * m2-DFT
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (n1 * m2 + n2 * m1) % (m1 * m2);
            
            DFT_data[index_m] = data_in[index_in];
        }    
    }
    
    for(n2=0;n2<m2;n2++)
    {
        //m1(= s_m1 * s_m2)-point cyc_PFA
        cyc_PFA2(data_tmp + m1 * n2, DFT_data + m1 * n2, s_m1, s_m2, s_inv1, s_inv2, prou_power(prou, m2, modular), modular);         
    }    
    
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_in = (k1 + k2 * m1) % (m1 * m2);
            
            DFT_data[index_m] = data_tmp[index_in];
        }    
    }
 
    {
        for(n1=0;n1<m1;n1++)
        {
            if(((n1 % s_m1) != 0) && ((n1 % s_m2) != 0)) //(gcd(n1, s_m1) != 1) and (gcd(n1, s_m2) != 1)
            {
                //m2-point cyc_DFT
                cyc_DFT(data_tmp + m2 * n1, DFT_data + m2 * n1, m2, prou_power(prou, m1, modular), modular);           
            }
        }
    }
             
    for(n2=0;n2<m2;n2++)
    {
        for(n1=0;n1<m1;n1++)
        {
            index_m = n1 + n2 * m1;
            k2 = index_m % m2;
            k1 = index_m / m2;
            index_out = (k1 * inv2 * m2 + k2 * inv1 * m1) % (m1 * m2);
            
            DFT_data[index_out] = data_tmp[index_m];
        }    
    }
}

// group are range, if 4-bit index, the group is 16
int LEGACY::Gray(int index,int group){
	int result_tmp;
	int weight_tmp;
	int bit_tmp;
	int index_tmp;
	int group_bit_size;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
	group_bit_size = (int) ceil(log2(group));
	bit_array.resize(group_bit_size);
	bit_array_tmp.resize(group_bit_size);
	
	index_tmp = index;
	for(int i = 0;i < group_bit_size ; i++){
        bit_tmp = index_tmp % 2;
        index_tmp = index_tmp >> 1;
        bit_array_tmp[i] = bit_tmp;		
	}
	
	result_tmp = 0;
	for(int i = 0;i < group_bit_size ;i++){
		if(i == (group_bit_size-1)) bit_array[i] = bit_array_tmp[i];
		else bit_array[i] = bit_array_tmp[i] ^ bit_array_tmp[i+1];
		
		if(bit_array[i] == 1) weight_tmp = 1 << i;
		else weight_tmp = 0;
		result_tmp = result_tmp + weight_tmp; 
	}
	
	return result_tmp;
}
//cyclic right shift
int LEGACY::RR(int BC, int shift_bit, int Bit_WIDTH){
	int    RR_out;
	std::vector<int> bit_array;
	bit_array.resize(Bit_WIDTH);
	//bit calculate
	for(int j=0; j < Bit_WIDTH;j++)
	{
		bit_array[j] = (BC >> j) & 1 ;
		//cout << bit_array[j] << endl ;
	} 		
	//cyclic right shift
	std::rotate(bit_array.begin(), bit_array.begin()+shift_bit, bit_array.end());
	RR_out = 0;
	for(int j=0; j < Bit_WIDTH;j++)
	{
		RR_out += bit_array[j] << j ;
	} 
	return RR_out;
}

int LEGACY::bit_reverse(int num, int Bit_WIDTH){
	int    output = 0;
	std::vector<int> bit_array, bit_array_tmp;
	bit_array.resize(Bit_WIDTH);
	bit_array_tmp.resize(Bit_WIDTH);	
	//bit calculate
	for(int j=0; j < Bit_WIDTH;j++)
	{
		bit_array[j] = (num >> j) & 1 ;
		//cout << bit_array[j] << endl ;
	} 		

	for(int j=0; j < Bit_WIDTH;j++)
	{
		bit_array_tmp[j] =  bit_array[Bit_WIDTH - j - 1];
		//cout << bit_array_tmp[j] << endl ;
	} 

	for(int j=0; j < Bit_WIDTH;j++)
	{
		output += bit_array_tmp[j] << j ;
		//cout << output << endl;
	} 
	//cout << output << endl;
	return output;
}


int LEGACY::Bit_convert(int addr){  //bit convert [0][1][2][3][4][5] --> [4][5][2][3][0][1]
	int    RR_out = 0;
	int Bit_WIDTH = 6;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;	
	bit_array.resize(Bit_WIDTH);
	bit_array_tmp.resize(Bit_WIDTH);	
	//bit calculate
	for(int j=0; j < Bit_WIDTH;j++)
	{
		bit_array[j] = (addr >> j) & 1 ;
		//cout << bit_array[j] << endl ;
	} 		
	bit_array_tmp[4] = bit_array[0];
	bit_array_tmp[5] = bit_array[1];
	bit_array_tmp[2] = bit_array[2];
	bit_array_tmp[3] = bit_array[3];
	bit_array_tmp[0] = bit_array[4];
	bit_array_tmp[1] = bit_array[5];	
	
	
	for(int j=0; j < Bit_WIDTH;j++)
	{
		RR_out += bit_array_tmp[j] << j ;
	} 
	return RR_out;
}

int LEGACY::unary_xor(int data_in, int Bit_WIDTH){
	int    xor_out;
	std::vector<int> bit_array;
	bit_array.resize(Bit_WIDTH);
	//bit calculate
	for(int j=0; j < Bit_WIDTH;j++)
	{
		bit_array[j] = (data_in >> j) & 1 ;
		//cout << bit_array[j] << endl ;
	} 		
	xor_out = 0;
	for(int j=0; j < Bit_WIDTH;j++)
	{
		xor_out += bit_array[j] ;
	} 
	xor_out %= 2 ;
	return xor_out;
}

void LEGACY::int2vec(int integer, int Bit_WIDTH, vector<int> &bit_array){
	bit_array.resize(Bit_WIDTH);
	//bit calculate
	for(int j=0; j < Bit_WIDTH;j++)
	{
		bit_array[j] = (integer >> j) & 1 ;
		//cout << bit_array[j] << endl ;
	} 	
}

int LEGACY::vec2int(vector<int> &bit_array, int Bit_WIDTH){
	

	//bit_array.resize(Bit_WIDTH);
	//bit calculate
	int integer = 0;
	for(int j=0; j < Bit_WIDTH;j++)
	{
		integer += bit_array[j] << j ;
		//cout << bit_array[j] << endl ;
	} 	
	return integer ;
}

void LEGACY::Radix_4_BU(vector<ZZ> &output, vector<ZZ> &input, ZZ tw_1_N, long N, ZZ modular)
{
	vector<ZZ> output_temp_s1(4), output_temp_s2(4);
	ZZ tw_1_4;
	PowerMod(tw_1_4, tw_1_N, N/4, modular);
	
	//stage 0
	AddMod(output_temp_s1[0],input[0],input[2],modular);  //a0 + a2
	SubMod(output_temp_s1[2],input[0],input[2],modular);  //a0 - a2
	AddMod(output_temp_s1[1],input[1],input[3],modular);  //a1 + a3
	SubMod(output_temp_s1[3],input[1],input[3],modular);  //a1 - a3
	MulMod(output_temp_s1[3],output_temp_s1[3],tw_1_4,modular); //(a1 - a3)*w
	//stage 1
	AddMod(output_temp_s2[0],output_temp_s1[0],output_temp_s1[1],modular);  //a0 + a2 + a1 + a3
	SubMod(output_temp_s2[1],output_temp_s1[0],output_temp_s1[1],modular);  //a0 + a2 - (a1 + a3)
	AddMod(output_temp_s2[2],output_temp_s1[2],output_temp_s1[3],modular);  //a0 - a2 + ((a1 - a3)*w)
	SubMod(output_temp_s2[3],output_temp_s1[2],output_temp_s1[3],modular);	//a0 - a2 - ((a1 - a3)*w)
	// bit-reverse
	output[0] = output_temp_s2[0];
	output[1] = output_temp_s2[2];	
	output[2] = output_temp_s2[1];	
	output[3] = output_temp_s2[3];	
}

void LEGACY::Relocation_4(vector<ZZ> &v0, vector<ZZ> &v1, vector<ZZ> &v2, vector<ZZ> &v3)
{
	vector<vector<ZZ>> Relocation_temp(4);
	for(int i = 0; i < 4; i++){	
		Relocation_temp[i].resize(4);
	}	
	for(int i = 0; i < 4; i++){	
		Relocation_temp[i][0] = v0[i];
		Relocation_temp[i][1] = v1[i];
		Relocation_temp[i][2] = v2[i];
		Relocation_temp[i][3] = v3[i];			
	}
	v0 = Relocation_temp[0];
	v1 = Relocation_temp[1];
	v2 = Relocation_temp[2];
	v3 = Relocation_temp[3];
	
}

void LEGACY::Relocation_2(vector<ZZ> &v0, vector<ZZ> &v1)
{
	
	//cout << "v0_in[0] = " << v0[0] <<endl;
	//cout << "v0_in[1] = " << v0[1] <<endl;		
	//cout << "v1_in[0] = " << v1[0] <<endl;
	//cout << "v1_in[1] = " << v1[1] <<endl;
	
	vector<vector<ZZ>> Relocation_temp(2);
	for(int i = 0; i < 2; i++){	
		Relocation_temp[i].resize(2);
	}	
	for(int i = 0; i < 2; i++){	
		Relocation_temp[i][0] = v0[i];
		Relocation_temp[i][1] = v1[i];
	}
	v0 = Relocation_temp[0];
	v1 = Relocation_temp[1];
	
	//cout << "v0_out[0] = " << v0[0] <<endl;
	//cout << "v0_out[1] = " << v0[1] <<endl;		
	//cout << "v1_out[0] = " << v1[0] <<endl;
	//cout << "v1_out[1] = " << v1[1] <<endl;	
	
}

void LEGACY::Radix_2_BU(vector<ZZ> &output, vector<ZZ> &input, ZZ modular){
	vector<ZZ> output_temp_s1(2);		
	
	//cout << "radix2_in[0] = " << input[0] <<endl;
	//cout << "radix2_in[1] = " << input[1] <<endl;	
	
	AddMod(output_temp_s1[0],input[0],input[1],modular);  //a0 + a1	
	SubMod(output_temp_s1[1],input[0],input[1],modular);  //a0 - a1		
	output[0] = output_temp_s1[0];
	output[1] = output_temp_s1[1];	
	
	//cout << "radix2_out[0] = " << output[0] <<endl;
	//cout << "radix2_out[1] = " << output[1] <<endl;	
}

ZZ LEGACY::expand_point_RA(ZZ m){
	ZZ m_prime;
	if( (m==3) || (m==5) || (m==17) || (m==257))
		m_prime = m-1;
	else{
		m_prime = find_m_prime(m);	
	}
	return m_prime;
}

void LEGACY::Config_PFA_Rader_FFT(vector<ZZ> &output, vector<ZZ> &input, ZZ m, ZZ modular){
	ZZ factor[20];
	long long 
	int cnt = Factorize(factor , m);
	//cout << cnt << endl;
	ZZ m_s[cnt]; // factor to several small primes 
	ZZ phi_m_s[cnt];
	for(int i = 0; i < cnt; i++){ // ususlly decompose to 3 factor m1 m2 m3
		m_s[i] = factor[i];
		phi_m_s[i] = m_s[i]-1;
		cout << "m" << i << " = "<< m_s[i] << endl;
	}
	//-----------------------------------------------------------
	// expand mi to mi'(power of 2)
	ZZ m_s_prime[cnt];
	for(int i = 0; i < cnt; i++){ // ususlly decompose to 3 factor m1 m2 m3
		m_s_prime[i] = expand_point_RA(m_s[i]);
		cout << "m" << i << "'= "<< m_s_prime[i] << endl;
	}	
	
	// counter for each stage
	ZZ stage_cnt_non_eliminate[cnt];	
	ZZ stage_cnt_eliminate[cnt];	
	for(int i = 0; i < cnt; i++){ // ususlly decompose to 3 factor m1 m2 m3
		stage_cnt_non_eliminate[i] = m/m_s[i];
		cout << "stage" << i << " = "<< stage_cnt_non_eliminate[i] << endl;
	}	
	
    mul(stage_cnt_eliminate[0],m_s[1],m_s[2]);
    mul(stage_cnt_eliminate[1],phi_m_s[0], m_s[2]);	
    mul(stage_cnt_eliminate[2],phi_m_s[0],phi_m_s[1]);		
	for(int i = 0; i < cnt; i++){ // ususlly decompose to 3 factor m1 m2 m3
		//cout << stage_cnt_eliminate[i] << endl;
	}		
	//------------------------------------------------------------
	vector<ZZ> Data_Mem(92837);
	vector<vector<ZZ>> FFT_1024_Mem(512);	
	for(int i = 0; i < 512; i++){	
		FFT_1024_Mem[i].resize(2);
	}	
	
}

void LEGACY::FFT_1024_radix2(vector<ZZ> &output, vector<ZZ> &input, int point, ZZ modular)
{
	int N = point ;
	int r = 2 ;
	int p,g,s;
	p = log(N)/log(r) ;
	g = N/(r*r) ;
	s = log2(r) ;
	int BC_WIDTH; 
	BC_WIDTH = (int)ceil(log2(N/r));	
	//main function
	int BC, MA ;
	vector<vector<ZZ>> Dual_port_mem_2w_512(512);
	vector<vector<ZZ>> input_buf(512);	
	int test_m = N;
	
	
	
	for(int i = 0; i < 512; i++){
		input_buf[i].resize(2);
		Dual_port_mem_2w_512[i].resize(2);
	}		
	
	//cout << " input = " << endl;	
	for(int i = 0; i < test_m/2 ; i++){
		for(int j = 0; j < 2; j++){
			Dual_port_mem_2w_512[i][j] = i + (test_m/2) * j + 1 ;
			input_buf[i][j] = Dual_port_mem_2w_512[i][j];
			//cout << Dual_port_mem_2w_512[i][j] << endl;
		}		
	}
	
	vector<vector<ZZ>> ROM_2w_512(512);
	vector<vector<ZZ>> ROM_2w_512_inv(512);	
	ZZ prou = find_prou(point, modular);
	//cout << "prou = " << prou << endl;
	
	
	for(int i = 0; i < 512; i++){
		ROM_2w_512[i].resize(2);
		ROM_2w_512_inv[i].resize(2);
	}	
	

	
	
	for(int i = 0; i < test_m/2 ; i++){
		for(int j = 0; j < 2; j++){
			if(j == 0){
				ROM_2w_512[i][j] = 1;
				ROM_2w_512_inv[i][j] = 1;
			}
			else {
				PowerMod(ROM_2w_512[i][j], prou, i, modular);
				PowerMod(ROM_2w_512_inv[i][j], prou, -i, modular);
			}
			//cout << ROM_2w_512_inv[i][j] << endl;
		}		
	}	
	
	
	
	int BC_tmp1, BC_tmp2;
	int MA_tmp1, MA_tmp2;
	int tw_idx;


	for (int t = 0; t < p; t++) //stage
	{
		//cout << "stage "<< t << endl ; 
		for(int i = 0; i < g; i++)  // relocation group
		{
			for(int j = 0; j < r; j++) // addr in group
			{
				BC = j*g + i ;
				MA = RR(BC, s*t, BC_WIDTH);
				//cout << "(BC, MA) = ";				
				//cout << "(" << BC << " , "<< MA << ")";	
				//cout << ") \n" ;
			//cout << Dual_port_mem_2w_512[MA][0] << endl;				
			//cout << Dual_port_mem_2w_512[MA][1] << endl;				
				Radix_2_BU(Dual_port_mem_2w_512[MA], Dual_port_mem_2w_512[MA], modular);
			//cout << Dual_port_mem_2w_512[MA][0] << endl;				
			//cout << Dual_port_mem_2w_512[MA][1] << endl;
			
				if(t == p-1){
					tw_idx = 0;
				}
				else {
					tw_idx = (MA % ((N/2)>>t)) << t  ;
				}
				
				//cout << tw_idx << endl;
				MulMod(Dual_port_mem_2w_512[MA][0], Dual_port_mem_2w_512[MA][0], ROM_2w_512[tw_idx][0], modular);
				MulMod(Dual_port_mem_2w_512[MA][1], Dual_port_mem_2w_512[MA][1], ROM_2w_512[tw_idx][1], modular);
				//cout << Dual_port_mem_2w_512[MA][0] << endl;				
				//cout << Dual_port_mem_2w_512[MA][1] << endl;
				//MulMod(Dual_port_mem_2w_512[MA][1], Dual_port_mem_2w_512[MA][1], ROM_2w_512[][1], modular);
			}
			if(t != p-1) {
				BC_tmp1 = 0*g + i;
				BC_tmp2 = 1*g + i;
				MA_tmp1 = RR(BC_tmp1, s*t, BC_WIDTH);
				MA_tmp2 = RR(BC_tmp2, s*t, BC_WIDTH);
				// cout << " g= " << g << endl;
				// cout << "BC_tmp1 = " << BC_tmp1 << endl;
				// cout << "BC_tmp2 = " << BC_tmp2 << endl;					
				// cout << "MA_tmp1 = " << MA_tmp1 << endl;
				// cout << "MA_tmp2 = " << MA_tmp2 << endl;			
				Relocation_2(Dual_port_mem_2w_512[MA_tmp1], Dual_port_mem_2w_512[MA_tmp2]);
			}
		}
	}

//cout << endl;


ZZ inv = find_inv((ZZ)N, modular) ;

//cout << " inv = " << inv << endl;

	for(int i = 0; i < test_m/2 ; i++){
		for(int j = 0; j < 2; j++){
			//Dual_port_mem_2w_512[i][j] = i+8*j;
			//cout << Dual_port_mem_2w_512[i][j] << endl;

		}		
	}



int MA_tmp;

	for (int t = 0; t < p; t++) //stage
	{
		//cout << "stage "<< t << endl ; 
		for(int i = 0; i < g; i++)  // relocation group
		{
			for(int j = 0; j < r; j++) // addr in group
			{
				BC = j*g + i ;
				MA_tmp = RR(BC, s*t, BC_WIDTH);
				MA = bit_reverse(MA_tmp, BC_WIDTH);
				//cout << " MA_tmp = " << MA_tmp << " MA =  " << MA << endl;
				//cout << "(BC, MA) = ";				
				//cout << "(" << BC << " , "<< MA << ")";	
				//cout << ") \n" ;
			//cout << Dual_port_mem_2w_512[MA][0] << endl;				
			//cout << Dual_port_mem_2w_512[MA][1] << endl;				
				Radix_2_BU(Dual_port_mem_2w_512[MA], Dual_port_mem_2w_512[MA], modular);
			//cout << Dual_port_mem_2w_512[MA][0] << endl;				
			//cout << Dual_port_mem_2w_512[MA][1] << endl;
			
				if(t == p-1){
					tw_idx = 0;
				}
				else {
					tw_idx = (MA_tmp % ((N/2)>>t)) << t  ;   ///!!!!! the order of tw is not change, only mapping the addr to inverse
				}
				
				//cout << tw_idx << endl;
				MulMod(Dual_port_mem_2w_512[MA][0], Dual_port_mem_2w_512[MA][0], ROM_2w_512_inv[tw_idx][0], modular);
				MulMod(Dual_port_mem_2w_512[MA][1], Dual_port_mem_2w_512[MA][1], ROM_2w_512_inv[tw_idx][1], modular);
				//cout << Dual_port_mem_2w_512[MA][0] << endl;				
				//cout << Dual_port_mem_2w_512[MA][1] << endl;
				//MulMod(Dual_port_mem_2w_512[MA][1], Dual_port_mem_2w_512[MA][1], ROM_2w_512[][1], modular);
			}
			if(t != p-1) {
				BC_tmp1 = 0*g + i;
				BC_tmp2 = 1*g + i;
				MA_tmp1 = bit_reverse(RR(BC_tmp1, s*t, BC_WIDTH),BC_WIDTH);
				MA_tmp2 = bit_reverse(RR(BC_tmp2, s*t, BC_WIDTH),BC_WIDTH);
				// cout << " g= " << g << endl;
				// cout << "BC_tmp1 = " << BC_tmp1 << endl;
				// cout << "BC_tmp2 = " << BC_tmp2 << endl;					
				// cout << "MA_tmp1 = " << MA_tmp1 << endl;
				// cout << "MA_tmp2 = " << MA_tmp2 << endl;			
				Relocation_2(Dual_port_mem_2w_512[MA_tmp1], Dual_port_mem_2w_512[MA_tmp2]);
			}
		}
	}

	//cout << endl;

	//cout << "output = " << endl;
	for(int i = 0; i < test_m/2 ; i++){
		for(int j = 0; j < 2; j++){
			MulMod(Dual_port_mem_2w_512[i][j], Dual_port_mem_2w_512[i][j], inv, modular);			
			//Dual_port_mem_2w_512[i][j] = i+8*j;

			//cout << Dual_port_mem_2w_512[i][j] << endl;
		}		
	}


//----------check correctness------//
	int k;
	for(int i = 0; i < test_m/2 ; i++){
		for(int j = 0; j < 2; j++){
			if( !(input_buf[i][j] == Dual_port_mem_2w_512[i][j])){
				cout << "error " << endl;
				break;
			}
			else k++;
		}	
	
		if(k == test_m)
			cout << "done " << endl;
	}
//-----------------------------------//
	


	
}



