#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "LEGACY.h"
#include <assert.h>
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
	
    for(i=0;i<m;i++)
    {
        DFT_data_tmp[i] = 0;
    }
    
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
    
    for(i=0;i<m;i++)
    {
        DFT_data[i] = DFT_data_tmp[i];
    }
    
	return 0;
}
/*
long long LEGACY::IDFT(long long *DFT_data, long long *data_in, long long m, long long prou, long long modular)
{
	long long prou_inv;
	long long n_inv;
	long long IDFT_tmp[n];
	
	prou_inv = find_inv(prou, modular);
	n_inv = find_inv(n, modular);	
	DFT(IFFT_tmp, data_in, n, prou_inv, modular);
	
	for (int i = 0; i< n ; i++){	
		IDFT_data[i] = ( n_inv * IFFT_tmp[i] ) % modular ;
	}	
	
	return 0;	
	
}
*/
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
        //m1-point DFT
        DFT(data_tmp + m1 * n2, DFT_data + m1 * n2, m1, prou_power(prou, m2, modular), modular);           
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
            //m2-point DFT
            DFT(data_tmp + m2 * n1, DFT_data + m2 * n1, m2, prou_power(prou, m1, modular), modular);           
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

long long LEGACY::PFA3(long long *DFT_data, long long *data_in, long long m1, long long m2, long long s_m1, long long s_m2, long long inv1, long long inv2, long long s_inv1, long long s_inv2, long long prou, long long modular) //primitive root of unity in m(= m1(= s_m1 * s_m2) * m2)-point PFA
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
        //m1(= s_m1 * s_m2)-point PFA
        PFA2(data_tmp + m1 * n2, DFT_data + m1 * n2, s_m1, s_m2, s_inv1, s_inv2, prou_power(prou, m2, modular), modular);         
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
    /*
    if((m2 == 17) || (m2 >= 29))
    {
        long long n_2_power;
        long long prom;         //primitive root modulo
        long long prom_temp;
        long long prou_temp;
        long long inv_temp;
        long long i;
        
        n_2_power = log2(2*m2-3) + 1;
        n_2_power = 1 << n_2_power;
    
        long long w_padd[n_2_power];
        long long w_FFT_out[n_2_power];
        
        long long d_padd[n_2_power];
        long long d_FFT_out[n_2_power];
        
        prom = find_prou(m2-1, m2);                 //primitive root modulo
        if(prom == 0)
        {
            printf("Error : can't find primitive root");
            return 0;
        }
        prou_temp = prou_power(prou, m1, modular);  //m2 primitive root of unity
        if(prou_temp == 0)
        {
            printf("Error : can't find primitive root of unity");
            return 0;
        }
        
        for(i = 0; i < n_2_power; i++)
        {
            w_padd[i] = 0;
        }
        
        for(i = 0; i < (m2 - 1); i++)               //Zero padding and cyclic
        {
            prom_temp = prou_power(prom, m2 - 1 - i, m2);    //in index
            w_padd[i] = prou_power(prou_temp, prom_temp, modular);
            if(i != 0)
            {
                w_padd[i + n_2_power - m2 + 1] = w_padd[i]; 
            }
        }
        
        prou_temp = find_prou(n_2_power, modular);  //n primitive root of unity
        if(prou_temp == 0)
        {
            printf("Error : can't find primitive root of unity");
            return 0;
        }
        
        //w FFT first
        FFT(w_FFT_out, w_padd, n_2_power, prou_temp, modular);
    
        //RA-FFT 
        for(n1=0;n1<m1;n1++)
        {
            for(i = 0; i < n_2_power; i++)
            {
                d_padd[i] = 0;
            }
            for(i = 0; i < (m2 - 1); i++)
            {
                d_padd[i] = DFT_data[m2 * n1 + prou_power(prom, i, m2)]; //Zero padding
            }
            //d FFT
            FFT(d_FFT_out, d_padd, n_2_power, prou_temp, modular);
               
            //point-wise mult
            for(i = 0; i < n_2_power; i++)
            {
                d_FFT_out[i] *= w_FFT_out[i];
                d_FFT_out[i] %= modular;
            }     
            
            //IFFT
            FFT(d_padd, d_FFT_out, n_2_power, find_inv(prou_temp, modular), modular);
            
            inv_temp = find_inv(n_2_power, modular);
            
            for(i = 0; i < (m2 - 1); i++)       //d IFFT
            {
                d_padd[i] *= inv_temp;
                d_padd[i] %= modular;
            } 
            
            //add d0 and sum
            data_tmp[m2 * n1] = DFT_data[m2 * n1 + m2 - 1];
            for(i = 0; i < (m2 - 1); i++)
            {
                data_tmp[m2 * n1 + prou_power(prom, m2 - 1 - i, m2)] = (d_padd[i] + DFT_data[m2 * n1]) % modular;
                data_tmp[m2 * n1] += DFT_data[m2 * n1 + i];
                data_tmp[m2 * n1] %= modular;
            } 
        }
    }
    else*/
    {
        for(n1=0;n1<m1;n1++)
        {
            //m2-point DFT
            DFT(data_tmp + m2 * n1, DFT_data + m2 * n1, m2, prou_power(prou, m1, modular), modular);           
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
    /*
    if((m2 == 17) || (m2 >= 29))
    {
        long long n_2_power;
        long long prom;         //primitive root modulo
        long long prom_temp;
        long long prou_temp;
        long long inv_temp;
        long long i;
        
        n_2_power = log2(2*m2-3) + 1;
        n_2_power = 1 << n_2_power;
    
        long long w_padd[n_2_power];
        long long w_FFT_out[n_2_power];
        
        long long d_padd[n_2_power];
        long long d_FFT_out[n_2_power];
        
        prom = find_prou(m2-1, m2);                 //primitive root modulo
        if(prom == 0)
        {
            printf("Error : can't find primitive root");
            return 0;
        }
        prou_temp = prou_power(prou, m1, modular);  //m2 primitive root of unity
        if(prou_temp == 0)
        {
            printf("Error : can't find primitive root of unity");
            return 0;
        }
        
        for(i = 0; i < n_2_power; i++)
        {
            w_padd[i] = 0;
        }
        
        for(i = 0; i < (m2 - 1); i++)               //Zero padding and cyclic
        {
            prom_temp = prou_power(prom, m2 - 1 - i, m2);    //in index
            w_padd[i] = prou_power(prou_temp, prom_temp, modular);
            if(i != 0)
            {
                w_padd[i + n_2_power - m2 + 1] = w_padd[i]; 
            }
        }
        
        prou_temp = find_prou(n_2_power, modular);  //n primitive root of unity
        if(prou_temp == 0)
        {
            printf("Error : can't find primitive root of unity");
            return 0;
        }
        
        //w FFT first
        FFT(w_FFT_out, w_padd, n_2_power, prou_temp, modular);
    
        //RA-FFT 
        for(n1=1;n1<m1;n1++)
        {
            for(i = 0; i < n_2_power; i++)
            {
                d_padd[i] = 0;
            }
            for(i = 0; i < (m2 - 1); i++)
            {
                d_padd[i] = DFT_data[m2 * n1 + prou_power(prom, i, m2)]; //Zero padding
            }
            //d FFT
            FFT(d_FFT_out, d_padd, n_2_power, prou_temp, modular);
               
            //point-wise mult
            for(i = 0; i < n_2_power; i++)
            {
                d_FFT_out[i] *= w_FFT_out[i];
                d_FFT_out[i] %= modular;
            }     
            
            //IFFT
            FFT(d_padd, d_FFT_out, n_2_power, find_inv(prou_temp, modular), modular);
            
            inv_temp = find_inv(n_2_power, modular);
            
            for(i = 0; i < (m2 - 1); i++)       //d IFFT
            {
                d_padd[i] *= inv_temp;
                d_padd[i] %= modular;
            } 
            
            //add d0
            for(i = 0; i < (m2 - 1); i++)
            {
                data_tmp[m2 * n1 + prou_power(prom, m2 - 1 - i, m2)] = (d_padd[i] + DFT_data[m2 * n1]) % modular;
            } 
        }
    }
    else*/
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
    /*
    if((m2 == 17) || (m2 >= 29))
    {
        long long n_2_power;
        long long prom;         //primitive root modulo
        long long prom_temp;
        long long prou_temp;
        long long inv_temp;
        long long i;
        
        n_2_power = log2(2*m2-3) + 1;
        n_2_power = 1 << n_2_power;
    
        long long w_padd[n_2_power];
        long long w_FFT_out[n_2_power];
        
        long long d_padd[n_2_power];
        long long d_FFT_out[n_2_power];
        
        prom = find_prou(m2-1, m2);                 //primitive root modulo
        if(prom == 0)
        {
            printf("Error : can't find primitive root");
            return 0;
        }
        prou_temp = prou_power(prou, m1, modular);  //m2 primitive root of unity
        if(prou_temp == 0)
        {
            printf("Error : can't find primitive root of unity");
            return 0;
        }
        
        for(i = 0; i < n_2_power; i++)
        {
            w_padd[i] = 0;
        }
        
        for(i = 0; i < (m2 - 1); i++)               //Zero padding and cyclic
        {
            prom_temp = prou_power(prom, m2 - 1 - i, m2);    //in index
            w_padd[i] = prou_power(prou_temp, prom_temp, modular);
            if(i != 0)
            {
                w_padd[i + n_2_power - m2 + 1] = w_padd[i]; 
            }
        }
        
        prou_temp = find_prou(n_2_power, modular);  //n primitive root of unity
        if(prou_temp == 0)
        {
            printf("Error : can't find primitive root of unity");
            return 0;
        }
        
        //w FFT first
        FFT(w_FFT_out, w_padd, n_2_power, prou_temp, modular);
    
        //RA-FFT 
        for(n1=0;n1<m1;n1++)
        {
            if(((n1 % s_m1) != 0) && ((n1 % s_m2) != 0)) //(gcd(n1, s_m1) != 1) and (gcd(n1, s_m2) != 1)
            {
                for(i = 0; i < n_2_power; i++)
                {
                    d_padd[i] = 0;
                }
                for(i = 0; i < (m2 - 1); i++)
                {
                    d_padd[i] = DFT_data[m2 * n1 + prou_power(prom, i, m2)]; //Zero padding
                }
                //d FFT
                FFT(d_FFT_out, d_padd, n_2_power, prou_temp, modular);
                   
                //point-wise mult
                for(i = 0; i < n_2_power; i++)
                {
                    d_FFT_out[i] *= w_FFT_out[i];
                    d_FFT_out[i] %= modular;
                }     
                
                //IFFT
                FFT(d_padd, d_FFT_out, n_2_power, find_inv(prou_temp, modular), modular);
                
                inv_temp = find_inv(n_2_power, modular);
                
                for(i = 0; i < (m2 - 1); i++)       //d IFFT
                {
                    d_padd[i] *= inv_temp;
                    d_padd[i] %= modular;
                } 
                
                //add d0
                for(i = 0; i < (m2 - 1); i++)
                {
                    data_tmp[m2 * n1 + prou_power(prom, m2 - 1 - i, m2)] = (d_padd[i] + DFT_data[m2 * n1]) % modular;
                } 
            }
        }
    }
    else*/
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
