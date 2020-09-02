#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <stdlib.h>
#include <iomanip>

#include "chi_SQ.h"

using namespace std;


int main() {
    
    double chi2[M];
    double n_i[M];
    Init(chi2,n_i);
    int i,j,k;
    double r;
    
    ofstream Numbers;
    Numbers.open("chi2.dat",ios::out);
    
    //N: numero di volte in cui genero 1000 numeri
    for(i=0; i<N;++i){
        Reset(n_i);   //Resetto i contatori per ogni sottointervallo di (0,1)
        
        // genero n_=10000 numeri random in (0,1) con distribuzione uniforme
        for(j=0; j<n_;++j){
            r = rnd.Rannyu();
            for(k=0;k<M;++k){
                if( (r<=(double)(k+1)/M)&&(r>(double)k/M) ){
                    n_i[k] = n_i[k] + 1.0;
                    cout << " " << n_i[k];
                }
            }
        }cout << endl;
        
        for(int k=0;k<M;++k){
            chi2[i]+=pow(n_i[k]-(double)(n_/M),2.)/(double)(n_/M);
        }
        Numbers << chi2[i] << endl;
    }
    
    Numbers.close();
    
    return 0;
}



void Init(double *X2,double *count){
    
    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
    for(int i=0;i<M;++i){
        X2[i]=0.0;
        count[i]=0.0;
    }
    
    return;
}



void Reset(double *count){
    for(int i=0;i<M;++i){
        count[i]=0.0;
    }
    
    return;
}
