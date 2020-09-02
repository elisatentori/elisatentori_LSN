#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <stdlib.h>
#include <iomanip>
#include "Options.h"

using namespace std;

int main(int argc, char ** argv) {

    if(argc != 2) {
        cerr << "Usage: enter temperature"<< endl;
        return -1;
    }
    int path_ = atof(argv[1]);

    Inizialise(path_);
        
    for(int iblk=1; iblk <= nblk; ++iblk) {
        Reset();
        for(int istep=1; istep <= nstep; ++istep) {
            Move();
            Accumulate();
        }
        Averages(iblk);
    }

    return 0;
}



void Inizialise(int path_){
        
    path = path_;

    P_glob_av = 0;
    P_glob_av2 = 0;
    C_glob_av = 0;
    C_glob_av2 = 0;
    
    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
}



void Move() {
    
    //Sampling directly the final asset price
    if(path==0) S_T = S_0 * exp( (r - sigma*sigma/2)*T + sigma * rnd.Gauss(0,1) * sqrt(T) );
    
    //sampling the discretized GBM(r,sigma^2) path of the asset price dividing [0,T] in 100 time intervals
    else{
        double S_new,S_old;
        S_new = S_0;
        
        for(int i=0;i<100;++i){
            S_old = S_new;
            S_new = S_old * exp( (r - sigma*sigma/2)*T/100 + sigma * rnd.Gauss(0,1) * sqrt(T/100) );
        
        }
        S_T = S_new;
    }
}



void Accumulate() {
    if(S_T-K>0) C_sum += exp(-r*T)*(S_T-K);
    else P_sum += exp(-r*T)*(K-S_T);
}



void Averages(int iblk){ //Print results for current block and global averages
    ofstream Put,Call;
    const int wd=15;
    
    cout << "Block number " << iblk << endl;
    
    P_sum /= nstep;
    C_sum /= nstep;
    
    if(path==0) Put.open("put_option_direct.dat",ios::app);
    else Put.open("put_option_discretized.dat",ios::app);
    P_glob_av  += P_sum;
    P_glob_av2 += P_sum * P_sum;
    P_err=Error(P_glob_av,P_glob_av2,iblk);
    Put << setw(wd) << iblk <<  setw(wd) << P_sum << setw(wd) << P_glob_av/(double)iblk << setw(wd) << P_err << endl;
    Put.close();
    
    if(path==0) Call.open("call_option_direct.dat",ios::app);
    else Call.open("call_option_discretized.dat",ios::app);
    C_glob_av  += C_sum;
    C_glob_av2 += C_sum * C_sum;
    C_err=Error(C_glob_av,C_glob_av2,iblk);
    Call << setw(wd) << iblk <<  setw(wd) << C_sum << setw(wd) << C_glob_av/(double)iblk << setw(wd) << C_err << endl;
    Call.close();

    cout << "----------------------------" << endl << endl;
}



void Reset() {  //Reset accumulators
    
    P_sum = 0;
    C_sum = 0;
}



double Error(double sum, double sum2, int iblk) {
    if(iblk==0) return 0;
    else {
        if(sum2/(double)iblk - pow(sum/(double)iblk,2)>0) return sqrt( (sum2/(double)iblk - pow(sum/(double)iblk,2)) / (double)iblk );
        else return sqrt( ( pow(sum/(double)iblk,2) - sum2/(double)iblk ) / (double)iblk );
    }
}
