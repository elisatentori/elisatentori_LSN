#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <stdlib.h>
#include <iomanip>
#include "Integral.h"

using namespace std;

int main() {
    
    SetSeed();
    
    for(int iblk=1; iblk <= nblk; ++iblk) {
        Reset(iblk);
        for(int istep=1; istep <= nstep; ++istep) {
            Accumulate();
        }
        Averages(iblk);
    }

    return 0;
}



void SetSeed(){
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


void Accumulate() {
    double r = rnd.Rannyu();
    r_sum += r;
    s_sum += (r-0.5)*(r-0.5);
}


void Averages(int iblk){ //Print results for current block and global averages
    ofstream Mean,Sigma;
    const int wd=15;
    
    cout << "Block number " << iblk << endl;
    
    r_sum /= nstep;
    s_sum /= nstep;
    
    Mean.open("integral_mean.dat",ios::app);
    r_glob_av  += r_sum;
    r_glob_av2 += r_sum * r_sum;
    r_err=Error(r_glob_av,r_glob_av2,iblk);
    Mean << setw(wd) << iblk <<  setw(wd) << r_sum << setw(wd) << r_glob_av/(double)iblk << setw(wd) << r_err << endl;
    Mean.close();
    
    Sigma.open("integral_sigma.dat",ios::app);
    s_glob_av  += s_sum;
    s_glob_av2 += s_sum * s_sum;
    s_err=Error(s_glob_av,s_glob_av2,iblk);
    Sigma << setw(wd) << iblk <<  setw(wd) << s_sum << setw(wd) << s_glob_av/(double)iblk << setw(wd) << s_err << endl;
    Sigma.close();

    cout << "----------------------------" << endl << endl;
}



void Reset(int iblk) {     //Reset block averages
    if(iblk == 1){
        r_glob_av = 0;
        r_glob_av2 = 0;
        s_glob_av = 0;
        s_glob_av2 = 0;
    }

    r_sum = 0;
    s_sum = 0;
}



double Error(double sum, double sum2, int iblk) {
    if(iblk==0) return 0;
    else {
        if(sum2/(double)iblk - pow(sum/(double)iblk,2)>0) return sqrt( (sum2/(double)iblk - pow(sum/(double)iblk,2)) / (double)iblk );
        else return sqrt( ( pow(sum/(double)iblk,2) - sum2/(double)iblk ) / (double)iblk );
    }
}
