#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_QUANTUM.h"

using namespace std;

int main() {

    Input();
    x = x_start;
    for(int iblk=1; iblk <= nblk; ++iblk) { //Simulation
        Reset(iblk);   //Reset block averages
        
        for(int istep=1; istep <= nstep; ++istep) {
            Metropolis();
            Measure();
            Accumulate(); //Update block averages
        }
        
        Averages(iblk);   //Print results for current block
    }
    return 0;
}



void Input(void) {
    
    ifstream ReadInput;

    cout << "Single quantum particle in a one dimensional (1D)         " << endl;
    cout << "Variational Monte Carlo simulation             " << endl << endl;
    cout << "Potential V(x) = x^4 - 5/2 x^2    " << endl << endl;
    cout << "Hamiltonian H = T + V(x)" << endl << endl;
    cout << "Trial wave function Psi(x) = exp(- (x - mean)^2 / sigma^2 + exp(- (x + mean)^2 / sigma^2" << endl << endl;

    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
  
    //Read input informations
    ReadInput.open("input.dat");

    ReadInput >> mean;
    cout << "Mean = " << mean << endl;

    ReadInput >> sigma;
    cout << "Sigma = " << sigma << endl;

    ReadInput >> x_start;
    cout << "Starting position = " << x_start << endl;
    
    ReadInput >> width;
    cout << "Range width = " << width << endl;

    ReadInput >> nblk;

    ReadInput >> nstep;  // nstep per blocco
    
    
    cout << "The program performs Metropolis moves with uniform transition probability" << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();

    //Prepare arrays for measurements
    ie = 0; //Potential energy
    n_props = 1; //Number of observables
    
    //measurement of p(x)
    ipofx = 1;
    nbins = 100;
    n_props = n_props + nbins;
    bin_size = width/(double)nbins;
}



void Metropolis(void) {
    double p, p_old, p_new;   // p = alfa
    double xold, xnew;

    xold = x;
    p_old = Psi2(xold,mean,sigma);
    
    xnew = rnd.Rannyu(xold-(mean+2*sigma),xold+(mean+2*sigma));
    //xnew = rnd.Rannyu(-mean-sigma/2,mean+sigma/2);
    p_new = Psi2(xnew,mean,sigma);

    p = p_new/p_old;
    if(p >= rnd.Rannyu()) {
        x = xnew;
        accepted = accepted + 1.0;
    }
    else{
        x=xold;
    }
    
    attempted = attempted + 1.0;

}



void Measure() {
    int bin;
    
    walker[ie] = ( -1*hbar/(2*m) * Der2Psi(x,mean,sigma) + Pot(x)*Psi(x,mean,sigma) ) / Psi(x,mean,sigma);
    
    //update of the histogram of p(x)
    bin = nbins/2 + (int)floorf( x / bin_size);
    walker[ipofx+bin] += 1.0;
    
    ofstream X;
    X.open("output.x.dat",ios::app);
    X << x << endl;
    X.close();
    
}



void Averages(int iblk) {  //Print results for current block

    double X_;
    ofstream Ene,Pofx,Pave;
    const int wd=14;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.dat",ios::app);
    Pofx.open("output.pofx.0",ios::app);
    Pave.open("output.pave.0",ios::app);
    
    //Potential energy per particle
    stima_ene = blk_av[ie]/blk_norm; //Potential energy
    glob_av[ie] += stima_ene;
    glob_av2[ie] += stima_ene*stima_ene;
    err_ene=Error(glob_av[ie],glob_av2[ie],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_ene << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_ene << endl;
    cout << endl << setw(wd) << iblk <<  setw(wd) << stima_ene << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_ene << endl;
    
    //p(x)
    for (int k=ipofx; k<ipofx+nbins; ++k) {
        Pofx << setw(wd) << walker[k] << endl;
        glob_av[k] += walker[k];
        glob_av2[k] += walker[k]*walker[k];
    }
    
    if (iblk==nblk){
        for(int k=ipofx; k<ipofx+nbins; ++k){
            err_pofx = Error(glob_av[k], glob_av2[k], iblk);
            X_ = -nbins/2*bin_size + (double)k * bin_size;
            Pave << setw(wd) << X_ << setw(wd) << glob_av[k]/(double)nblk/width/(double)nbins << setw(wd) << err_pofx/width/(double)nbins << endl;
        }
    }
    
    cout << "----------------------------" << endl << endl;

    Ene.close();
    Pofx.close();
    Pave.close();

}



//  *  *  *    Psi(x) and derivatives, p(x), V(x)    *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

double Psi(double xx, double mean, double sigma) {
    return ( exp(-1*(xx-mean)*(xx-mean)/(2*sigma*sigma)) + exp(-1*(xx+mean)*(xx+mean)/(2*sigma*sigma)) );
}



double Psi2(double xx, double mean, double sigma) {
    //return Psi(xx,mean,sigma)*Psi(xx,mean,sigma);
    double psi;
    psi = exp(-1*(xx-mean)*(xx-mean)/(2*sigma*sigma)) + exp(-1*(xx+mean)*(xx+mean)/(2*sigma*sigma));
    return(psi*psi);
}



double DerPsi(double xx, double mean, double sigma) {
    double p;
    p = (xx-mean) * exp( -1*(xx-mean)*(xx-mean)/(2*sigma*sigma) )  +  (xx+mean) * exp( -1*(xx+mean)*(xx+mean)/(2*sigma*sigma) );
    return -1/(sigma*sigma) * p;
}



double Der2Psi(double xx, double mean, double sigma) {
    double p2;
    p2 = -Psi(xx,mean,sigma) + pow(( (xx-mean)/sigma ),2.) * exp(-1*(xx-mean)*(xx-mean)/(2*sigma*sigma)) + pow( ((xx+mean)/sigma) ,2.) * exp(-1*(xx+mean)*(xx+mean)/(2*sigma*sigma));
    return ( 1/(sigma*sigma) * p2 );
}



double Pot(double xx) {
    return xx*xx*(xx*xx-5/2);
}

//  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *



void Reset(int iblk) { //Reset block averages
    //reset the hystogram of p(x)
    for (int k=ipofx; k<ipofx+nbins; ++k) walker[k]=0.0;
    
    if(iblk == 1) {
        for(int i=0; i<n_props; ++i) {
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }

    for(int i=0; i<n_props; ++i) {
        blk_av[i] = 0;
    }
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}



void Accumulate(void) { //Update block averages

    for(int i=0; i<n_props; ++i) {
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;

}



double Error(double sum, double sum2, int iblk) {
    if( iblk == 1 ) return 0.0;
    else {
        if(sum2/(double)iblk - pow(sum/(double)iblk,2)>0) return sqrt( (sum2/(double)iblk - pow(sum/(double)iblk,2)) / (double)(iblk-1) );
        else return sqrt( ( pow(sum/(double)iblk,2) - sum2/(double)iblk ) / (double)(iblk-1) );
    }

}
