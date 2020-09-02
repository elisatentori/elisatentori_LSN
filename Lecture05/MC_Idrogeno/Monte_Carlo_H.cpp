#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_H.h"

using namespace std;

int main() {

    Input();
    Equilibration();  //comment if you don't want equilibration
    
    for(int iblk=1; iblk <= nblk; ++iblk) { //Simulation
        Reset(iblk);   //Reset block averages and starting point
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

    cout << endl << "Sampling expectation values for the radius <r> of an electron  " << endl;
    cout << "in Hidrogen atom, for wave functions " << endl;
    cout << "Psi_100(x) = a_0^-3/2 1/sqrt(pi) exp(- r/a_0)" << endl;
    cout << "Psi_210(x) = a_0^-5/2 1/sqrt(32 pi) exp(- r/2a_0)" << endl << endl;
    cout << "where a_0 = 0.0529 nm is the Bohr radius." << endl << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "The program uses a_0 units   " << endl << endl;

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

    ReadInput >> x_start;
    ReadInput >> y_start;
    ReadInput >> z_start;
    
    cout << "Starting position = " << x_start << "  " << y_start << "  " << z_start << endl;
    
    ReadInput >> orbital;
    if(orbital==1){
        cout << "Sampling Psi_100 " << endl;   //se orbital==1 1s, se orbital==0 2p
    } else {
        cout << "Sampling Psi_210 " << endl;
    }

    ReadInput >> Tdistr;        //se Tdistr==1 UNIFORM, se no GAUSS
    if(Tdistr==1){
        cout << "The program performs Metropolis moves with uniform transition probability T(x|y) " << endl;
    } else {
        cout << "The program performs Metropolis moves with gaussian transition probability T(x|y) " << endl;
    }
    
    ReadInput >> x_width;       // acceptance 50% empirical rule
    ReadInput >> y_width;
    ReadInput >> z_width;
    
    //starting point
    x = x_start;
    y = y_start;
    z = z_start;
    
    ReadInput >> nblk;

    ReadInput >> nstep;  // nstep for each block
        
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();

    //Prepare arrays for measurements
    ir = 0;
    n_props = 1; //Number of observables

}


void Equilibration(void){
    const int wd=14;
    
    ofstream R;
    R.open("equilibrate.xyz.dat",ios::app);
    for(int i=0;i<=200;++i){
        Metropolis();
        R << setw(wd) << x << setw(wd) << y << setw(wd) << z << setw(wd) << sqrt(x*x+y*y+z*z) << endl;
    }
    R.close();
    return;
}

void Metropolis(void) {   // i numeri sono generati nel metropolis in unità di a_0
    double p, p_old, p_new;
    double xold,xnew,yold,ynew,zold,znew;

    xold = x;
    yold = y;
    zold = z;
    p_old = Psi2(xold,yold,zold);
    
    if(Tdistr==1){ // T uniform
        xnew = rnd.Rannyu(xold-x_width,xold+x_width);
        ynew = rnd.Rannyu(yold-y_width,yold+y_width);
        znew = rnd.Rannyu(zold-z_width,zold+z_width);
    } else {  // T gaussian
        xnew = rnd.Gauss(xold,x_width);
        ynew = rnd.Gauss(yold,y_width);
        znew = rnd.Gauss(zold,z_width);
    }
    p_new = Psi2(xnew,ynew,znew);

    p = p_new/p_old;
    if(p >= rnd.Rannyu()) {
        x = xnew;
        y = ynew;
        z = znew;
        accepted = accepted + 1.0;
    }
    else{
        x = xold;
        y = yold;
        z = zold;
    }
    
    attempted = attempted + 1.0;
}


//Psi2 = probability density
double Psi2(double xx, double yy, double zz) {
    if(orbital==1)
        return exp(-2*sqrt(xx*xx+yy*yy+zz*zz));
    else
        return exp(-1*sqrt(xx*xx+yy*yy+zz*zz))*zz*zz;
}



void Measure() {
    const int wd=14;
    
    walker[ir] = sqrt(x*x+y*y+z*z);
    
    ofstream XYZ;
    XYZ.open("output.xyz.dat",ios::app);
    XYZ << setw(wd) << x << setw(wd) << y << setw(wd) << z << setw(wd) << sqrt(x*x+y*y+z*z) << endl;
    XYZ.close();
}



void Reset(int iblk) { //Reset block averages
    
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



void Averages(int iblk) {  //Print results for current block

    ofstream R;
    const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    R.open("output.distance.dat",ios::app);
    
    //Potential energy per particle
    stima_R = blk_av[ir]/blk_norm; //Potential energy
    glob_av[ir] += stima_R;
    glob_av2[ir] += stima_R*stima_R;
    err_R=Error(glob_av[ir],glob_av2[ir],iblk);
    R << setw(wd) << iblk <<  setw(wd) << stima_R << setw(wd) << glob_av[ir]/(double)iblk << setw(wd) << err_R << endl;
    cout << endl << setw(wd) << iblk <<  setw(wd) << stima_R << setw(wd) << glob_av[ir]/(double)iblk << setw(wd) << err_R << endl;
    
    cout << "----------------------------" << endl << endl;

    R.close();
}



double Error(double sum, double sum2, int iblk) {
    if( iblk == 1 ) return 0.0;
    else return sqrt( (sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1) );
}
