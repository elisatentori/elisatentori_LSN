/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <stdlib.h>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(int argc, char ** argv) {

    if(argc != 2) { cerr << "Usage: enter temperature"<< endl;
        return -1;}
    double temperature = atof(argv[1]);
    
    Input_Restart(temperature);  //Inizialization startin from config equilibrated

    //Equilibrazione
    for(int istep=1;istep<=300;++istep){   // faccio equilibrare il sistema per 300 steps
        Move();
    }
    Measure();
    cout << "Energy after equilibration = " << walker[iu]/(double)nspin << endl;

    
    //Simulazione
    for(int iblk=1; iblk <= nblk; ++iblk) { //Simulation
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep) {
            Move();
            Measure();
            Accumulate(); //Update block averages
        }
        Averages(iblk);   //Print results for current block
    }
    
    ConfFinal(); //Write final configuration
    PrintFinal(nblk);  // Print (append) the final averages of U,M,H,C in a file
    return 0;
}



void Input_Restart(double t) {
    ifstream ReadInput,ReadConf;

    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;

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

    temp=t;
    beta = 1.0/temp;
    cout << "Temperature = " << temp << endl;

    ReadInput >> nspin;
    cout << "Number of spins = " << nspin << endl;

    ReadInput >> J;
    cout << "Exchange interaction = " << J << endl;

    ReadInput >> h;
    cout << "External field = " << h << endl << endl;
    
    ReadInput >> metro; // if=1 Metropolis else Gibbs
    ReadInput >> nblk;
    ReadInput >> nstep;

    if(metro==1) cout << "The program perform Metropolis moves" << endl;
    else cout << "The program perform Gibbs moves" << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();


    //Prepare arrays for measurements
    iu = 0; //Energy
    iu2 = 1; //Energy square
    im = 2; //Magnetization
    im2 = 3; //Sum square of all the spins of a configuration
    ic = 4; //Heat capacity
    ix= 5; //Magnetic susceptibility
 
    n_props = 6; //Number of observables

    // initial configuration from file
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<nspin; ++i){
        ReadConf >> s[i];
    }
    ReadConf.close();
  
    //Evaluate energy etc. of the initial configuration
    Measure();

    //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}



void Move() {
    // accepted=0;
    int o;
    double r, p, sm;
    double energy_old, energy_new;
    double energy_up,energy_down;

    for(int i=0; i<nspin; ++i) {
        o = (int)(rnd.Rannyu()*nspin);  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
        attempted +=1;
        
        if(metro==1){       //Metropolis
            sm=s[o];
            energy_old = Boltzmann(sm,o);
            energy_new = Boltzmann(sm*(-1),o);
            if(energy_new-energy_old<0) {
                s[o]*=(-1);     //if E_new<E_old -> flip   (we definitely accept the move: A=1)
                accepted+=1;
            }
            else {
                p = exp((energy_old-energy_new)/temp);   //  p = exp{ -beta(E_n-E_o) }      acceptance
                r = rnd.Rannyu();
                if(r<=p) {
                    s[o]*=(-1); //if r<A -> flip (else don't filp)
                    accepted+=1;
                }
            }
        }
        else {      //Gibbs sampling
            r = rnd.Rannyu();
            energy_down = 2/temp * (J*(s[Pbc(o-1)]+s[Pbc(o+1)])+h);
            if(r<=(1/(1+exp(energy_down)))){
                s[o]=-1;
                accepted+=1;
            } else {
                s[o]=+1;
            }
        }
    }
}



double Boltzmann(int sm, int ip) {
    double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
    return ene;
}



void Measure(){  // calcolare u e u2 per calore specifico
    //int bin;
    double u = 0.0, m = 0.0;

    //cycle over spins
    for (int i=0; i<nspin; ++i) {
        u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
        m += s[i];
    }
    walker[iu] = u;   // energia della configurazione
    walker[iu2] = u*u;  // serve per C
    walker[im] = m;  //sum of spins
    walker[im2] = m*m; //sum square of spins
}



void Accumulate(void){ //Update block averages
    for(int i=0; i<n_props-2; ++i) {
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}



void Averages(int iblk){ //Print results for current block
    ofstream Ene, Heat, Mag, Chi;
    const int wd=15;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output_ene.dat",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy per particle
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();
    
    Heat.open("output_heat.dat",ios::app);
    stima_c = ( blk_av[iu2]/blk_norm - pow(blk_av[iu]/blk_norm,2.) )/pow(temp,2.)/(double)nspin; //Heat Capacity per particle
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();
    
    Mag.open("output_mag.dat",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();
    
    Chi.open("output_chi.dat",ios::app);
    if(h==0) stima_x = (blk_av[im2]/blk_norm)/temp/(double)nspin;                           //Magnetic Susceptibility
    else stima_x = (blk_av[im2]/blk_norm-pow(blk_av[im]/blk_norm,2.))/temp/(double)nspin;
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();

    // INCLUDE YOUR CODE HERE

    cout << "----------------------------" << endl << endl;
}




void PrintFinal(int nblk) {
    const int wd=15;
    if(h==0){
        ofstream Ene,Heat,Chi;
        
        Ene.open("final_ene_T.dat",ios::app);
        Ene << setw(wd) << temp << setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
        Ene.close();
        
        Heat.open("final_heat_T.dat",ios::app);
        Heat << setw(wd) << temp << setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
        Heat.close();
        
        Chi.open("final_chi_T.dat",ios::app);
        Chi << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x << endl;
        Chi.close();
    } else {
        ofstream Mag;
        Mag.open("final_mag_T.dat",ios::app);
        Mag << setw(wd) << temp << setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
        Mag.close();
    }
}



void Reset(int iblk) { //Reset block averages
    if(iblk == 1){
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



void ConfFinal(void) {
    ofstream WriteConf;

    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");
    for (int i=0; i<nspin; ++i) {
        WriteConf << s[i] << endl;
    }
    WriteConf.close();

    rnd.SaveSeed();
}



int Pbc(int i){  //Algorithm for periodic boundary conditions
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}



double Error(double sum, double sum2, int iblk) {
    if(iblk==0) return 0;
    else {
        if(sum2/(double)iblk - pow(sum/(double)iblk,2)>0) return sqrt( (sum2/(double)iblk - pow(sum/(double)iblk,2)) / (double)iblk );
        else return sqrt( ( pow(sum/(double)iblk,2) - sum2/(double)iblk ) / (double)iblk );
    }
}



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
