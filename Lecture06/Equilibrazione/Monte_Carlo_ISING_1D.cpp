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
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main() {
    
    //Input(); //Inizialization
    Input_Restart();  //Inizialization startin from config equilibrated

    for(int istep=1; istep <= nstep; ++istep) {
        Move(metro);
        Measure();
    }
    
    ConfFinal(); //Write final configuration

    return 0;
}




void Input(void) {
    ifstream ReadInput;
    
    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl << endl;
    cout << "Equilibration of the sistem" << endl << endl;

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

    ReadInput >> temp;
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

    //initial configuration
    for (int i=0; i<nspin; ++i) {
        if(rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
    }
  
    //Evaluate energy etc. of the initial configuration
    Measure();
}



void Input_Restart(void) {
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

    ReadInput >> temp;
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

    // initial configuration from file
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<nspin; ++i){
        ReadConf >> s[i];
    }
    ReadConf.close();
  
    //Evaluate energy etc. of the initial configuration
    Measure();
}



void Move(int metro) {
    int o;
    double r, A, energy_old, energy_new, sm;

    for(int i=0; i<nspin; ++i) {
        o = (int)(rnd.Rannyu()*nspin);  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)

        if(metro==1){       //Metropolis
            sm=s[o];
            energy_old = Boltzmann(sm,o);
            energy_new = Boltzmann(sm*(-1),o);
            if(energy_new-energy_old<0) s[o]*=(-1);     //if E_new<E_old -> flip
            else {
                A = pow(M_E,(energy_old-energy_new)/temp);  //acceptance
                r = rnd.Rannyu();
                if(r<A) s[o]*=(-1); //if r<A -> flip (else don't filp)
            }
        }
        else {      //Gibbs sampling
            return;// INCLUDE YOUR CODE HERE
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
    
    ofstream Ene, Mag;

    Ene.open("output.ene.0",ios::app);
    Ene << u/(double)nspin << endl;
    Ene.close();
    
    if(h!=0) {
        Mag.open("output.mag.0",ios::app);
        Mag << m << endl;
        Mag.close();
    }
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



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
