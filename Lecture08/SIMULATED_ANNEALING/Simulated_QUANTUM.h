/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __1DP__
#define __1DP__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, ie, ipofx;
double bin_size,nbins,sd;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_ene,err_ene;
double accepted, attempted;
double err_pofx;

//configuration
double x;

//quantum state
int npart;
double sigma,mean;
const int hbar=1;
const int m=1;

// simulation
int nstep, nblk;

// simulated annealing
int isimul=0;
int print=0;
double E_T;
double M_accepted, M_attempted;
double T,dT,d_opt;


//functions
void Input(void);
void Annealing(void);
void Print(void);
void Metropolis(void);
void Simulation(void);
void Measure(void);
void Istogram(void);
void Averages(int);
void Reset(int);
void Accumulate(void);
double Psi(double, double, double);
double Psi2(double, double, double);
double DerPsi(double, double, double);
double Der2Psi(double, double, double);
double Pot(double);
double Error(double,double,int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
