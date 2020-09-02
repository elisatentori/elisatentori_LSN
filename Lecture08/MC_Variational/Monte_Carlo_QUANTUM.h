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
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_ene,err_ene;
double err_pofx;

//configuration
double x;

//quantum state
int npart;
double mean,sigma,x_start,width;
const int hbar=1;
const int m=1;

// simulation
int nstep, nblk;


//functions
void Input(void);
void Accumulate(void);
void Metropolis(void);
void Measure(void);
void Averages(int);
double Psi(double, double, double);
double Psi2(double, double, double);
double Der2Psi(double, double, double);
double Pot(double);
void Reset(int);
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
