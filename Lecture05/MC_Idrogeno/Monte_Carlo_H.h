#ifndef __1DP__
#define __1DP__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=100;
int n_props, ir;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_R,err_R;

//configuration
double x,y,z;

//quantum state
double x_start,y_start,z_start;
double x_width,y_width,z_width;
int orbital;
const double a_0 = 0.0529;

// simulation
int nstep, nblk;
int Tdistr;

//pigreco
//const double pi=3.1415927;

//functions
void Input(void);
void Equilibration(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Metropolis(void);
void Measure(void);
double Psi2(double, double, double);
double Error(double,double,int);

#endif
