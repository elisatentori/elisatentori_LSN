#ifndef __BUFFON__
#define __BUFFON__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const double L=0.25;
const double d=1.;
int n_hits,n_throws;
double pi;

// averages
double blk_av,accepted,attempted;
double glob_av,glob_av2;
double err_p;

// simulation
int nstep=100000;
int nblk=100;

//functions
void SetSeed(void);
void Throw(void);
void Measure(void);
void Averages(int);
void Reset(int);
double Error(double,double,int);

#endif
