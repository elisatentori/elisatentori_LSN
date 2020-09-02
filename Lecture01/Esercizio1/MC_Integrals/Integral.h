#ifndef __INTEGRAL__
#define __INTEGRAL__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;


// averages
double r_sum,s_sum;
double r_blk_av,s_blk_av;
double r_glob_av,r_glob_av2;
double s_glob_av,s_glob_av2;
double r_err,s_err;


// simulation
int nstep=100000;
int nblk=100;

//functions
void SetSeed(void);
void Accumulate(void);
void Averages(int);
void Reset(int);
double Error(double,double,int);

#endif
