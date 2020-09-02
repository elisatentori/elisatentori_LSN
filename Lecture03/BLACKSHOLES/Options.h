#ifndef __INTEGRAL__
#define __INTEGRAL__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

// parameters
int path;
double S_T;
double S_0=100.;  //asset price at t=0
double K = 100.;    //strike price
double T = 1.;      //delivery time
double r = 0.1;     //risk-free interest rate
double sigma = 0.25;//volatility

// averages for PUT and CALL
double P_sum,C_sum;
double P_glob_av,P_glob_av2;
double C_glob_av,C_glob_av2;
double C_err,P_err;


// simulation
int nstep=10000;
int nblk=100;

//functions
void Inizialise(int);
void Move(void);
void Accumulate(void);
void Averages(int);
void Reset(void);
double Error(double,double,int);

#endif
