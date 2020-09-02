#ifndef __CHI_SQ__
#define __CHI_SQ__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int N=100;   // numero di volte che faccio il test
const int M=100;   // numero di intervalli
const int n_=1000;  // numero di step


//functions
void Init(double *,double *);
void Reset(double *);

#endif
