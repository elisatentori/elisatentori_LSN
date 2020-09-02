#ifndef __INTEGRAL__
#define __INTEGRAL__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;


// averages
double RW[3];
int avanzamento[6]={1,-1,1,-1,1,-1};  //avanti o indietro dir x, avanti e indietro dir y, avanti e indietro dir z
int direz[6]={0,0,1,1,2,2};  // x, y,z
double glob_av[100],glob_av2[100];
double err=0;


// simulation
int nstep=100;
int nblk=10000;
int metrica;

//functions
void Inizialise(int);
void Move(void);
void Accumulate(int);
void Averages(void);
void Reset(void);

#endif
