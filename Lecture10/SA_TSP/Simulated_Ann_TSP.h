/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __GA__
#define __GA__

#include "Popolazione.h"
#include "Cromosoma.h"
#include "Random.h"

//Random numbers
int seed[4];
Random rnd;

//parameters, observables
int const mcity=32;
double city_x[mcity], city_y[mcity];

double m_mut, m_cross;
double T, dT;

int nstep, ncity, npop;
double attempted=0.0,accepted=0.0;


//functions
void Input(void);


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
