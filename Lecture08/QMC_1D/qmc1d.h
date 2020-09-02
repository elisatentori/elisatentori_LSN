/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __QMC__
#define __QMC__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

/*
This is the definition of the physical constants. The constants have been
put to 1 in order to work with more comfortable numbers (i.e. the potential
energy of the harmonic oscillator on its ground state is 0.25hbar^2*w)
*/

double hbar = 1;
double boltzmann = 1;
double mass = 1;

double mu = 0.81;   //0.664
double sigma = 0.65;   //0.95


/****** GLOBAL VARIABLES ************/
/*
PIGS is only a flag variable that is used to determine whether you are running
a zero temperature or a finite temperature Path Integral.
lambda is hbar*hbar/2m, a constant widely used in Path Integral Theory.
dtau is the timestep, that is, the "small imaginary-time" by which the total
propagation time is divided by the Path Integral. Remember that known
propagators are "only" approximations of the true propagator that are valid for
small imaginary-times.
*/

double lambda;
double dtau;
int PIGS;
double alpha;
/*
The following declarations are the variables used by QMC1D. Don't worry, they
are self explaining if you roughly know how a PIMC works.
*/

int timeslices, brownianBridgeReconstructions, brownianBridgeAttempts, brownianMotionReconstructions;
int MCSTEPS, equilibration, blocks, histogram_bins;
int timeslices_averages_start, timeslices_averages_end;
double temperature, imaginaryTimePropagation, delta_variational, delta_translation;
double histogram_start, histogram_end;

int acceptedTranslations, acceptedVariational, acceptedBB, acceptedBM;
int totalTranslations, totalVariational, totalBB, totalBM;

double* positions;
double* potential_energy;
double* potential_energy_accumulator;
double* potential_energy_square_accumulator;
                                                                                                                 
double* kinetic_energy;
double* kinetic_energy_accumulator;
double* kinetic_energy_square_accumulator;
                                                                                                                 
double* positions_histogram;
double* positions_histogram_accumulator;
double* positions_histogram_square_accumulator;


void readInput();   // reads input from the file "input.dat"
void deleteMemory(); // handles the dynamic allocation of memory
void initialize();  // initializes the variables
void consoleOutput(); // writes the output on the screen
                                                                                                                 
                                                                                                                 
double potential_density_matrix(double val, double val_next);
double u_prime(double x, int m);
double u_sec(double x, int m);

/*
potential_density_matrix returns only the potential part of the correlation between two adjacent timeslices.
*/

double external_potential(double);  // this is the external potential definition
double external_potential_prime(double); // ...and here goes its first derivative
double external_potential_second(double); // ... and its second derivative 

/*
The derivatives are necessary for the evaluation of the kinetic estimator, because it contains
the laplacian operator ! 
*/                                                                                                    
                                                                                                                 
void translation(); // performs a rigid translation
void brownianBridge();  // reconstructs a segment of the polymer with a free particle propagation. 
void brownianMotion(int);  // reconstructs a segment at the extremities of the polymer with a free particle propagation. 
                                                                                                                 
double variationalWaveFunction(double);  
/*variationalWaveFunction is the variational wave function that is
projected in a PIGS simulation.
*/
double variationalWaveFunction_second(double);
double variationalLocalEnergy(double val);
/*
as for the potential, you have to specify its first and second derivative for the evaluation
of the kinetic local energy.
*/                                                                                                                

int index_mask(int); 
/* index_mask is just a compatibility function that takes into account whether the
polymer is open (PIGS) or closed in periodic boundary contitions (PIMC-ring polymer).
*/

void upgradeAverages(); // at every MCSTEP accumulates the estimators values.

void upgradeHistogram(); // fills the histogram of positions foreach MCSTEP
void endBlock(); // finalizes the averages at the end of each block

double kineticEstimator(double,double);  // evaluates the kinetic energy along the polymer
void finalizePotentialEstimator();
void finalizeKineticEstimator();
void finalizeHistogram();
/*
The last three functions are called at the end of the simulation, basically they average over each
block and evaluate the error on the block average. This is an application of the central limit
theorem.
*/
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
