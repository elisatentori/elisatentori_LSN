/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;      // dichiara la classe
   int seed[4];     // vettore di 4 seed
   int p1, p2;
   ifstream Primes("Primes.dms");   // legge i primi due numeri da file
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;

   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");       // legge i 4 semi da file
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    double x,y;
    const int wd=15;
    
    ofstream fileout,fileout2;
    fileout.open("config.square");
    if (fileout.is_open()){
        for(int i=0; i<32; i++){
            x = rnd.Rannyu(-1,1);
            y = rnd.Rannyu(-1,1);
            fileout << x << setw(wd) << y << endl;
        }
    } else cout << "PROBLEM: Unable to open random.out" << endl;
    fileout.close();
    
    fileout2.open("config.circle");
    if (fileout2.is_open()){
        for(int i=0; i<32; i++){
            x = rnd.Rannyu(-1,1);
            if(i%2) y = sqrt(1-x*x);
            else y = -sqrt(1-x*x);
            fileout2 << x << setw(wd) << y << endl;
        }
    } else cout << "PROBLEM: Unable to open random.out" << endl;
    fileout2.close();

    rnd.SaveSeed();      // scrive su file i 4 semi
    return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
