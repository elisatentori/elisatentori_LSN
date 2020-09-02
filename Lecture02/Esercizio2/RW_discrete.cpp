#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <stdlib.h>
#include <iomanip>
#include "RW_discrete.h"

using namespace std;

int main(int argc, char ** argv) {

    if(argc != 2) {
        cerr << "Usage: enter temperature"<< endl;
        return -1;
    }
    int metr = atof(argv[1]);
    
    Inizialise(metr);
    
    for(int iblk=0; iblk < nblk; ++iblk) {
        Reset();
        for(int istep=0; istep < nstep; ++istep) {
            Move();
            Accumulate(istep);
        }
    }
    Averages();

    return 0;
}


void Inizialise(int metr){
    
    metrica=metr;
    
    //inizializza il vettore distanza
    for (int i=0;i<3;++i){
        RW[i] = 0.0;
    }
    
    //inizializza le variabili globali
    for (int i=0;i<nstep;++i){
        glob_av[i] = 0;
        glob_av2[i] = 0;
    }
    
    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
}


void Move() {
    
    if(metrica==0){
        //avanzamento=[1,-1,1,-1,1,-1]  //avanti o indietro nelle 3 dir
        //direz=[0,0,1,1,2,2]  // direzione x, y o z

        double r = rnd.Rannyu();
        int direzione=0,verso=0;
        
        //Ad ogni step genero uniformemente un numero r in (0,1).
        //Per generare un RW discreto in un reticolo 3D, divido l'intervallo (0,1) in
        //6 sottointervalli di lunghezza equivalente.
        //Se r appartiene a [0,1/6), il camminatore si sposta di +1 in direzione x;
        //se r appartiene a [1/6,2/6), il camminatore si sposta di -1 in direzione x;
        //se r appartiene a [2/6,3/6), il camminatore si sposta di +1 in direzione y;
        //e così via.
        //Quindi il camminatore può muoversi equiprobabilmente in ognuna delle direzioni,
        //sia avanti che indietro.
        for (int j=0; j<6; ++j){
            if( (r>(double)j/6)&(r<(double)(j+1)/6) ){
                verso = avanzamento[j];
                direzione = direz[j];
            }
        }
        
        RW[direzione]+=verso;   //incremento una delle 3 coordinate del vettore posizione di +/-1
    }
    
    else {
        double theta = rnd.Rannyu(0,M_PI);
        double phi = rnd.Rannyu(0,2*M_PI);
        RW[0]+=sin(theta)*cos(phi);    // registro in RW l'avanzamento della posizione
        RW[1]+=sin(theta)*sin(phi);    // passando da sferiche a cartesiane
        RW[2]+=cos(theta);
    }
}


void Accumulate(int istep) {
    
    glob_av[istep] += sqrt(RW[0]*RW[0]+RW[1]*RW[1]+RW[2]*RW[2]);
    glob_av2[istep] += RW[0]*RW[0]+RW[1]*RW[1]+RW[2]*RW[2];
    
}


void Averages(){ //Print results for current block
    ofstream Pos;
    if(metrica==0) Pos.open("RW_discrete.dat",ios::out);
    else Pos.open("RW_continuum.dat",ios::out);
    
    const int wd=15;
    
    for(int istep=0;istep<nstep;++istep){
        glob_av[istep] /= nblk;
        glob_av2[istep] /= nblk;
        err = glob_av2[istep]/nblk - (glob_av[istep]/nblk)*(glob_av[istep]/nblk);
        cout << endl << glob_av[istep] << setw(wd) << glob_av2[istep] << setw(wd) << sqrt(err) << endl;  // err
        Pos << setw(wd) << istep <<  setw(wd) << sqrt(glob_av2[istep]) << setw(wd) << sqrt(err) << endl;
        err=0;
    }
    
    Pos.close();
}



void Reset() {     //Reset block averages
    for (int i=0;i<3;++i){
        RW[i] = 0.0;
    }
}
