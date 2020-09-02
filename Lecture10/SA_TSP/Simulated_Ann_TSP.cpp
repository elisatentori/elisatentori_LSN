#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <string.h>

#include "Cromosoma.h"
#include "Popolazione.h"
#include "Random.h"
#include "Simulated_Ann_TSP.h"


using namespace std;


int main() {
    
    Input();  //carico la configurazione di punti
    
    //creazione una classe popolazione e una classe cromosoma
    Popolazione P(nstep, ncity, npop, m_mut, m_cross);
    Cromosoma C(nstep, ncity, npop, city_x,city_y);
    
    C.Add_Cities();   // associa a ogni punto(x,y) un numero da 1 a 32 e crea il primo cromosoma
    P.AddCromo(C);  // aggiunge il cromosoma alla popolazione
    P.Add_Permutations(C,rnd);   // crea npop permutazioni del primo cromosoma e le aggiuge alla popolazione
    P.BubbleSort();  // ordina la popolazione
    
    // l e m servono per memorizzare gli indici dei due cromosomi su cui fare mutazioni, crossover e metropolis
    int l,m;
    cout << endl << "------------------------------------------" << endl << endl;
        
    while(T>=dT){
        
        accepted=0, attempted=0;
        
        for(int i=1;i<=nstep;++i){
                    
            cout << " step " << i << "    T " << T << endl << endl;
            
            //Selezione di 2 cromosomi con prob p(L)=1-L_cromo/L_tot
            l=P.Selector2(rnd,-1);
            m=P.Selector2(rnd,l);
            
            //simulated annealing
            P.Mutate(rnd,l,m);
            P.Crossover(rnd);  // muto e faccio crossover
            
            //Metropolis
            double p = P.Boltz_weigh(npop-2,npop-1,npop,npop+1,T);  //confronto le L dei due figli con le due della coda
            //double p = P.Boltz_weigh(l,m,npop,npop+1,T);
            if(p >= rnd.Rannyu()) {      // se la mossa è accettata sostituisco i due
                                         // figli alla coda della popolazione
                P.Sostit_tail(npop-2,npop-1,npop,npop+1);
                //P.Sostit_tail(l,m,npop,npop+1);
                P.EraseLast2();
                
                accepted = accepted + 1.0;
            } else {                    // altrimenti cancello le due mutazioni appena fatte e tengo
                                        // la popolazione così com'è
                P.EraseLast2();
            }
            attempted = attempted + 1.0;
            
            cout << endl << " Acceptance rate:  " << accepted/attempted << endl << endl;
            
            
            P.BubbleSort();   // ordina i cromosomi in funzione di L, da L_min a L_max
            
            //P.Print_Pop();
            
            P.Best_L();       // stampa su file L del best path a ogni step
            P.Mean_L();       // stampa su file la media di L e il suo errore in funzione dello step
            cout << endl << "------------------------------------------" << endl << endl;
        }
        
        T -= dT;

    }
    
    P.Best_Path();         // stampa su file le coordinate del best path finale
    
    return 0;
}




void Input(void) {
    
    ifstream ReadInput,ReadConf;

    cout << endl << " Problema del commesso viaggiatore per N città" << endl << endl;
    cout << " Il programma ottimizza la funzione costo  L = |c1-c2| + ... + |cN-1 - cN| + |cN - c1| " << endl;
    cout << " con un Algoritmo Genetico." << endl << endl;
    cout << " Evoluzione di una popolazione di \"npop\" cromosomi;" << endl;
    cout << " ogni cromosoma corrisponde a un percorso del commesso viaggiatore." << endl << endl;
    cout << " Il programma utilizza operatori di Permutazione, Mutazione e Crossover. " << endl << endl;

    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
  
    //Read input informations
    ReadInput.open("input.dat");

    ReadInput >> ncity;
    cout << "Numero di città = " << ncity << endl;
    
    ReadInput >> npop;
    cout << "Numero di cromosomi per generazione = " << npop << endl;
    
    ReadInput >> nstep;
    cout << "Numero di steps = " << nstep << endl;
    
    ReadInput >> m_mut;
    cout << "Parametro probabilità di mutazione = " << m_mut << endl;
    
    ReadInput >> m_cross;
    cout << "Parametro probabilità di crossover = " << m_cross << endl;
    
    ReadInput >> T;
    cout << "Temperatura di annealing = " << T << endl;
    
    ReadInput >> dT;
    cout << "dT di annealing = " << dT << endl;
    
    // initial configuration from file
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<ncity; ++i){
        ReadConf >> city_x[i] >> city_y[i];
    }
    ReadConf.close();
    
}
