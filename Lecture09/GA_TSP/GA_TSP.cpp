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
#include "GA_TSP.h"


using namespace std;


int main() {
    
    Input();  //carico la configurazione di punti
    
    //creazione una classe popolazione e una classe cromosoma
    Popolazione P(nstep, ncity, npop, m_mut, m_cross);
    Cromosoma C(nstep, ncity, npop,city_x,city_y);
    
    //i cromosomi sono vettori di numeri interi tra 1 e 32, ogni numero è associato univocamente
    //alla posizione di una città
    
    C.Add_Cities();   // associa a ogni punto(x,y) un numero da 1 a 32 e crea il primo cromosoma
    P.AddCromo(C);  // aggiunge il cromosoma alla popolazione
    P.Add_Permutations(C,rnd);   // crea npop permutazioni del primo cromosoma e le aggiuge alla popolazione
    P.BubbleSort();  // ordina la popolazione
    
    // l e m servono per memorizzare gli indici dei due cromosomi su cui fare mutazioni e crossover
    int l,m;
    cout << endl << "------------------------------------------" << endl << endl;
    
    //simulazione
    for(int i=1;i<=nstep;++i){
        
        cout << " step " << i << endl << endl;
        
        //Selezione di 2 cromosomi con prob r^p   -  non p(L)=1-L_cromo/L_tot
        l=P.Selector2(rnd,-1);
        m=P.Selector2(rnd,l);
        
        //Creazione di due figli identici ai genitori in coda alla popolazione. Li muto con probabilità m
        P.Mutate(rnd,l,m);
        
        //Crossover dei figli con probabilità c
        P.Crossover(rnd);
        
        P.BubbleSort();   // Sorting dei cromosomi in funzione di L, in ordine crescente
        P.EraseLast2();   // Cancella i 2 cromosomi con la peggior L
        
        //P.Print_Pop();
        
        P.Best_L();       // Stampa su file della L del miglior path
        P.Mean_L();       // Stampa su file della media degli L della prima metà ordinata della popolazione
        cout << endl << "------------------------------------------" << endl << endl;
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
    
    // initial configuration from file
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<ncity; ++i){
        ReadConf >> city_x[i] >> city_y[i];
    }
    ReadConf.close();
    
    return;
}
