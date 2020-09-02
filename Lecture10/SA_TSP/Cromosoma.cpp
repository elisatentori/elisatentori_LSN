#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <stdlib.h>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <string.h>

#include <algorithm>
#include <ctime>
#include <chrono>
#include <random>

#include "Cromosoma.h"
#include "Random.h"

using namespace std;


Cromosoma :: Cromosoma(int _nstep,int _ncity,int _npop,double *_x,double *_y) {
    nstep = _nstep;
    ncity = _ncity;
    npop = _npop;
    for (int i=0;i<_ncity;++i){
        x[i]=_x[i];
        y[i]=_y[i];
    }
    for (int i=0;i<_ncity;++i){
        cout <<x[i] << "  "  <<y[i] << endl;
    }    
}

Cromosoma :: Cromosoma() {}

Cromosoma :: ~Cromosoma() {}


//Aggiunge una città a un cromosoma
void Cromosoma :: Add_Cities() {
    for(int i=0;i<ncity;++i) {
        cromo.push_back(i+1);
    }
    return;
}


//Stampa un cromosoma
void Cromosoma :: PrintCromo() {
    cout << endl << endl;
    for (int i=0;i<ncity;++i){
        cout << cromo[i] << " ";
    }
    return;
}


//Funzione costo
double Cromosoma :: L() {

    double sum=0;
    for (int i=0; i<ncity-1; ++i){
        sum+=sqrt( pow(y[cromo[i]-1] - y[cromo[i+1]-1],2.) + pow(x[cromo[i]-1] - x[cromo[i+1]-1],2.) );
    }
    sum+=sqrt( pow(y[cromo[ncity-1]-1] - y[cromo[0]-1],2.) + pow(x[cromo[ncity-1]-1] - x[cromo[0]-1],2.) );

    return sum;
}


//   *  *  *  Mutazioni   *  *  *

// Permutazione di un cromosoma
void Cromosoma :: Permuta(Random *rand){
    int j;
    for(int i=(ncity-1);i>1;--i){
        j=int(rand->Rannyu(1,i));
        //swap
        cromo[i] += cromo[j];
        cromo[j] = cromo[i] - cromo[j];
        cromo[i] -= cromo[j];
    }
    return;
}



//Mutazione: permutazione di m città contigue con altre m città contigue
void Cromosoma :: Group_Permutation(Random *rand){
    int m = int(rand->Rannyu(1,ncity/2));   // m: numero dei due gruppi di città contigue
    int d = int(rand->Rannyu(1,ncity-2*m));  // d: distanza tra un gruppo di città e l'altro
    
    int i1 = int(rand->Rannyu(1,ncity-2*m-d)); // i1: città da cui parte il primo gruppo
    int i2 = i1 + m + d;  // i2: città da cui parte il secondo gruppo
    swap_ranges(cromo.begin()+i1,cromo.begin()+i1+m,cromo.begin()+i2);
}



//Mutazione: inversione di m città contigue
void Cromosoma :: Inversion(Random *rand){
    int m = int(rand->Rannyu(1,ncity));   //numero delle città da invertire
    int i = int(rand->Rannyu(1,ncity-m)); //città da cui partire
    reverse(cromo.begin()+i,cromo.begin()+i+m);
}



//Mutazione: shift di n città a destra di m posizioni
void Cromosoma:: Shift(Random *rand){
    int m = int(rand->Rannyu(2,ncity-2));  //numero di città da shiftare
    int n = int(rand->Rannyu(2,ncity-m));  //da dove parte il blocco di città
    int h = int(rand->Rannyu(1,ncity-3));  //di quante posizioni voglio shiftare il blocco
    vector <int> c(ncity,0);
    
    //Copio in tre vettori la parte iniziale, il blocco da shiftare e la coda del cromo
    int init[n-1],block[m],tail[ncity-m-n];
    int k=1;
    for(int i=0;i<n-1;++i) {  //il pezzo iniziale arriva fino alla posiz n-1
        init[i]=cromo[k];
        k++;
    }
    for(int i=0;i<m;++i){  //il blocco da spostare va dalla pos n alla pos n+m-1
        block[i]=cromo[k];
        k++;
    }
    for(int i=0;i<ncity-m-n;++i) { //la coda va dalla pos n+m alla pos ncity-1
        tail[i]=cromo[k];
        k++;
    }
    
    //Shifto il blocco
    int j,i;
    k=0;
    for(i=n+h;i<n+h+m;++i){
        if(i<32) {
            c[i]=block[k];
            j=i;
            k++;
        } else{
            j=i-31;
            c[j]=block[k];
            k++;
        }
    }
    int l;
    k=0;
    j++;
    for(i=j;i<j+(ncity-m-n);++i){
        if(i<32) {
            c[i]=tail[k];
            l=i;
            k++;
        } else{
            l=i-31;
            c[l]=tail[k];
            k++;
        }
    }
    int p;
    k=0;
    l++;
    for(i=l;i<l+(n-1);++i){
        if(i<32) {
            c[i]=init[k];
            k++;
        } else{
            p=i-31;
            c[p]=init[k];
            k++;
        }
    }
    c[0]=1;
    
    cromo=c;
    return;
}


// Mutazione: scambia due posizioni del cromosoma
void Cromosoma:: Swap2(Random *rand){
    int M;
    int L = int(rand->Rannyu(1,ncity));
    do{
        M = int(rand->Rannyu(1,ncity));
    } while(M==L);

    //swap
    cromo[L] += cromo[M];
    cromo[M] = cromo[L] - cromo[M];
    cromo[L] -= cromo[M];
    return;
}



//  *  *  *   Metodi per il Crossover   *  *  *

//Restituisce il numero intero (corrispondente a una città) contenuto alla posizione i-esima del cromosoma
int Cromosoma :: GiveMeCity(int i){
    int daaai;
    daaai = cromo[i];  // AAA: passaggi inutili, ma il mio compilatore si prende male se no.
    return daaai;
}



//GiveMeOrder of 2 elements
int Cromosoma :: GiveMeOrder(int city_ult,int city_penult){
    int a,b;
    for(int i=1; i<ncity; ++i){
        if(cromo[i]==city_penult) a=i;
        if(cromo[i]==city_ult) b=i;
    }
    if(a>b) return 1;      // se return 1 -> l'ordine è invertito e in crossover scambio la coda dell'altro cromosoma
    else return 0;
}



//scrive gli n elementi della coda del cromosoma secondo come in a
void Cromosoma :: ChangeOrder(int *a,int n){
    int k=0;
    for(int i=n;i<ncity;++i){
        cromo[i]=a[k];
        k++;
    }
}




//stampa in un file le coordinate delle città
void Cromosoma :: PrintPositions(){
    const int wd=15;
    ofstream Pos;
    Pos.open("best_path.dat",ios::app);
    for (int i=0;i<ncity;i++){
        Pos << setw(wd) << x[cromo[i]-1] << setw(wd) << y[cromo[i]-1] << endl;
    }
    Pos.close();
}


//stampa in un file le coordinate delle città
void Cromosoma :: PrintPositions2(){

    const int wd=15;
    ofstream Pos;
    Pos.open("initial_path.dat",ios::app);
    for (int i=0;i<ncity;i++){
        Pos << setw(wd) << x[cromo[i]-1] << setw(wd) << y[cromo[i]-1] << endl;
    }
    Pos.close();
}

