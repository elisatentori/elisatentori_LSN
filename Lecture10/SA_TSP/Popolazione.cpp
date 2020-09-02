#include "Popolazione.h"
#include "Cromosoma.h"
#include "Random.h"

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <string.h>

using namespace std;


Popolazione::Popolazione(int _nstep,int _ncity,int _npop,double _m_mut,double _m_cross) {
    nstep = _nstep;
    ncity = _ncity;
    npop = _npop;
    m_mut = _m_mut;
    m_cross = _m_cross;
}

Popolazione :: Popolazione() {}

Popolazione :: ~Popolazione() {}



//aggiugne un cromosoma alla coda della popolazione
void Popolazione :: AddCromo(Cromosoma chr) {
	pop.push_back(chr);
    return;
}



// gli passi il primo cromo, ti fa npop permutazioni e le carica nella popolazione
void Popolazione ::Add_Permutations(Cromosoma chr, Random &rand){
    _rnd=&rand;
    for(int i=0;i<npop-1;i++){
        chr.Permuta(_rnd);
        pop.push_back(chr);
    }
    return;
}



//ordina i cromosomi della popolazione
void Popolazione :: BubbleSort(){
    int i;
    Cromosoma temp;
    bool swapp = true;
    while(swapp){
        swapp = false;
        for(i=0;i<pop.size()-1;i++){
            if(pop[i].L() > pop[i+1].L()){
                
                temp=pop[i];
                pop[i]=pop[i+1];
                pop[i+1]=temp;
                
                swapp=true;
            }
        }
    }
    
    return;
}


    

// Muta (eventualmente più volte il cromosoma i-esimo della popolazione
void Popolazione :: Mutate(Random &rand,int i,int j){
    Cromosoma son1,son2;
    
    son1=pop[i];
    son2=pop[j];

    _rnd=&rand;
    
    if(_rnd->Rannyu()<=m_mut) son1.Group_Permutation(_rnd);
    if(_rnd->Rannyu()<=m_mut) son1.Inversion(_rnd);
    if(_rnd->Rannyu()<=m_mut) son1.Swap2(_rnd);
    if(_rnd->Rannyu()<=m_mut) son1.Shift(_rnd);

    if(_rnd->Rannyu()<=m_mut) son2.Group_Permutation(_rnd);
    if(_rnd->Rannyu()<=m_mut) son2.Inversion(_rnd);
    if(_rnd->Rannyu()<=m_mut) son2.Swap2(_rnd);
    if(_rnd->Rannyu()<=m_mut) son2.Shift(_rnd);

    pop.push_back(son1);
    pop.push_back(son2);
    
    return;
}



// Crossover
void Popolazione :: Crossover(Random &rand) {
    _rnd=&rand;
    // Genitore 1 - Genitore 2
    int parent[ncity],parent2[ncity];
    int k,l;
    
    //copio gli interi cromosomi che faccio prima
    for(int i=0;i<ncity;++i){
        parent[i] = pop[npop].GiveMeCity(i);
        parent2[i] = pop[npop+1].GiveMeCity(i);
    }
    
    if(_rnd->Rannyu()<=m_cross){
        int n = (int)(_rnd->Rannyu(1,ncity-1));  // dove tagliare la coda
        int tail[ncity-n],tail2[ncity-n],tail_ordered[ncity-n],tail2_ordered[ncity-n];
        
        //copio le code del primo e del secondo
        k=0;
        for(int i=n;i<ncity;++i){
            tail[k] = parent[i];
            tail2[k] = parent2[i];
            k++;
        }
     
        //ordina la coda del primo in base all'ordine delle posizioni nel secondo e viceversa
        k=0;
        for(int i=0; i<ncity; ++i) {
            for(int j=0;j<ncity-n;j++){
                if(parent2[i]==tail[j]){
                    tail_ordered[k] = parent2[i];
                    k++;
                }
            }
        }
     
        l=0;
        for(int i=0; i<ncity; ++i) {
            for(int j=0;j<ncity-n;j++){
                if(parent[i]==tail2[j]){
                    tail2_ordered[l] = parent[i];
                    l++;
                }
            }
        }
        
        pop[npop].ChangeOrder(tail_ordered,n);      //riscrive le code ordinate rispetto all'altro
        pop[npop+1].ChangeOrder(tail2_ordered,n);
    }
    
    return;

}


int Popolazione :: Boltz_weigh(int old1,int old2,int new1,int new2, double T){
    return exp( -1/T  *  (pop[new1].L()+pop[new2].L()-pop[old1].L()-pop[old2].L())  );
}


void Popolazione :: Sostit_tail(int old1,int old2,int new1,int new2){
    pop[old1]=pop[new1];
    pop[old2]=pop[new2];
}


//seleziona un cromosoma i secondo la probabilità p(L[i]) = L[i]/L_tot
int Popolazione :: Selector(Random &rand,int exclude){
    _rnd=&rand;
    int i=0;
    int k=0;
    double cost,r;
    double Ltot = 0;
    
    //somma delle funzioni costo della popolazione
    while(i<pop.size()){
        Ltot+=pop[i].L();
        ++i;
    }
    
    //scelgo un cromosoma con indice != exclude
    i=exclude;
    while(k==0){
        while(i==exclude){
            r = _rnd->Rannyu();
            i = int(r * pop.size());
        }
        cost = pop[i].L();
        
        if(_rnd->Rannyu() <= 1-cost/Ltot){
            k+=1;
        }
    }
    return i;
}


//seleziona un cromosoma i secondo la probabilità p(i) = int (ncity r^2)
int Popolazione :: Selector2(Random &rand, int exclude){
    _rnd=&rand;
    int i;
    double r=_rnd->Rannyu();
    i=int(ncity*pow(r,2));
    while(i==exclude){
        r=_rnd->Rannyu();
        i=int(ncity*pow(r,2));
    }
    return i;
}



void Popolazione :: Init_Path(){
    pop[0].PrintPositions2();
}


//Stampa le coordinate del percorso associato a un certo cromosoma
void Popolazione :: Best_Path(){
    pop[0].PrintPositions();
}


//stampa il percorso migliore della popolazione (ordinata in ordine di L crescente)
void Popolazione :: Best_L(){
    ofstream BestL;
    BestL.open("best_L.dat",ios::app);
    BestL << pop[0].L() << endl;
    BestL.close();
    
    return;
}


// stampa su file la media di L sui migliori npop/2 cromosomi
void Popolazione :: Mean_L(){
    ofstream Mean;
    double stima=0;
    
    //medio L sui migliori npop/2 cromosomi
    int i=0;
    while(i<pop.size()/2){
        stima += pop[i].L();
        ++i;
    }
    stima /= (pop.size()/2);
    
    Mean.open("mean_L.dat",ios::app);
    Mean << stima << endl;
    Mean.close();
    
    return;
}


//stampa un singolo vettore di città singolo (cromosoma)
void Popolazione :: Print_Single(int i){
    pop[i].PrintCromo();
    return;
}


//stampa a video la popolazione
void Popolazione ::Print_Pop(){
    int i=0;
    while(i<pop.size()){
        pop[i].PrintCromo();
        cout << "   L " << i+1 << "  : " << pop[i].L();
        ++i;
    }
    cout << endl;
    return;
}



//cancella gli ultimi 2 cromo
void Popolazione :: EraseLast2(void){
    pop.pop_back();
    pop.pop_back();
}


//cancella l'ultimo cromo
void Popolazione :: EraseLast(void){
    pop.pop_back();
}
