#ifndef _Cromosoma_h_
#define _Cromosoma_h_

#include <vector>

#include "Random.h"


using namespace std;



class Cromosoma {
    
	public:
        Cromosoma(int _nstep,int _ncity,int _npop,double *_x,double *_y);
        Cromosoma();
        ~Cromosoma();
    
        // methods
        void Add_Cities();
        void PrintCromo();
        double L();
        void Permuta(Random *rand);
        void PrintPositions();
        void PrintPositions2();
    
        //mutazioni
        //void Shift(Random *rand);
        void Group_Permutation(Random *rand);
        void Inversion(Random *rand);
        void Shift(Random *rand);
        void Swap2(Random *rand);
        
    
        //crossover
        int GiveMeCity(int i);
        int GiveMeOrder(int city_ult,int city_penult);
        void CopyTail(int n,int *a);
        void CopyOrder(int* _tail,int n,int *a);
        void ChangeOrder(int *a,int n);
    
    private:
        vector<int> cromo;
        double x[32],y[32];
        Random *_rnd;
        int nstep, ncity, npop;    
    
};

#endif
