#ifndef _Popolazione_h_
#define _Popolazione_h_

#include "Cromosoma.h"
#include "Random.h"
#include <vector>


class Popolazione  {
    
    private:
            
        Random *_rnd;
        int nstep, ncity, npop;
        double m_mut, m_cross;
        vector <Cromosoma> pop;
    
    public:
		Popolazione(int _nstep,int _ncity,int _npop,double _m_mut,double _m_cross);
        Popolazione();
		~Popolazione();
        
		//Methods
		void AddCromo(Cromosoma chr);
        void Add_Permutations(Cromosoma chr, Random &rand);
		void BubbleSort();
    
        int Selector(Random &rand,int exclude);
        int Selector2(Random &rand, int exclude);
        void Mutate(Random &rand, int i,int j);
        void Crossover(Random &rand);
        
        void Best_Path();
        void Init_Path();
        void Best_L();
        void Mean_L();
        void Print_Single(int i);
        void Print_Pop();
        void EraseLast2(void);
        
};

#endif
