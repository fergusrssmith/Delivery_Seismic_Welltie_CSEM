// File  : EVOLVER.HPP
// Author: Mikal Keenan <keenan@cis,uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>


#ifndef EVOLVERHPP
#define EVOLVERHPP


#include "ctrl.hpp"
#include "vec.hpp"


class    Evolver:   public Controls
{         public:   Evolver (void);
                    Evolver (Controls&);
                    Evolver (int, char**);
                    Evolver (ifstream&);
                    Evolver (char*);
                   ~Evolver (void);

            bool    Init (Controls&);
            bool    Init (int, char**);
            bool    Init (ifstream&);
            bool    Init (char*);
            bool    Viable (void);
            bool    EnableWrite (void);
            void    Configure (Controls&);
            bool    EvolveInterval (void);
         Evolver&   Evolve (void);
         Evolver&   Report (void);
         Evolver&   Record (void);
         Evolver&   Trace (void);

            bool    Successful (void);
            bool    Stagnant (void);
            bool    Exhausted (void);
            bool    OutOfTime (void);

            bool    NotRefresh (void);

            long    GetEvaluation (void);
            long    GetGeneration (void);
            long    GetRedundancy (void);
          double    GetMinimum (void);
	       double    GetPrevious (void);
          Vector*   GetBest (void);
             int    IndexBest (void);

   static double    MinFConstant, MaxFConstant;
   static double    MinCrossover, MaxCrossover;
      static int    MinCount, MaxCount;
     static long    MinGenerations, MaxGenerations;
     static long    MinStasis, MaxStasis;
     static long    MinTime, MaxTime;
     static long    MinRefresh, MaxRefresh;

         private:
         
    typedef bool   (Evolver::*pTest) (void);
           pTest    Terminal;

         Evolver&   Zero (void);
         Evolver&   Void (void);
            bool    Build (void);
            bool    Member (int, int, ...);
            void    Select (int&, int&, int&, int&, int&, int&);
            bool    ExhaustedSuccessfulStagnantOutOfTime (void);
            bool    ExhaustedSuccessfulOutOfTime (void);
            bool    ExhaustedStagnantOutOfTime (void);
            bool    SuccessfulStagnantOutOfTime (void);
            bool    ExhaustedSuccessfulStagnant (void);
            bool    ExhaustedSuccessful (void);
            bool    ExhaustedStagnant (void);
            bool    SuccessfulStagnant (void);
            bool    SuccessfulOutOfTime (void);
            bool    StagnantOutOfTime (void);
            bool    ExhaustedOutOfTime (void);
            bool    NoTerminalState (void);

          double    BestValue;
          double    Previous;
             int    BestIndex;
          Vector   *P0, *P1;
          Vector   *G0, *G1;
          Vector    BestYet;
          Vector    GenBest;
          Vector    Temp;
        ofstream    OutStream;
            long    Evaluation;
            long    Generation;
            long    Redundancy;
};

#endif
