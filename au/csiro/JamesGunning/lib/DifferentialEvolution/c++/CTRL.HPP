// File  : CTRL.HPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>


#ifndef CTRLHPP
#define CTRLHPP


#include <fstream.h>


#include "msg.hpp"
#include "state.hpp"   // Includes varargs


class   Controls:  public State
{         public:  Controls (void);
                   Controls (Controls&);
                   Controls (int, char**);
                   Controls (ifstream&);
                   Controls (char*);

        Controls&  operator = (Controls&);
    virtual bool   Init (Controls&);
    virtual bool   Init (int, char**);
    virtual bool   Init (ifstream&);
    virtual bool   Init (char*);
        unsigned   Command (ifstream&);
            char*  Token (ifstream&);
            void   Ignore (ifstream&);
             int   Compare (const char*, const char*);
    virtual bool   Viable (void);

            char*  GetOutPath (char*);
             int   GetStrategy (void);         // Choice of strategy
             int   GetCount (void);            // Size of population
             int   GetVariables (void);        // Number of variables
            long   GetRefresh (void);          // Console refresh rate
            long   GetGenerations (void);      // Stopping criterion
            long   GetStasis (void);           // Stopping criterion
            long   GetTime (void);             // Stopping criterion
          double   GetMinimum (void);          // Stopping criterion
          double   GetVarMaximum (void);       // Upper variable limit
          double   GetVarMinimum (void);       // Lower variable limit
          double   GetFConstant (void);        // F constant (?)
          double   GetCrossover (void);        // Crossover rate
            long   GetSeed (void);             // RNG seed

            enum { AccessBit,
                   DoExpBest1, DoExpRand1, DoExpRandToBest1,
                   DoExpBest2, DoExpRand2, DoBinBest1, DoBinRand1,
                   DoBinRandToBest1, DoBinBest2, DoBinRand2, StrategyMax,
                   DoParameters, DoMembers, DoOutput, DoRefresh,
                   DoVariables, DoVarMinimum, DoVarMaximum, VarRange,
                   DoGenerations, DoMinimum, DoStasis, DoTime,
                   DoStrategy, StrategyToN, DoFConstant, DoCrossover, DoSeed,
                   NoConfig, NoCommand, NoMemory, NoPath,
                   FileExists, FileWriteOK, FileOpen, FileRead,
                   DoComment, EndOfFile,
                   Usage, Unknown
                 };

            enum { Silent, Interact, PathLength = 65, LineSize = 255 };
      static int   MinArgC;

       protected:
       
virtual Controls&  Zero (void);
            bool   SetError (unsigned);

     static char   Comment;

            char   Buffer [LineSize];
            char   OutPath [PathLength];
             int   Strategy;              // Choice of strategy
             int   Count;                 // Size of population
             int   Variables;             // Number of variables
            long   Refresh;               // Console refresh rate
            long   Generations;           // Stopping criterion
            long   Stasis;                // Stopping criterion
            long   Seconds;               // Stopping criterion
          double   Minimum;               // Stopping criterion
          double   VarMaximum;            // Upper variable limit
          double   VarMinimum;            // Lower variable limit
          double   FConstant;             // F constant (?)
          double   Crossover;             // Crossover rate
            long   Seed;                  // RNG seed
};

struct    Symbol
{       unsigned   Value;
            char*  Text;
        operator   unsigned (void);
        operator   char* (void);
};


#endif
