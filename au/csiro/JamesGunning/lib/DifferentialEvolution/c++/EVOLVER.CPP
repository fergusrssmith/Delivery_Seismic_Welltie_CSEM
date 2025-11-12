// File  : EVOLVER.CPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>


#include <stdio.h>
#include <float.h>


#include "evolver.hpp"
#include "timer.hpp"
#include "rng.hpp"
#include "err.hpp"


   int Evolver::MinCount       = 4;
   int Evolver::MaxCount       = 500;

double Evolver::MinFConstant   = 0.5;
double Evolver::MaxFConstant   = 1.5;      // Fix: Realistic?

double Evolver::MinCrossover   = 0.0;
double Evolver::MaxCrossover   = 1.0;

  long Evolver::MinRefresh     = 1;
  long Evolver::MaxRefresh     = 500;      // Fix: Realistic?

  long Evolver::MinGenerations = 100;
  long Evolver::MaxGenerations = 100000;

  long Evolver::MinStasis      = 10;
  long Evolver::MaxStasis      = 1000;     // Fix: Realistic?
  
  long Evolver::MinTime        = 10;
  long Evolver::MaxTime        = 100000;   // Fix: Realistic?


Evolver::Evolver (void)
{ Zero();
}

Evolver::Evolver (Controls& C)
{ Zero().Init (C);
}

Evolver::Evolver (int argc, char* argv[])
{ Zero().Init (argc, argv);
}

Evolver::Evolver (char* Path)
{ Zero().Init (Path);
}

Evolver::Evolver (ifstream& In)
{ Zero().Init (In);
}

Evolver::~Evolver (void)
{ Void();
}

Evolver& Evolver::Zero (void)
{ Controls::Zero();
  Terminal = &Evolver::NoTerminalState;
  BestValue = DBL_MAX;
  Previous = BestValue;
  BestIndex = 0;
  Evaluation = 0;
  Generation = 0;
  Redundancy = 0;
  P0 = 0;
  P1 = 0;
  G0 = 0;
  G1 = 0;
  return *this;
}

Evolver& Evolver::Void (void)
{ if (P0)
    delete [] P0;
  if (P1)
    delete [] P1;
  return Zero();
}

bool Evolver::Init (Controls& C)
{ if (Test (AccessBit))
    Void();
    
  if (C.Viable())
  { Configure (C);

    if (EnableWrite())
    { Random.Init (Seed);
      if (Build())
      { BestYet.Copy (P0[IndexBest()]);              // Best ever
        GenBest.Copy (BestYet);                      // Best in generation
        G0 = P0, G1 = P1;                            // Generation t & t+1
        return 1;
      }
    }
  }
  return 0;
}

bool Evolver::Init (int argc, char* argv[])
{ return Init (Controls (argc, argv));
}

bool Evolver::Init (char* Path)
{ return Init (Controls (Path));
}

bool Evolver::Init (ifstream& In)
{ return Init (Controls (In));
}

bool Evolver::Build (void)
{ int N = Count;
  if (N > 0)
  { if (! ((P0 = new Vector [N]) && (P1 = new Vector [N])))
      Exception.Post (NoMemory, "Evolver::Build");
    else
    { Count = 0;
      if (! (BestYet.Init (Variables) && GenBest.Init (Variables)
      && Temp.Init (Variables)))
        Exception.Post (NoMemory, "Evolver::Build");
      else
      { while (Count < N)
          if (P0[Count].Init (Variables, VarMinimum, VarMaximum)
           && P1[Count].Init (Variables))
            ++Count;
          else
            break;
        if (Count == N)
          return 1;
      }
    }
    Void();
  }
  return 0;
}

bool Evolver::Viable (void)
{ return Controls::Viable() && Count && Test (DoOutput);
}

void Evolver::Configure (Controls& C)
{ Controls::Init (C);
  int Exhaust = Test (DoGenerations);
  int Succeed = Test (DoMinimum);
  int Stagnate = Test (DoStasis);
  int TimeOut = Test (DoTime);

  if (Exhaust && Succeed && Stagnate && TimeOut)
    Terminal = &Evolver::ExhaustedSuccessfulStagnantOutOfTime;
  else if (Exhaust && Succeed && Stagnate)
    Terminal = &Evolver::ExhaustedSuccessfulStagnant;
  else if (Exhaust && Succeed && TimeOut)
    Terminal = &Evolver::ExhaustedSuccessfulOutOfTime;
  else if (Exhaust && Stagnate && TimeOut)
    Terminal = &Evolver::ExhaustedStagnantOutOfTime;
  else if (Succeed && Stagnate && TimeOut)
    Terminal = &Evolver::SuccessfulStagnantOutOfTime;
  else if (Exhaust && Succeed)
    Terminal = &Evolver::ExhaustedSuccessful;
  else if (Exhaust && Stagnate)
    Terminal = &Evolver::ExhaustedStagnant;
  else if (Succeed && Stagnate)
    Terminal = &Evolver::SuccessfulStagnant;
  else if (Exhaust && TimeOut)
    Terminal = &Evolver::ExhaustedOutOfTime;
  else if (Succeed && TimeOut)
    Terminal = &Evolver::SuccessfulOutOfTime;
  else if (Stagnate && TimeOut)
    Terminal = &Evolver::StagnantOutOfTime;
  else if (Exhaust)
    Terminal = &Evolver::Exhausted;
  else if (Succeed)
    Terminal = &Evolver::Successful;
  else if (Stagnate)
    Terminal = &Evolver::Stagnant;
  else if (TimeOut)
    Terminal = &Evolver::OutOfTime;
  else
    Terminal = &Evolver::NoTerminalState;
}

bool Evolver::Member (int I, int N, ...)
{ va_list ArgPtr;
  va_start (ArgPtr, N);

  bool True = 0;
  while (N--)
  { int X = va_arg (ArgPtr, int);
    if (I == X)
    { True = 1;
      break;
    }
  }

  va_end (ArgPtr);

  return True;
}

void Evolver::Select (int& I, int& I1, int& I2, int& I3, int& I4, int& I5)
// Pick random Vector.
{ do I1 = Random.Next (Count);
  while (I1 == I);

  do I2 = Random.Next (Count);
  while ((I2 == I) || (I2 == I1));

  do I3 = Random.Next (Count);
  while (Member (I3, 3, I, I1, I2));

  do I4 = Random.Next (Count);
  while (Member (I4, 4, I, I1, I2, I3));

  do I5 = Random.Next (Count);
  while (Member (I5, 5, I, I1, I2, I3, I4));
}

Evolver& Evolver::Evolve (void)
{ BestIndex = 0;

  for (int I = 0;  I < Count;  I++)                  // Loop through pool
  { Temp.Copy (G0[I]);

    int I1, I2, I3, I4, I5;                             // Random indices
    Select (I, I1, I2, I3, I4, I5);

    switch (Strategy)
//  EXPONENTIAL CROSSOVER
	 { case DoExpBest1:
        Temp.DEBest1Exp (GenBest, G0[I2], G0[I3], FConstant, Crossover);
	     break;
      case DoExpRand1:
        Temp.DERand1Exp (G0[I1], G0[I2], G0[I3], FConstant, Crossover);
	     break;
      case DoExpRandToBest1:
        Temp.DERandToBest1Exp (GenBest, G0[I1], G0[I2], FConstant, Crossover);
	     break;
      case DoExpBest2:
        Temp.DEBest2Exp (GenBest, G0[I1], G0[I2], G0[I3], G0[I4], FConstant, Crossover);
        break;
	   case DoExpRand2:
        Temp.DERand2Exp (G0[I5], G0[I1], G0[I2], G0[I3], G0[I4], FConstant, Crossover);
        break;
//  BINOMIAL CROSSOVER with same strategies
      case DoBinBest1:
        Temp.DEBest1Bin (GenBest, G0[I2], G0[I3], FConstant, Crossover);
	     break;
      case DoBinRand1:
        Temp.DERand1Bin (G0[I1], G0[I2], G0[I3], FConstant, Crossover);
	     break;
      case DoBinRandToBest1:
        Temp.DERandToBest1Bin (GenBest, G0[I1], G0[I2], FConstant, Crossover);
	     break;
      case DoBinBest2:
        Temp.DEBest2Bin (GenBest, G0[I1], G0[I2], G0[I3], G0[I4], FConstant, Crossover);
        break;
	   default:
        Temp.DERand2Bin (G0[I5], G0[I1], G0[I2], G0[I3], G0[I4], FConstant, Crossover);
        break;
    }

    double Test = Temp.Evaluate();          // Trial mutation now in Temp
    Evaluation++;

    if (Test <= P0[I])                      // Better value?
    { G1[I].Copy (Temp);                    // Propagate new
      if (Test < BestValue)                 // If new low...
      { BestYet.Copy (Temp);                //  Update best
        BestValue = Test;                   //  and reset...
        BestIndex = I;
      }
    }
    else                                    // Otherwise,
      G1[I].Copy (G0[I]);                   //  Propagate old
  }

  GenBest.Copy (BestYet);                   // Store current best

  Vector* GX = G0;                          // Swap pool arrays
  G0 = G1;                                  // New generation becomes old
  G1 = GX;

  BestValue == Previous?  ++Redundancy : (Redundancy = 0);

  Generation++;

  return Trace();
}

bool Evolver::EvolveInterval (void)
{ do Evolve();
  while (NotRefresh());
  return ! (this->*Terminal)();
}

long Evolver::GetEvaluation (void)
{ return Evaluation;
}

long Evolver::GetGeneration (void)
{ return Generation;
}

long Evolver::GetRedundancy (void)
{ return Redundancy;
}

Vector* Evolver::GetBest (void)
{ return &BestYet;
}

int Evolver::IndexBest (void)
{ BestValue = P0[BestIndex = 0];
  for (int I = 1;  I < Count;  I++)
  { double X = P0[I];
    if (X < BestValue)
    { BestValue = X;
      BestIndex = I;
    }
  }
  return BestIndex;
}

bool Evolver::NotRefresh (void)
{ return Generation % Refresh;                   // Test refresh interval
}

bool Evolver::NoTerminalState (void)
{ return 0;
}

bool Evolver::Exhausted (void)
{ return Generation >= Generations;
}

bool Evolver::Successful (void)
{ return BestValue <= Minimum;
}

bool Evolver::Stagnant (void)
{ return Redundancy >= Stasis;
}

bool Evolver::OutOfTime (void)
{ return Time.Elapsed();
}

bool Evolver::ExhaustedOutOfTime (void)
{ return Exhausted() || OutOfTime();
}

bool Evolver::SuccessfulOutOfTime (void)
{ return Successful() || OutOfTime();
}

bool Evolver::StagnantOutOfTime (void)
{ return Stagnant() || OutOfTime();
}

bool Evolver::ExhaustedSuccessful (void)
{ return Exhausted() || Successful();
}

bool Evolver::ExhaustedStagnant (void)
{ return Exhausted() || Stagnant();
}

bool Evolver::SuccessfulStagnant (void)
{ return Successful() || Stagnant();
}

bool Evolver::ExhaustedSuccessfulStagnant (void)
{ return ExhaustedSuccessful() || Stagnant();
}

bool Evolver::ExhaustedSuccessfulOutOfTime (void)
{ return ExhaustedSuccessful() || OutOfTime();
}

bool Evolver::ExhaustedStagnantOutOfTime (void)
{ return ExhaustedStagnant() || OutOfTime();
}

bool Evolver::SuccessfulStagnantOutOfTime (void)
{ return SuccessfulStagnant() || OutOfTime();
}

bool Evolver::ExhaustedSuccessfulStagnantOutOfTime (void)
{ return ExhaustedSuccessful() || StagnantOutOfTime();
}

Evolver& Evolver::Report (void)
{ cout << "\n\n\n Best evaluation so far: " << BestValue;
  BestYet.Report();
  char S [81];
  sprintf (S, " Generation: %-8ld   Evaluations: %ld", Generation, Evaluation);
  cout << "\n\n" << S;
  sprintf (S, "   Strategy: %-8d       Members: %d", Strategy, Count);
  cout << "\n" << S;
  cout << "      FConstant: " << FConstant << "   Crossover: " << Crossover;
  return *this;
}

bool Evolver::EnableWrite (void)
{ bool OK = Exception.Test (FileWriteOK, OutPath);

  if (! OK)
  { cout << "\nOK to overwrite " << OutPath << " (y/n)? ";
    char C;  cin >> C;
    if (OK = (C == 'y' || C == 'Y'))      // Fix: Get new file name here?
      Exception.Remove (Exception.Last());
  }

  if (OK)
  { OutStream.open (OutPath);
    if (! OutStream)
    { Exception.Post (FileOpen, OutPath);
      OK = 0;
    }
  }

  Set (DoOutput, OK);

  return OK;
}

Evolver& Evolver::Record (void)
{ OutStream << "\n\n\n Best evaluation so far: " << BestValue;
  BestYet.Report (OutStream);
  char S [81];
  sprintf (S, " Generation: %-8ld   Evaluations: %ld", Generation, Evaluation);
  OutStream << "\n\n" << S;
  sprintf (S, "   Strategy: %-8d       Members: %d", Strategy, Count);
  OutStream << "\n" << S;
  OutStream << "      FConstant: " << FConstant << "   Crossover: " << Crossover;
  return *this; //Fix: OutStream success?
}

Evolver& Evolver::Trace (void)
{ OutStream << '\n' << Evaluation << "   " << BestValue;
  return *this;
}


