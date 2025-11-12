// File  : CTRL.CPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>


#include <ctype.h>


#include "ctrl.hpp"
#include "evolver.hpp"
#include "err.hpp"


char Controls::Comment = '#';
 int Controls::MinArgC = 2;

 
Symbol::operator unsigned (void)
{ return Value;
}

Symbol::operator char* (void)
{ return Text;
}

Symbol SymbolTable[] =
{// In alphanumeric (text) order
  { Controls::DoComment,        "#"              },
  { Controls::DoBinBest1,       "BINBEST1"       },
  { Controls::DoBinBest2,       "BINBEST2"       },
  { Controls::DoBinRand1,       "BINRAND1"       },
  { Controls::DoBinRand2,       "BINRAND2"       },
  { Controls::DoBinRandToBest1, "BINRANDTOBEST1" },
  { Controls::DoCrossover,      "CROSSOVER"      },
  { Controls::DoExpBest1,       "EXPBEST1"       },
  { Controls::DoExpBest2,       "EXPBEST2"       },
  { Controls::DoExpRand1,       "EXPRAND1"       },
  { Controls::DoExpRand2,       "EXPRAND2"       },
  { Controls::DoExpRandToBest1, "EXPRANDTOBEST1" },
  { Controls::DoFConstant,      "FCONSTANT"      },
  { Controls::DoGenerations,    "GENERATIONS"    },
  { Controls::DoMembers,        "MEMBERS"        },
  { Controls::DoMinimum,        "MINIMUM"        },
  { Controls::DoOutput,         "OUTPUT"         },
  { Controls::DoParameters,     "PARAMETERS"     },
  { Controls::DoRefresh,        "REFRESH"        },
  { Controls::DoSeed,           "SEED"           },
  { Controls::DoStasis,         "STASIS"         },
  { Controls::DoStrategy,       "STRATEGY"       },
  { Controls::DoTime,           "TIME"           },
  { Controls::DoVariables,      "VARIABLES"      },
  { Controls::DoVarMaximum,     "VARMAXIMUM"     },
  { Controls::DoVarMinimum,     "VARMINIMUM"     }
};

const Symbols = sizeof (SymbolTable) / sizeof (Symbol);

Controls::Controls (void)
{ Zero();
}

Controls::Controls (Controls& C)
{ Zero().Init (C);
}

Controls::Controls (int argc, char* argv[])
{ Zero().Init (argc, argv);
}

Controls::Controls (char* Path)
{ Zero().Init (Path);
}

Controls::Controls (ifstream& In)
{ Zero().Init (In);
}

Controls& Controls::Zero (void)
{ Strategy = 0;
  Minimum = 0;
  Generations = 0;
  Stasis = 0;
  Seconds = 0;
  Refresh = 0;
  Count = 0;
  Variables = 0;
  VarMaximum = 0;
  VarMinimum = 0;
  FConstant = 0;
  Crossover = 0;
  Seed = 0;
 *OutPath = 0;
  Reset();
  return *this;
}

Controls& Controls::operator = (Controls& C)
{ Init (C);
  return *this;
}

bool Controls::Init (Controls& C)
{ if (C.Viable())
  { *((State*) this) = (State) C;
    Strategy = C.Strategy;
    Minimum = C.Minimum;
    Generations = C.Generations;
    Stasis = C.Stasis;
    Seconds = C.Seconds;
    Refresh = C.Refresh;
    Count = C.Count;
    Variables = C.Variables;
    VarMaximum = C.VarMaximum;
    VarMinimum = C.VarMinimum;
    FConstant = C.FConstant;
    Crossover = C.Crossover;
    Seed = C.Seed;
    strcpy (OutPath, C.OutPath);
    return 1;
  }
  return 0;
}

bool Controls::Init (int argc, char* argv[])
{ if (argc < MinArgC)
    SetError (Usage);
  else if (Init (argv[1]))
  { if (argc == 3)
      strcpy (OutPath, argv[2]);
    return 1;
  }
  return 0;
}

bool Controls::Init (char* Path)
{ ifstream Config (Path);
  if (Config)
  { Set (AccessBit);
    return Init (Config);
  }
  Exception.Post (FileOpen, Path);
  Reset (AccessBit);
  return 0;                                 // Configuration is necessary
}

bool Controls::Init (ifstream& In)
{ if (In)
  { unsigned State = Command (In);
    while (State != EndOfFile && State != NoCommand)
    { switch (State)
      { case DoComment    :  Ignore (In);              break;
        case DoParameters :  Init (Token (In));        break;
        case DoOutput     :  In >> OutPath;            break;
        case DoStrategy   :  Strategy = Command (In);  break;
        case DoMinimum    :  In >> Minimum;
                             Set (DoMinimum);          break;
        case DoGenerations:  In >> Generations;
                             Set (DoGenerations);      break;
        case DoStasis     :  In >> Stasis;
                             Set (DoStasis);           break;
        case DoTime       :  In >> Seconds;
                             Set (DoTime);             break;
        case DoRefresh    :  In >> Refresh;            break;
        case DoMembers    :  In >> Count;              break;
        case DoVariables  :  In >> Variables;          break;
        case DoVarMinimum :  In >> VarMinimum;         break;
        case DoVarMaximum :  In >> VarMaximum;         break;
        case DoFConstant  :  In >> FConstant;          break;
        case DoCrossover  :  In >> Crossover;          break;
        case DoSeed       :  In >> Seed;               break;
        case NoCommand    :  Zero();                   break;
      }
      State = Command (In);
    }
    return State != NoCommand;
  }
  return 0;
}

int Controls::Compare (const char *S1, const char *S2)
{ char X1[LineSize], *C1 = X1;

  for (strcpy (X1, S1);  *C1 && *S2;  ++C1, ++S2)
    if ((*C1 = toupper (*C1)) > *S2)
      return 1;
    else if (*C1 < *S2)
      return -1;

  if ((*C1 = toupper (*C1)) == *S2)
    return 0;
  else if (! *C1 || (*C1 < *S2))
    return -1;

  return 1;
}

unsigned Controls::Command (ifstream& In)
{ if (Token (In))
  { int Found = 0, Low = 0, High = Symbols - 1;
    int Seeking = 1;

    Symbol* Entry;
    while (Seeking && (Low <= High))
    { int Median = (Low + High) >> 1;
      Entry = SymbolTable + Median;
      if ((Seeking = Compare (Buffer, *Entry)) < 0)
        High = Median - 1;
      else if (Seeking)
        Low = Median + 1;
      else
        Found = 1;
    }

    if (Found)
      return *Entry;
    else if (*Buffer == Comment)
      return DoComment;
    else if (! In.eof())
      return NoCommand;
  }
  return EndOfFile;
}

char* Controls::Token (ifstream& In)
{ if (! In.eof())
  { In >> Buffer;
    return Buffer;
  }
  *Buffer = 0;
  return 0;
}

void Controls::Ignore (ifstream& In)
{ In.getline (Buffer, LineSize);
}

template <class T> bool OutOfRange (T& X, T& Min, T& Max)
{ return X < Min || X > Max;
}

bool Controls::SetError (unsigned Type)
{ Exception.Post (Type);
  Reset (Type);
  return 0;
}

bool Controls::Viable (void)
{ bool OK = 0;

  if (! Test (AccessBit))
    Exception.Post (NoConfig);
  else
  { OK = 1;

    if (! *OutPath)
      OK = SetError (DoOutput);

    if (OutOfRange (Count, Evolver::MinCount, Evolver::MaxCount))
      OK = SetError (DoMembers);

    if (OutOfRange (Variables, Vector::MinCount, Vector::MaxCount))
      SetError (DoVariables);

    if (VarMinimum >= VarMaximum)
      OK = SetError (VarRange);
      
    if (OutOfRange (Strategy, 0, (int)StrategyMax))
      OK = SetError (DoStrategy);
    else switch (Strategy)
    { case DoExpBest1      :
      case DoExpRandToBest1:
      case DoBinBest1      :
      case DoBinRandToBest1: if (Count >= 4)
                               break;
      case DoExpRand1      :
      case DoBinRand1      : if (Count >= 5)
                               break;
      case DoExpBest2      :
      case DoExpRand2      :
      case DoBinBest2      :
      case DoBinRand2      : if (Count >= 6)
                               break;
      default              : OK = SetError (StrategyToN);
    }

    if (OutOfRange (FConstant, Evolver::MinFConstant, Evolver::MaxFConstant))
      OK = SetError (DoFConstant);

    if (OutOfRange (Crossover, Evolver::MinCrossover, Evolver::MaxCrossover))
      OK = SetError (DoCrossover);

    if (Test (DoGenerations)
    && OutOfRange (Generations, Evolver::MinGenerations, Evolver::MaxGenerations))
      OK = SetError (DoGenerations);

    if (Test (DoStasis) && OutOfRange (Stasis, Evolver::MinStasis, Evolver::MaxStasis))
      OK = SetError (DoStasis);

    if (Test (DoTime) && OutOfRange (Seconds, Evolver::MinTime, Evolver::MaxTime))
      OK = SetError (DoTime);

    if (OutOfRange (Refresh, Evolver::MinRefresh, Evolver::MaxRefresh))
      OK = SetError (DoRefresh);
  }
  Set (AccessBit, OK);
  return OK;
}

char* Controls::GetOutPath (char* Path)
{ if (Path)
    return strcpy (Path, OutPath);
  return OutPath;
}

int Controls::GetStrategy (void)
{ return Strategy;
}

int Controls::GetCount (void)
{ return Count;
}

int Controls::GetVariables (void)
{ return Variables;
}

long Controls::GetRefresh (void)
{ return Refresh;
}

long Controls::GetGenerations (void)
{ return Generations;
}

long Controls::GetStasis (void)
{ return Stasis;
}

long Controls::GetTime (void)
{ return Seconds;
}

double Controls::GetMinimum (void)
{ return Minimum;
}

double Controls::GetVarMaximum (void)
{ return VarMaximum;
}

double Controls::GetVarMinimum (void)
{ return VarMinimum;
}

double Controls::GetFConstant (void)
{ return FConstant;
}

double Controls::GetCrossover (void)
{ return Crossover;
}

long Controls::GetSeed (void)
{ return Seed;
}


