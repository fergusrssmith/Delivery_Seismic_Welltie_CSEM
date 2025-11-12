// File  : ERR.CPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>


#include <stdio.h>
#include <io.h>


#include "err.hpp"
#include "ctrl.hpp"


Errors Exception;

Symbol ErrorTable [] =
{// In Numeric Order
  { Controls::DoMembers,     "Population size out of range."        },
  { Controls::DoOutput,      "Can't write to output file."          },
  { Controls::DoRefresh,     "Refresh rate is out of range."        },
  { Controls::DoVariables,   "Number of variables is out of range." },
  { Controls::VarRange,      "VarMinimum >= VarMaximum."            },
  { Controls::DoGenerations, "Generation limit is out of range."    },
  { Controls::DoMinimum,     "Target minimum is out of range."      },
  { Controls::DoStasis,      "Stasis limit is out of range."        },
  { Controls::DoTime,        "Time limit is out of range."          },
  { Controls::DoStrategy,    "Strategy out of range."               },
  { Controls::StrategyToN,   "Bad match: Strategy to N."            },
  { Controls::DoFConstant,   "FConstant is out of range."           },
  { Controls::DoCrossover,   "Crossover rate is out of range."      },
  { Controls::NoConfig,      "No configuration."                    },
  { Controls::NoCommand,     "Unknown command: %s."                 },
  { Controls::NoMemory,      "%s: Out of memory."                   },
  { Controls::NoPath,        "Missing file path."                   },
  { Controls::FileExists,    "'%s' does not exist."                 },
  { Controls::FileWriteOK,   "'%s' already exists."                 },
  { Controls::FileOpen,      "Failed opening '%s'."                 },
  { Controls::FileRead,      "Failed reading '%s'."                 },
  { Controls::Usage,         "Usage: DE ConfigFile [OutPathFile]"   },
  { Controls::Unknown,       "Unknown error."                       }
};

const TableSize = sizeof (ErrorTable) / sizeof (Symbol);


Errors::Errors (void)
{ //nada mucho
}

Errors::~Errors (void)
{ //nada mucho
}


bool Errors::Occurred (void)
{ return GetCount();
}

char* Errors::Text (unsigned I)
{ Symbol* Entry;
  int Found = 0, Low = 0, High = TableSize - 1;
  while (! Found && (Low <= High))
  { int Median = (Low + High) >> 1;
    unsigned EntryValue = *(Entry = ErrorTable + Median);
    if (I < EntryValue)
      High = Median - 1;
    else if (I > EntryValue)
      Low = Median + 1;
    else
      Found = 1;
  }
  return Found? *Entry : ErrorTable[TableSize-1];
}

bool Errors::Post (unsigned Type)
{ return Add (Text (Type));
}

bool Errors::Post (unsigned Type, char* S)
{ sprintf (Buffer, Text (Type), S);
  return Add (Buffer);
}

bool Errors::Post (char* S)
{ return Add (S);
}

bool Errors::Test (unsigned Type, char* S)
{ bool Result = 1;

  switch (Type)
  { case Controls::FileExists : Result = ! access (S, 0);
                                break;
    case Controls::FileWriteOK: Result = access (S, 0);
                                break;
//  Gotta plan for the future...
            default : break;
  }

  if (! Result)
    Post (Type, S);

  return Result;
}

bool Errors::Test (unsigned Type)
{ return Test (Type, 0);
}

