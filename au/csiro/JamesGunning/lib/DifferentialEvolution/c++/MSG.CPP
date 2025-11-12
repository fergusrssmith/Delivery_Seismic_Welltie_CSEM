// File  : MSG.CPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>


#include "msg.hpp"


Messages::Messages (void)
{ Zero();
}

Messages::~Messages (void)
{ Void();
}

bool Messages::Save (ofstream& F)
{ if (Reset(), F << "\nMessages")
  { while (Next())
      if (! (F << '\n' << *Current()))
        break;
    F << endl;
    return ! Cursor && F.good();
  }
  return 0;
}

bool Messages::Add (const string& J)
{ string* X = new string (J);
  if (! X)
    cerr << "\nMessages::Add(&): NoMemory";
  else if (! Add (X))
    delete X;
  return X != 0;
}

bool Messages::Add (string* J)
{ if (J)
  { STRING* X = new STRING;
    if (! X)
      cerr << "\nMessages::Add(*): NoMemory";
    else
    { X->Data = J;
      if (Begin)
      { X->Prev = End;
        X->Next = Begin;
        Begin->Prev = X;
        End->Next = X;
      }
      else
      { Begin = X;
        Begin->Prev = X;
        Begin->Next = X;
      }
      End = X;
      return ++Count;
    }
  }
  return 0;
}

bool Messages::Viable (void)
{ return Begin != 0;
}

int Messages::GetCount (void)
{ return Count;
}

Messages& Messages::Reset (void)
{ Cursor = 0;
  return *this;
}

string* Messages::First (void)
{ return (Cursor = Begin)?  Cursor->Data : 0;
}

string* Messages::Last (void)
{ return (Cursor = End)?  Cursor->Data : 0;
}

string* Messages::Next (void)
{ if (Cursor = Cursor?  (Cursor != End?  Cursor->Next : 0) : Begin)
    return Cursor->Data;
  return 0;
}

string* Messages::operator [] (unsigned I)
//Indexed from 1, not 0
{ if (I && (I <= Count))
  { Reset();
    while (I && Next())
      --I;
    if (Cursor)
      return Cursor->Data;
  }
  return 0;
}

string* Messages::Current (void)
{ return Cursor?  Cursor->Data : 0;
}

string* Messages::Current (string* X)
{ Reset();
  while (Next())
    if (Cursor && (Cursor->Data == X))
      return Cursor->Data;
  return 0;
}

string* Messages::Remove (string* X)
{ return Current (X)? Remove() : NULL;
}

string* Messages::Remove (void)
{ STRING* X = Cursor;
  string* J = Cursor->Data;

  if (Cursor)
  { if (Cursor->Next != Cursor)
    { Cursor->Prev->Next = Cursor->Next;
      Cursor->Next->Prev = Cursor->Prev;
      if ((Cursor = Cursor->Next) == Begin)
        Cursor = 0;
      if (X == Begin)
        Begin = Begin->Next;
      End = Begin->Prev;
      Count--;
    }
    else
      Zero();
    delete X;
    return J;
  }
  return 0;
}

Messages& Messages::Zero (void)
{ Count = 0;
  Cursor = 0;
  Begin = 0;
  End = 0;
  return *this;
}

void Messages::Void (void)
{ if (Begin)
  { End->Next = 0;
    while (Begin)
    { Cursor = Begin;
      Begin = Begin->Next;
      if (Cursor->Data)
        delete Cursor->Data;
      delete Cursor;
    }
  }
  Zero();
}


