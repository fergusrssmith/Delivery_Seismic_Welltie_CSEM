// File..........:  State.CPP
// Content.......:  Implementation of State class
// Author........:  Mikal Keenan <keenan@cis.uab.edu>
//                     <http://www.cis.uab.edu/info/grads/mbk>
// Update........:  03/24/95


#include <limits.h>


#include "state.hpp"


unsigned State::Bytes = sizeof (StateData);
unsigned State::Bits = State::Bytes * CHAR_BIT;
unsigned State::BitLimit = Bits - 1;


State::State (void)
{ Reset();
}


State& State::operator = (const State& S)
{ Data = S.Data;
  return *this;
}


inline State::StateData State::MaskBit (unsigned Offset)
{ return (StateData) 1 << (Offset & BitLimit);
}


bool State::Set (unsigned Offset, int On)
{ return On?  Set (Offset) : Reset (Offset);
}


bool State::Set (unsigned Offset)
{ return (Offset < Bits) && (Data |= MaskBit (Offset));
}


bool State::Reset (unsigned Offset)
{ return (Offset < Bits) && (Data &= ~MaskBit (Offset));
}


bool State::Test (unsigned Offset)
{ return (Offset < Bits) && (Data & MaskBit (Offset));
}


bool State::SetAnd (unsigned N, ...)
{ bool And = 1;
  va_list V;
  va_start (V, N);
  for (int Index;  N && (Index = va_arg (V, int));  N--)
    if (! Set (Index))
      And = 0;
  va_end (V);
  return And;
}


bool State::TestAnd (unsigned N, ...)
{ bool And = 1;
  va_list V;
  va_start (V, N);
  for (int Index;  N && (Index = va_arg (V, int));  N--)
    if (! Test (Index))
      And = 0;
  va_end (V);
  return And;
}


bool State::TestOr (unsigned N, ...)
{ bool Or = 0;
  va_list V;
  va_start (V, N);
  for (int Index;  ! Or && N && (Index = va_arg (V, int));  N--)
    if (Test (Index))
      Or = 1;
  va_end (V);
  return Or;
}


bool State::ResetAnd (unsigned N, ...)
{ bool And = 1;
  va_list V;
  va_start (V, N);
  for (int Index;  N && (Index = va_arg (V, int));  N--)
    if (! Reset (Index))
      And = 0;
  va_end (V);
  return And;
}


State& State::Reset (void)
{ Data = 0;
  return *this;
}

State& State::Set (void)
{ unsigned char* Q = (unsigned char*) &Data;
  unsigned char* P = Q + Bytes;
  while (--P >= Q)
    *P = UCHAR_MAX;
  return *this;
}
