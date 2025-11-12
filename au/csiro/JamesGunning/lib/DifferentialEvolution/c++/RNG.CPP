// File  : RNG.CPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>

#include <time.h>

#include "rng.hpp"

RNG Random;

long RNG::IM1  = 2147483563;
long RNG::IMM1 = IM1 - 1;
long RNG::IM2  = 2147483399;
long RNG::IA1  = 40014;
long RNG::IA2  = 40692;
long RNG::IQ1  = 53668;
long RNG::IQ2  = 52774;
long RNG::IR1  = 12211;
long RNG::IR2  = 3791;

double RNG::NDiv = 1 + IMM1 / Count;
double RNG::AM   = 1.0 / IM1;
double RNG::Eps  = 1.2e-7;
double RNG::RnMx = 1.0 - Eps;

RNG::RNG (void)
{ Init();
}

RNG::RNG (long& SeedValue)
{ Init (SeedValue);
}

RNG& RNG::Init (void)
{ Init (time (0));
  return *this;
}

RNG& RNG::Init (long& SeedValue)
{ if (! SeedValue)
    return Init();
  Z = (Seed = SeedValue > 0? SeedValue : 1);
  for (long I = Count + 7;  I >= 0;  I--)
  { long J = Seed / IQ1;
    if ((Seed = IA1 * (Seed - J * IQ1) - J * IR1) < 0)
      Seed += IM1;
    if (I < Count)
      X[I] = Seed;
  }
  Y = X[0];
  return *this;
}

double RNG::Next (void)
{ long I = Seed / IQ1;
  if ((Seed = IA1 * (Seed - I * IQ1) - I * IR1) < 0)
    Seed += IM1;

  I = Z / IQ2;
  if ((Z = IA2 * (Z - I * IQ2) - I * IR2) < 0)
    Z += IM2;

  Y = X[I = Y / NDiv] - Z;
  X[I] = Seed;
  if (Y < 1)
    Y += IMM1;

  double Temp = AM * Y;

  return Temp > RnMx?  RnMx : Temp;
}

double RNG::Next (double& Min, double& Max)
{ return Min + Next() * (Max - Min);
}

int RNG::Next (int& Range)
{ return Next() * Range;
}

RNG::operator double (void)
{ return Next();
}

