// File  : VEC.CPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>


#include <mem.h>


#include "vec.hpp"
#include "rng.hpp"
#include "ctrl.hpp"
#include "err.hpp"


int Vector::MinCount = 1;
int Vector::MaxCount = 35;


Vector::Vector (void)
{ Zero();
}

Vector::Vector (Vector& M)
{ Zero() = M;
}

Vector::Vector (int& NParameters)
{ Zero().Init (NParameters);
}

Vector::Vector (int& NParameters, double& Minimum, double& Maximum)
{ Zero().Init (NParameters, Minimum, Maximum);
}

Vector::~Vector (void)
{ Void();
}

Vector& Vector::Zero (void)
{ Count = 0;
  Value = 0.0;
  X = 0;
  return *this;
}

Vector& Vector::Void (void)
{ if (X)
    delete [] X;
  return Zero();
}

bool Vector::Init (int& NParameters)
{ Void();

  if (NParameters && (NParameters <= MaxCount))
    if (! (X = new double [NParameters]))
      Exception.Post (Controls::NoMemory, "Vector::Init (int&)");
    else
      Count = NParameters;

  return X && Count;
}

bool Vector::Init (int& NParameters, double& Minimum, double& Maximum)
{ Void();

  if (NParameters && (NParameters <= MaxCount))
    if (X = new double [NParameters])
    { Count = NParameters;
      Init (Minimum, Maximum);
    }
    else
      Exception.Post (Controls::NoMemory, "Vector::Init (int&, dbl&, dbl&)");

  return X && Count;
}

Vector& Vector::Init (double& Minimum, double& Maximum)
// Right now this part is kept fairly simple and just generates
// random numbers in the range [Minimum, Maximum]. You might
// want to extend the init part such that you can initialize
// each parameter separately.
{ double* P = X;
  for (int I = Count;  I--;  ++P)
    *P = Random.Next (Minimum, Maximum);
  Evaluate();
  return *this;
}

Vector& Vector::operator = (Vector& V)
{ if (Init (V.Count))
    Copy (V);
  return *this;
}

Vector::operator double (void)
{ return Value;
}

double& Vector::operator [] (int Index)
{ return *(X + (Index % Count));
}

Vector& Vector::Copy (Vector& V)
{ memcpy (X, V.X, sizeof (double) * Count);
  Value = V.Value;
  return *this;
}

void Vector::Report (void)
{ for (int I = 0;  I < Count;  I++)
    cerr << "\n Best[" << I << "] = " << X[I];
}

ofstream& Vector::Report (ofstream& Out)
{ for (int I = 0;  I < Count;  I++)
    Out << "\n Best[" << I << "] = " << X[I];
  return Out;
}

//  STRATEGIES
//  Naming convention is DE/x/y/z
//  DE:  Differential Evolution
//  x :  String denoting the vector to be perturbed
//  y :  Number of difference vectors taken for perturbation of x
//  z :  Crossover method, exp = exponential, bin = binomial
//       Some simple rules worth following...
//   1)  FConstant is usually between 0.5 and 1 (in rare cases > 1)
//   2)  Crossover is in (0..1), 0.0, 0.3, 0.7 and 1.0 are worth trying first.
//   3)  To start off NMembers = 10*NParameters is a reasonable choice.
//       Increase NMembers if misconvergence  happens.                                                                     */
//   4)  If you increase NMembers, FConstant usually has to be decreased
//   5)  When the DE/BestYet... schemes fail DE/Rand... usually works
//       and vice versa

void Vector::DEBest1Exp (Vector& Best, Vector& V1, Vector& V2,
double& FConstant, double& Crossover)
// Strategy DE0 (not in our paper)
// Our oldest strategy but still not bad. However, we have found
// several optimization problems where misconvergence occurs.
{ double *XB = Best.X, *X1 = V1.X, *X2 = V2.X;
  int I = Random.Next (Count), Counter = 0;
  do    X[I] = XB[I] + FConstant * (X1[I] - X2[I]);
  while (++I %= Count, (Crossover >= Random) && (++Counter < Count));
}

void Vector::DERand1Exp (Vector& V1, Vector& V2, Vector& V3,
double& FConstant, double& Crossover)
// Strategy DE1 in the techreport
// This is one of my favourite strategies. It works especially well when
// "GenBest"-schemes experience misconvergence. Try, e.g., FConstant=0.7
// and Crossover=0.5 as a first guess.
{ double *X1 = V1.X, *X2 = V2.X, *X3 = V3.X;
  int I = Random.Next (Count), Counter = 0;
  do    X[I] = X1[I] + FConstant * (X2[I] - X3[I]);
  while (++I %= Count, (Crossover >= Random) && (++Counter < Count));
}

void Vector::DERandToBest1Exp (Vector& Best, Vector& V1,
Vector& V2, double& FConstant, double& Crossover)
// Similar to DE2 but generally better
// This strategy seems to be one of the best. Try FConstant=0.85 and
// Crossover=1.  If you get misconvergence try increasing NMembers. If it
// doesn't help you should play around with all three control variables.
{ double *XB = Best.X, *X1 = V1.X, *X2 = V2.X;
  int I = Random.Next (Count), Counter = 0;
  do    X[I] += FConstant * ((XB[I] - X[I]) + (X1[I] - X2[I]));
  while (++I %= Count, (Crossover >= Random) && (++Counter < Count));
}

void Vector::DEBest2Exp (Vector& Best, Vector& V1, Vector& V2,
Vector& V3, Vector& V4, double& FConstant, double& Crossover)
// Another powerful strategy worth trying
{ double *XB = Best.X, *X1 = V1.X, *X2 = V2.X, *X3 = V3.X, *X4 = V4.X;
  int I = Random.Next (Count), Counter = 0;
  do    X[I] = XB[I] + FConstant * (X1[I] + X2[I] - X3[I] - X4[I]);
  while (++I %= Count, (Crossover >= Random) && (++Counter < Count));
}

void Vector::DERand2Exp (Vector& R5, Vector& V1, Vector& V2, Vector& V3,
Vector& V4, double& FConstant, double& Crossover)
// Seems to be a robust optimizer for many functions
{ double *X5 = R5.X, *X1 = V1.X, *X2 = V2.X, *X3 = V3.X, *X4 = V4.X;
  int I = Random.Next (Count), Counter = 0;
  do    X[I] = X5[I] + FConstant * (X1[I] + X2[I] - X3[I] - X4[I]);
  while (++I %= Count, (Crossover >= Random) && (++Counter < Count));
}

void Vector::DEBest1Bin (Vector& Best, Vector& V1, Vector& V2,
double& FConstant, double& Crossover)
// Strategy DE0 (not in our paper) with binomial crossover
// Our oldest Strategy but still not bad. However, we have found
// several optimization problems where misconvergence occurs.
{ double *XB = Best.X, *X1 = V1.X, *X2 = V2.X;
  for (int Counter = 0, I = Random.Next (Count);  Counter < Count; ++I %= Count)
	 if ((Crossover >= Random) || (++Counter == Count))
	   X[I] = XB[I] + FConstant * (X1[I] - X2[I]);
}

void Vector::DERand1Bin (Vector& V1, Vector& V2, Vector& V3,
double& FConstant, double& Crossover)
// Strategy DE1 in the techreport with binomial crossover
// This is one of my favourite strategies. It works especially well
// when "GenBest"-schemes experience misconvergence. Try e.g. FConstant=0.7
// and Crossover=0.5 as a first guess.
{ double *X1 = V1.X, *X2 = V2.X, *X3 = V3.X;
  for (int Counter = 0, I = Random.Next (Count);  Counter < Count; ++I %= Count)
	 if ((Crossover >= Random) || (++Counter == Count))        // Change at least
	   X[I] = X1[I] + FConstant * (X2[I] - X3[I]);             // one value
}

void Vector::DERandToBest1Bin (Vector& Best, Vector& V1,
Vector& V2, double& FConstant, double& Crossover)
// Similar to DE2 (with binomial crossover) but generally better
// This strategy seems to be one of the best. Try FConstant=0.85 and
// Crossover=1.  If you get misconvergence try increasing NMembers. If it
// doesn't help you should play around with all three control variables.
{ double *XB = Best.X, *X1 = V1.X, *X2 = V2.X;
  for (int Counter = 0, I = Random.Next (Count);  Counter < Count; ++I %= Count)
	 if ((Crossover >= Random) || (++Counter == Count))        // Change at least
	   X[I] += FConstant * ((XB[I] - X[I]) + (X1[I]- X2[I]));  // one value
}

void Vector::DEBest2Bin (Vector& Best, Vector& V1, Vector& V2,
Vector& V3, Vector& V4, double& FConstant, double& Crossover)
// Like DEBest2Exp but with binomial crossover
{ double *XB = Best.X, *X1 = V1.X, *X2 = V2.X, *X3 = V3.X, *X4 = V4.X;
  for (int Counter = 0, I = Random.Next (Count);  Counter < Count; ++I %= Count)
    if ((Crossover >= Random) || (++Counter == Count))        // Change at least
	   X[I] = XB[I] + FConstant * (X1[I] + X2[I] - X3[I] - X4[I]);  // one value
}

void Vector::DERand2Bin (Vector& R5, Vector& V1, Vector& V2, Vector& V3,
Vector& V4, double& FConstant, double& Crossover)
// Like DERand2Exp but with binomial crossover
{ double *X5 = R5.X, *X1 = V1.X, *X2 = V2.X, *X3 = V3.X, *X4 = V4.X;
  for (int Counter = 0, I = Random.Next (Count);  Counter < Count; ++I %= Count)
    if ((Crossover >= Random) || (++Counter == Count))        // Change at least
	   X[I] = X5[I] + FConstant * (X1[I] + X2[I] - X3[I] - X4[I]);  // one value
}
