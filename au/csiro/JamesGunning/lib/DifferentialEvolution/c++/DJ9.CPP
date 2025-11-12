// File  : DJ9.CPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>

#include "vec.hpp"

double Vector::Evaluate (void)
//Polynomial fitting problem
{  Value = 0.0;

   int N, I, M = 60;

   double xx = -1, px;
   double dx = 2.0 / (double) M;
   double* P;

   for (N = M+1;  N--;  xx += dx)
   {  for (P = X, px = *P++, I = Count;  --I;  P++)
	     px = *P + xx * px;
      if (px < -1 || px > 1)
        Value += (1 - px) * (1 - px);
   }

   for (P = X, px = *P++, I = Count;  --I;  P++)
     px = *P + 1.2 * px;
   if ((px -= 72.661) < 0)
     Value += px * px;

   for (P = X, px = *P++, I = Count;  --I;  P++)
     px = *P - 1.2 * px;
   if ((px -= 72.661) < 0)
     Value += px * px;
     
   return Value;
}
