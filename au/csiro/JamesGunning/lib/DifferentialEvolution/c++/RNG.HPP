// File  : RNG.HPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>


#ifndef RNGHPP
#define RNGHPP


class           RNG
// HISTORY    : random_uniform
// DESCRIPTION: Generates an equally distributed random number
//                  in the interval [0, 1].
// AUTHOR     : Knuth?  See Numerical Recipes in C, Press et al
//              Cambridge University Press, 1992.
//                  OO class by Mikal Keenan <keenan@cis.uab.edu>
//                                         <kalki@compuserve.com>
{      public:  RNG (void);
                RNG (long&);

          RNG&  Init (void);
          RNG&  Init (long&);
       double   Next (void);
       double   Next (double&, double&);
          int   Next (int&);
                operator double (void);

  static long   IM1, IMM1, IM2;
  static long   IA1, IA2;
  static long   IQ1, IQ2;
  static long   IR1, IR2;

static double   NDiv, AM;
static double   Eps, RnMx;

         enum { Count = 32 };

      private:

         long   Seed;
         long   X[Count];
         long   Y, Z;
};


extern RNG Random;


#endif
