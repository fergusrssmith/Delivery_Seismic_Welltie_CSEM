// File  : VEC.HPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>


#ifndef VECHPP
#define VECHPP


#include <fstream.h>


class   Vector
{       public:   Vector (void);
                  Vector (Vector&);
                  Vector (int&);
                  Vector (int&, double&, double&);
                 ~Vector (void);

          bool    Init (int&);
        Vector&   Init (double&, double&);
          bool    Init (int&, double&, double&);
        Vector&   Copy (Vector&);
        Vector&   operator = (Vector&);

        double    Evaluate (void);
          void    Report (void);
      ofstream&   Report (ofstream&);


           void   DEBest1Exp (Vector&, Vector&, Vector&, double&, double&);
           void   DERand1Exp (Vector&, Vector&, Vector&, double&, double&);
           void   DERandToBest1Exp (Vector&, Vector&, Vector&, double&, double&);
           void   DEBest2Exp (Vector&, Vector&, Vector&, Vector&, Vector&, double&, double&);
           void   DERand2Exp (Vector&, Vector&, Vector&, Vector&, Vector&, double&, double&);

           void   DEBest1Bin (Vector&, Vector&, Vector&, double&, double&);
           void   DERand1Bin (Vector&, Vector&, Vector&, double&, double&);
           void   DERandToBest1Bin (Vector&, Vector&, Vector&, double&, double&);
           void   DEBest2Bin (Vector&, Vector&, Vector&, Vector&, Vector&, double&, double&);
           void   DERand2Bin (Vector&, Vector&, Vector&, Vector&, Vector&, double&, double&);

        double&   operator [] (int);
                  operator double (void);

    static int    MinCount, MaxCount;

       private:
        Vector&   Zero (void);
        Vector&   Void (void);
        
        double*   X;
      unsigned    Count;
        double    Value;
};

#endif
