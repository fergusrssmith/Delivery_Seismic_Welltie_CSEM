// File..........:  State.hpp
// Content.......:  Header for State class
// Author........:  Mikal Keenan <keenan@cis.uab.edu>
// Update........:  03/24/95


#ifndef  StateHPP
#define  StateHPP


#include <stdarg.h>



class      State
{         public:  State (void);

           State&  operator = (const State&);

            bool   Set (unsigned, int);
            bool   Set (unsigned);
            bool   Reset (unsigned);
            bool   Test (unsigned);
            bool   SetAnd (unsigned, ...);
            bool   TestAnd (unsigned, ...);
            bool   TestOr (unsigned, ...);
            bool   ResetAnd (unsigned, ...);
           State&  Reset (void);
           State&  Set (void);

         typedef   unsigned long StateData;

         private:
inline StateData   MaskBit (unsigned);
 static unsigned   BitLimit;
 static unsigned   Bits;
 static unsigned   Bytes;
       StateData   Data;
};


#endif
