// File..........:  Timer.H
// Content.......:  Header for Timer class
// Author........:  Mikal Keenan <keenan@cis.uab.edu>
//                   <http://cis.uab.edu/info/grads/mbk
// Update........:  02/01/96


#ifndef TIMERHPP
#define TIMERHPP


#include <time.h>


class       Timer
{          public:  Timer (void);
                    Timer (unsigned long);

             void   Start (void);
             void   Start (const unsigned long);
             void   Restart (void);
             void   Current (unsigned&, unsigned&, unsigned&);
             char*  Current (void);
              int   Remains (void);
              int   Elapsed (void);
    unsigned long   Stop (void);
    unsigned long   Final (void);

        protected:

           time_t   Interim;
           time_t   Interval;
           time_t   ClockStart;
           time_t   ClockStop;
             char   Clock[9];
};


extern Timer Time;


#endif
