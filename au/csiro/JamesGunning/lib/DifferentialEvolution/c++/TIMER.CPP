// File..........:  Timer.C
// Content.......:  Time management
// Author........:  Mikal Keenan <keenan@cis.uab.edu>
//                   <http://www.cis.uab.edu/info/grads/mbk
// Update........:  03/24/95


#include <stdio.h>
#include <time.h>
#include <math.h>


#include "timer.hpp"


Timer Time;

Timer::Timer (void)
{ ClockStart = 0;
  ClockStop = 0;
  Interval = 0;
  Interim = 0;
}


Timer::Timer (unsigned long Seconds)
{ ClockStart = 0;
  ClockStop = 0;
  Interim = 0;
  Start (Seconds);
}


void Timer::Start (void)
{ ClockStop = time (&ClockStart);
}


void Timer::Start (const unsigned long Time)
{ Interval = Time;
  Start();
}


void Timer::Restart (void)
{ ClockStart += (time (NULL) - ClockStop);
}


void Timer::Current (unsigned& H, unsigned& M, unsigned& S)
{ time_t Lapse = Interim % 86400L;
  H =  Lapse / 3600;
  M = (Lapse % 3600L) / 60;
  S =  Lapse % 60;
}


char* Timer::Current (void)
{ unsigned H, M, S;
  Current (H, M, S);
  sprintf (Clock, "%02d:%02d:%02d", H, M, S);
  return Clock;
}


int Timer::Remains (void)
{ return Stop() < Interval;
}


int Timer::Elapsed (void)
{ return Stop() >= Interval;
}


unsigned long Timer::Stop (void)
{ return Interim = (time_t) difftime (time (&ClockStop), ClockStart);
}


unsigned long Timer::Final (void)
{ return Interim = (time_t) difftime (ClockStop, ClockStart);
}


