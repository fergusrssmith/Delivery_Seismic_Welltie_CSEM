// File   : DE.CPP
// Author : Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>
//           P.O. Box 19065
//           Birmingham, AL 35209
// Notes  : Based on DE.C by Dr. Rainer Storn and Kenneth Price
//          Retained contents of most original comments for reference
// Version: 1.00 <3.7>
// Date   : 09/05/96
// Update : 09/05/96 Redesigned using OT, all data/ops to objects/methods
//          09/05/96 Changed all float vars to double, faster solutions?
//          MM/DD/YY See Modified Original Header, end section


#include "de37.hpp"
#include "rng.hpp"
#include "timer.hpp"
#include "evolver.hpp"
#include "err.hpp"


/* Modified Original Header Information                              <MK>

 D I F F E R E N T I A L     E V O L U T I O N

 Program: de.c
 Version: 3.6     

 Authors: Dr. Rainer Storn
          c/o ICSI, 1947 Center Street, Suite 600
          Berkeley, CA 94707
          Tel   : 510-642-4274 (extension 192)
          Fax   : 510-643-7684
          E-mail: storn@icsi.berkeley.edu
          WWW   : http://http.icsi.berkeley.edu/~storn/
          on leave from
          Siemens AG, ZFE T SN 2, Otto-Hahn Ring 6
          D-81739 Muenchen, Germany
          Tel   : 636-40502
          Fax   : 636-44577
          E-mail: rainer.storn@zfe.siemens.de

          Kenneth Price
          836 Owl Circle
          Vacaville, CA 95687
          E-mail: kprice@solano.community.net

 This program implements some variants of Differential Evolution (DE)
 as described in part in the techreport tr-95-012.ps of ICSI.
 Available via ftp.icsi.berkeley.edu/pub/techreports/1995/tr-95-012.ps.Z
 or via WWW: http://http.icsi.berkeley.edu/~storn/litera.html*
 A more extended version of tr-95-012.ps has been submitted for
 publication in the Journal of Evolutionary Computation.

 You may use this program for any purpose, give it to any person or
 change it according to your needs as long as you are referring to
 Rainer Storn and Ken Price as the originators of the DE idea <and MK
 for the C++ implementation>.

 If you have questions concerning DE feel free to contact us.
 We also will be happy to know about your experiences with DE
 and your suggestions for improvement.

  #  Ver   Date    Request    Modification                         Author
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  1  3.1  5/18/95     -       Strategy DE/Rand-to-Best/1           Storn
                              included
  1  3.2  6/06/95  C.Fleiner  Change loops into memcpy             Storn
  2  3.2  6/06/95     -       Update comments                      Storn

  1  3.3  6/15/95  K.Price    Included strategy DE/Best/2          Storn
  2  3.3  6/16/95     -       Comments and beautifying             Storn
  3  3.3  7/13/95     -       Upper and lower bound for            Storn
                              initialization
  1  3.4  2/12/96     -       Increased printout prec.             Storn

  1  3.5  5/28/96     -       Revisited strategies                 Storn
  2  3.5  5/28/96     -       Included strategy DE/Rand/2          Storn

  1  3.6  8/06/96  K.Price    Added Binomial Crossover             Storn

  1  3.7  9/05/96     -       C++ redesign/reimplementation        Keenan
                              as Controls, RNG, Vectors, Evolver.
                              Converted to double precision
  2  3.7  9/06/96     -       Minimum halting criterion            Keenan
  3  3.7  9/13/96     -       Halting criteria Time, Stasis,       Keenan
                              Generations & all combinations
                              (incl. Minimum)
                      -       Integrated helper classes, e.g.,
                              Messages (usese STL string), Errors,
                              Timer, Evolver derived from Controls
                              derived from State
*/


main (int argc, char *argv[])
// Arguments...
//   Argument count, argument strings
// Usage...  DE <ConfigFilePath> [OutputFilePath]
//   ConfigFile: Viable inputs according to classes Controls and Evolver.
//   Can contain OutputFilePath.
// Role(s)...
//   DE driver with periodic (Refresh) console report (currently managed by
//   Evolver but should be independent of/external to component classes...
//   FIX: Should use access methods?).
//   If no errors occur before or during evolution, results go to
//   OutputFilePath (also managed by Evolver).
{ int Outcome = 0;

  Evolver DE (argc, argv);
  if (DE.Viable())
  { if (DE.Test (Controls::DoTime))
      Time.Start (DE.GetTime());
    while (DE.EvolveInterval() && ! kbhit())
    { DE.Report();
#     ifdef CONIO
        cout << "\n\nPress any key to abort...";
#     endif
    }
    DE.Report().Record();
  }

  if (Exception.Occurred())
  { for (string* E = Exception.First();  E;  E = Exception.Next())
      cout << '\n' << *E;
    Outcome = 1;
  }
# ifdef CONIO
    cout << "\n\nDone.  Press any key...\n";
    while (! kbhit());
# endif
  return Outcome;
}





