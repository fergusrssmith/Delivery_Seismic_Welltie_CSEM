// File  : MSG.HPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>


#ifndef MESSAGESHPP
#define MESSAGESHPP


#include <string>                                                 // STL
#include <fstream.h>


class   Messages
{         public:  Messages (void);
                  ~Messages (void);

            bool   Viable (void);
             int   GetCount (void);
             
            bool   Add (const string&);
        Messages&  Reset (void);
          string*  First (void);
          string*  Next (void);
          string*  Last (void);
          string*  operator [] (unsigned);
          string*  Current (void);
          string*  Current (string*);
          string*  Remove (string*);
          string*  Remove (void);
            bool   Save (ofstream&);

         private:

            bool   Add (string*);
        Messages&  Zero (void);
            void   Void (void);

        unsigned   Count;
          struct   STRING
                   { string*  Data;
                     STRING*  Prev;
                     STRING*  Next;
                   };
          STRING*  Cursor;
          STRING*  Begin;
          STRING*  End;
};


#endif
