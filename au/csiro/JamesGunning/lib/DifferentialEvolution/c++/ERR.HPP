// File  : ERR.HPP
// Author: Mikal Keenan <keenan@cis.uab.edu>
//           <http://cis.uab.edu/info/grads/mbk> <kalki@compuserve.com>


#ifndef ERRORSHPP
#define ERRORSHPP


#include "msg.hpp"


class     Errors:  public Messages
{         public:  Errors (void);
                  ~Errors (void);

            bool   Post (unsigned);
            bool   Post (unsigned, char*);
            bool   Post (char*);
            bool   Test (unsigned, char*);
            bool   Test (unsigned);
            bool   Occurred (void);

         private:

            char*  Text (unsigned);

            enum { BufferSize = 255 };
            char   Buffer[BufferSize];
};

extern    Errors   Exception;

#endif
