#include "pre8.c"

void r1024(register real *a)
{
  rpass(a,d1024,128);
  r512(a);
  c256((complex *)(a + 512));
}
