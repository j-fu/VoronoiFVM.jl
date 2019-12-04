/* This file is not part of the Triangle distribution*/

#include <stdio.h>
#define REAL double


/* Function type for triunsuitable func*/
typedef int (*unsuitable_func)(REAL org_x, REAL org_y, REAL dest_x, REAL dest_y, REAL apex_x, REAL apex_y, REAL area);

/* Trivial default triunsuitable func */
int trivial_triunsuitable(REAL org_x, REAL org_y, REAL dest_x, REAL dest_y, REAL apex_x, REAL apex_y, REAL area)
{
  return 0;
}

/* Static variable containing the actual triunsuitable function.
   This makes things non-threadsafe. An alternative would consist
   of patching Triangle and adding this variable to the triangulateio struct.
*/
static unsuitable_func jl_unsuitable_func=trivial_triunsuitable;


/* Set triunusuitable function */
void triunsuitable_callback(unsuitable_func f)
{
  jl_unsuitable_func=f;
}

/* Function called by Triangle (when compiled with -DEXTERNAL_TEST */
int triunsuitable(REAL* triorg, REAL* tridest, REAL* triapex, REAL area)
{
  return jl_unsuitable_func(triorg[0], triorg[1], tridest[0], tridest[1], triapex[0], triapex[1], area);
}

