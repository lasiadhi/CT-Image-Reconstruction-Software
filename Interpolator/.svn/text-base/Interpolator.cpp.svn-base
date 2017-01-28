#include "Interpolator.H"
#include "Interps.H"
#include <stddef.h>

using namespace std;

Interpolator * Interpolator::interpFactory(int a_interp)
{
  if(a_interp==1) return new linearInterpolator;
  else if(a_interp==2) return new nearestNeighborInterpolator;
  else return NULL;
}
