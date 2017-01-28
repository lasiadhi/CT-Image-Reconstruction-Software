#include "Filter.H"
#include "filters.H"
#include <stddef.h>
using namespace std;

Filter* Filter::filterFactory(int a_filter, int a_m, double a_d)
{
  if(a_filter==0)
   { 
    Filter* temp = new SheppLoganFilter(a_m, a_d);
    return temp;
   }
  else if(a_filter==1) 
    {
     Filter* temp = new RamLakFilter(a_m, a_d);
     return temp;
    }
  else if(a_filter==3) 
    {
     Filter* temp = new LowPassCosineFilter(a_m, a_d);
     return temp;
    }
  else 
    return NULL;
}
