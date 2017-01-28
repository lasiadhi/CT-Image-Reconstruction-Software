#ifndef _RamLakFILTER_H
#define _RamLakFILTER_H
#include <cmath>
#include "Filter.H"
using namespace std;

class RamLakFilter: public Filter
{
public:
  RamLakFilter(int a_size, double a_L): Filter(a_size, a_L) { }; //constructor of the filter initializing m_filter to zeros.
 void inverseFFT()
  {
    for(int j=-m_size; j <= m_size; j++)
      { 
	if( (j % 2) > 0) // odd j wil have nonzero values
        {
          int index = m_size - j;
          m_filter[index] = -2.0*m_L*m_L/(M_PI*M_PI*M_PI*j*j);
        }
      }
  }
  void initial(){};
};
#endif
