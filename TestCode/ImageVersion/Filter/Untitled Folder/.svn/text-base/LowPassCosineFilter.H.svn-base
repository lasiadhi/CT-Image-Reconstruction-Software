#ifndef _LOWPASSCOSINEFILTER_H
#define _LOWPASSCOSINEFILTER_H
#include "Filter.H"

class LowPassCosineFilter: public Filter
{
public:
  LowPassCosineFilter(int a_size, double a_d): Filter(a_size, a_d) { }; //constructor of the filter initializing m_filter to zeros.
  void inverseFFT(); 
  void initial() ;
};
#endif
