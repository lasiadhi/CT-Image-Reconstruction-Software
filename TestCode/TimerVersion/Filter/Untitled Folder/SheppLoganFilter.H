#ifndef _SHEPPLOGANFILTER_H
#define _SHEPPLOGANFILTER_H
#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include "Filter.H"
using namespace std;

class SheppLoganFilter: public Filter
{
public:
  SheppLoganFilter(int a_size, double a_L): Filter(a_size, a_L) { }; //constructor of the filter initializing m_filter to zeros.
  void inverseFFT()
  {
    for(int j=-m_size; j <= m_size; j++)
      {
        int index = m_size - j;
        m_filter[index] = 4.0*m_L*m_L/(M_PI*M_PI*M_PI*(1-4*j*j));
      }
  }
  void initial(){};
};
#endif
