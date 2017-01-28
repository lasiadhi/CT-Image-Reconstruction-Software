#ifndef _FILTERS_H
#define _FILTERS_H
#include "Filter.H"
#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
using namespace std;

class LowPassCosineFilter: public Filter
{
public:
  LowPassCosineFilter(int a_size, double a_d): Filter(a_size, a_d) { }; //constructor of the filter initializing m_filter to zeros.
  void inverseFFT()
  {
    for(int j=-m_size; j <= m_size; j++)
      {
        int index = m_size - j;
        double factor = 2.0*m_L*m_L/(M_PI*M_PI*M_PI);
        double dem = 1.0-4.0*j*j;
        double term1 = M_PI*cos(M_PI*j)/dem;
        double term2 = 2.0*(4.0*j*j+1.0)/(dem * dem);
        m_filter[index] = factor * (term1 - term2);
      }
  } 
  void initial(){};
};

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
        else if (j==0)
        {
          int index = m_size;
          m_filter[index] = m_L*m_L/(M_PI*M_PI);
        } 
      }
  }
  void initial(){};
};

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
