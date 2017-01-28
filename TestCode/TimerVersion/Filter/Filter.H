#ifndef _FILTER_H_
#define _FILTER_H_
#include <vector>
using namespace std;


class Filter
{
public:
  Filter(int a_m, double a_d){ m_size=2*a_m; m_L = 1.0/(2.*a_d); m_filter.resize(2*m_size+1); } //constructor of the filter initializing m_filter to zeros.
  // a_m = m where t has 2m+1 values
  virtual ~Filter(){};
  virtual void inverseFFT()=0; 
  virtual void initial() = 0;
  static Filter* filterFactory(int a_filter, int a_m, double a_d);

  vector<double> m_filter; // convolution need this data
  int m_size; // the size of filter = 2*m_size +1 = 4*a_size +1
  int m_L;
};


#endif
