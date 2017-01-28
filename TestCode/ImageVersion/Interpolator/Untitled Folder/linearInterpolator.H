#ifndef _LINEARINTERPOLATOR_H_
#define _LINEARINTERPOLATOR_H_
#include <cmath>
#include "Interpolator.H"

using namespace std;

class linearInterpolator:public Interpolator
{
public:
  linearInterpolator(){};
  virtual double Interpolation(double* a_radonptr, vector<double> & a_tArray, double a_t, int thetaIndex, int thetaSize) const
  {
    double d = a_tArray[0]-a_tArray[1];
    int M = (a_tArray.size()-1)/2;
    int tlowIndex = floor(a_t/d);
    int LowIndex = M - tlowIndex;
    int HighIndex = M - tlowIndex -1;
    int linearLowIndex = thetaIndex + thetaSize*(LowIndex);
    if(a_tArray[LowIndex]==a_t)
      return a_radonptr[linearLowIndex];
    int linearHighIndex = thetaIndex + thetaSize*(HighIndex);
    double slope = (a_radonptr[linearHighIndex]-a_radonptr[linearLowIndex] )/(a_tArray[HighIndex]-a_tArray[LowIndex]);
    double Delta = a_t - a_tArray[LowIndex];
    return slope * Delta + a_radonptr[linearLowIndex];
  }
};  
#endif
