#ifndef _INTERPS_H_
#define _INTERPS_H_
#include <cmath>
#include "Interpolator.H"

using namespace std;

class linearInterpolator:public Interpolator
{
public:
  double Interpolation(double * a_radonptr, vector<double> & a_tArray, double a_t, int thetaIndex, int thetaSize) const
 {
  double a_d = a_tArray[0]-a_tArray[1];
  int m = (a_tArray.size()-1)/2;

  double tot = 0;
  for (int m1 = -m ; m1 <= m ; m1++)
    {
      int linearIndex = thetaIndex + thetaSize*(m-m1);
      double temp = tent((a_t/a_d)-m1);
      tot += a_radonptr[linearIndex] * temp;
    }
  return tot;
 }

  double tent(double a_x) const
  {
    if(abs(a_x) <= 1)
      {
        return 1-abs(a_x);
      }
    else
      {
        return 0;
      }
  }

};  

class nearestNeighborInterpolator: public Interpolator
{
public:
  double Interpolation(double * a_radonptr, vector<double> & a_tArray, double a_t, int thetaIndex, int thetaSize) const {};
};

#endif
