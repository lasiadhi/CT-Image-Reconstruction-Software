#include <cmath>
#include "Interps.H"
#include "CH_Timer.H"
using namespace std;

inline double linearInterpolator::Interpolation(double * a_radonptr, vector<double> & a_tArray, 
                                         double a_t, int thetaIndex, int thetaSize) const
 {
  CH_TIMERS("Linear Interpolation");
  double a_d = a_tArray[0]-a_tArray[1];
  int m = (a_tArray.size()-1)/2;

  double tot = 0;
  for (int m1 = -m ; m1 <= m ; m1++)
    {
      double temp = tent((a_t/a_d)-m1);
      int linearIndex = thetaIndex + thetaSize*(m-m1);
      tot += a_radonptr[linearIndex] * temp;
    }
  return tot;
 }

inline double linearInterpolator::tent(double a_x) const
  {
    //CH_TIMERS("tent");
    if(abs(a_x) <= 1)
      {
        return 1-abs(a_x);
      }
    else
      {
        return 0;
      }
  }

inline double nearestNeighborInterpolator::Interpolation(double * a_radonptr, vector<double> & a_tArray,
                                                  double a_t, int thetaIndex, int thetaSize) const
 {
  CH_TIMERS("Nearest Neighbor Interpolation");
  double a_d = a_tArray[0]-a_tArray[1];
  int m = (a_tArray.size()-1)/2;

  double tot = 0;
  for (int m1 = -m ; m1 <= m ; m1++)
    {
      int linearIndex = thetaIndex + thetaSize*(m-m1);
      double temp = squareHalf((a_t/a_d)-m1);
      tot += a_radonptr[linearIndex] * temp;
    }
  return tot;
 }

inline double nearestNeighborInterpolator::squareHalf(double a_x) const
  {
    //CH_TIMERS("squareHalf");
    if(abs(a_x) < 0.5 || a_x == 0.5)
      {
        return 1.0;
      }
    else
      {
        return 0.0;
      }
  }
