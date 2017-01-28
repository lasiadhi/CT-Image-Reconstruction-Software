#ifndef _INTERPS_H_
#define _INTERPS_H_
#include <cmath>
#include "Interpolator.H"

using namespace std;

class linearInterpolator:public Interpolator
{
public:
  double Interpolation(double * a_radonptr, vector<double> & a_tArray, 
                       double a_t, int thetaIndex, int thetaSize) const;
  double tent(double a_x) const;

};  

class nearestNeighborInterpolator: public Interpolator
{
public:
  double Interpolation(double * a_radonptr, vector<double> & a_tArray, 
                       double a_t, int thetaIndex, int thetaSize) const;
  double squareHalf(double a_x) const;
};

#endif
