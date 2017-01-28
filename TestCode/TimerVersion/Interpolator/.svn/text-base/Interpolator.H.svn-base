#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_
#include <vector>
#include <stddef.h>
using namespace std;

class Interpolator
{
public:
  Interpolator(){};
  virtual ~Interpolator(){};
  virtual double Interpolation(double * a_radonptr, vector<double> & a_tArray, double a_t, int thetaIndex, int thetaSize) const =0;
  static Interpolator * interpFactory(int a_interp);

};  
#endif
