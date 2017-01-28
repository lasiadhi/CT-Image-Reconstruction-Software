#include "ImageTest.H"
#include <iostream>
#include "RectMDArray.H"
#include "PPoint.H"
#include <cmath>


using namespace std;

ImageTest::ImageTest(double a_d, int a_thetaSize)
{ 
  m_M = 1./a_d; 
  m_d = a_d; 
  m_thetaSize = a_thetaSize;  
}

double ImageTest::testRadon(double t)
{
  
  if (abs(t) <= 0.25)
          return (1/2.)*sqrt(9-16*t*t) - (3/4.)*sqrt(1-4*t*t)+(1/8.)*sqrt(1-16*t*t);
  else
    if (abs(t) > 0.25 && abs(t) <= 0.5)
            return (1/2.)*sqrt(9-16*t*t) - (3/4.)*sqrt(1-4*t*t);
    else
      if (abs(t) > 0.5 && abs(t) <= 0.75)
            return (1/2.)*sqrt(9-16*t*t);
      else
	if (abs(t) > 0.75)
            return 0;
       
}

RectMDArray<double > ImageTest::randonTransform()
{
  int a_N = m_thetaSize;
  double a_d = m_d;
  int beams = (2/a_d)+1;
  double tval[beams];

  for (int i = 0; i < beams ; i++)
    {
      tval[i] = -1 + (i* a_d);
      //cout << t[i]<<endl;
    }

  int lo[] = {0,0};
  int hi[] = {a_N-1,beams-1};

  PPoint lowCorner(lo);
  PPoint highCorner(hi);
  Box box(lo,hi);
  RectMDArray<double > sinogram(box);

  double* sinptr = sinogram.getPointer();

  for (int theta = 0 ; theta < a_N ; theta++)
    {
      for (int t = 0 ; t < beams ; t++)
	{
	  int linearIndex = theta + (a_N * ((beams-1)-t));
	  sinptr[linearIndex] = testRadon(tval[t]);
	  //cout<< linearIndex<<endl;
	}
    }

  return sinogram;

}

