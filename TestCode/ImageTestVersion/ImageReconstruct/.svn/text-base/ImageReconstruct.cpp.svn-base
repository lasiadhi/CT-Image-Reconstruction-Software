#include "ImageReconstruct.H"
#include "Interpolator.H"
#include "ImageTest.H"
//#include "Image.H"
#include "Filter.H"
#include "RectMDArray.H"
#include "PPoint.H"
#include <string>
#include <cstdlib>
#include <cmath>
using namespace std;
//ImageReconstruct::ImageReconstruct(Image & a_image, int a_filter, int a_interp, int a_pixel)
ImageReconstruct::ImageReconstruct(ImageTest & a_image, int a_filter, int a_interp, int a_pixel)
{
  m_array = a_image.randonTransform();// record the result of Randon transformation
  m_M = a_image.m_M;
  m_d = a_image.m_d;
  m_thetaSize = a_image.m_thetaSize;
  m_pixel = a_pixel;
  m_h = 2.0/(m_pixel-1);
  int l[2]={0, 0};
  int h[2]={m_pixel-1, m_pixel-1};
  Box bx = Box(PPoint(l), PPoint(h));
  RectMDArray<double> matrix(bx);
  m_matrix = matrix;

  m_filterptr = Filter::filterFactory(a_filter,m_M, m_d);
  m_interpptr = Interpolator::interpFactory(a_interp);
}
ImageReconstruct::~ImageReconstruct()
{
  delete m_filterptr;
  delete m_interpptr;
}

void ImageReconstruct::construct()
{
  convolution();
  backProjection();
}

void ImageReconstruct::convolution()
{
  m_filterptr->inverseFFT();
  Box bx = m_array.getBox();
  RectMDArray<double> conv(bx);
  conv.setVal(0.0);
  double * radonptr = m_array.getPointer();
  double * convptr = conv.getPointer();
  for(int thetaIndex=0; thetaIndex<m_thetaSize; thetaIndex++)
    {
      for(int m1=-m_M; m1<=m_M; m1++)
        {
           int mIndex=m_M-m1;
           int linearIndexConv = thetaIndex +mIndex*m_thetaSize;
           double convM1Theta = 0.0;  
           for(int j = -m_M; j<=m_M; j++)
             {
               int tIndex = m_M - j; // tIndex= when t=-1
               int filterIndex = 2*m_M -(m1-j);
               int linearIndexRadon = thetaIndex + tIndex*m_thetaSize;
               convM1Theta += m_filterptr->m_filter[filterIndex]*radonptr[linearIndexRadon];
             }
           convptr[linearIndexConv] = convM1Theta;
        }
    }
  m_array = conv;
}

void ImageReconstruct::backProjection()
{
  vector<double> tValues;
  tValues.resize(2*m_M+1);
  for(int j=-m_M; j <= m_M; j++)
    tValues[m_M-j]=j*m_d;
  
  m_matrix.setVal(0.0); // initial the output imageMatrix
  double * radonptr = m_array.getPointer();
  double * matrixptr = m_matrix.getPointer();
  for(int y= 0; y<m_pixel; y++)
    {
      for(int x= 0; x<m_pixel; x++)
        {
          double xi=-1+x*m_h, yi=1-y*m_h;
          double sum=0.0;
          for (int thetaIndex = 0; thetaIndex< m_thetaSize; thetaIndex++)
            {
              double angle = thetaIndex*M_PI/m_thetaSize;
              double t = xi*cos(angle)+yi*sin(angle);
              sum += m_interpptr->Interpolation(radonptr, tValues, t, thetaIndex, m_thetaSize);
            }
          int linearIndex = x + y*m_pixel;
          matrixptr[linearIndex] = sum /m_thetaSize;
         }
    }
}
