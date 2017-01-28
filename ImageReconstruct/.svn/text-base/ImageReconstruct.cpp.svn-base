#include "Image.H"
#include "CH_Timer.H"
#include "ImageReconstruct.H"
#include "Interpolator.H"
#include "Filter.H"
#include "RectMDArray.H"
#include "PPoint.H"
#include <string>
#include <cstdlib>
#include <cmath>
using namespace std;

ImageReconstruct::ImageReconstruct(Image & a_image, int a_filter, int a_interp, int a_pixel)
{
  m_array = a_image.getRadon();// record the result of Randon transformation
  m_M = a_image.getM();
  m_d = a_image.getd();
  m_thetaSize = a_image.getThetaSize();
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
  CH_TIMERS("Convolution");
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
  CH_TIMERS("Back Projection");
  //CH_TIMER("Summation over all angles", t1);
  CH_TIMER("Linear Interpolation", t4);
  CH_TIMER("Linear Index", t2);
  CH_TIMER("calculation of t", t3);
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
          //CH_START(t1);
          for (int thetaIndex = 0; thetaIndex< m_thetaSize; thetaIndex++)
            {
              CH_START(t3);
              double angle = thetaIndex*M_PI/m_thetaSize;
              double t = xi*cos(angle)+yi*sin(angle);
              CH_STOP(t3);
              CH_START(t4);
              sum += m_interpptr->Interpolation(radonptr, tValues, t, thetaIndex, m_thetaSize);
              CH_STOP(t4);
            }
          //CH_STOP(t1);
          CH_START(t2);
          int linearIndex = x + y*m_pixel;
          matrixptr[linearIndex] = sum /m_thetaSize;
          CH_STOP(t2);
	  // Progress bar
	  cout << "\r" << ceil((double)linearIndex/(m_pixel*m_pixel)*100.0) << "% completed.       " << flush;
         }
    }
}

void ImageReconstruct::getBounds(RectMDArray<double>& a_array,double & a_min, double & a_max)
{
  double min = 1e10;
  double max = -1e10;
  Box bx = a_array.getBox();
  double val;
  for (PPoint pt=bx.getLowCorner(); bx.notDone(pt); bx.increment(pt))
    {
      val = a_array[pt];
      if(val < min) min = val;
      if(val > max) max = val;
    }
  a_min = min;
  a_max = max;
}

void ImageReconstruct::WriteArray(std::string a_fileName, RectMDArray<double> a_Array)
{
  double max, min;
  // Find max and min to scale the image back to 0-255
  ImageReconstruct::getBounds(a_Array,min,max);
  double sgn = min > 0.0 ? -1.0 : 1.0;
  double minp = 0.0;
  double maxp = max + sgn*abs(min);
  double mult = 255.0/abs(maxp);
  std::ofstream myfile(a_fileName);
  if (myfile.is_open())
    {
      Box bx = a_Array.getBox();
      myfile << "P2" << std::endl;
      myfile << bx.getHighCorner()[0]+1;
      myfile << ' ';
      myfile << bx.getHighCorner()[1]+1 << std::endl;
      myfile << 255 << std::endl;
      int i = 0;
      int mx = bx.getHighCorner()[0];
      int my = bx.getHighCorner()[1];
      int val;
      double *work = a_Array.getPointer();
      for(int iy=my; iy>=0; --iy)
      	{
      	  for(int ix=mx; ix>=0; --ix)
      	    {
      	      int idx = ix + (mx+1)*iy;
	      double *loc = work + idx;
      	      val = (int) ((*(loc) + sgn*abs(min))*mult);
      	      myfile << val << ' ';
      	      ++i;
      	      if(i%mx == 0) myfile << '\n';
      	    }
      	}
      myfile.close();
    }
  else
    {
      std::cout << "Error opening file\n";
      abort();
    }
}

  void ImageReconstruct::WriteImage(std::string a_fileName) {WriteArray(a_fileName, m_matrix); }
  void ImageReconstruct::WriteConvolution(std::string a_fileName) {WriteArray(a_fileName, m_array); }
