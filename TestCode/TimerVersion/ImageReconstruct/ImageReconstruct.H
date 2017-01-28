#ifndef _IMAGERECONSTRUCT_H_
#define _IMAGERECONSTRUCT_H_
//#include "CH_Timer.H"
#include "Interpolator.H"
#include "Filter.H"
#include "Image.H"
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

class ImageReconstruct
{
public:
  RectMDArray<double> m_matrix;
  ImageReconstruct(Image & a_image, int a_filter, int a_interp, int a_pixel);
  ~ImageReconstruct();
  void construct();
  void WriteImage(std::string a_fileName);
  void WriteConvolution(std::string a_fileName);
  void WriteArray(std::string a_fileName, RectMDArray<double> a_Array);
  void getBounds(RectMDArray<double>& a_array,double & a_min, double & a_max);
 
private:
  void convolution();
  //Interpolates the value of the function at the point (xcos(\theta)+ysin(\theta), \theta) based on the values of the function at all points of (t, \theta) with fixed theta 
  void backProjection(); //update m_array using the convolution.

  RectMDArray<double> m_array;
  double m_d; // Delt_t separation between neighbor beams
  int m_M;
  int m_thetaSize;
  int m_pixel;
  double m_h; // the Delta for pixel
  Filter * m_filterptr;
  Interpolator * m_interpptr;
};

#endif
