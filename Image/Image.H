#ifndef _IMAGE_H_
#define _IMAGE_H_
#include "RectMDArray.H"
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
// #include "CH_Timer.H"
using namespace std;
class Image
{
public:
  ///
  Image(std::string a_fileName, double a_d);
  void radonTransform();
  //RectMDArray<double>& randonTransform();
  void WriteImage(std::string a_fileName);
  void getBounds(RectMDArray<double>& a_array,double & a_min, double & a_max);
 //  // void compareImages(Mat & a_matrix);
  int getXsize() {return m_xsize;}
  int getYsize() {return m_ysize;}
  int getM() {return m_M;}
  int getThetaSize() {return m_thetaSize;}
  double getd() {return m_d;}
  double gethalfLen() {return m_rc;}
  RectMDArray<double> getRadon() {return m_work;}
  std::string m_fileName; // Record the name of phantom image file.
 private:
  RectMDArray<double> m_original;
  RectMDArray<double> m_work;
  int m_xsize; // Record the size in the x direction of the matrix data.
  int m_ysize; // Record the sizx in the y direction of the matrix data.
  int m_thetaSize;// Record the number of rotations for the image.
  int m_M; // Record the half number of beams except zero position beam.
  double m_d; // Record the beam spacing in units
  double m_dp; // Record the beam spacing in pixels
  double m_rc; // Record the center of the image along diagonal
};
#endif
