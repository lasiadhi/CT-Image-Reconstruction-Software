#include "opencv2/highgui/highgui.hpp"
#include <opencv2/opencv.hpp>
#include <iostream>
#include "RectMDArray.H"
#include "PPoint.H"
#include <cmath>
#include "ImageTest.H"
#include "ImageReconstruct.H"
#include "Filter.H"
#include "Interpolator.H"

using namespace cv;
using namespace std;


void printRectMDarray(RectMDArray<double > a_recMD);
Mat rectToMat(RectMDArray<double > a_recMD);


int main( int argc, const char** argv )
{


  double d = 0.05;
  int N = 30;
  double di= 0.01;
  int pixel = 2.0/di+1;

  ImageTest imagT(d, N);

  double L = 1.0/(2*d);
  int filterInd = 0;
  int InterpInd = 0;
  ImageReconstruct imageRecon(imagT, filterInd, InterpInd, pixel);
  imageRecon.convolution();
  //printRectMDarray(imageRecon.m_array);
  imageRecon.backProjection();
  Mat reImage = rectToMat(imageRecon.m_matrix);
  //Mat reImage = backProjection(0.01, N, conv,d);
  //writeImg("reconstructed.pgm", reImage);

  
  namedWindow("MyWindow", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
  imshow("MyWindow", reImage); //display the image which is stored in the 'img' in the "MyWindow" window
  waitKey(0); //wait infinite time for a keypress
  destroyWindow("MyWindow"); //destroy the window with the name, "MyWindow"

  return 0;
}


void printRectMDarray(RectMDArray<double > a_recMD)
{
   Box bx = a_recMD.getBox();
   double* recMDptr = a_recMD.getPointer();
   int x1 = bx.getHighCorner()[0];
   int x2 = bx.getHighCorner()[1];

   for (int y = x2;y >= 0; y--)
     {
       for (int x = 0; x < x1+1; x++)
	 {
	   int kkl = x + y*(x1+1);
	   cout << recMDptr[kkl] << "    ";	   
	 }
       cout<<endl;
     }
   
}

Mat rectToMat(RectMDArray<double > a_recMD)
{
   Box bx = a_recMD.getBox();
   double* recMDptr = a_recMD.getPointer();
   int x1 = bx.getHighCorner()[0];
   int y1 = bx.getHighCorner()[1];

   Mat mat = cv::Mat::zeros(y1+1, x1+1, CV_BGR2GRAY);

   for (int y = 0;y < y1+1; y++)
     {
       for (int x = 0; x < x1+1; x++)
	 {
	   int kkl = x + y*(x1+1);
	   mat.at<double>(x,y) =  recMDptr[kkl];
	 }
     }
   return mat;
}



