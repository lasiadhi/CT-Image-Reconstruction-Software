#ifndef _ERROR_H_
#define _ERROR_H_
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/opencv.hpp>
#include<RectMDArray.H>
#include<random>
#include<cassert>
class Error
{
public:
  Error();
  Error(cv::Mat& a_img);
  Error(RectMDArray<int>& a_true);
  RectMDArray<int> addNoise(double a_p, double a_scale); // Bernoulli / poisson params
  double computeError(RectMDArray<int>& a_recon);
  void printMat(); 		// Debugging
private:
  RectMDArray<int> m_true;
};

Error::Error()
{
}

Error::Error(RectMDArray<int>& a_true)
{
  m_true = a_true;
}

void Error::printMat()
{
  using namespace std;
  cout << "bah" << endl;
  Box bx = m_true.getBox();
  int nx = bx.getHighCorner()[0];
  int i=0;
  for(PPoint pt=bx.getLowCorner(); bx.notDone(pt); bx.increment(pt))
    {
      if (i%nx == 0) cout << endl;
      cout << m_true[pt] << "    ";
      ++i;
    }
}

Error::Error(cv::Mat& a_img)
{
  using namespace std;
  unsigned char *input = (unsigned char*)(a_img.data);
  int nx = a_img.rows;
  int ny = a_img.cols;
  int low[2] = {0,0};
  int high[2] = {nx,ny};
  PPoint l = PPoint(low);
  PPoint h = PPoint(high);
  Box bx = Box(l,h);
  RectMDArray<int> img(bx);
  PPoint pt;
  cv::Scalar sc;
  for(int y=0; y<ny; ++y)
    {
      for(int x=0; x<nx; ++x)
  	{
  	  int kkl =  x+y*ny;
  	  pt = bx.getPoint(kkl);
  	  sc = a_img.at<uchar>(y,x);
  	  img[pt] = (int) sc.val[0];
  	}
    }
  m_true = img;
}
RectMDArray<int> Error::addNoise(double a_p, double scale)
{
  using namespace std;
  Box bx = m_true.getBox();
  random_device rd;
  default_random_engine bern_gen(rd());
  default_random_engine pois_gen(rd());
  bernoulli_distribution bern_dist(a_p);
  poisson_distribution<int> pois_dist(scale);
  int mult;
  RectMDArray<int> noisy(bx);
  noisy.setVal(0);
  for (PPoint pt=bx.getLowCorner(); bx.notDone(pt); bx.increment(pt))
    {
      noisy[pt] = m_true[pt];
      if(bern_dist(bern_gen))
	{
	  mult = bern_dist(bern_gen) ? 1 : -1;
	  int bah = pois_dist(pois_gen);
	  cout << "mult " << mult <<  " bah " << bah << endl;
	  noisy[pt] += mult*pois_dist(pois_gen);
	}
    }
  // Fix the max and min limits
  for (PPoint pt=bx.getLowCorner(); bx.notDone(pt); bx.increment(pt))
    {
      noisy[pt] = (noisy[pt] > 255) ? 255 : noisy[pt];
      noisy[pt] = (noisy[pt] < 0) ? 0 : noisy[pt];
    }
  return noisy;
}

double Error::computeError(RectMDArray<int>& a_recon)
{
  Box bxt = m_true.getBox();
  Box bxr = a_recon.getBox();
  assert(bxt == bxr);
  int err = 0;
  for(PPoint pt = bxr.getLowCorner(); bxr.notDone(pt); bxr.increment(pt))
    {
      err += (m_true[pt] - a_recon[pt])*(m_true[pt] - a_recon[pt]);
    }
  return (double) err;
}

#endif
