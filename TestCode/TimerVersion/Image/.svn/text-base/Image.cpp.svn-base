#include "Image.H"
#include "RectMDArray.H"
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

void Image::getBounds(RectMDArray<double>& a_array,double & a_min, double & a_max)
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

void Image::WriteImage(std::string a_fileName)
{
  double max, min;
  // Find max and min to scale the image back to 0-255
  Image::getBounds(m_work,min,max);
  double sgn = min > 0.0 ? -1.0 : 1.0;
  double minp = 0.0;
  double maxp = max + sgn*abs(min);
  double mult = 255.0/abs(maxp);
  std::ofstream myfile(a_fileName);
  if (myfile.is_open())
    {
      Box bx = m_work.getBox();
      myfile << "P2" << std::endl;
      myfile << bx.getHighCorner()[0]+1;
      myfile << ' ';
      myfile << bx.getHighCorner()[1]+1 << std::endl;
      myfile << 255 << std::endl;
      int i = 0;
      int mx = bx.getHighCorner()[0];
      int my = bx.getHighCorner()[1];
      int val;
      double *work = m_work.getPointer();
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

Image::Image(std::string a_fileName, double a_d)
{
  int xs, ys, max;
  std::string line;
  m_fileName = a_fileName;
  std::ifstream myfile;
  myfile.open(a_fileName);
  if (myfile.is_open())
    {
      getline(myfile,line);
      myfile >> xs;
      myfile >> ys;
      myfile >> max;
      m_xsize = xs;
      m_ysize = ys;
      //m_d = a_d;
      //int m = (int)ceil(sqrt(8)/m_d);
      //m_M = (m % 2 == 0) ? m/2 : (m-1)/2;
      int l[2] = {0,0};
      int h[2] = {xs-1,ys-1};
      Box bx = Box(PPoint(l),PPoint(h));
      RectMDArray<double> orig(bx);
      double dmax = (double) max;
      int val;
      for (int i=0; i<xs; ++i)
	{
	  for (int j=0; j<ys; ++j)
	    {
	      int ind = j + i*ys;
	      myfile >> val;
	      orig[ind] = ((double)val)/dmax;
	    }
	}
      myfile.close();
      m_original = orig;

      m_thetaSize = 180; // set no of scans to 180

      int rhomax = ceil(sqrt((xs*xs) + (ys*ys))); // total pixels along diagonal of the image
      
      double optd = (double)sqrt(8)/rhomax;
      m_d = a_d;
      if (a_d < optd)
	{
	  m_d = optd;
	}
      
      int dp = ceil((m_d/sqrt(8))*rhomax);  // pixel offset
      //m_dp  = (dp == 0) ? 1 : dp;  // set the minimum pixel offset to 1
      m_dp = dp;
      
      int rhomaxt = ceil((double)rhomax/m_dp);
      if (rhomaxt % 2 == 0) // make no of beams to odd 
	{
	  rhomaxt += 1;
	  rhomax += dp;
	}
      m_M = (rhomaxt-1)/2; // no of beams in each side
      m_rc = round(rhomax/2.); // center of the image along diagonal

    }
  else
    {
      std::cout << "Error opening file while reading image\n";
      abort();
    }
}

void Image::radonTransform()
{
  //center of the image
  int m = floor((m_xsize)/2);
  int n = floor((m_ysize)/2);
 
  int rc = m_rc;
  int mt =  m_thetaSize;
  int rhomaxt = 2*m_M + 1;

  // construct the sinogram
  int lo[] = {0,0};
  int hi[] = {mt-1,2*m_M};
  PPoint lowCorner(lo);
  PPoint highCorner(hi);
  Box box(lo,hi);
  RectMDArray<double > sinogram(box);
  sinogram.setVal(0.);

  for (int t = 1; t<= 45; t++) //below 45 degrees,use y as variable
    {
      double costheta = cos(t*M_PI/180.);
      double sintheta = sin(t*M_PI/180.);
      double a = -costheta/sintheta; //y=ax+b

      for (int r = 1; r <= rhomaxt; r++)
  	{
  	  int rho = (1+m_dp*(r-1)) - rc;
  	  double b = rho/sintheta; //y=ax+b
  	  int ymax = min((int)round(-a*m + b),n-1);
          int ymin = max((int)round(a*m + b),-n);

  	  for (int y = ymin; y <= ymax; y++)
  	    {
  	      double x = (y-b)/a;
  	      int xfloor = floor(x); //The integer part of x
  	      double xup = x-xfloor; //The decimals of x
  	      double xlow = 1-xup;
  	      x = xfloor;
  	      x = max(x,(double)-m);
  	      x = min(x,(double)m-2);

  	      int rpt[] = {(int)mt-t,(int)rhomaxt-r}; PPoint respt(rpt);
  	      int pt1[] = {(int)x+m+1,(int)y+n+1}; PPoint f1pt(pt1);
  	      int pt2[] = {(int)x+m+2,(int)y+n+1}; PPoint f2pt(pt2);

  	      //cout<<pt2[0]<<":"<<pt2[1]<<endl;
  	      sinogram[respt] += xlow * m_original[f1pt];
  	      sinogram[respt] += xup * m_original[f2pt];
  	    }
  	}
    }


  for (int t = 46; t<= 90; t++) //below 46-90 degrees
    {
      double costheta = cos(t*M_PI/180.);
      double sintheta = sin(t*M_PI/180.);
      double a = -costheta/sintheta; //y=ax+b
      for (int r = 1; r <= rhomaxt; r++)
  	{
  	  int rho = (1+m_dp*(r-1)) - rc;
  	  double b = rho/sintheta; //y=ax+b
  	  double xmax = min(round((-n-b)/a),(double)m-1);
          double xmin = max(round((n-b)/a),(double)-m);

  	  for (double x = xmin; x <= xmax; x++)
  	    {
  	      double y = a*x + b;
  	      int yfloor = floor(y); //The integer part of y
  	      double yup = y-yfloor; //The decimals of y
  	      double ylow = 1-yup;
  	      y = yfloor;
  	      y = max(y,(double)-n);
  	      y = min(y,(double)n-2);

  	      int rpt[] = {(int)mt-t,(int)rhomaxt-r}; PPoint respt(rpt);
  	      int pt1[] = {(int)x+m+1,(int)y+n+1}; PPoint f1pt(pt1);
  	      int pt2[] = { (int)x+m+1,(int)y+n+2}; PPoint f2pt(pt2);

  	      //cout<<pt2[0]<<":"<<pt2[1]<<endl;
  	      sinogram[respt] += ylow * m_original[f1pt];
  	      sinogram[respt] += yup * m_original[f2pt];
  	    }
  	}
    }


  for (int t = 91; t<= 135; t++) //below 91-135 degrees
    {
      double costheta = cos(t*M_PI/180.);
      double sintheta = sin(t*M_PI/180.);
      double a = -costheta/sintheta; //y=ax+b
      for (int r = 1; r <= rhomaxt; r++)
  	{
  	  int rho = (1+m_dp*(r-1)) - rc;
  	  double b = rho/sintheta; //y=ax+b
  	  double xmax = min(round((n-b)/a),(double)m-1);
          double xmin = max(round((-n-b)/a),(double)-m);
  	  for (double x = xmin; x <= xmax; x++)
  	    {
  	      double y = a*x + b;
  	      int yfloor = floor(y); //The integer part of y
  	      double yup = y-yfloor; //The decimals of y
  	      double ylow = 1-yup;
  	      y = yfloor;
  	      y = max(y,(double)-n);
  	      y = min(y,(double)n-2);

  	      int rpt[] = {(int)mt-t,(int)rhomaxt-r}; PPoint respt(rpt);
  	      int pt1[] = {(int)x+m+1,(int)y+n+1}; PPoint f1pt(pt1);
  	      int pt2[] = {(int)x+m+1,(int)y+n+2}; PPoint f2pt(pt2);

	      //cout<<pt2[0]<<":"<<pt2[1]<<endl;
  	      sinogram[respt] += ylow * m_original[f1pt];
  	      sinogram[respt] += yup * m_original[f2pt];
  	    }
  	}
    }


  for (int t = 136; t<= 179; t++) //below 135-179 degrees
    {
      double costheta = cos(t*M_PI/180.);
      double sintheta = sin(t*M_PI/180.);
      double a = -costheta/sintheta; //y=ax+b
      for (int r = 1; r <= rhomaxt; r++)
  	{
  	  int rho = (1+m_dp*(r-1)) - rc;
  	  double b = rho/sintheta; //y=ax+b
  	  int ymax = min((int)round(a*m+b),n-1);
          int ymin = max((int)round(-a*m+b),-n);
  	  for (int y = ymin; y <= ymax; y++)
  	    {
  	      double x = (y-b)/a;
  	      int xfloor = floor(x); //The integer part of x
  	      double xup = x-xfloor; //The decimals of x
  	      double xlow = 1-xup;
  	      x = xfloor;
  	      x = max(x,(double)-m);
  	      x = min(x,(double)m-2);

  	      int rpt[] = {(int)mt-t,(int)rhomaxt-r}; PPoint respt(rpt);
  	      int pt1[] = {(int)x+m+1,(int)y+n+1}; PPoint f1pt(pt1);
  	      int pt2[] = {(int)x+m+2,(int)y+n+1}; PPoint f2pt(pt2);

	      //cout<<rpt[1]<<":"<<rpt[0]<<endl;
  	       sinogram[respt] += xlow * m_original[f1pt];
  	       sinogram[respt] += xup * m_original[f2pt];
  	    }
  	}
   }

  //printRectMDarray(sinogram);
  m_work = sinogram;
  //return m_work;
}

