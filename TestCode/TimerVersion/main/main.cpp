#include <iostream>
#include <string>
#include "RectMDArray.H"
#include "CH_Timer.H"
#include "PPoint.H"
#include "Image.H"
#include "ImageReconstruct.H"
#include "Filter.H"
#include "Interpolator.H"

using namespace std;


int main( int argc, const char** argv )
{
  CH_TIMERS("main");
  CH_TIMER("Radon Transformation", t1);
  CH_TIMER("Image ReconStruction", t2);
  // Name of the file to be reconstructed.
  string filename;
  // Name of output file
  string outfile1, outfile2, outfile3;

  cout << "Please input the image file number to be constructed: "<<endl;
  cout << " 1 - Shepp_Logan image with size 256 x 256\n";
  cout << " 2 - Shepp_Logan image with size 200 x 200\n";
  cout << " 3 - Shepp_Logan image with size 100 x 100\n";
  cout << "Your choice: ";
  int no;
  cin >> no;

  switch (no)
    {
    case 1: 
      {
      cout << "\nThe minimum beam spacing for 256 x 256 image is 0.008 (with 363 X ray beams)\n";
      filename = "./phantoms/shepp256.pgm";
      outfile1 = "./RadonTransform256.pgm";
      outfile2 = "./Convolution256.pgm";
      outfile3 = "./ReconstructedImage256.pgm";
      break;
      }
    case 2: 
      {
	cout << "\nThe minimum beam spacing for 200x200 image is 0.01 (with 283 X ray beams)\n";
	filename = "./phantoms/shepp200.pgm";
        outfile1 = "./RadonTransform200.pgm";
        outfile2 = "./Convolution200.pgm";
        outfile3 = "./ReconstructedImage200.pgm";
	break;
      }
    case 3: 
      {
	cout << "\nThe minimum beam spacing for 100x100 image is 0.02 (with 143 X ray beams)\n";
	filename = "./phantoms/shepp100.pgm";
        outfile1 = "./RadonTransform100.pgm";
        outfile2 = "./Convolution100.pgm";
        outfile3 = "./ReconstructedImage100.pgm";
	break;
      }
    }
  cout << "\nPlease input your beam spacing (any value smaller than the minimum spacing will be treated as minimum): ";
  double d;
  cin>> d;
  
  // create Image Class object
  Image img(filename, d);
  
  cout << " \nPlease input the type of the filter: \n";
  cout << " 1 : Shepp-Logan filter\n";
  cout << " 2 : Ram-Lak filter \n";
  cout << " 3 : Low-pass cosine filter \n";
  cout << "Your choice: ";
  int filterInd;
  cin >> filterInd;

  cout << " \nPlease input the type of the Interporlator: \n";
  cout << " 1 : Linear Interpolator\n";
  cout << " 2 : Nearest neighbor Interpolator \n";
  cout << "Your choice: ";
  int InterpInd;
  cin >> InterpInd;

  cout << "\nNow we are using "<<2*img.getM()+1<<" X-ray beams for reconstruction process."<<endl;
  cout << "Please wait ...now image is processing.\n";

  int pixelsX = img.getXsize();
  
  CH_START(t1);
  img.radonTransform();
  CH_STOP(t1);
  img.WriteImage(outfile1);


  ImageReconstruct imageRecon(img, filterInd, InterpInd, pixelsX);
  CH_START(t2);
  imageRecon.construct();
  CH_STOP(t2);
  imageRecon.WriteConvolution(outfile2);
  imageRecon.WriteImage(outfile3);
  std::cout << "\nDone writing image... Check main folder for reconstructed image.\n";

  return 0;
}

