#include <iostream>
#include <string>
#include "RectMDArray.H"
#include "PPoint.H"
#include "Image.H"
#include "ImageReconstruct.H"
#include "Filter.H"
#include "Interpolator.H"

using namespace std;


int main( int argc, const char** argv )
{
  // Name of the file to be reconstructed.
  string filename;

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
      cout << "\nYou can reconstruct the image with use of 363 beams if you set beam space < 0.008, otherwise it will use lower number of beams.\n";
      filename = "./phantoms/shepp256.pgm";
      break;
      }
    case 2: 
      {
	cout << "\nYou can reconstruct the image with use of 283 beams if you set beam space < 0.01, otherwise it will use lower number of beams.\n";
	filename = "./phantoms/shepp200.pgm";
	break;
      }
    case 3: 
      {
	cout << "\nYou can reconstruct the image with use of 143 beams if you set beam space < 0.02, otherwise it will use lower number of beams.\n";
	filename = "./phantoms/shepp100.pgm";
	break;
      }
    }
  cout << "Please input the beam spacing: ";
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
  
 
  img.randonTransform();
  std::string outfile1 = "./radon Transform.pgm";
  img.WriteImage(outfile1);

  ImageReconstruct imageRecon(img, filterInd, InterpInd, pixelsX);
  imageRecon.construct();
 
  std::string outfile2 = "./Convolution.pgm";
  imageRecon.WriteConvolution(outfile2);

  std::string outfile3 = "./reconstructedImage.pgm";
  imageRecon.WriteImage(outfile3);
  std::cout << "\nDone writing image... Check main folder for reconstructed image.\n";

  return 0;
}

