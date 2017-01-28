## Discrete Image Reconstruction using Parallel Beam Geometry
# A Software for computerized axial tomography (CAT or CT) based on parallel beam geometry implemented using C++. 

# How to run: 

The main.cpp is in the main folder, run "make all" will automatically generate the main.exe which the excutable file for our project. 

The main.exe will take the following input information:

1) ask the user to choose the filename for the image. There are three images to choose and they have different dimensions.

2) ask the user to input the beam spacing to do the radon transformation. There is a minimum spacing due to pixels of the image. Giving the beam spacing equals to choose how many beams to do the radon transformation. The smaller the beam spacing the better resolution for the image back projection. (e.g: 0.005)

3) ask the user to choose the filter and Interpolator.

4) main.exe won't ask the user to input the grid spacing for reconstructed image because we already choose the same pixel size as the original image.

5) At end of the execution, the pgm files for the results of Radon transform, convolution and backprojection will be generated so that they could be compared with the original figures (see Phantoms folder) and the results of other programs.

Several different test cases can be found in the Chapter 4 of Final Report. 

Run "make clean" will clear *.o *.exe *.d *.pgm from all subfolders.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Folders and their contents:
  main: main.cpp and makefiles
  Image: Image Class and its implementation.
  ImageReconstruct: ImageReconstruct Class and its implementation.
  Filter: Filter Class and its implementation.
  Interpolator: Interpolator Class and its implementation.
  RectArray: RectMDArray class
  timer: CH_TIMER class
  FinalReportLatex: Latex files for the report
  TestCode: Older version of the codes
  main\phantom : original images to be reconstructed. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

We used UBUNTU 12.04 and g++ compiler 4.6.4 to compile and test this c++ project.

Thank you.

Best,
Dan Hu / Lasith Adhikari / Garnet Vaz
(Dec, 2013)
