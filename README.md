This program generates a 3D surface mesh from a segmentation mask provided in the nifti format.

Compilation for Linux:

The libraries VTK, QT5 and CGAL and Eigen 3 must be installed to compile this program.

Once dependencies are installed, simply execute in the directory of the repository:
  cmake .
  make
  
For windows and MacOs:

No compilation was tested with these platforms, however, the compilation should be possible using cmake.
