# About

This program creates a 3D surface mesh from a segmentation mask provided in NIfTI format.

# Compilation

## Required 3rd party libraries

The libraries below are required to be compiled and installed on the system.

- VTK [BSD License](https://vtk.org/about/#license)
- QT5 [GPL3 and LGPL2.1 license](https://www.qt.io/qt-licensing)
- CGAL [GPL3 and LGPL3 license](https://www.cgal.org/license.html)
- Eigen3 [MPL2 license](https://www.mozilla.org/en-US/MPL/2.0/)

## For Linux:

Simply execute in the directory of the repository:  
`cmake .`  
`make` 

This was tested with VTK 9, QT5 5.15.3, CGAL 5.4.1 and Eigen3 3.4.0 on Ubuntu 22.04.5 LTS.
  
## For windows and MacOs:

Compilation on these platforms has not been tested, but it should be achievable using CMake.

# Usage

When executing the "mesh_from_nifti" binary, a prompt will request a CSV file containing a list of files to process along with configuration parameters for mesh extraction. An example file is provided as "example.csv".

The file is a 4-column CSV containing the paths to the NIfTI files to be processed, followed by the pixel coordinates of a point near the center of the region to be converted into a surface mesh. The structure is as follows:

'[path to nifti file 1];[x coordinate];[y coordinate];[z coordinate]  
[path to nifti file 2];[x coordinate];[y coordinate];[z coordinate]  
...`

Subsequently, the program will prompt for the following inputs:

- The desired edge size in millimeters (mm).
- The number of smoothing iterations.
- The percentage of edges to retain during the simplification process.

Once these parameters are provided, the mesh extraction is performed, and the resulting meshes are saved in an "STL_files" folder within the same directory as the NIfTI files.

# Example

The example is generated from data provided in https://github.com/neurolabusc/niivue-images
