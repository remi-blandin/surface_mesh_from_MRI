# About

This program creates a 3D surface mesh from a segmentation mask provided in NIfTI format.

# Compilation

## For Linux:

The libraries VTK, QT5 and CGAL and Eigen 3 must be installed to compile this program.

Once dependencies are installed, simply execute in the directory of the repository:  
`cmake .`  
`make` 
  
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
