This program generates a 3D surface mesh from a segmentation mask provided in the nifti format.

# Compilation

## For Linux:

The libraries VTK, QT5 and CGAL and Eigen 3 must be installed to compile this program.

Once dependencies are installed, simply execute in the directory of the repository:  
`cmake .`  
`make` 
  
## For windows and MacOs:

No compilation was tested with these platforms, however, the compilation should be possible using cmake.

# Usage

When runing the binary "mesh_from_nifti", a prompt requires to provide a file.
This file must be generated according to the data which are processed.
It is a 4 columns `csv`file constaining the paths to the nifti files to process followed by the coordinates
in pixel unit of a point situated in the vicinity of the area to be converted to a surface mesh.
It is organized as follows:

[path to nifti file 1];[x coordinate];[y coordinate];[z coordinate]  
[path to nifti file 2];[x coordinate];[y coordinate];[z coordinate]  
...`

After that:
- the desired edge size is asked in mm,
- the number of smoothing iterations is asked,
- the percentage of edges to keep in the simplification process is asked

Finally the mesh extraction is done and the resulting meshes are saved in a folder named STL in the same directory aas the nifti files.
