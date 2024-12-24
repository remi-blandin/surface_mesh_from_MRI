#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/read_vtk_image_data.h>
#include <fstream>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>

// CGAL for mesh simplification
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

// CGAL for mesh smoothing
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

// VTK includes
#include <vtkNew.h>
#include <vtkNIFTIImageReader.h>

// QT includes
#include <QApplication>
#include <QFileDialog>
#include <QDialog>
#include <QInputDialog>
#include <QLabel>
#include <QPushButton>
#include <QVBoxLayout>

// STL includes
#include <iostream>
#include <fstream>
#include <filesystem>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;
typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;

typedef CGAL::Simple_cartesian<double>               Kernel;
typedef Kernel::Point_3                              Point_3;
typedef CGAL::Surface_mesh<Point_3>                  Surface_mesh;

// Typedef for mesh smoothing
typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char *argv[]) {

  QApplication app(argc, argv);

  std::string directory_stl_files("/STL_files/");
  std::filesystem::path stl_file_path;
  std::string extension;

  //*********************************************
  // Load the mesh extraction parameter file
  //*********************************************
  
  // create the file dialog
  QString fileName = QFileDialog::getOpenFileName(nullptr, "Open File", "All Files (*.*)");
  
  std::ifstream inputFile(fileName.toStdString());
  std::string line, item;
  char separator(';');
  int cnt;
  double coord[3];
  std::vector<std::filesystem::path> filePathList;
  std::vector<Point_3> sphereCenterList;
  
  if (!inputFile.is_open())
  {
    std::cout << "Cannot open " << fileName.toStdString() << endl;
    return(false);
  }
  else if (inputFile.peek() != std::ifstream::traits_type::eof())
  {
    while (getline(inputFile, line))
    {
      std::stringstream coordLine(line);
      cnt = 0;
      for (int i(0); i < 3; i++){coord[i] = 0.;}
      
      while (getline(coordLine, item, separator) && (cnt < 4))
      {
        if (item.size() > 0)
        {
          if (cnt == 0)
          {
            filePathList.push_back(item);
            cnt ++;
          }
          else
          {
            try
            {
              coord[cnt-1] = std::stod(item);
              cnt++;
            }
            catch (const std::invalid_argument& ia) {
              std::cout << ia.what() << endl;
            }
            

          }
        }
      }
      sphereCenterList.push_back(Point_3(coord[0], coord[1], coord[2]));
    }
  }
  
  //*********************************************
  // Request meshing parameters
  //*********************************************

  // Prompt the user for a number using QInputDialog and get the result
  bool ok;
  double max_dim = QInputDialog::getDouble(nullptr, "Maximal dimension", "Enter mesh maximal dimension (mm)", 1., 0., 100., 1, &ok);

	unsigned int nb_iterations = QInputDialog::getInt(nullptr, "Number of smoothing iterations", "Enter the number of smoothing iterations", 6, 0., 50, 1, &ok);
    
  double stop_ratio = QInputDialog::getDouble(nullptr, "Percentage simplification", "Enter the percentage of edges to remove", 10., 0., 100., 1, &ok);
  stop_ratio /= 100.;
  
  if (ok){
  
    for (int img(0); img < filePathList.size(); img ++)
    {
      
      //*********************************************
      // Load Nifti file
      //*********************************************
      
      std::cout << "\n**********************************************************" << endl;
      std::cout << "Processing file " << img + 1 << "/" << filePathList.size() << " :\n" 
      << filePathList[img] << std::endl;
      
      // load NifTi image with VTK
      vtkNew<vtkNIFTIImageReader> reader;
      reader->SetFileName(filePathList[img].c_str()); 
      reader->Update();
      auto vtk_image = reader->GetOutput();
      Gray_level_image image(CGAL::IO::read_vtk_image_data(vtk_image), 2.9f);
      std::cout << "\nImage successfully imported" << endl;
      
      // print the dimensions of the image
      std::cout << "xdim: " << image.xdim() 
      << " ydim: " << image.ydim()
      << " zdim: " << image.zdim() << std::endl;

      //*********************************************
      // Create surface mesh
      //*********************************************
    
      std::cout << "\nMesh will be generated with maximal dimension " << max_dim 
        << " mm" << std::endl; 

      Tr tr;            // 3D-Delaunay triangulation
      C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

      // Carefully chosen bounding sphere: the center must be inside the
      // surface defined by 'image' and the radius must be high enough so that
      // the sphere actually bounds the whole image.
      GT::Point_3 bounding_sphere_center(
        sphereCenterList[img].x(),
        sphereCenterList[img].y(),
        sphereCenterList[img].z());

      double radius_mm = 260.;
      
      GT::FT bounding_sphere_squared_radius = radius_mm * radius_mm;
      GT::Sphere_3 bounding_sphere(bounding_sphere_center,
                                       bounding_sphere_squared_radius);
                                   
      // definition of the surface, with 10^-5 as relative precision
      Surface_3 surface(image, bounding_sphere, 1e-5);

      // defining meshing criteria
      CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,
                                                         max_dim,
                                                         max_dim);

      // meshing surface, with the "manifold without boundary" algorithm
      CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
      
    //*********************************************
    // Export generated mesh
    //*********************************************
    
      std::ofstream out("out.off");
      CGAL::output_surface_facets_to_off (out, c2t3);
      std::cout << "Final number of points: " << tr.number_of_vertices() << endl;
      
    //*********************************************
    // Smooth mesh
    //*********************************************
    
      const std::string filename_smooth = (argc > 1) ? argv[1] : CGAL::data_file_path("out.off");

      Mesh mesh;
      if(!PMP::IO::read_polygon_mesh(filename_smooth, mesh))
      {
        std::cerr << "Invalid input." << std::endl;
        return 1;
      }
      
      const double time = (argc > 3) ? std::atof(argv[3]) : 0.05;
    
      std::set<Mesh::Vertex_index> constrained_vertices;
      for(Mesh::Vertex_index v : vertices(mesh))
      {
        if(is_border(v, mesh))
          constrained_vertices.insert(v);
      }
      
      std::cout << "\nMesh will be smoothed with " << nb_iterations << " iterations of time step "
      << time << endl;
      std::cout << "Mesh contains: " << constrained_vertices.size() << " border vertices" << std::endl;
      CGAL::Boolean_property_map<std::set<Mesh::Vertex_index> > vcmap(constrained_vertices);
      
      PMP::smooth_shape(mesh, time, PMP::parameters::number_of_iterations(nb_iterations)
                                                    .vertex_is_constrained_map(vcmap));
      std::cout << "Smoothing successful" << std::endl;
      
      CGAL::IO::write_polygon_mesh((argc > 3) ? argv[3] : "out_smooth.off", mesh, CGAL::parameters::stream_precision(17));
      
    //*********************************************
    // Simplify mesh
    //*********************************************
    
      SMS::Count_ratio_stop_predicate<Surface_mesh> stop(stop_ratio);
      
      std::cout << "\nMesh will be simplified keeping " << 100. * stop_ratio << "% of the edges" << endl;
      
      int r = SMS::edge_collapse(mesh, stop);
    
      std::cout << "Mesh simplification successful:" << r << " edges removed, " << mesh.number_of_edges() << " final edges." << endl;
      
      CGAL::IO::write_polygon_mesh((argc > 3) ? argv[3] : "out_simp.off", mesh, CGAL::parameters::stream_precision(17));

      // check if the STL files directory exists and create it if not
      stl_file_path = filePathList[img].parent_path().string() + directory_stl_files;
      if (!std::filesystem::is_directory(stl_file_path)) 
      {
        std::filesystem::create_directory(stl_file_path);
      }
      
      // check if the file have already been processed before
      cnt = 0;
      for (int f(0); f < img; f++)
      {
        if (filePathList[f] == filePathList[img]) {cnt ++;}
      }
      
      // Add the filename with the correct extension
      stl_file_path = stl_file_path.string() + filePathList[img].filename().string();
      stl_file_path.replace_extension();
      stl_file_path.replace_extension();
      if (cnt > 0)
      {
        extension = "_" + std::to_string(cnt) + ".stl";
      }
      else
      {
        extension = ".stl";
      }
      stl_file_path = stl_file_path.string() + extension;
      
      // Export mesh in STL format
      CGAL::IO::write_STL(stl_file_path.string(), mesh, CGAL::parameters::stream_precision(17));
    }
  
  }
  else
  {
    std::cout << "Abort" << std::endl;
  }
}
