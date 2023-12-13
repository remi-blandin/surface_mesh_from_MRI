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

namespace SMS = CGAL::Surface_mesh_simplification;

int main(int argc, char *argv[]) {

  //*********************************************
  // Load Nifti file
  //*********************************************
  
  // create the file dialog
  QApplication app(argc, argv);
  QString fileName = QFileDialog::getOpenFileName(nullptr, "Open File", "./images", "All Files (*.*)");
  std::cout << "Processing file:\n" << fileName.toStdString() << std::endl;
  
  // load NifTi image with VTK
  vtkNew<vtkNIFTIImageReader> reader;
  reader->SetFileName(fileName.toStdString().c_str()); 
//  reader->SetFileName("images/tense-a-envelop2.nii.gz"); // tense-a-envelop
  reader->Update();
  auto vtk_image = reader->GetOutput();
  Gray_level_image image(CGAL::IO::read_vtk_image_data(vtk_image), 2.9f);
  
  // print the dimensions of the image
  std::cout << "xdim: " << image.xdim() 
  << " ydim: " << image.ydim()
  << " zdim: " << image.zdim() << std::endl;
  
    // print the dimensions of the image
  std::cout << "vx: " << image.vx() 
  << " vy: " << image.vy()
  << " vz: " << image.vz() << std::endl;
  
      // print the dimensions of the image
  std::cout << "tx: " << image.tx() 
  << " ty: " << image.ty()
  << " tz: " << image.tz() << std::endl;
  
//  // print the values of the image
//  for (int i(0); i < image.xdim(); i++)
//  {
//    for (int j(0); j < image.ydim(); j++)
//    {
//      for (int k(0); k < image.zdim(); k++)
//      {
//        std::cout << image.value(i, j, k) << std::endl;
//      }
//    }
//  }

  //*********************************************
  // Request meshing parameters
  //*********************************************

    // Prompt the user for a number using QInputDialog and get the result
    bool ok;
    double max_dim = QInputDialog::getDouble(nullptr, "Maximal dimension", "Enter mesh maximal dimension (mm)", 5., 0., 100., 1, &ok);
    
    double center_x = QInputDialog::getDouble(nullptr, "Center X", "Enter x coordinate of the center (mm)", double(image.xdim()) / 2., -1000., 1000., 1, &ok);
    
    double center_y = QInputDialog::getDouble(nullptr, "Center Y", "Enter y coordinate of the center (mm)", double(image.ydim()) / 2., -1000., 1000., 1, &ok);
    
    double center_z = QInputDialog::getDouble(nullptr, "Center Z", "Enter z coordinate of the center (mm)", double(image.zdim()) / 2., -1000., 1000., 1, &ok);
    
    double stop_ratio = QInputDialog::getDouble(nullptr, "Percentage simplification", "Enter the percentage of edges to remove", 10., 0., 100., 1, &ok);
    stop_ratio /= 100.;


  //*********************************************
  // Create surface mesh
  //*********************************************
  
  if (ok){
  
    std::cout << "Mesh will be generated with maximal dimension " << max_dim 
      << " mm" << std::endl; 

    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

    // Carefully chosen bounding sphere: the center must be inside the
    // surface defined by 'image' and the radius must be high enough so that
    // the sphere actually bounds the whole image.
    GT::Point_3 bounding_sphere_center(
      center_x,
      center_y,
      center_z);
//    GT::Point_3 bounding_sphere_center(double(image.xdim()) / 2.,
//      double(image.ydim()) / 2.,
//      double(image.zdim()) / 2.);
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
    std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";
    
  //*********************************************
  // Simplify mesh
  //*********************************************
  
    SMS::Count_ratio_stop_predicate<Surface_mesh> stop(stop_ratio);
  
    Surface_mesh surface_mesh;
    const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("out.off");
    std::ifstream is(filename);
    if(!is || !(is >> surface_mesh))
    {
      std::cerr << "Failed to read input mesh: " << filename << std::endl;
      return EXIT_FAILURE;
    }
    if(!CGAL::is_triangle_mesh(surface_mesh))
    {
      std::cerr << "Input geometry is not triangulated." << std::endl;
      return EXIT_FAILURE;
    }

    
    int r = SMS::edge_collapse(surface_mesh, stop);
  
    std::cout << "\nFinished!\n" << r << " edges removed.\n" << surface_mesh.number_of_edges() << " final edges.\n";
    
    CGAL::IO::write_polygon_mesh((argc > 3) ? argv[3] : "out_simp.off", surface_mesh, CGAL::parameters::stream_precision(17));
  
  }
  else
  {
    std::cout << "Abort" << std::endl;
  }
}
