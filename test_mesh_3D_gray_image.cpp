#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/IO/read_vtk_image_data.h>

// for weighted image
#include <CGAL/Mesh_3/generate_label_weights.h>
#include <CGAL/IO/File_binary_mesh_3.h>
#include <CGAL/tags.h>



// VTK includes
#include <vtkNew.h>
#include <vtkNIFTIImageReader.h>

// STL includes
#include <filesystem>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char* argv[])
{
  /// [Loads image]
	const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("tense-i.nii.gz");
  vtkNew<vtkNIFTIImageReader> reader;
      reader->SetFileName(fname.c_str()); 
      reader->Update();
      auto vtk_image = reader->GetOutput();
      CGAL::Image_3 image(CGAL::IO::read_vtk_image_data(vtk_image));
      std::cout << "\nImage successfully imported" << endl;
  /// [Loads image]

  /// [Domain creation]
  Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image);
  /// [Domain creation]

  // Mesh criteria
  const double facet_angle_val = (argc > 2) ? std::atof(argv[2]) : 30;
  const double facet_size_val = (argc > 3) ? std::atof(argv[3]) : 3;
  const double facet_distance_val = (argc > 4) ? std::atof(argv[4]) : 4;
  const double cell_radius_edge_ratio_val = (argc > 5) ? std::atof(argv[5]) : 3;
  const double cell_size_val = (argc > 6) ? std::atof(argv[6]) : 8;
  
  cout << "facet_angle: " << facet_angle_val << endl;
  cout << "facet_size: " << facet_size_val << endl;
  cout << "facet_distance: " << facet_distance_val << endl;
  cout << "cell_radius_edge_ratio: " << cell_radius_edge_ratio_val << endl;
  cout << "cell_size: " << cell_size_val << endl;
  Mesh_criteria criteria(facet_angle=facet_angle_val, facet_size=facet_size_val, 
  facet_distance=facet_distance_val, cell_radius_edge_ratio=cell_radius_edge_ratio_val, 
  cell_size=cell_size_val);  
//  Mesh_criteria criteria(facet_angle, facet_size, 
//  facet_distance, cell_radius_edge_ratio, cell_size);

  /// [Meshing]
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
  /// [Meshing]

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);
  medit_file.close();
  
  std::ofstream surface_mesh_file("out.off");
  std::string dir_output_surfaces("output_surfaces");
  
  // check if the output surface directory exists, and create it if necessary
  if (!std::filesystem::is_directory(dir_output_surfaces)) 
  {
    std::filesystem::create_directory(dir_output_surfaces);
  }
  
  std::string file_name;
  
  std::vector<std::string> name_list{
  "Fat",
  "Fat",
  "Skull",
  "Eyes",
  "Upper_teeth",
  "Air",
  "Velum",
  "Mucosa",
  "Cerebrospinal_fluid",
  "Grey_matter",
  "White_matter",
  "Mandible",
  "Lower_teeth",
  "Tongue",
  "Skin",
  "Additional_tongue"
  };
  
  // extract all the surfaces
  for (int i(1); i < 16; i++)
  {
		file_name = dir_output_surfaces + "/" + name_list[i] + "_" + std::to_string(i) + ".off";
		cout << file_name << endl;
		surface_mesh_file.open(file_name);
	  c3t3.output_boundary_to_off(surface_mesh_file, i);
	  surface_mesh_file.close();
  }
  
  //*****************************************
  // Make image with weights
  //*****************************************
  
  const float sigma = 10.f; //(std::max)(image.vx(), (std::max)(image.vy(), image.vz()));
  CGAL::Image_3 img_weights = 
  	CGAL::Mesh_3::generate_label_weights(image, sigma);

	domain = Mesh_domain::create_labeled_image_mesh_domain(image,
                                                    weights = img_weights,
                                                    relative_error_bound = 1e-6);
	
	c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
	
	medit_file.open("out_w.mesh");

  return 0;
}
