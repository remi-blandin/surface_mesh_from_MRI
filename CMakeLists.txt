# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.8...3.12 FATAL_ERROR)
project(mesh_from_nifti)

include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

find_package(VTK
  COMPONENTS
    CommonCore
    CommonDataModel
    CommonTransforms
    FiltersCore
    FiltersGeneral
    FiltersGeometry
    FiltersModeling
    IOImage
    IOXML
    ImagingCore
    ImagingStatistics
    InteractionStyle
    RenderingCore
    RenderingVolume
    RenderingOpenGL2
    RenderingVolumeOpenGL2
  OPTIONAL_COMPONENTS
    TestingCore
    TestingRendering)
if (NOT VTK_FOUND)
  message("Skipping example: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif ()

# Find the Qt5 package
find_package(Qt5 COMPONENTS Widgets REQUIRED)

find_package(CGAL REQUIRED COMPONENTS ImageIO)

# include for local package
find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)

add_executable(mesh_from_nifti MACOSX_BUNDLE mesh_from_nifti.cpp)
add_executable(test_mesh_3D_gray_image MACOSX_BUNDLE test_mesh_3D_gray_image.cpp)
# create_single_source_cgal_program("mesh_from_nifti.cpp")
target_link_libraries(mesh_from_nifti 
	PRIVATE ${VTK_LIBRARIES}
	Qt5::Widgets
	CGAL::CGAL
	CGAL::Eigen3_support
	)
target_link_libraries(test_mesh_3D_gray_image 
	PRIVATE ${VTK_LIBRARIES}
	Qt5::Widgets
	CGAL::CGAL
	CGAL::Eigen3_support
	)
