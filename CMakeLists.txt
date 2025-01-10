# ****************************************************************************
# Author: Rémi Blandin
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# ****************************************************************************

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
# create_single_source_cgal_program("mesh_from_nifti.cpp")
target_link_libraries(mesh_from_nifti 
	PRIVATE ${VTK_LIBRARIES}
	Qt5::Widgets
	CGAL::CGAL
	CGAL::Eigen3_support
	)
