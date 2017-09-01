// Copyright (c) 2015-2017   The University of Melbourne.
// All rights reserved.
//
// This file is part of SPNC: Sphere Packing Network Construction
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s)     : Joost van der Linden <joosthvanderlinden@gmail.com>
//
// Please cite the following paper if you use this code:
//
//          J.H. van der Linden, A. Sufian, G. Narsilio, A.R. Russell, A. 
//          Tordesillas, A Computational Geometry Approach to Pore Network 
//          Construction for Granular Packings (2017)

#ifndef _pack_h_included_
#define _pack_h_included_

// -------------------------------------------------------------------------------- INCLUDES
// -----------------------------------------------------------------------------------------
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

#include <cassert>
#include <fstream>
#include <chrono>

#include "globals.h"
#include "configure.h"
#include "porenet.h"
#include "contactnet.h"
#include "disjointset.h"
#include "triangulate.h"
#include "rangetree.h"

// --------------------------------------------------------------------------------- CLASSES
// -----------------------------------------------------------------------------------------

// This class is the main interface between main and all other classes. It stores the
// various representations of the assembly of particles, including the particle data, pore-
// and contact-network, Delaunay triangulation, range tree and disjoint-set. All public
// functions in this class are called from main.
class Pack
{
private:
    // Representations for the packing, see class implementations for details.
    Particles               particles_;
    PoreNet                 pore_network_;
    ContactNet              contact_network_;
    Triangulate             triangulation_;
    RangeTree               range_tree_;
    DisjointSet             disjoint_set_;
 
    // Value maps for cell and facet volumes and areas. These parameters are computed
    // by using functions and properties of multiple representations (shown above), hence
    // they are stored here.
    glob::CellValueMap      cell_void_volume_;
    glob::CellValueMap      cell_solid_volume_;
    glob::CellValueMap      cell_surface_area_;
    glob::FacetValueMap     facet_void_area_;
    glob::FacetValueMap     facet_solid_area_;
    
    // Computes the cell surface area (area of the particles intersecting the cell adjacent
    // to the void space) analytically, using the solid angle.
    void AnalyticalCellSurfaceArea(glob::CellHandle cell);
    
    // Computes the void volume in a cell analytically, also using the solid angle. This
    // approach will fail if the particles overlap too much, in which case the user should
    // use VoxelizedCellVolume().
    void AnalyticalCellVolume(glob::CellHandle cell);
    
    // Computes the void volume in a cell using voxels. The lower_corner and upper_corner
    // define the vicinity of the cell to be voxelized, while the cell_vertices and
    // nearby_vertices should contain the vertices corresponding to the particles that
    // intersect the cell.
    void VoxelizedCellVolume(glob::CellHandle cell, glob::Point3D lower_corner,
                                                    glob::Point3D upper_corner,
                                                    glob::VertexVector cell_vertices,
                                                    glob::VertexVector nearby_vertices);
    
public:
    
    // Read the input data, containing (at least) the particle coordinates and radii.
    // See README for details on the input data format.
    void ReadPack(const std::string inputfile);
    
    // Initialization functions, interfacing to the corresponding classes.
    void InitializeTriangulation();
    void InitializeRangeTree();
    void InitializeDisjointSet();
    
    // Interfaces to the Triangulate class, to determine which cells and facets are within
    // the limits of the sample.
    void                    FindInteriorCellsAndFacets();
    
    glob::ConstFacetIter    InteriorFacetsBegin();
    glob::ConstFacetIter    InteriorFacetsEnd();
    int                     NumInteriorFacets();
    
    // Interfaces to the DisjointSet class, to merge the cells on either side of the facet.
    void MergeFacetCells(glob::Facet f);
    
    // Interfaces to the Triangulate class, and returns the vertices at the three corners
    // of the facet.
    glob::VertexVector GetFacetVertices(glob::Facet f);
    
    // Uses the RangeTree class to perform a window query around the facet vertices, to
    // determine which facets are nearby.
    glob::VertexVector FindNearbyVertices(glob::VertexVector facet_vertices);
    
    // Sets the facet void area and solid area values in the corresponding maps.
    void SetFacetAreas(glob::Facet it, double void_area, double solid_area);
    
    // Performs checks to decide if the cell volume can be computed analytically, or
    // if it should be computed using voxels. Always computes the surface area analyticaly.
    void CellParameters(Configure config);
    
    // Interfaces to the network classes.
    void AddContactNetworkNodes();
    void AddContactNetworkEdges();
    void AddPoreNetworkNodes();
    void AddPoreNetworkEdges();
    
    // Writes output to file.
    void WriteTriangulation(std::string path);
    void WritePoreNetwork(std::string path);
    void WriteContactNetwork(std::string path);
    void WriteCellProperties(std::string path);
    void WriteFacetProperties(std::string path);
};


#endif