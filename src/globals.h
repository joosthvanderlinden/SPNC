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

#ifndef _globals_h_included_
#define _globals_h_included_

// ------------------------------------------------------------------------- GLOBAL INCLUDES
// -----------------------------------------------------------------------------------------
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/squared_distance_3.h>

// ------------------------------------------------------------------------ GLOBAL NAMESPACE
// -----------------------------------------------------------------------------------------

// Namespace for CGAL types and containers that are used throughout the program.
namespace glob
{
    // Kernel for computations. Inexact constructions implies a double is used for storage.
    // See also http://www.cgal.org/FAQ.html#predicates_vs_constructions .
    typedef CGAL::Exact_predicates_inexact_constructions_kernel EpicKernel;
    
    typedef EpicKernel::Point_3                                 Point3D;
    typedef EpicKernel::Plane_3                                 Plane3D;
    
    // Stores custom variables associated with a vertex in the triangulation
    struct VertexInfo
    {
        double      r;              // radius of sphere
        Point3D     rotated_point;  // point coordinates after rotation, see Transform class.
        std::string location;       // location of the vertex, see Triangulate::SetVertexLocations()
        int         original_id;    // original ID of the particle
    };

    // Typedefs for the Delaunay triangulation. https://doc.cgal.org/latest/Triangulation_3/
    typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, EpicKernel> VertexBase;
    typedef CGAL::Triangulation_data_structure_3<VertexBase>    TriangulationDataStructure;

    typedef TriangulationDataStructure::Vertex_handle           VertexHandle;
    typedef TriangulationDataStructure::Cell_handle             CellHandle;
    typedef TriangulationDataStructure::Facet                   Facet;
    
    typedef CGAL::Delaunay_triangulation_3<EpicKernel, TriangulationDataStructure> DelaunayTriangulation;
    typedef DelaunayTriangulation::Finite_vertices_iterator     TriangulationVertexIter;
    
    // Commonly used containers for vertices, cells and facets.
    typedef std::vector<VertexHandle>                           VertexVector;
    typedef std::vector<Facet>                                  FacetVector;
    typedef std::vector<CellHandle>                             CellVector;
    
    typedef VertexVector::const_iterator                        ConstVertexIter;
    typedef CellVector::const_iterator                          ConstCellIter;
    typedef FacetVector::const_iterator                         ConstFacetIter;
    
    // CGAL defines a facet by a CellHandle + an index (0-3). This was found to be too abstract
    // so instead, a pair is defined for the two cell handles on either side of the facet.
    typedef std::pair<CellHandle,CellHandle>                    FacetCells;
    
    typedef std::map<CellHandle,double>                         CellValueMap;
    typedef std::map<FacetCells,double>                         FacetValueMap;
    
    typedef std::vector< std::pair<Point3D,VertexInfo> >        ParticleVector;
    typedef ParticleVector::const_iterator                      ParticleVectorIter;
}

#endif