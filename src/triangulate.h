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
// Author(s)     : Joost van der Linden <joostv@student.unimelb.edu.au>
//
// Please cite the following paper if you use this code:
//
//          J.H. van der Linden, A. Sufian, G. Narsilio, A.R. Russell, A. Tordesillas,
//          (2016), Delaunay-based pore network construction for granular packings,
//          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#ifndef _triangulate_h_included_
#define _triangulate_h_included_

// -------------------------------------------------------------------------------- INCLUDES
// -----------------------------------------------------------------------------------------
#include <CGAL/Mesh_3/dihedral_angle_3.h>

#include "globals.h"
#include "particles.h"

// --------------------------------------------------------------------------------- CLASSES
// -----------------------------------------------------------------------------------------

// Stores the Delaunay triangulation, and related functions. For more details see:
// https://doc.cgal.org/latest/Triangulation_3/.
class Triangulate
{
private:
	// A tetrahedra on the hull of the triangulation is adjacent to an infinite cell in 
	// CGAL. These infinite cells and corresponding infinite facets should not be considered, 
	// hence the finite cell/facet iterators are used.
    typedef glob::DelaunayTriangulation::Finite_cells_iterator    TriangulationCellIter;
    typedef glob::DelaunayTriangulation::Finite_facets_iterator   TriangulationFacetIter;
    
    glob::DelaunayTriangulation triangulation_;

    // The Delaunay triangulation suffers from highly elongated tetrahedra near the hull of
    // the triangulation. To deal with this, "interior" cells and facets are identified.
    glob::CellVector 	interior_cells_;
    glob::FacetVector   interior_facets_;
    
public:

    glob::TriangulationVertexIter   FiniteVerticesBegin();
    glob::TriangulationVertexIter   FiniteVerticesEnd();
    glob::ConstCellIter             InteriorCellsBegin();
    glob::ConstCellIter             InteriorCellsEnd();
    glob::ConstFacetIter            InteriorFacetsBegin();
    glob::ConstFacetIter            InteriorFacetsEnd();
    
    int                     		NumFiniteVertices();
    int                     		NumInteriorCells();
    int                     		NumInteriorFacets();

	// Performs the triangulation in CGAL, using the particles from the Particles class.
    void Initialize(glob::ParticleVectorIter particles_begin,
                    glob::ParticleVectorIter particles_end);
    
    // Computes the solid angle at point p[i], where p contains the 4 vertices of
    // a tetrahedron. See also: https://en.wikipedia.org/wiki/Solid_angle#Tetrahedron.
    double SolidAngle(glob::Point3D p[4], int i);

    // Tests if the point p is inside the tetrahedron corresponding to the cell.
    bool InsideTetrahedron(glob::Point3D& p, glob::CellHandle cell);

    // Tests if the spherical particle at the nearby_vertex is (perpendicular) less than 
    // it's radius away from the plane that the facet resides in.
    bool SphereFacetIntersect(glob::VertexHandle nearby_vertex, 
    						  glob::VertexVector facet_vertices);
    
    // Determine for each vertex if the corresponding point is located at the top of the
    // packing, in the middle, or at the bottom.
    void SetVertexLocations(Particles& particles);
    
    // The Delaunay triangulation will form highly elongated tetrahedra on the hull of the
    // triangulation, which are undesirable for the construction of the pore network.
    // To this end, 'interior' cells and facets are identified, away from the hull.
    void FindInteriorCells(Particles& particles);
    void FindInteriorFacets(Particles& particles);
    bool CellInsideLimits(glob::CellHandle c, Particles& particles);
    bool FacetInsideLimits(glob::Facet f, Particles& particles);

    glob::VertexVector GetFacetVertices(glob::Facet it);
};


#endif