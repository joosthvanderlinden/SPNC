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

#ifndef _rangetree_h_included_
#define _rangetree_h_included_

// -------------------------------------------------------------------------------- INCLUDES
// -----------------------------------------------------------------------------------------
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>

#include <cassert>
#include <fstream>

#include "globals.h"

// --------------------------------------------------------------------------------- CLASSES
// -----------------------------------------------------------------------------------------

// A range tree segments a collection of points in a tree structure, for fast point
// retrieval. This range tree stores the vertices of the triangulation, as a pair of 
// coordinates and the handle to the vertex (for retrieval of other information, such as the
// radius of the corresponding particle). 
// For more details, see http://doc.cgal.org/latest/SearchStructures/index.html
class RangeTree
{
private:
	// Define the range tree type as a pair of coordinates and the vertex handle.
    typedef CGAL::Range_tree_map_traits_3<glob::EpicKernel,
                                          glob::VertexHandle>       TreeClass;
    typedef CGAL::Range_tree_3<TreeClass>                           RangeTreeType;
    typedef TreeClass::Key                                          TreeEntry;
    typedef TreeClass::Interval                                     TreeInterval;

    RangeTreeType range_tree_;

    // The lower corner and upper corner are the respective minimum coordinates (i.e. closest
    // to the origin) and maximum coordinates (i.e. furthest from the origin) of the 3D, 
    // rectangular box-shaped search window.
    glob::Point3D lower_corner_;
    glob::Point3D upper_corner_;
    
    // Computes corners are defined above, with an additional offset if needed.
    void ComputeLowerCorner(glob::VertexVector& window_vertices, double offset);
    void ComputeUpperCorner(glob::VertexVector& window_vertices, double offset);
    
public:
	glob::Point3D WindowLowerCorner();
    glob::Point3D WindowUpperCorner();

	// Initializes the range tree by providing iterators to the beginning and the end of the
	// container containing the triangulation vertices.
    void Initialize(glob::TriangulationVertexIter vertex_it_begin,
                    glob::TriangulationVertexIter vertex_it_end);
    
    // Performs the "window query" by returning all vertices in the range tree that fall into
    // the search window. Here, the search window is defined by the lower_corner_ and 
    // upper_corner_ (with offset, if needed) of the vertices in window_input.
    glob::VertexVector WindowQuery(glob::VertexVector& window_input, double offset);
};


#endif