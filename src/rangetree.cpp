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

#include "rangetree.h"

void RangeTree::Initialize(glob::TriangulationVertexIter vertex_it_begin,
                           glob::TriangulationVertexIter vertex_it_end)
{
    std::vector<TreeEntry> tree_entries(std::distance(vertex_it_begin,
                                                      vertex_it_end));
    
    glob::TriangulationVertexIter vertex_it = vertex_it_begin;
    int i = 0;
    for ( ; vertex_it != vertex_it_end; ++vertex_it, ++i)
    {
        tree_entries[i] = TreeEntry(vertex_it->point(), vertex_it);
    }
    
    // Construct the range with CGAL
    range_tree_.make_tree(tree_entries.begin(), tree_entries.end());
}


void RangeTree::ComputeLowerCorner(glob::VertexVector& window_vertices, double offset)
{
    // Set initial minimum
    glob::ConstVertexIter vertex_it = window_vertices.begin();
    glob::Point3D p = (*vertex_it)->point();
    double min_x    = p.x();
    double min_y    = p.y();
    double min_z    = p.z();
    ++vertex_it;
    
    // Check remaining vertices
    for (; vertex_it != window_vertices.end(); ++vertex_it)
    {
        p = (*vertex_it)->point();
        
        if( p.x() < min_x )
            min_x = p.x();
        
        if( p.y() < min_y )
            min_y = p.y();
        
        if( p.z() < min_z )
            min_z = p.z();
    }
    
    lower_corner_ = glob::Point3D(min_x - offset,
                                  min_y - offset,
                                  min_z - offset);
}


void RangeTree::ComputeUpperCorner(glob::VertexVector& window_vertices, double offset)
{
    // Set initial maximum
    glob::ConstVertexIter vertex_it = window_vertices.begin();
    glob::Point3D p = (*vertex_it)->point();
    double max_x    = p.x();
    double max_y    = p.y();
    double max_z    = p.z();
    ++vertex_it;
    
    // Check remaining vertices
    for(;vertex_it != window_vertices.end(); ++vertex_it)
    {
        p = (*vertex_it)->point();
        
        if( p.x() > max_x )
            max_x = p.x();
        
        if( p.y() > max_y )
            max_y = p.y();
        
        if( p.z() > max_z )
            max_z = p.z();
    }
    
    upper_corner_ = glob::Point3D (max_x + offset,
                                   max_y + offset,
                                   max_z + offset);
}


glob::VertexVector RangeTree::WindowQuery(glob::VertexVector& window_input, double offset)
{
    // Minimum coordinates - offset
    ComputeLowerCorner(window_input, offset);
    
    // Maximum coordinates + offset
    ComputeUpperCorner(window_input, offset);
    
    // According to this http://doc.cgal.org/latest/SearchStructures/index.html
    // TreeInterval should allow for 2 Point3D's to be passed, but it only seems
    // to take 2 TreeEntry's. TreeEntry takes a vertex handle as the 2nd argument
    // but this is not used.
    glob::ConstVertexIter vertex_it = window_input.begin();
    TreeInterval window(TreeInterval(TreeEntry(lower_corner_,(*vertex_it)),
                                     TreeEntry(upper_corner_,(*vertex_it))));
    
    // window_query() takes a search window and a back_inserter()
    std::vector<TreeEntry> window_output;
    range_tree_.window_query(window, std::back_inserter(window_output));
    
    // Some vertices in window_input will also be in window_output. These vertices
    // should not be included in nearby_vertices, so including a check. Note that
    // vertex_it is already set to window_input.begin() above.
    glob::VertexVector nearby_vertices;
    for(std::vector<TreeEntry>::iterator entry_it = window_output.begin();
        entry_it != window_output.end(); ++entry_it)
    {
        if (std::find(window_input.begin(), window_input.end(),
                      entry_it->second) == window_input.end())
        {
            nearby_vertices.push_back(entry_it->second);
        }
    }
    
    return nearby_vertices;
}


glob::Point3D RangeTree::WindowLowerCorner()
{
    return lower_corner_;
}


glob::Point3D RangeTree::WindowUpperCorner()
{
    return upper_corner_;
}
