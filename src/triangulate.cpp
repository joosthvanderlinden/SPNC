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

#include "triangulate.h"

void Triangulate::Initialize(glob::ParticleVectorIter particles_begin,
                             glob::ParticleVectorIter particles_end)
{
    // Perform triangulation with CGAL
    triangulation_.insert(particles_begin, particles_end);
    
    assert(triangulation_.is_valid());
}


bool Triangulate::SphereFacetIntersect(glob::VertexHandle vertex, glob::VertexVector facet_vertices)
{
    glob::Plane3D facet_plane = glob::Plane3D(facet_vertices[0]->point(),
                                              facet_vertices[1]->point(),
                                              facet_vertices[2]->point());
    
    double distance = std::sqrt(CGAL::to_double(CGAL::squared_distance(vertex->point(), facet_plane)));
    
    return ( distance < vertex->info().r );
}


bool Triangulate::InsideTetrahedron(glob::Point3D& p, glob::CellHandle cell)
{
    // loc, li, lj are not used
    glob::DelaunayTriangulation::Locate_type loc;
    int li, lj;
    if (triangulation_.side_of_cell(p, cell, loc, li, lj) == CGAL::ON_BOUNDED_SIDE)
    {
        return true;
    }
    return false;
}

double Triangulate::SolidAngle(glob::Point3D p[4], int i)
{
    // Set indices for vertex rotation (=/= i)
    int idxs[3];
    for(int j = 0; j < 3; ++j)
    {
        idxs[j] = (i + j + 1) % 4;
    }
    
    // The dihedral angle is the angle between two intersecting planes
    double dihedral_angle[3];
    for(int j = 0; j < 3; ++j)
    {
        // CGAL::Mesh_3::dihedral_angle returns in degrees, hence /180*pi to get radians
        dihedral_angle[j] = std::abs(CGAL::Mesh_3::dihedral_angle(p[i], p[idxs[(j+1)%3]],
                                                                        p[idxs[(j+2)%3]],
                                                                        p[idxs[(j+3)%3]]) / 180 * CGAL_PI);
    }
    
    return dihedral_angle[0] + dihedral_angle[1] + dihedral_angle[2] - CGAL_PI;
}


bool Triangulate::CellInsideLimits(glob::CellHandle c, Particles& particles)
{
    if (triangulation_.is_infinite(c))
    {
        return false;
    }
    
    for(int i = 0; i < 4; ++i)
    {
        glob::Facet f = std::make_pair(c, i);
        if (FacetInsideLimits(f, particles))
        {
            return true;
        }
    }
    
    // If none of the vertices are inside the limits, return false
    return false;
}


bool Triangulate::FacetInsideLimits(glob::Facet f, Particles& particles)
{
    glob::VertexVector facet_vertices = GetFacetVertices(f);
    
    // If one of the vertices on the facet is outside the limits,
    // the whole facet is deemed outside the limits.
    for (glob::ConstVertexIter vertex_it = facet_vertices.begin();
         vertex_it != facet_vertices.end(); vertex_it++)
    {
        if (particles.OutsideLimits((*vertex_it)->point(), 0.0))
        {
            return false;
        }
    }
    
    // Also check if either of the adjacent cells are infinite
    glob::CellHandle c1 = f.first;
    glob::CellHandle c2 = f.first->neighbor(f.second);
    if (triangulation_.is_infinite(c1) || triangulation_.is_infinite(c2))
    {
        return false;
    }
    return true;
}


void Triangulate::SetVertexLocations(Particles& particles)
{
    // This is an inelegant implementation but it does the job and is
    // (currently) not a bottleneck, so leaving as is.
    
    for(TriangulationCellIter cell_it = triangulation_.finite_cells_begin();
        cell_it != triangulation_.finite_cells_end(); ++cell_it)
    {
        bool cell_in_vertical_limits = true;
        bool cell_in_horizontal_limits = true;
        std::vector<std::string> locations(4);
        for(int i = 0; i < 4; ++i)
        {
            glob::VertexHandle v = cell_it->vertex(i);
            
            // Check if cell is inside the limits of the subsample. For this check,
            // distinguish between vertical limits (y) and horizontal limits (x and z).
            // For the cell to be inside the limits, every vertex needs to be within the
            // limits. If the cell is (1) inside the horizontal limits, and (2) outside the
            // vertical limits, i.e. the cell is at the top (bottom), then we want
            // all vertices for this cell to be marked as top (bottom). By doing so, we get
            // a truthful division of top/middle/bottom vertices for all the cells in
            // interior_cells_ (which is going to be used later to decide if pores are top/middle/bottom).
            std::string v_location = particles.VerticalLocation(v->point().y(), 0.0);
            std::string h_location = particles.HorizontalLocation(v->point().x(), v->point().z(), 0.0);
            
            if ((v_location == "top") || (v_location == "bottom"))
            {
                cell_in_vertical_limits = false;
            }
            
            if (h_location == "outside")
            {
                cell_in_horizontal_limits = false;
            }
            
            // Now get the vertical location of the sphere (with radius r) at this vertex.
            // If y + r > y_max or y - r < y_min, then the sphere is top or bottom, respectively.
            v_location = particles.VerticalLocation(v->point().y(), v->info().r);
            
            // If the vertex location has already been assigned "top" or "bottom", then
            // ignore the location from particles.VerticalLocation() and use existing
            // location instead.
            if ((v->info().location == "top") || (v->info().location == "bottom"))
            {
                locations[i] = v->info().location;
            }
            else
            {
                locations[i] = v_location;
            }
        }
        
        // If the cell is outside the vertical limits, inside the horizontal limits,
        // and one of the vertices is "top", make all vertices "top".
        if ((!cell_in_vertical_limits) && cell_in_horizontal_limits &&
            (std::find(locations.begin(), locations.end(), "top") != locations.end()))
        {
            std::fill(locations.begin(), locations.end(), "top");
        }
        
        // Same for "bottom". Assumes the sample is large enough such that no cells
        // are both top and bottom.
        if ((!cell_in_vertical_limits) && cell_in_horizontal_limits &&
            (std::find(locations.begin(), locations.end(), "bottom") != locations.end()))
        {
            std::fill(locations.begin(), locations.end(), "bottom");
        }
        
        // Finally, assign the locations
        for(int i = 0; i < 4; ++i)
        {
            glob::VertexHandle v = cell_it->vertex(i);
            v->info().location = locations[i];
        }
    }
}


void Triangulate::FindInteriorCells(Particles& particles)
{
    for(TriangulationCellIter cell_it = triangulation_.finite_cells_begin();
        cell_it != triangulation_.finite_cells_end(); ++cell_it)
    {
        if (CellInsideLimits(cell_it, particles))
        {
            interior_cells_.push_back(cell_it);
        }
    }
}


void Triangulate::FindInteriorFacets(Particles& particles)
{
    for(TriangulationFacetIter facet_it = triangulation_.finite_facets_begin();
        facet_it != triangulation_.finite_facets_end(); ++facet_it)
    {
        if (FacetInsideLimits(*facet_it, particles))
        {
            interior_facets_.push_back(*facet_it);
        }
    }
}


glob::VertexVector Triangulate::GetFacetVertices(glob::Facet f)
{
    glob::VertexVector facet_vertices(3);
    
    // Facet is defined as a pair(cell handle, idx), where
    // idx is the index of the vertex opposite to the facet
    for(int i = 1; i < 4; ++i)
    {
        facet_vertices[i-1] = f.first->vertex((f.second+i)%4);
    }
    
    return facet_vertices;
}

glob::TriangulationVertexIter Triangulate::FiniteVerticesBegin()
{
    return triangulation_.finite_vertices_begin();
}

glob::TriangulationVertexIter Triangulate::FiniteVerticesEnd()
{
    return triangulation_.finite_vertices_end();
}

glob::ConstCellIter Triangulate::InteriorCellsBegin()
{
    return interior_cells_.begin();
}

glob::ConstCellIter Triangulate::InteriorCellsEnd()
{
    return interior_cells_.end();
}

glob::ConstFacetIter Triangulate::InteriorFacetsBegin()
{
    return interior_facets_.begin();
}

glob::ConstFacetIter Triangulate::InteriorFacetsEnd()
{
    return interior_facets_.end();
}

int Triangulate::NumFiniteVertices()
{
    return triangulation_.number_of_vertices();
}

int Triangulate::NumInteriorFacets()
{
    return interior_facets_.size();
}

int Triangulate::NumInteriorCells()
{
    return interior_cells_.size();
}

