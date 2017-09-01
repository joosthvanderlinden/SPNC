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

#include "porenet.h"

int PoreNet::NumNodes()
{
    return nodes_.size();
}


int PoreNet::NumEdges()
{
    return edges_.size();
}


PoreNet::NodeIter PoreNet::FindNode(glob::CellHandle root_cell)
{
    NodeIter node_end = nodes_.end();
    
    for (NodeIter node_it = nodes_.begin(); node_it != node_end; ++node_it)
    {
        // Nodes in the pore network are identified by the cell handle
        // of the "root cell" (root node in the disjoint-set).
        if (node_it->root_cell == root_cell)
        {
            return node_it;
        }
    }
    
    return node_end;
}


void PoreNet::AddNode(glob::CellHandle cell, glob::CellHandle root_cell)
{
    NodeIter node_it = FindNode(root_cell);
    NodeIter missing = nodes_.end();
    
    if (node_it == missing)
    {
        Node node;
        node.root_cell = root_cell;
        node.cells.push_back(cell);
        nodes_.push_back(node);
        
        cell_to_node_idx_[cell] = nodes_.size() - 1;
    }
    else
    {
        node_it->cells.push_back(cell);
        
        cell_to_node_idx_[cell] = node_it - nodes_.begin();
    }
}


PoreNet::EdgeIter PoreNet::FindEdge(std::pair<int,int> edge_node_indices)
{
    // FindEdge expects the pair<a,b> of node indices to be ordered, i.e. a < b
    // (can use std::minmax). This prevents duplication of edges, i.e. a->b and b->a.
    std::map<std::pair<int, int>, int>::iterator find_it = edge_to_edge_idx.find(edge_node_indices);
    
    if (find_it == edge_to_edge_idx.end())
    {
        return edges_.end();
    }
    else
    {
        // In this case, the edge already exists in edges_ and we return an iterator to
        // that edge (find_it->second is the index of our edge in the edges_ vector) by
        // doing some arithmetic with STL iterators.
        return edges_.begin() + find_it->second;
    }
}

bool PoreNet::CellEdgeEquals(glob::FacetCells facet_cells,
                             glob::CellHandle cell,
                             glob::CellHandle neighbor_cell)
{
    return (facet_cells.first  == cell &&
            facet_cells.second == neighbor_cell) ||
           (facet_cells.first  == neighbor_cell &&
            facet_cells.second == cell);
}


PoreNet::CellEdgeIter PoreNet::FindCellEdge(EdgeIter edge_it,
                                            glob::CellHandle cell,
                                            glob::CellHandle neighbor_cell)
{
    CellEdgeIter cell_edge_end = edge_it->cell_edges.end();
    
    for (CellEdgeIter cell_edge_it = edge_it->cell_edges.begin(); cell_edge_it != cell_edge_end; ++cell_edge_it)
    {
        if (CellEdgeEquals(*cell_edge_it, cell, neighbor_cell))
        {
            return cell_edge_it;
        }
    }
    return cell_edge_end;
}


void PoreNet::AddEdge(glob::CellHandle cell, glob::CellHandle neighbor_cell)
{
    // nodes_ vector has been filled, so every cell will have a nodes_ vector index
    Node* node           = &nodes_[cell_to_node_idx_[cell]];
    Node* neighbor_node  = &nodes_[cell_to_node_idx_[neighbor_cell]];
    
    std::pair<int,int> edge_node_indices = std::minmax(cell_to_node_idx_[cell],
                                                       cell_to_node_idx_[neighbor_cell]);
    
    EdgeIter edge_it = FindEdge(edge_node_indices);
 
    if (edge_it == edges_.end())
    {
        Edge e;
        e.void_area  = 0.0;
        e.solid_area = 0.0;
        e.edge_nodes = std::make_pair(node, neighbor_node);
        e.cell_edges.push_back(std::make_pair(cell, neighbor_cell));
        edges_.push_back(e);
        
        edge_to_edge_idx[edge_node_indices] = edges_.size() - 1;
    }
    else
    {
        // Assuming the network edge already exists, there are two remaining possibilities:
        // #1. No entry in cell_edges exists yet for this pair of cell and neighbor_cell
        // #2. The entry does exist (possibly in reversed order)
        // In case #1, add the edge to cell_edges of the existing network edge
        // In case #2, do nothing
        CellEdgeIter cell_edge_it = FindCellEdge(edge_it, cell, neighbor_cell);
        if (cell_edge_it == edge_it->cell_edges.end())
        {
            edge_it->cell_edges.push_back(std::make_pair(cell, neighbor_cell));
        }
    }
}


PoreNet::ConstNodeIter PoreNet::NodesBegin()
{
    return nodes_.begin();
}


PoreNet::ConstNodeIter PoreNet::NodesEnd()
{
    return nodes_.end();
}


PoreNet::ConstEdgeIter PoreNet::EdgesBegin()
{
    return edges_.begin();
}


PoreNet::ConstEdgeIter PoreNet::EdgesEnd()
{
    return edges_.end();
}


void PoreNet::AggregateNodeCells(glob::CellValueMap& cell_void_volume,
                                 glob::CellValueMap& cell_solid_volume,
                                 glob::CellValueMap& cell_surface_area)
{
    for (NodeIter node_it = nodes_.begin(); node_it != nodes_.end(); ++node_it)
    {
        glob::VertexVector pore_vertices;
        glob::FacetVector  pore_facets;
        node_it->void_volume  = 0.0;
        node_it->solid_volume = 0.0;
        node_it->surface      = 0.0;
        
        for(glob::ConstCellIter cell_it = node_it->cells.begin();
            cell_it != node_it->cells.end(); ++cell_it)
        {
            node_it->void_volume  += cell_void_volume[*cell_it];
            node_it->solid_volume += cell_solid_volume[*cell_it];
            node_it->surface      += cell_surface_area[*cell_it];
            
            // Add (unique) vertices and facets of every cell to pore_vertices
            for(int i = 0; i < 4; ++i)
            {
                glob::VertexHandle v = (*cell_it)->vertex(i);
                glob::Facet        f = std::make_pair((*cell_it), i);
                
                if (std::find(pore_vertices.begin(),
                              pore_vertices.end(), v) == pore_vertices.end())
                {
                    pore_vertices.push_back(v);
                }
                
                if (std::find(pore_facets.begin(),
                              pore_facets.end(), f) == pore_facets.end())
                {
                    pore_facets.push_back(f);
                }
            }
        }
        
        node_it->center   = PoreCentroid(pore_vertices);
        node_it->location = PoreLocation(pore_facets);
    }
}


// Compute void/solid area, conductance permeability and equivalent radius for each edge in the
// pore network. Conductance permeability and equivalent radius are computed through the weighted
// arithmetic mean (for throats in parallel) and the harmonic mean (for the pore-throat-pore
// series). For full details, refer to the publication shown at the top of this file. Notation
// below follows the formula in that paper as much as possible.
void PoreNet::CalcEdgeWeights(glob::FacetValueMap& facet_void_area,
                              glob::FacetValueMap& facet_solid_area)
{
    for (EdgeIter edge_it = edges_.begin();
         edge_it != edges_.end(); ++edge_it)
    {
        Node* p1 = edge_it->edge_nodes.first;
        Node* p2 = edge_it->edge_nodes.second;
        
        double t_conductance       = 0.0;
        double t_permeability      = 0.0;
        double t_equivalent_length = 0.0;
        double t_equivalent_radius = 0.0; // stores r_eqv^4
        
        for (CellEdgeIter cell_edge_it = edge_it->cell_edges.begin();
             cell_edge_it != edge_it->cell_edges.end(); ++cell_edge_it)
        {
            glob::FacetCells facet_cells = std::make_pair(cell_edge_it->first, cell_edge_it->second);
            
            // Assume the throat is a circle with the same area as the void area on the facet
            double t_radius      = CGAL::sqrt(facet_void_area[facet_cells] / CGAL_PI);
            double cell_distance = CellDistance(facet_cells);
            
            edge_it->void_area  += facet_void_area[facet_cells];
            edge_it->solid_area += facet_solid_area[facet_cells];
            
            t_conductance       += Conductance(t_radius, cell_distance) * facet_void_area[facet_cells];
            t_permeability      += Permeability(t_radius) * facet_void_area[facet_cells];
            t_equivalent_length += cell_distance * facet_void_area[facet_cells];
            t_equivalent_radius += pow(t_radius, 4);
        }
        
        t_conductance       = t_conductance / edge_it->void_area;
        t_permeability      = t_permeability / edge_it->void_area;
        t_equivalent_length = t_equivalent_length / edge_it->void_area;
        
        // First assume the pores are spheres, compute the radii
        double p1_radius = pow( 3 * p1->void_volume / (4 * CGAL_PI), 1/3.);
        double p2_radius = pow( 3 * p2->void_volume / (4 * CGAL_PI), 1/3.);

        // Now assume the pores are cylinders, compute the lengths.
        double p1_length = p1->void_volume / (CGAL_PI * pow(p1_radius, 2));
        double p2_length = p2->void_volume / (CGAL_PI * pow(p2_radius, 2));
        
        // 1. Divide by pore length by 2 because otherwise every pore would be double-counted.
        // 2. Subtract half the throat length, otherwise the pore tubes would be too long.
        //p1_length = 0.5 * p1_length - 0.5 * t_equivalent_length;
        //p2_length = 0.5 * p2_length - 0.5 * t_equivalent_length;
        
        // Finally, get the conductance and permeability using the Hagen-Poiseuille equation
        // and the harmonic mean
        double p1_conductance  = Conductance(p1_radius, p1_length);
        double p2_conductance  = Conductance(p2_radius, p2_length);
        double p1_permeability = Permeability(p1_radius);
        double p2_permeability = Permeability(p2_radius);
        
        edge_it->equivalent_length = p1_length + t_equivalent_length + p2_length;
        
        edge_it->conductance = edge_it->equivalent_length /
                               ( (p1_length / p1_conductance) +
                                 (t_equivalent_length / t_conductance) +
                                 (p2_length / p2_conductance) );
        
        edge_it->permeability = edge_it->equivalent_length /
                                ( (p1_length / p1_permeability) +
                                  (t_equivalent_length / t_permeability) +
                                  (p2_length / p2_permeability) );
        
        // The equivalent_radius is the radius of the equivalent cylinder, replacing
        // the three cylinders representing the pore-throat(s)-pore series.
        edge_it->equivalent_radius = pow(edge_it->equivalent_length /
                                         ( (p1_length / pow(p1_radius, 4)) +
                                           (t_equivalent_length / t_equivalent_radius) +
                                           (p2_length / pow(p2_radius, 4)) ), 1/4.);
        
        // Use Hagen-Poisuille to get the corresponding conductance. For the length, use the
        // actual distance between the pore centroids, rather than the equivalent_length.
        double pore_distance = PointDistance(p1->center, p2->center);
        edge_it->equivalent_conductance = Conductance(edge_it->equivalent_radius, pore_distance);
    }
}


double PoreNet::Conductance(double r, double L)
{
    double viscosity = 0.001002; // at room temperature, 293 K
    return (CGAL_PI * std::pow(r, 4)) / (8 * L * viscosity);
}

double PoreNet::Permeability(double r)
{
    return std::pow(r, 2) / 8 ;
}

double PoreNet::CellDistance(glob::FacetCells facet_cells)
{
    glob::Point3D c1_center = CellCentroid(facet_cells.first);
    glob::Point3D c2_center = CellCentroid(facet_cells.second);
    
    return PointDistance(c1_center, c2_center);
}

double PoreNet::PointDistance(glob::Point3D p1, glob::Point3D p2)
{
    return CGAL::sqrt(CGAL::squared_distance(p1, p2));
}

glob::Point3D PoreNet::CellCentroid(glob::CellHandle c)
{
    // The average coordinates give the centroid of the tetrahedron
    double x_avg = 0.0;
    double y_avg = 0.0;
    double z_avg = 0.0;
    
    for(int i = 0; i < 4; ++i)
    {
        x_avg += c->vertex(i)->point().x();
        y_avg += c->vertex(i)->point().y();
        z_avg += c->vertex(i)->point().z();
    }
    
    x_avg = x_avg / 4;
    y_avg = y_avg / 4;
    z_avg = z_avg / 4;
    
    return glob::Point3D(x_avg, y_avg, z_avg);
}

glob::Point3D PoreNet::PoreCentroid(glob::VertexVector vertices)
{
    double x_avg = 0.0;
    double y_avg = 0.0;
    double z_avg = 0.0;
    
    for(glob::ConstVertexIter vertex_it = vertices.begin();
        vertex_it != vertices.end(); ++vertex_it)
    {
        x_avg += (*vertex_it)->point().x();
        y_avg += (*vertex_it)->point().y();
        z_avg += (*vertex_it)->point().z();
    }
    
    x_avg = x_avg / vertices.size();
    y_avg = y_avg / vertices.size();
    z_avg = z_avg / vertices.size();
    
    return glob::Point3D(x_avg, y_avg, z_avg);
}

std::string PoreNet::PoreLocation(glob::FacetVector facets)
{
    // Assume the pore is located in the middle and check
    // for anything different
    for(glob::ConstFacetIter facet_it = facets.begin();
        facet_it != facets.end(); ++facet_it)
    {
        // For the first vertex on the facet, get the location
        int i = 1;
        glob::VertexHandle vertex = facet_it->first->vertex((facet_it->second+i) % 4);
        std::string first_location = vertex->info().location;
        
        // Check if remaining vertices on the facet have the same location
        bool same_location = true;
        for(int i = 2; i < 4; ++i)
        {
            glob::VertexHandle vertex = facet_it->first->vertex((facet_it->second+i) % 4);
            std::string other_location = vertex->info().location;
            
            if (first_location != other_location)
                same_location = false;
        }
        
        // If all all vertices have the same location, and this location is
        // not "middle", then return the this location
        if( (first_location != "middle") && same_location )
        {
            return first_location;
        }
    }
    
    return "middle";
}
