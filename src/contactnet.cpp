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

#include "contactnet.h"

int ContactNet::NumNodes()
{
    return nodes_.size();
}


int ContactNet::NumEdges()
{
    return edges_.size();
}


ContactNet::NodeIter ContactNet::FindNode(glob::VertexHandle vertex)
{
    NodeIter node_end = nodes_.end();
    
    for (NodeIter node_it = nodes_.begin(); node_it != node_end; ++node_it)
    {
        // Nodes in the contact network are identified through the
        // corresponding vertex of the particle centroids
        if (node_it->vertex == vertex)
        {
            return node_it;
        }
    }
    
    return node_end;
}


void ContactNet::AddNode(glob::VertexHandle vertex,
                         double volume,
                         double area)
{
    NodeIter node_it = FindNode(vertex);
    NodeIter missing = nodes_.end();
    
    if( node_it == missing )
    {
        Node node;
        node.vertex = vertex;
        
        // Volume of the intersection with neighbouring
        // particles will be removed in AddEdge()
        node.solid_volume = volume;
        node.void_volume  = 0.0; // currently not used
        node.surface_area = area;
        node.location = vertex->info().location;
        
        nodes_.push_back(node);
        
        vertex_to_node_idx_[vertex] = nodes_.size() - 1;
    }
    else
    {
        std::cout << "Duplicate contact network node found." << std::endl;
    }
    
    
}


ContactNet::EdgeIter ContactNet::FindEdge(std::pair<int,int> edge_node_indices)
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


void ContactNet::AddEdge(glob::VertexHandle vertex,
                         glob::VertexHandle neighbor_vertex,
                         double contact_area,
                         double contact_cap_height,
                         double contact_cap_area,
                         double contact_cap_volume)
{
    // nodes_ vector has already been filled, so it's safe to use a pointer
    Node* node           = &nodes_[vertex_to_node_idx_[vertex]];
    Node* neighbor_node  = &nodes_[vertex_to_node_idx_[neighbor_vertex]];
    
    // If looping over all vertices (as is done in pack.cpp) to add the contact network
    // edges, every edge will be added twice (a<->b and b<->a). That's why it is checked
    // here if the edge has already been added.
    std::pair<int,int> edge_node_indices = std::minmax(vertex_to_node_idx_[vertex],
                                                       vertex_to_node_idx_[neighbor_vertex]);
    EdgeIter edge_it = FindEdge(edge_node_indices);
    
    if (edge_it == edges_.end())
    {
        Edge e;
        e.contact_area      = contact_area;
        e.edge_nodes        = std::make_pair(node, neighbor_node);
        e.void_area         = 0.0; // not used
        e.conductance       = 0.0; // not used
        e.conductivity      = 0.0; // not used
        edges_.push_back(e);
        
        edge_to_edge_idx[edge_node_indices] = edges_.size() - 1;
    }
    
    node->surface_area -= contact_cap_area;
    node->solid_volume -= contact_cap_volume;
    
    // Check for negative volume/area due to rounding
    // (theoretically possible, unlikely in practice)
    if( node->solid_volume < 0 )
    {
        node->solid_volume = 0.0;
    }
    if( node->surface_area < 0 )
    {
        node->surface_area = 0.0;
    }
}


ContactNet::ConstNodeIter ContactNet::NodesBegin()
{
    return nodes_.begin();
}


ContactNet::ConstNodeIter ContactNet::NodesEnd()
{
    return nodes_.end();
}

ContactNet::ConstEdgeIter ContactNet::EdgesBegin()
{
    return edges_.begin();
}


ContactNet::ConstEdgeIter ContactNet::EdgesEnd()
{
    return edges_.end();
}