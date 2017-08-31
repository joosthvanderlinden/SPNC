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
// Author:  Joost van der Linden <joostv@student.unimelb.edu.au>
//
// Please cite the following paper if you use this code:
//
//          J.H. van der Linden, A. Sufian, G. Narsilio, A.R. Russell, A. Tordesillas,
//          (2016), Delaunay-based pore network construction for granular packings,
//          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


#ifndef _contactnet_h_included_
#define _contactnet_h_included_

// -------------------------------------------------------------------------------- INCLUDES
// -----------------------------------------------------------------------------------------
#include "globals.h"
#include "particles.h"

#include <vector>
#include <string>

// ----------------------------------------------------------------------------------- CLASS
// -----------------------------------------------------------------------------------------

// Contains contents of the contact network, defined by nodes (particle centroids) and edges
// (if particles are in contact), and related functions.
class ContactNet
{
private:
    
    struct Node; // Forward declaration
    
    // Two nodes are connected by an edge if the corresponding particles touch or overlap.
    struct Edge
    {
        std::pair <Node*, Node*>    edge_nodes;
        double                      contact_area;
        double                      void_area;      // Currently not used
        double                      conductance;    // Currently not used
        double                      conductivity;   // Currently not used
    };
    
    // Every particle in the assembly is assigned a node in the contact network.
    struct Node
    {
        std::vector <Edge*>         neighbors;
        double                      solid_volume;   // Volume of the particle
        double                      void_volume;    // Currently not used
        
        // Handle of the vertex in the triangulation (= centroid of the particle).
        // This handle is the unique identifier used
        glob::VertexHandle          vertex;
        
        // Surface area adjecent to the void space
        double                      surface_area;
        
        // Location of the node in the packing. Can be "top" (above limits), "bottom"
        // (below limits) or "middle" (anything else), where the limits are defined in
        // the Particles class.
        std::string                 location;
    };

    typedef std::vector<Node>       NodeVector;
    typedef std::vector<Edge>       EdgeVector;
    typedef NodeVector::iterator    NodeIter;
    typedef EdgeVector::iterator    EdgeIter;
    NodeVector                      nodes_;
    EdgeVector                      edges_;

    // Convenience map from any vertex to index of the node in the nodes_ vector that
    // the vertex belongs to. Using integer indices because vectors resize() after which
    // pointers are invalidated.
    std::map<glob::VertexHandle,int> vertex_to_node_idx_;
    
    // Convenience map from any edge (defined by indices of adjacent nodes in the nodes_
    // vector) to the index of the edge in the edges_ vector. This map is used for
    // (O(log(n))) edge find operations, which is much faster than if we were to do find
    // operations on the edges_ vector.
    std::map<std::pair<int, int>, int> edge_to_edge_idx;
    
    // Finds a node in the contact network using the node's identifier (VertexHandle).
    NodeIter FindNode(glob::VertexHandle vertex);
    
    // Finds an edge by providing the indices of the nodes in the nodes_ vector of the two
    // adjacent nodes. Uses the edge_to_edge_idx map to find the edge.
    EdgeIter FindEdge(std::pair<int,int> edge_node_indices);
    
public:
    typedef NodeVector::const_iterator      ConstNodeIter;
    typedef EdgeVector::const_iterator      ConstEdgeIter;
    int                                     NumNodes();
    int                                     NumEdges();
    ConstNodeIter                           NodesBegin();
    ConstNodeIter                           NodesEnd();
    ConstEdgeIter                           EdgesBegin();
    ConstEdgeIter                           EdgesEnd();
    
    // Adds a node to the contact network. Assumes the location has been set in vertex->info().
    void AddNode(glob::VertexHandle vertex, double volume, double area);
    
    // Adds an edge to the contact network, defined by the vertices of the two contacting particles.
    void AddEdge(glob::VertexHandle vertex, glob::VertexHandle neighbor_vertex,
                 double contact_area,
                 double contact_cap_height,
                 double contact_cap_area,
                 double contact_cap_volume);
};

#endif