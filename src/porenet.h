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

#ifndef _porenet_h_included_
#define _porenet_h_included_

// -------------------------------------------------------------------------------- INCLUDES
// -----------------------------------------------------------------------------------------
#include "globals.h"
#include "particles.h"

#include <vector>
#include <string>

// --------------------------------------------------------------------------------- CLASSES
// -----------------------------------------------------------------------------------------

// Contains contents of the pore network, defined by nodes (pores, i.e. encapsulated void
// spaces) and edges (if pores are connected by a throat), and related functions.
class PoreNet
{
private:
    struct Node; // Forward declaration

    // Two nodes are connected by an edge if the nodes (pores) are connected by a throat.
    struct Edge
    {
        std::pair<Node*,Node*> edge_nodes;
        
        // List of all pairs of adjacent cells (stored in FacetCells type) that make up 
        // this edge, i.e. that connect the two nodes (pores) in edge_nodes.
        std::vector<glob::FacetCells> cell_edges;

        // Summed void/solid area of facets corresponding to the adjacent cells in cell_edges.
        double void_area;
        double solid_area;

        // Conductance of the throat (see CalcEdgeWeights() for details)
        double conductance;
        
        // Permeability of the throat (see CalcEdgeWeights() for details)
        double permeability;
        
        // Equivalent radius of the throat (see CalcEdgeWeights() for details)
        double equivalent_radius;
        
        // Equivalent length of the cylinder representing the pore-throat(s)-pore series
        double equivalent_length;
        
        // Equivalent conductance of the cylinder representing the pore-throat(s)-pore series
        double equivalent_conductance;
    };
    
    // Every set of cells, as defined by the disjoint-set, is assigned a node in the pore
    // network and represents a pore in the physical sense.
    struct Node
    {
        std::vector<Edge*> neighbors;
        glob::CellVector cells;

        // Cell at the root node of the disjoint-set, acting as an identifier for this node.
        glob::CellHandle root_cell;

        // Summed void/solid volume of the cells
        double void_volume;
        double solid_volume;

        // Summed surface area (adjacent to the void space) in the cells
        double surface;

        // Location (top, middle or bottom) of the pore (see PoreLocation() for details)
        std::string location;

        // Average coordinates of the cell centroids, defining the pore centroid.
        glob::Point3D center;
    };
    

    typedef std::vector<Node>                               NodeVector;
    typedef std::vector<Edge>                               EdgeVector;
    typedef NodeVector::iterator                            NodeIter;
    typedef EdgeVector::iterator                            EdgeIter;
    typedef std::vector<glob::FacetCells>::const_iterator   CellEdgeIter;
    typedef std::vector<glob::CellHandle>::iterator         NodeCellIter;
    NodeVector                              nodes_;
    EdgeVector                              edges_;

    // Convenience map from any vertex to index of the node in the nodes_ vector that
    // the vertex belongs to. Using integer indices because vectors resize() after which
    // pointers are invalidated.

    // Convenience map from any cell to index of the node in the nodes_ vector that the 
    // cell belongs to. Using integer indices because vectors resize() after which 
    // pointers are invalidated.
    std::map<glob::CellHandle,int> cell_to_node_idx_;      
    
    // Convenience map from any edge (defined by indices of adjacent nodes in the nodes_
    // vector) to the index of the edge in the edges_ vector. This map is used for
    // (O(log(n))) edge find operations, which is much faster than if we were to do find
    // operations on the edges_ vector.


    // Convenience map from any edge (defined by indices of adjacent nodes in the nodes_ 
    // vector) to the index of the edge in the edges_ vector. This map is used for 
    // (O(log(n))) edge find operations, which is much faster than if we were to do find 
    // operations on the edges_ vector.
    std::map<std::pair<int, int>, int> edge_to_edge_idx;      
    
    // Finds a node in the contact network using the node's identifier (root_cell).
    NodeIter FindNode(glob::CellHandle root_cell);

    // Finds an edge by providing the indices of the nodes in the nodes_ vector of the two
    // adjacent nodes. Uses the edge_to_edge_idx map to find the edge.
    EdgeIter FindEdge(std::pair<int,int> edge_node_indices);
    
    // Determines if the two cells in facet_cells are equal to cell and neighbor_cell
    bool CellEdgeEquals(glob::FacetCells facet_cells, glob::CellHandle cell,
                                                      glob::CellHandle neighbor_cell);
    
    // Finds an edge in edge_it->cells using the provided cell and neighbor_cell.
    CellEdgeIter FindCellEdge(EdgeIter edge_it, glob::CellHandle cell,
                                                glob::CellHandle neighbor_cell);
    
    // Implements the Hagen-Poiseuille equation to compute conductance
    double Conductance(double r, double L);
    
    // Implements the Hagen-Poiseuille equation to compute permeability
    double Permeability(double r);

    // Computes the distance between the centroids of the cells in facet_cells.
    double CellDistance(glob::FacetCells facet_cells);
    
    // Computes the distance between the points p1 and p2
    double PointDistance(glob::Point3D p1, glob::Point3D p2);
    
    // Computes the centroid of c by averaging the coordinates of the vertices.
    glob::Point3D CellCentroid(glob::CellHandle c);

    // Computes the centroid of the pore, by averaging the coordinates of the vertices.
    glob::Point3D PoreCentroid(glob::VertexVector vertices);

    // Returns the location of the pore (top, middle or bottom). The criterion for
    // top (or bottom) is that at least one of the facets is top (or bottom), which,
    // in turn, requires all vertices on the facet to be top (or bottom) - see
    // SetVertexLocations() in the Triangulation class.
    std::string PoreLocation(glob::FacetVector facets);
    
public:
    typedef std::vector <Node>::const_iterator  ConstNodeIter;
    typedef std::vector <Edge>::const_iterator  ConstEdgeIter;
    int                                         NumNodes();
    int                                         NumEdges();
    ConstNodeIter                               NodesBegin();
    ConstNodeIter                               NodesEnd();
    ConstEdgeIter                               EdgesBegin();
    ConstEdgeIter                               EdgesEnd();

    // Adds a node to the pore network. Checks if a node (identified by root_cell)
    // already exists, and if so, adds cell to the cells vector of that node.
    // Otherwise, a new node is created.
    void AddNode(glob::CellHandle cell, glob::CellHandle root_cell);
    
    // Adds an edge to the pore network, defined by a cell and its neighbor. If the
    // edge already exists but this connection between cell and neighbor_cell is not
    // contained yet in the cell_edges of this edge, then it's added.
    void AddEdge(glob::CellHandle cell, glob::CellHandle neighbor_cell);
    
    // Assuming the cells vector of every node has been filled, this function loops
    // over those cells and computes the summed void volume, solid volume and surface
    // area of the node (pore), as well as the pore centroid and location.
    void AggregateNodeCells(glob::CellValueMap& cell_void_volume, 
                            glob::CellValueMap& cell_solid_volume,
                            glob::CellValueMap& cell_surface_area);
    
    // Assuming cell_edges vector of every edge has been filled, this function loops over
    // those pairs of cells and computes the void area, solid area, conductance and permeability.
    void CalcEdgeWeights(glob::FacetValueMap& facet_void_area, 
                         glob::FacetValueMap& facet_solid_area);
};

#endif