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

#ifndef _disjointset_h_included_
#define _disjointset_h_included_

// -------------------------------------------------------------------------------- INCLUDES
// -----------------------------------------------------------------------------------------
#include "globals.h"

#include <vector>
#include <cassert>

// --------------------------------------------------------------------------------- CLASSES
// -----------------------------------------------------------------------------------------

// Implementation of a disjoint-set data structure, based on code by Emil Stefanov (original
// link is offline, but see also https://en.wikipedia.org/wiki/Disjoint-set_data_structure ).
// Note: CGAL also has a union_find (=disjoint-set) class but implementing it here provides
// better insight into how it works.
class DisjointSet
{
private:
    
    // Each disjoint-set is defined as a tree structure (a set), consisting of nodes (not to
    // be confused with pore/contact network nodes). Unless a node is a root node, the node
    // has a parent node. Nodes are identified by a CellHandle (from the triangulation).
    struct Node
    {
        // Roughly represents the maximum height of the node in the tree.
        int                 rank;
        
        // The identifier of this node
        glob::CellHandle    cell;
        
        // The parent node of this node
        Node*               parent;
    };
    
    // The num_sets_ sets (trees) together make up all elements (cells).
    int num_sets_;
    
    // Map from identifier (CellHandle) to a node in disjoint-set data structure
    std::map<glob::CellHandle,Node*> nodes_;
    
public:

    ~DisjointSet();
    
    // Initializes a disjoint-set with the interior cells, as defined in pack.cpp.
    void Initialize(glob::ConstCellIter interior_cells_begin,
                    glob::ConstCellIter interior_cells_end);
    
    // Find the root node for cell, by moving up the tree until the root node is reached.
    Node* FindRoot(glob::CellHandle c);
    
    // Merge two sets into one, by making one of the nodes the parent of the other node.
    void Merge(glob::CellHandle c1, glob::CellHandle c2);
    
    // Flatten the tree, by attaching every node directly to its root node. Also called
    // "path compression".
    void FlattenTree();
    
    // Returns the parent cell of a cell
    glob::CellHandle GetParentCell(glob::CellHandle c);
};

#endif