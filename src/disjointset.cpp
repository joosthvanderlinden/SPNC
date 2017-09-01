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

#include "disjointset.h"

DisjointSet::~DisjointSet()
{
    std::map<glob::CellHandle,Node*>::iterator it;
    
    for(it = nodes_.begin(); it != nodes_.end(); ++it)
        delete it->second;
    
    nodes_.clear();
    num_sets_ = 0;
}

void DisjointSet::Initialize(glob::ConstCellIter interior_cells_begin,
                             glob::ConstCellIter interior_cells_end)
{
    for(glob::ConstCellIter cell_it = interior_cells_begin;
        cell_it != interior_cells_end; ++cell_it)
    {
        Node* n          = new Node();
        n->parent        = NULL;
        n->cell          = *cell_it;
        n->rank          = 0;
        nodes_[*cell_it] = n;
    }
    
    num_sets_ = std::distance(interior_cells_begin, interior_cells_end);
}


DisjointSet::Node* DisjointSet::FindRoot(glob::CellHandle c)
{
    // Check that node exists
    Node* node = nodes_[c];
    assert(node != NULL);
    
    // Find the root node of the set that c belongs to
    while (node->parent != NULL)
    {
        node = node->parent;
    }
    Node* root_node = node;
    
    // Walk from the node that c belongs to, to the root note, and make the traversed nodes the direct
    // children of the root node (path compression). This optimizes the tree for future FindRoot invokations.
    node = nodes_[c];
    while (node != root_node)
    {
        Node* next   = node->parent;
        node->parent = root_node;
        node         = next;
    }
    
    return root_node;
}

void DisjointSet::Merge(glob::CellHandle c1, glob::CellHandle c2)
{
    Node* c1_root = FindRoot(c1);
    Node* c2_root = FindRoot(c2);
    
    if(c1_root == c2_root)
        return; // already merged

    // Determine which node representing a set has a higher rank. The node with the higher rank is
    // likely to have a bigger subtree so in order to better balance the tree representing the
    // union, the node with the higher rank is made the parent of the one with the lower rank and
    // not the other way around.
    if(c1_root->rank > c2_root->rank)
    {
        c2_root->parent = c1_root;
    }
    else if(c1_root->rank < c2_root->rank)
    {
        c1_root->parent = c2_root;
    }
    else
    {
        c2_root->parent = c1_root;
        ++c1_root->rank; // update rank
    }
    
    --num_sets_;
}


void DisjointSet::FlattenTree()
{
    // Attach every node in a set directly to it's root node. This improves efficiency.
    Node* temp;
    for(std::map<glob::CellHandle, Node*>::const_iterator it = nodes_.begin();
        it != nodes_.end(); ++it)
    {
        temp = FindRoot(it->first);
    }
}

glob::CellHandle DisjointSet::GetParentCell(glob::CellHandle c)
{
    Node* node = nodes_[c];
    // If node has no parent then it is a root node
    if( node->parent != NULL)
        return node->parent->cell;
    else
        return c;
}


