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

#include "pack.h"

void Pack::ReadPack(const std::string inputfile)
{
    std::ifstream file;
    file.open(inputfile);
    
    if(!file.good())
    {
        std::cout << "File not found." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // First 4 lines are assumed to be comments
    std::string temp;
    for(int i=0; i<4; ++i)
    {
        getline(file,temp);
    }
    
    // Read simulation parameter data
    particles_.ReadData(file);
    
    // Read subsample limits
    // (token-based parsing doesn't recognize newlines)
    getline(file,temp); getline(file,temp);
    particles_.ReadLimits(file);
    getline(file,temp);
    
    // Read particle data
    particles_.ReadParticles(file);
    
    file.close();
    return;
}

void Pack::InitializeTriangulation()
{
    triangulation_.Initialize(particles_.ParticlesBegin(),
                              particles_.ParticlesEnd());
    
    triangulation_.SetVertexLocations(particles_);
}


void Pack::InitializeRangeTree()
{
    range_tree_.Initialize(triangulation_.FiniteVerticesBegin(),
                           triangulation_.FiniteVerticesEnd());
}


void Pack::InitializeDisjointSet()
{
    disjoint_set_.Initialize(triangulation_.InteriorCellsBegin(),
                             triangulation_.InteriorCellsEnd());
}

void Pack::FindInteriorCellsAndFacets()
{
    triangulation_.FindInteriorCells(particles_);
    triangulation_.FindInteriorFacets(particles_);
}


glob::ConstFacetIter Pack::InteriorFacetsEnd()
{
    return triangulation_.InteriorFacetsEnd();
}

int Pack::NumInteriorFacets()
{
    return triangulation_.NumInteriorFacets();
}

void Pack::MergeFacetCells(glob::Facet f)
{
    glob::CellHandle c1 = f.first;
    glob::CellHandle c2 = f.first->neighbor(f.second);
    
    disjoint_set_.Merge(c1, c2);
}

glob::VertexVector Pack::GetFacetVertices(glob::Facet f)
{
    return triangulation_.GetFacetVertices(f);
}

glob::ConstFacetIter Pack::InteriorFacetsBegin()
{
    return triangulation_.InteriorFacetsBegin();
}

glob::VertexVector Pack::FindNearbyVertices(glob::VertexVector facet_vertices)
{
    // Increase the search window by the maximum particle radius
    double window_offset = particles_.MaxParticleRadius();
    glob::VertexVector window_vertices = range_tree_.WindowQuery(facet_vertices,
                                                                 window_offset);
    
    glob::VertexVector nearby_vertices;
    for(glob::ConstVertexIter vertex_it = window_vertices.begin();
        vertex_it != window_vertices.end(); ++vertex_it)
    {
        // Check if sphere at vertex_it is less than it's radius
        // away from the plane that the facet resides in.
        if (triangulation_.SphereFacetIntersect(*vertex_it, facet_vertices))
        {
            nearby_vertices.push_back(*vertex_it);
        }
    }
    
    return nearby_vertices;
}

void Pack::SetFacetAreas(glob::Facet f, double void_area, double solid_area)
{
    // Add the areas to the map twice, for the pairs (c1,c2) and (c2,c1)
    // Not very efficient in terms of storage, but it gets the job done
    glob::CellHandle c1 = f.first;
    glob::CellHandle c2 = f.first->neighbor(f.second);
    
    facet_void_area_[std::make_pair(c1, c2)] = void_area;
    facet_void_area_[std::make_pair(c2, c1)] = void_area;
    
    facet_solid_area_[std::make_pair(c1, c2)] = solid_area;
    facet_solid_area_[std::make_pair(c2, c1)] = solid_area;
}

void Pack::CellParameters(Configure config)
{
    std::chrono::time_point<std::chrono::system_clock> start;
    start = std::chrono::system_clock::now();
    
    int count     = 0;
    int num_cells = triangulation_.NumInteriorCells();
    
    for(glob::ConstCellIter cell_it = triangulation_.InteriorCellsBegin();
        cell_it != triangulation_.InteriorCellsEnd(); ++cell_it)
    {
        std::cout << "\r" << "Computing cell volume and area: " << ++count << "/"
                  << num_cells << std::flush;
        
        // Currently, the cell surface area can only be computed analytically and
        // does not account for overlapping particles.
        AnalyticalCellSurfaceArea(*cell_it);
        
        // Next, retrieve the four vertices for this cell, and do a window query.
        glob::VertexVector cell_vertices(4);
        for(int i = 0; i < 4; ++i)
        {
            cell_vertices[i] = (*cell_it)->vertex(i);
        }
        
        // Increase the search window by the maximum particle radius (which is the
        // furthest away an intersecting particle could possibly be).
        double offset = particles_.MaxParticleRadius();
        
        // lower_corner and upper_corner are calculated in the WindowQuery().
        // Retrieve them for the voxelization of the vicinity of the cell.
        glob::VertexVector nearby_vertices = range_tree_.WindowQuery(cell_vertices, offset);
        
        // If no other vertices are inside the window, and the overlap between
        // particles is non-existent or negligible (e.g. hard spheres) then the
        // volume can be computed analytically. This is often not the case. Comment
        // out the nearby_vertices.empty() below to significantly speed up the
        // cell volume calculation.
        bool no_particle_overlap = config.GetNoParticleOverlap();
        
        if (nearby_vertices.empty() && no_particle_overlap)
        {
            AnalyticalCellVolume(*cell_it);
        }
        else
        {
            // lower_corner and upper_corner have already been calculated
            // in the WindowQuery(), so retrieve them.
            glob::Point3D lower_corner = range_tree_.WindowLowerCorner();
            glob::Point3D upper_corner = range_tree_.WindowUpperCorner();
            
            VoxelizedCellVolume(*cell_it, lower_corner, upper_corner,
                                cell_vertices, nearby_vertices);
        }
    }
    
    std::chrono::duration<double> duration = (std::chrono::system_clock::now() - start);
    std::cout << " done in " << duration.count() << " seconds." << std::endl;
    
}

void Pack::AnalyticalCellSurfaceArea(glob::CellHandle cell)
{
    // NOTE: this function only works well if
    //       - overlap between spheres is negligible
    //       - spheres located at tetrahedron vertices are the only spheres intersecting the tetrahedron
    //       - triangles are not flat (ie spheres at a vertex do not stick out of tetrahedron on the other side)
    
    // A cell has 4 vertices, retrieved through vertex(i)
    glob::Point3D   p[4];
    double          r[4];
    for(int i = 0; i < 4; ++i)
    {
        p[i] = cell->vertex(i)->point();
        r[i] = cell->vertex(i)->info().r;
    }
    
    double surface_area = 0.0;
    for(int i = 0; i < 4; ++i)
    {
        // Compute the angle subtended by tetrahedron at vertex i
        double solid_angle = triangulation_.SolidAngle(p, i);
        
        // Add the surface area for each spherical excess
        surface_area += solid_angle * pow(r[i],2);
    }
    
    cell_surface_area_[cell] = surface_area;
}


void Pack::AnalyticalCellVolume(glob::CellHandle cell)
{
    // NOTE: this function only works well if
    //       - overlap between spheres is negligible
    //       - spheres located at tetrahedron vertices are the only spheres intersecting the tetrahedron
    //       - triangles are not flat (ie spheres at a vertex do not stick out of tetrahedron on the other side)
    // if the latter case occurs, void_volume will sometimes be negative
    //
    // This method also accounts for the top/bottom boundaries.
    
    // A cell has 4 vertices, retrieved through vertex(i)
    glob::Point3D   p[4];
    double          r[4];
    for(int i = 0; i < 4; ++i)
    {
        p[i] = cell->vertex(i)->point();
        r[i] = cell->vertex(i)->info().r;
    }

    // Void volume = tetrahedron volume - solid volume
    double tetrahedron_volume  = CGAL::volume(p[0],p[1],p[2],p[3]);
    double void_volume = tetrahedron_volume;
    
    for(int i = 0; i < 4; ++i)
    {
        // Compute the angle subtended by tetrahedron at vertex i
        double solid_angle = triangulation_.SolidAngle(p, i);

        // Subtract the spherical excess volume, where
        // volume = [solid angle] * r^2 / [sphere area] * [sphere volume] = solid_angle / 3 * r^3
        void_volume -= solid_angle / 3 * pow(r[i],3);
    }
    
    // If the volume is negative, there is too much overlap to use the analytical
    // approach and the voxelized approach should be used instead.
    if( void_volume < 0 )
    {
        std::cout << std::endl << "Negative cell volume detected (probably due to particle overlap or"
                  << std::endl << " elongated tetrahedra). Switch to voxelized approach." << std::endl;
        exit(1);
    }
    
    cell_void_volume_[cell]  = void_volume;
    cell_solid_volume_[cell] = tetrahedron_volume - void_volume;
}


void Pack::VoxelizedCellVolume(glob::CellHandle cell,
                               glob::Point3D lower_corner,
                               glob::Point3D upper_corner,
                               glob::VertexVector cell_vertices,
                               glob::VertexVector nearby_vertices)
{
    // After voxelizing the vicinity of the cell, check if
    // - Voxel is inside tetrahedron
    // - Voxel is inside spheres residing at vertices of tetrahedron
    // - Voxel is inside other spheres that intersect the tetrahedron
    // - Voxel is inside subsample limits
    
    double Dx = particles_.ParticleVoxelSize(10);
    double Dy = Dx;
    double Dz = Dx;
    double voxel_volume = Dx*Dy*Dz;

    int Nx = std::round( (upper_corner.x() - lower_corner.x()) / Dx);
    int Ny = std::round( (upper_corner.y() - lower_corner.y()) / Dy);
    int Nz = std::round( (upper_corner.z() - lower_corner.z()) / Dz);
    
    // Using variable length arrays, even though it's not recommended
    // (speed gain over using std::vectors is significant)
    double* x = new double[Nx];
    for(int i=0; i<Nx; ++i)
    {
        x[i] = lower_corner.x() + Dx*i;
    }
    
    double* y = new double[Ny];
    for(int i=0; i<Ny; ++i)
    {
        y[i] = lower_corner.y() + Dy*i;
    }
    
    double* z = new double[Nz];
    for(int i=0; i<Nz; ++i)
    {
        z[i] = lower_corner.z() + Dz*i;
    }
    
    double void_volume  = 0.0;
    double solid_volume = 0.0;
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<Nz; ++k)
            {
                // query_point is the center of the voxel
                glob::Point3D query_point = glob::Point3D(x[i] + Dx/2, y[j] + Dy/2, z[k] + Dz/2);
                
                if (particles_.InsideLimits(query_point) &&
                    triangulation_.InsideTetrahedron(query_point, cell))
                {
                    bool in_void = true;
                    
                    // Test if query_point is inside any of the spheres at the 4 cell vertices
                    for(glob::ConstVertexIter cell_v_it  = cell_vertices.begin();
                        cell_v_it != cell_vertices.end(); ++cell_v_it)
                    {
                        if (particles_.WithinSphere(query_point,
                                                    (*cell_v_it)->point(),
                                                    (*cell_v_it)->info().r))
                            in_void = false;
                    }
                    
                    // If in_void is still true then test if query_point is inside any of the nearby spheres
                    if (in_void)
                    {
                        // Loop over nearby vertices
                        for(glob::ConstVertexIter near_v_it  = nearby_vertices.begin();
                            near_v_it != nearby_vertices.end(); ++near_v_it)
                        {
                            if (particles_.WithinSphere(query_point,
                                                        (*near_v_it)->point(),
                                                        (*near_v_it)->info().r))
                                in_void = false;
                        }
                    }
                
                    if (in_void)
                        void_volume  += voxel_volume;
                    else
                        solid_volume += voxel_volume;
                }
            }
        }
    }

    // Cleaning up
    delete [] x;
    delete [] y;
    delete [] z;
    
    cell_void_volume_[cell]  = void_volume;
    cell_solid_volume_[cell] = solid_volume;
}

void Pack::AddContactNetworkNodes()
{
    int num_vertices = triangulation_.NumFiniteVertices();
    int counter      = 0;
    for (glob::TriangulationVertexIter v_it  = triangulation_.FiniteVerticesBegin();
         v_it != triangulation_.FiniteVerticesEnd(); ++v_it)
    {
        std::cout << "\r" << "Adding contact network nodes: "
                  << ++counter << " / " << num_vertices << std::flush;
        
        double volume = particles_.SphereVolume(v_it);
        double area   = particles_.SphereArea(v_it);
        contact_network_.AddNode(v_it, volume, area);
    }
    std::cout << " done." << std::endl;
}

void Pack::AddContactNetworkEdges()
{
    int num_nodes = contact_network_.NumNodes();
    int counter   = 0;
    for(ContactNet::ConstNodeIter node_it = contact_network_.NodesBegin();
        node_it != contact_network_.NodesEnd(); ++node_it)
    {
        std::cout << "\r" << "Adding contact network edges for each node: "
                  << ++counter << " / " << num_nodes << std::flush;
        
        glob::VertexHandle vertex = node_it->vertex; // for clarity
        glob::VertexVector window_vertices(1, vertex); // WindowQuery() takes a VertexVector
        
        // Increase the search window by max particle radius + v's radius
        double window_offset = particles_.MaxParticleRadius() + vertex->info().r;
        glob::VertexVector nearby_vertices = range_tree_.WindowQuery(window_vertices,
                                                                     window_offset);
        
        for(glob::ConstVertexIter vertex_it  = nearby_vertices.begin();
            vertex_it != nearby_vertices.end(); ++vertex_it)
        {
            glob::VertexHandle nearby_vertex = *vertex_it; // for clarity
            
            // Check if sphere for vertex intersects with any of the
            // spheres corresponding to the vertices in nearby_vertices
            double extended_radius = nearby_vertex->info().r + vertex->info().r;
            if (particles_.WithinSphere(vertex->point(),
                                        nearby_vertex->point(),
                                        extended_radius))
            {
                double contact_radius     = particles_.SphereContactRadius(vertex, nearby_vertex);
                double contact_area       = CGAL_PI * pow(contact_radius, 2);
                double contact_cap_height = particles_.SphereCapHeight(vertex, nearby_vertex);
                double contact_cap_area   = particles_.SphereCapArea(vertex->info().r, contact_cap_height);
                double contact_cap_volume = particles_.SphereCapVolume(vertex->info().r, contact_cap_height);
                
                contact_network_.AddEdge(node_it->vertex, nearby_vertex,
                                         contact_area,
                                         contact_cap_height,
                                         contact_cap_area,
                                         contact_cap_volume);
            }
        }
    }
    std::cout << " done." << std::endl;
}


void Pack::AddPoreNetworkNodes()
{
    // Flatten the tree in the disjoint-set for faster retrieval
    disjoint_set_.FlattenTree();
    
    // Add each interior cell to a pore network node, where the node
    // identifier is the root cell in the disjoint-set
    int num_cells = triangulation_.NumInteriorCells();
    int counter   = 0;
    for(glob::ConstCellIter cell_it = triangulation_.InteriorCellsBegin();
        cell_it != triangulation_.InteriorCellsEnd(); ++cell_it)
    {
        std::cout << "\r" << "Adding pore network nodes: "
                  << ++counter << " / " << num_cells << std::flush;
        
        glob::CellHandle root_cell = disjoint_set_.GetParentCell(*cell_it);
        pore_network_.AddNode(*cell_it, root_cell);
    }
    
    // Aggregate the solid volume, void volume and surface area
    // (over all cells) for every pore node
    pore_network_.AggregateNodeCells(cell_void_volume_,
                                     cell_solid_volume_,
                                     cell_surface_area_);
    
    std::cout << " done." << std::endl;
}


void Pack::AddPoreNetworkEdges()
{
    // Note: cell edge = connection between two tetrahedra
    //       pore network edge = connection between two (merged sets of) tetrahedra
    // A pore network edge is not created unless a cell edge exists between the two pores
    
    // Loop over all pores
    int num_nodes = pore_network_.NumNodes();
    int counter   = 0;
    for(PoreNet::ConstNodeIter node_it = pore_network_.NodesBegin();
        node_it != pore_network_.NodesEnd(); ++node_it)
    {
        std::cout << "\r" << "Adding pore network edges for each node: "
                  << ++counter << " / " << num_nodes << std::flush;
        
        // Loop over all cells in this pore
        for(glob::ConstCellIter cell_begin = node_it->cells.begin(),
                                cell_end   = node_it->cells.end(),
            cell_it = cell_begin; cell_it != cell_end; ++cell_it)
        {
            glob::CellHandle cell = (*cell_it);
            
            // Loop over all neighbors of this cell
            for(int i = 0; i < 4; ++i)
            {
                glob::CellHandle neighbor_cell = cell->neighbor(i);
                
                // Check if neighbor is inside limits
                //     & if facet between cell and neighbor is inside limits
                //     & if facet between cell and neighbor has void_area > 0
                //     & if neighbor is in another disjoint set
                if( (triangulation_.CellInsideLimits(neighbor_cell, particles_)) &&
                    (triangulation_.FacetInsideLimits(std::make_pair(*cell_it, i), particles_)) &&
                    (facet_void_area_[std::make_pair((*cell_it), neighbor_cell)] > 0) &&
                    (std::find(cell_begin, cell_end, neighbor_cell) == cell_end) )
                {
                    pore_network_.AddEdge(*cell_it, neighbor_cell);
                }
            }
        }
    }
    
    // Having determined the facets corresponding to each edge, calculate the weights
    pore_network_.CalcEdgeWeights(facet_void_area_,
                                  facet_solid_area_);
    
    std::cout << " done." << std::endl;
}


void Pack::WriteTriangulation(std::string path)
{
    // Initialize output
    std::ofstream file(path.append("/networks/delaunay_triangulation.txt"),std::ios::out);
    file << "# Triangulation output format: " << std::endl;
    file << "# N_vertices; [N_vertices x 4] array of coordinates & original particle ID; "
         <<   "N_cells; [N_cells x 4] array of vertex indices; "
         <<   "N_pores; [N_pores x .] list of cells in a pore; " << std::endl;
    
    // Write vertices while saving vertex row number in vertex_to_row
    std::map<glob::VertexHandle,int> vertex_to_row;
    glob::TriangulationVertexIter vertex_it  = triangulation_.FiniteVerticesBegin();
    int i = 1;
    file << triangulation_.NumFiniteVertices() << std::endl;
    for ( ; vertex_it != triangulation_.FiniteVerticesEnd(); ++vertex_it, ++i)
    {
        file << *vertex_it << " " << vertex_it->info().original_id << std::endl;
        vertex_to_row[vertex_it] = i;
    }
    
    // Write cells while saving cell row number in cell_to_row
    std::map<glob::CellHandle,int> cell_to_row;
    glob::ConstCellIter cell_it = triangulation_.InteriorCellsBegin();
    file << triangulation_.NumInteriorCells() << std::endl;
    for(i = 1; cell_it != triangulation_.InteriorCellsEnd(); ++cell_it, ++i)
    {
        cell_to_row[*cell_it] = i;
        for(int j = 0; j < 4; ++j)
        {
            file << vertex_to_row[(*cell_it)->vertex(j)] << "\t";
        }
        file << std::endl;
    }
    
    // Write pores
    file << pore_network_.NumNodes() << std::endl;
    for(PoreNet::ConstNodeIter pore_node_it = pore_network_.NodesBegin();
        pore_node_it != pore_network_.NodesEnd(); ++pore_node_it)
    {
        for(glob::ConstCellIter pore_cell_it  = pore_node_it->cells.begin();
            pore_cell_it != pore_node_it->cells.end(); ++pore_cell_it)
        {
            file << cell_to_row[*pore_cell_it] << " ";
        }
        file << std::endl;
    }
    
    file.close();
}


void Pack::WritePoreNetwork(std::string path)
{
    int num_node_metrics = 5;
    int num_edge_metrics = 8;
    
    std::ofstream file(path.append("/networks/pore_network.txt"),std::ios::out);
    file << "# Pore network output format: " << std::endl;
    file << "num_nodes & num_node_metrics ; names of columns ; [num_nodes x num_node_metrics] array; "
         << "num_edges & num_edge_metrics ; names of columns ; [num_edges x num_edge_metrics] array" << std::endl;
    
    // Write pore nodes while saving line number
    int i = 1;
    PoreNet::ConstNodeIter node_it = pore_network_.NodesBegin();
    std::map<glob::CellHandle,int> pore_to_row;
    file << pore_network_.NumNodes() << " " << num_node_metrics << std::endl;
    file << "x y z void_volume solid_volume surface location n_cells" << std::endl;
    for(; node_it != pore_network_.NodesEnd(); ++node_it, ++i)
    {
        pore_to_row[node_it->root_cell] = i;
        file << node_it->center       << "\t" << node_it->void_volume  << "\t"
             << node_it->solid_volume << "\t" << node_it->surface      << "\t"
             << node_it->location     << "\t" << node_it->cells.size() << std::endl;
    }
    
    // Write pore edges while referring to line number of nodes above
    file << pore_network_.NumEdges() << " " << num_edge_metrics << std::endl;
    file << "from_node to_node void_area solid_area conductance permeability equivalent_radius equivalent_length equivalent_conductance n_cell_edges" << std::endl;
    for(PoreNet::ConstEdgeIter edge_it = pore_network_.EdgesBegin();
        edge_it != pore_network_.EdgesEnd(); ++edge_it)
    {
        file <<         pore_to_row[edge_it->edge_nodes.first->root_cell]
             << "\t" << pore_to_row[edge_it->edge_nodes.second->root_cell]
             << "\t" << edge_it->void_area
             << "\t" << edge_it->solid_area
             << "\t" << edge_it->conductance
             << "\t" << edge_it->permeability
             << "\t" << edge_it->equivalent_radius
             << "\t" << edge_it->equivalent_length
             << "\t" << edge_it->equivalent_conductance
             << "\t" << edge_it->cell_edges.size() << std::endl;
    }
    
    file.close();
}


void Pack::WriteContactNetwork(std::string path)
{
    int num_node_metrics = 5;
    int num_edge_metrics = 4;
    
    std::ofstream file(path.append("/networks/contact_network.txt"),std::ios::out);
    file << "# Contact network output format: " << std::endl;
    file << "num_nodes & num_node_metrics ; names of columns ; [num_nodes x num_node_metrics] array; "
         << "num_edges & num_edge_metrics ; names of columns ; [num_edges x num_edge_metrics] array" << std::endl;
    
    // Write contact nodes while saving line number
    int i = 1;
    ContactNet::ConstNodeIter node_it = contact_network_.NodesBegin();
    std::map<glob::VertexHandle,int> vertex_to_row;
    file << contact_network_.NumNodes() << " " << num_node_metrics << std::endl;
    file << "x y z void_volume solid_volume surface location original_id" << std::endl;
    for(; node_it != contact_network_.NodesEnd(); ++node_it, ++i)
    {
        vertex_to_row[node_it->vertex] = i;
        file << node_it->vertex->point() << "\t" << node_it->void_volume
             << "\t" << node_it->solid_volume
             << "\t" << node_it->surface_area
             << "\t" << node_it->location
             << "\t" << node_it->vertex->info().original_id << std::endl;
    }
    
    // Write contact network edges while referring to line number of nodes above
    file << contact_network_.NumEdges() << " " << num_edge_metrics << std::endl;
    file << "from_node to_node void_area solid_area conductance conductivity" << std::endl;
    for(ContactNet::ConstEdgeIter edge_it  = contact_network_.EdgesBegin();
        edge_it != contact_network_.EdgesEnd(); ++edge_it)
    {
        file <<         vertex_to_row[edge_it->edge_nodes.first->vertex]
             << "\t" << vertex_to_row[edge_it->edge_nodes.second->vertex]
             << "\t" << edge_it->void_area
             << "\t" << edge_it->contact_area
             << "\t" << edge_it->conductance
             << "\t" << edge_it->conductivity << std::endl;
    }
    
    file.close();
}


void Pack::WriteCellProperties(std::string path)
{
    int N_metrics = 2;
    
    // Open output file
    std::ofstream file(path.append("/networks/delaunay_cell_properties.txt"), std::ios::out);
    
    // Write header
    file << "# Cell body metrics output format: " << std::endl;
    file << "# N_cells N_metrics; [N_cells x 3] array of cell void volumes & cell solid volume & cell surface areas"  << std::endl;
    
    // Loop over entries and write
    file << triangulation_.NumInteriorCells() << " " << N_metrics << std::endl;
    file << "void_volume solid_volume surface" << std::endl;
    for(glob::ConstCellIter cell_it = triangulation_.InteriorCellsBegin();
        cell_it != triangulation_.InteriorCellsEnd(); ++cell_it)
    {
        file << cell_void_volume_[*cell_it] << "\t"
             << cell_solid_volume_[*cell_it] << "\t"
             << cell_surface_area_[*cell_it] << std::endl;
    }
    
    // Close output file
    file.close();
}


void Pack::WriteFacetProperties(std::string path)
{
    int N_metrics = 2;
    
    // Open output file
    std::ofstream file(path.append("/networks/delaunay_facet_properties.txt"),std::ios::out);
    
    // Write header
    file << "# Throat metrics output format: " << std::endl;
    file << "# N_faces N_metrics; [N_faces x 2] array of throat void area & throat solid area "  << std::endl;
    
    // Loop over facets and write
    file << triangulation_.NumInteriorFacets() << " " << N_metrics << std::endl;
    file << "void_area solid_area" << std::endl;
    for(glob::ConstFacetIter facet_it = triangulation_.InteriorFacetsBegin();
        facet_it != triangulation_.InteriorFacetsEnd(); ++facet_it)
    {
        // Check if cell is interior
        if (triangulation_.FacetInsideLimits(*facet_it, particles_))
        {
            file << facet_void_area_[std::make_pair(facet_it->first,
                                                    facet_it->first->neighbor(facet_it->second))] << "\t"
                 << facet_solid_area_[std::make_pair(facet_it->first,
                                                     facet_it->first->neighbor(facet_it->second))] << std::endl;
        }
    }
    
    // Close output file
    file.close();
}
