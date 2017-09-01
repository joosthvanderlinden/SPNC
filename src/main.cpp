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

#include "configure.h"
#include "pack.h"
#include "transform.h"
#include "intersect.h"

int main(int argc, char **argv)
{
    // ------------------------------------------------------------- 1. INPUT SETUP
    // ----------------------------------------------------------------------------
    
    // 1.1 Read in the command-line arguments
    Configure config;
    config.ReadArguments(argc, argv);
    double porosity_threshold = config.GetPorosityThreshold();
    
    std::chrono::time_point<std::chrono::system_clock> start;
    std::vector<std::string> path_list = config.GetPathList();
    for(std::vector<std::string>::const_iterator file_it = path_list.begin();
        file_it != path_list.end(); ++file_it)
    {
        // ---------------------------------------------------------- 2. INITIALIZATION
        // ----------------------------------------------------------------------------
        start = std::chrono::system_clock::now();
        std::cout << "Initialising..." << std::flush;
        
        // 2.1 Create packing in which most data is stored
        Pack packing;
        
        // 2.2 Read input data from text file
        packing.ReadPack(*file_it + "/DEM/subsample.txt");
        
        // 2.3 Perform Delaunay triangulation
        packing.InitializeTriangulation();
        
        // 2.4 Compute range tree of the packing
        packing.InitializeRangeTree();
        
        // 2.5 Determine interior cells and facets
        packing.FindInteriorCellsAndFacets();
        
        // 2.6 Initialize disjoint-set data structure
        packing.InitializeDisjointSet();
        
        std::cout << " done." << std::endl;
        // ------------------------------------------------ 3. CONTRUCT CONTACT NETWORK
        // ----------------------------------------------------------------------------

        // 3.1 Add nodes to the contact network
        packing.AddContactNetworkNodes();
        
        // 3.2 Add edges to the contact network, and compute edges weights
        packing.AddContactNetworkEdges();

        // ------------------------------------------------- 4. COMPUTE CELL PARAMETERS
        // ----------------------------------------------------------------------------
        
        // 4.1 Compute the void volume and surface area in each tetrahedron
        packing.CellParameters(config);
        
        // ------------------------------------------------------------------ MAIN LOOP
        // ----------------------------------------------------------------------------
        int num_facets = packing.NumInteriorFacets();
        int counter    = 0;
        for(glob::ConstFacetIter facet_it = packing.InteriorFacetsBegin();
            facet_it != packing.InteriorFacetsEnd(); ++facet_it)
        {
            std::cout << "\r" << "Checking facet merge criterion: "
                      << ++counter << " / " << num_facets << std::flush;
            
            // Retrieve vertices of the facet and vertices nearby the facet
            glob::VertexVector facet_vertices  = packing.GetFacetVertices(*facet_it);
            glob::VertexVector nearby_vertices = packing.FindNearbyVertices(facet_vertices);
            
            // --------------------------------------------------------- 5. ROTATE VERTICES
            // ----------------------------------------------------------------------------
            
            // 5.1 Set the facet plane and horizontal plane to rotate to
            Transform transformation(facet_vertices);
            
            // 5.2 Compute the rotation matrix needed to rotate vertices
            transformation.RotationMatrix();
        
            // 5.3 Rotate the facet vertices and nearby vertices to a horizontal plane
            transformation.Rotate(facet_vertices);
            transformation.Rotate(nearby_vertices);
            
            // --------------------------------------------------------- 6. INTERSECT FACET
            // ----------------------------------------------------------------------------
            
            // 6.1 Initialize the intersection by providing the facet vertices
            Intersect intersection(facet_vertices);
            
            // 6.2 Set vertices to be considered for intersection with the facet
            intersection.AddNearbyVertices(nearby_vertices);
        
            // 6.3 Intersect the circles with the facet
            intersection.FacetIntersection();
            
            // --------------------------------------------------- 7. CHECK MERGE CRITERION
            // ----------------------------------------------------------------------------
            
            // 7.1 Retrieve the resulting void area and solid area
            double solid_area = intersection.SolidArea();
            double void_area  = intersection.VoidArea();
            
            // 7.2 Store the void and solid area for later use
            packing.SetFacetAreas(*facet_it, void_area, solid_area);
            
            // 7.3 If the porosity on the facet is larger than the threshold,
            //     merge the two adjacent cells in the distjoint-set
            if (void_area / (void_area + solid_area) > porosity_threshold)
            {
                packing.MergeFacetCells(*facet_it);
            }
        }
        
        std::cout << " done." << std::endl;
        // -------------------------------------------------- 8. CONSTRUCT PORE NETWORK
        // ----------------------------------------------------------------------------

        // 8.1 Use the disjoint-set to define the pore network nodes
        packing.AddPoreNetworkNodes();

        // 8.2 Add pore network edges and compute edge weights
        packing.AddPoreNetworkEdges();
        
        std::chrono::duration<double> duration = (std::chrono::system_clock::now() - start);
        std::cout << "Total construction time: " << duration.count() << " seconds." << std::endl;
        // --------------------------------------------------- 9. WRITE RESULTS TO FILE
        // ----------------------------------------------------------------------------
        std::cout << "Writing output files..." << std::flush;
        
        // 9.1 Write Delaunay triangulation (and pores) to output file
        packing.WriteTriangulation(*file_it);
        
        // 9.2 Write pore network to output file (txt and vtk)
        packing.WritePoreNetwork(*file_it);

        // 9.3 Write contact network to output file
        packing.WriteContactNetwork(*file_it);
        
        // 9.4 Write cell properties to output file
        packing.WriteCellProperties(*file_it);
        
        // 9.5 Write facet properties to output file
        packing.WriteFacetProperties(*file_it);
        
        std::cout << " done." << std::endl;
    }
    
    return 0;
}



