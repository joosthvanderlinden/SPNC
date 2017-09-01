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

#ifndef _intersect_h_included_
#define _intersect_h_included_

// -------------------------------------------------------------------------------- INCLUDES
// -----------------------------------------------------------------------------------------
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_traits_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Polygon_2.h>

#include <math.h>

#include "globals.h"

// --------------------------------------------------------------------------------- CLASSES
// -----------------------------------------------------------------------------------------

// Contains the types and functions to perform an intersection of a facet (triangle) and
// three or more circles (stemming from spherical particles that intersect the facet).
// Main "output" of the class is the void area and solid area of the intersection.
class Intersect
{
private:
    // Types for exact constructions and polygons (to do the intersection analytically).
    // See also http://www.cgal.org/FAQ.html#predicates_vs_constructions.
    typedef glob::EpicKernel::Point_2                               EpicPoint2D;
    typedef CGAL::Exact_predicates_exact_constructions_kernel       EpecKernel;
    typedef EpecKernel::Point_2                                     EpecPoint2D;
    typedef EpecKernel::Circle_2                                    EpecCircle2D;
    typedef EpecKernel::Vector_2                                    EpecVector2D;
    typedef CGAL::Polygon_with_holes_2<EpecKernel>                  PolygonWithHoles2D;
    typedef CGAL::Polygon_set_2<EpecKernel>                         PolygonSet2D;
    typedef CGAL::Polygon_2<EpecKernel>                             Polygon2D;
    typedef CGAL::Cartesian_converter<glob::EpicKernel,EpecKernel>  KernelConverter;
    
    // Defines another vertex container for convenience
    struct IntersectionVertex
    {
        // Radius and centroid of spherical particle
        double          r;
        glob::Point3D   p_3D;
        
        // Centroid of the particle after it's converted to 2D
        EpecPoint2D     p_2D;
        
        // Circle of the particle after it's intersected with the plane that the facet
        // resides in.
        EpecCircle2D    circle;
        
        // Polygon approximation of the circle.
        Polygon2D       polygon;
    };
    
    typedef std::vector<IntersectionVertex>     IntersectionVertexVector;
    typedef IntersectionVertexVector::iterator  IntersectionVertexIter;
    
    // Facet vertices are the three vertices located at the corners of the facet.
    IntersectionVertexVector facet_vertices_;
    
    // Nearby vertices are vertices *not in facet_vertices) whose corresponding particles
    // also intersect the facet.
    IntersectionVertexVector nearby_vertices_;
    
    // Polygon representation of the facet triangle
    Polygon2D facet_polygon_;
    
    // Y-coordinate of the plane that the facet resides in. It is assumed the facet has
    // been rotated to align with the horizontal plane (see transform.cpp).
    double horizontal_plane_height_;
    
    // The void area of the intersection is the area not intersected by the particles
    // corresponding to the vertices in facet_vertices_ and nearby_vertices_.
    double intersection_void_area_;
    double intersection_solid_area_;
    
    // Convenience function to turn a VertexHandle into an IntersectionVertex.
    IntersectionVertex GetIntersectionVertex(glob::VertexHandle vertex);
    
    // Convert from the EPIC (inexact constructions) to EPEC (exact constructions) kernel.
    EpecPoint2D Epic3DToEpec2D(glob::Point3D p);
    
    // Approximate a circle using a polygon
    Polygon2D CircleToPolygon(const EpecCircle2D& circle);
    
    // Analytical intersection of the circles and the triangle, using basic geometry.
    void IntersectAnalytically();
    
    // Semi-analytical intersection of the triangle and the polygon-approximations of the circles.
    void IntersectSemiAnalytically();
    
    // Computes the area of the lens created by intersecting circles c1 and c2
    double SphereIntersectionLensArea(EpecCircle2D c1, EpecCircle2D c2);
    
    // Computes the circles obtained by intersecting the particle corresponding to each vertex
    // with the horizontal plane (defined by horizontal_plane_height_).
    void PlaneIntersectionCircles(IntersectionVertexVector& vertices);
    
    // Calls PlaneIntersectionCircles() but subsequently converts the circles to polygons.
    void PlaneIntersectionPolygons(IntersectionVertexVector& vertices);
    
    // Accepts circle_union by reference, and adds circles in the IntersectionVertexVector
    // to the union.
    void CircleUnion(IntersectionVertexVector& vertices, PolygonSet2D& circle_union);
    
public:
    // Constructor; takes the 3 vertices located at the corners of the facet and sets
    // horizontal_plane_height_, facet_vertices_ and facet_polygon_.
    Intersect(glob::VertexVector& facet_vertices);
    
    // Appends any nearby_vertices to nearby_vertices_ in the IntersectionVertex format.
    void AddNearbyVertices(glob::VertexVector& nearby_vertices);
    
    // Performs the main checks to decide if the intersection can be done analytically (which is
    // much faster) or should be done semi-analytically.
    void FacetIntersection();
    
    double VoidArea();
    double SolidArea();
    
};

#endif