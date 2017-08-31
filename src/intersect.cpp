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

#include "intersect.h"

Intersect::Intersect(glob::VertexVector& facet_vertices)
{
    // y-coordinate is the same for all rotated points
    horizontal_plane_height_ = facet_vertices[0]->info().rotated_point.y();
    
    facet_vertices_.resize(3);
    for(int i = 0; i < 3; ++i)
    {
        facet_polygon_.push_back(Epic3DToEpec2D(facet_vertices[i]->info().rotated_point));
        facet_vertices_[i] = GetIntersectionVertex(facet_vertices[i]);
    }
    
    // Check orientation; if it's not counter-clockwise, reverse
    if (facet_polygon_.orientation() == CGAL::CLOCKWISE)
        facet_polygon_.reverse_orientation();
}

Intersect::EpecPoint2D Intersect::Epic3DToEpec2D(glob::Point3D p)
{
    // To do the intersection, we need the epec kernel instead of
    // epic kernel (epec = exact predicates, exact constructions)
    // See: http://www.cgal.org/FAQ.html#predicates_vs_constructions
    KernelConverter from_epic_to_epec;
    
    // Convert 3D epic point to 2D epic point
    EpicPoint2D epic_point(p.x(), p.z());
    
    // Convert 2D epic point to 2D epec (CGAL function)
    return from_epic_to_epec(epic_point);
}

void Intersect::AddNearbyVertices(glob::VertexVector& nearby_vertices)
{
    for (glob::ConstVertexIter vertex_it = nearby_vertices.begin();
         vertex_it != nearby_vertices.end(); ++vertex_it)
    {
        IntersectionVertex v = GetIntersectionVertex(*vertex_it);
        nearby_vertices_.push_back(v);
    }
}

Intersect::IntersectionVertex Intersect::GetIntersectionVertex(glob::VertexHandle vertex)
{
    IntersectionVertex v;
    v.r     = vertex->info().r;
    v.p_3D  = vertex->info().rotated_point;
    v.p_2D  = Epic3DToEpec2D(v.p_3D);
    
    return v;
}


void Intersect::PlaneIntersectionCircles(IntersectionVertexVector& vertices)
{
    for(IntersectionVertexIter vertex_it = vertices.begin();
        vertex_it != vertices.end(); ++vertex_it)
    {
        // Compute distance from vertex to horizontal plane
        // If vertex_it is a facet vertex, then vertex_height = 0
        double vertex_height = vertex_it->p_3D.y() - horizontal_plane_height_;
        double sphere_radius = vertex_it->r;
        
        // Compute the radius of circle that is the result of intersecting
        // the sphere and the horizontal plane using Pythagoras
        double squared_circle_radius = std::pow(sphere_radius,2) -
                                       std::pow(vertex_height,2);
        
        vertex_it->circle = EpecCircle2D(vertex_it->p_2D, squared_circle_radius);
    }
}

void Intersect::PlaneIntersectionPolygons(IntersectionVertexVector& vertices)
{
    PlaneIntersectionCircles(vertices);
    
    for(IntersectionVertexIter vertex_it = vertices.begin();
        vertex_it != vertices.end(); ++vertex_it)
    {
        vertex_it->polygon = CircleToPolygon(vertex_it->circle);
    }
}

Intersect::Polygon2D Intersect::CircleToPolygon(const EpecCircle2D& circle)
{
    Polygon2D           polygon;
    
    // Approximating a circle with 64-sided polygon
    double              num_points = 64;
    
    double              r   = std::sqrt(CGAL::to_double(circle.squared_radius()));
    double              c_x = CGAL::to_double(circle.center().x());
    double              c_y = CGAL::to_double(circle.center().y());
    
    for(int i=0; i<num_points; ++i)
    {
        // Rotating counter-clockwise across the circle
        double x = c_x + r * cos( (i / num_points) * 2 * M_PI );
        double y = c_y + r * sin( (i / num_points) * 2 * M_PI );
        
        polygon.push_back(EpecPoint2D(x,y));
    }
    
    return polygon;
}


void Intersect::FacetIntersection()
{
    // If nearby_vertices_ is empty then there are no particles intersecting the facet
    // other than the three particles located at the vertices of the facet itself.
    // In this case, compute the intersection analytically.
    if (nearby_vertices_.empty())
    {
        PlaneIntersectionCircles(facet_vertices_);
        IntersectAnalytically();
    }
    // Else, use a polygon approximation of the intersecting circles and compute the
    // intersection semi-analytically. This function accounts for most of the
    // execution time of the program.
    else
    {
        PlaneIntersectionPolygons(facet_vertices_);
        PlaneIntersectionPolygons(nearby_vertices_);
        
        IntersectSemiAnalytically();
    }
}


void Intersect::IntersectAnalytically()
{
    intersection_solid_area_ = 0.0;
    for(int i = 0; i < 3; ++i)
    {
        EpecPoint2D p1 = facet_vertices_[(i+1)%3].circle.center();
        EpecPoint2D p2 = facet_vertices_[(i+2)%3].circle.center();
        EpecPoint2D p3 = facet_vertices_[(i+3)%3].circle.center();
        
        EpecVector2D v1 = p1-p2;
        v1 = v1 / std::sqrt(CGAL::to_double(v1 * v1));
        EpecVector2D v2 = p3-p2;
        v2 = v2 / std::sqrt(CGAL::to_double(v2 * v2));
        
        double angle = std::acos(CGAL::to_double(v1 * v2)) * 180 / CGAL_PI;
        
        // solid circle section area = angle / 360 * circle area
        intersection_solid_area_ += (angle / 360) * CGAL_PI * CGAL::to_double(facet_vertices_[(i+2)%3].circle.squared_radius());
        
        // Reduce the solid area by 1/2 of the area of the lens created by the overlap between the circles at p1 and p2,
        // otherwise this area is double counted.
        intersection_solid_area_ -= 0.5 * SphereIntersectionLensArea(facet_vertices_[(i+1)%3].circle,
                                                                     facet_vertices_[(i+2)%3].circle);
    }
    double facet_area = CGAL::to_double(facet_polygon_.area());
    
    // If particles overlap a lot, then the solid area may be larger than the facet area.
    // In that case, assume the void area is 0.
    if (intersection_solid_area_ > facet_area)
        intersection_solid_area_ = facet_area;
    
    // Just in case
    if (intersection_solid_area_ < 0)
        intersection_solid_area_ = 0.0;
    
    intersection_void_area_ = facet_area - intersection_solid_area_;
}



double Intersect::SphereIntersectionLensArea(EpecCircle2D c1, EpecCircle2D c2)
{
    // Adapted from Scott Marley's answer at
    // http://stackoverflow.com/questions/4247889/area-of-intersection-between-two-circles
    
    // Using to_double to get double, so not an exact answer but good enough
    double c1_r_squared = CGAL::to_double(c1.squared_radius());
    double c1_r         = CGAL::sqrt(c1_r_squared);
    double c2_r_squared = CGAL::to_double(c2.squared_radius());
    double c2_r         = CGAL::sqrt(c2_r_squared);
    double distance     = CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(c1.center(), c2.center())));
    
    // Circles do not overlap
    if (distance > c1_r + c2_r)
    {
        return 0;
    }
    
    // Circle c2 is completely inside circle c1
    else if (distance <= std::abs(c1_r - c2_r) && c1_r >= c2_r)
    {
        // Return area of circle c2
        return M_PI * c2_r_squared;
    }
    
    // Circle c1 is completely inside circle c2
    else if (distance <= std::abs(c1_r - c2_r) && c1_r < c2_r)
    {
        // Return area of circle0
        return M_PI * c1_r_squared;
    }
    
    // Circles partially overlap and create a "lens"
    else
    {
        double phi   = (std::acos((c1_r_squared + (distance * distance) - c2_r_squared) / (2 * c1_r * distance))) * 2;
        double theta = (std::acos((c2_r_squared + (distance * distance) - c1_r_squared) / (2 * c2_r * distance))) * 2;
        double area1 = 0.5 * theta * c2_r_squared - 0.5 * c2_r_squared * std::sin(theta);
        double area2 = 0.5 * phi * c1_r_squared - 0.5 * c1_r_squared * std::sin(phi);
        
        return area1 + area2;
    }
}


void Intersect::IntersectSemiAnalytically()
{
    // Take the union of all intersecting polygon circles. If no polygons in nearby_
    // vertices intersect the facet, circle_union will be empty and we can resort to
    // the analytical solution and be done.
    PolygonSet2D circle_union;
    CircleUnion(nearby_vertices_, circle_union);
    
    if (circle_union.is_empty())
    {
        IntersectAnalytically(); return;
    }
    // Otherwise, also add the facet vertices to the union and continue with the
    // semi-analytical calculation.
    else
    {
        CircleUnion(facet_vertices_, circle_union);
    }
    
    // Intersect the union of circle polygons with the facet with a CGAL function
    circle_union.intersection(facet_polygon_);
    
    // Retrieve the resulting polygon (which may have holes)
    std::vector<PolygonWithHoles2D> intersection_polygon;
    circle_union.polygons_with_holes(back_inserter(intersection_polygon));
    
    // Unclear if there are any scenarios in which the intersection would consist of more
    // than one polygon, but using an iterator just in case.
    intersection_solid_area_ = 0.0;
    for(std::vector<PolygonWithHoles2D>::const_iterator poly_it = intersection_polygon.begin();
        poly_it != intersection_polygon.end(); ++poly_it)
    {
        intersection_solid_area_ += std::abs(CGAL::to_double(poly_it->outer_boundary().area()));
        
        // Subtract the area of any holes (if present)
        for(PolygonWithHoles2D::Hole_const_iterator hole_it = poly_it->holes_begin();
            hole_it != poly_it->holes_end(); ++hole_it)
        {
            intersection_solid_area_ -= std::abs(CGAL::to_double(hole_it->area()));
        }
    }
    
    double facet_area = CGAL::to_double(facet_polygon_.area());
    intersection_void_area_ = facet_area - intersection_solid_area_;
}

void Intersect::CircleUnion(IntersectionVertexVector& vertices, PolygonSet2D& circle_union)
{
    // Add a circle polygon of a vertex to the union only
    // if it intersects the facet
    for(IntersectionVertexIter vertex_it = vertices.begin();
        vertex_it != vertices.end(); ++vertex_it)
    {
        if ((CGAL::do_intersect(vertex_it->polygon, facet_polygon_)))
        {
            circle_union.join(vertex_it->polygon);
        }
    }
}

double Intersect::SolidArea()
{
    return intersection_solid_area_;
}

double Intersect::VoidArea()
{
    return intersection_void_area_;
}


