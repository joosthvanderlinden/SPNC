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

#include "transform.h"

Transform::Transform(glob::VertexVector& facet_vertices)
{
    facet_plane_      = glob::Plane3D(facet_vertices[0]->point(),
                                      facet_vertices[1]->point(),
                                      facet_vertices[2]->point());
    
    horizontal_plane_ = glob::Plane3D(glob::Point3D(1,0,0),
                                      glob::Point3D(0,0,1),
                                      glob::Point3D(1,0,1));
}


Transform::Vector3D Transform::normalize(const Vector3D v)
{
    double len = sqrt(v.squared_length());
    
    if (len == 0.0)
        return v;
    else
        return v / len;
}


Transform::Vector3D Transform::RotationAxis()
{
    // intersection() is a CGAL function and returns a CGAL object
    CGAL::Object obj = intersection(facet_plane_,horizontal_plane_);
    
    // If succesful, return rotation axis (as a vector)
    // If unsuccesful, then planes are parallel
    Line3D line1;
    if (assign(line1, obj))
       return normalize(line1.to_vector());
    else
       return Vector3D(1.0,0.0,0.0);
    
}


double Transform::RotationAngle()
{
    Vector3D facet_ortho      = facet_plane_.orthogonal_vector();
    Vector3D horizontal_ortho = horizontal_plane_.orthogonal_vector();
    
    double dot_product = facet_ortho * horizontal_ortho;
    double lengths     = std::sqrt(facet_ortho.squared_length() * horizontal_ortho.squared_length());
    
    if (lengths)
        return std::acos(dot_product/lengths);
    else
        return 0.0;
}


void Transform::RotationMatrix()
{
	// Derived from https://www.fastgraph.com/makegames/3drotation/.
    Vector3D axis = RotationAxis();
    double angle  = RotationAngle();
    
    double c = std::cos(angle);
    double s = std::sin(angle);
    double t = 1.0 - c;
    
    double m00 = c + axis.x() * axis.x() * t;
    double m11 = c + axis.y() * axis.y() * t;
    double m22 = c + axis.z() * axis.z() * t;
    
    double m10 = axis.x() * axis.y() * t + axis.z() * s;
    double m01 = axis.x() * axis.y() * t - axis.z() * s;
    
    double m20 = axis.x() * axis.z() * t - axis.y() * s;
    double m02 = axis.x() * axis.z() * t + axis.y() * s;
    
    double m21 = axis.y() * axis.z() * t + axis.x() * s;
    double m12 = axis.y() * axis.z() * t - axis.x() * s;
    
    rotation_matrix_ = Transform3D(m00, m01, m02, m10, m11, m12, m20, m21, m22);
}


void Transform::Rotate(glob::VertexVector& vertices)
{
    // Using CGAL's transform() function
    for(glob::VertexVector::iterator vertex_it = vertices.begin();
        vertex_it != vertices.end(); ++vertex_it)
    {
        (*vertex_it)->info().rotated_point = (*vertex_it)->point().transform(rotation_matrix_);
    }
}


