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

#ifndef _transform_h_included_
#define _transform_h_included_

// -------------------------------------------------------------------------------- INCLUDES
// -----------------------------------------------------------------------------------------
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Object.h>

#include "globals.h"

// --------------------------------------------------------------------------------- CLASSES
// -----------------------------------------------------------------------------------------

// Class to perform a rotation of a collection of vertices. This class is used to take the
// 3 vertices located at the facet, and any vertices whose corresponding particles may 
// intersect the facet, and rotate all of them such that the three facet vertices align with
// a horizontal plane.
class Transform
{
private:
    typedef CGAL::Aff_transformation_3<glob::EpicKernel> Transform3D;
    typedef glob::EpicKernel::Line_3                     Line3D;
    typedef glob::EpicKernel::Vector_3                   Vector3D;
    
    // The plane that the three facet vertices reside in.
    glob::Plane3D facet_plane_;

    // A simple horizontal plane defined by the points (1,0,0), (0,0,1) and (1,0,1).
    glob::Plane3D horizontal_plane_;

    // The 3 by 3 matrix to rotate a vector with (x,y,z) coordinates.
    Transform3D rotation_matrix_;

    // Normalizes a vector v by dividing by its length.
    Vector3D normalize(const Vector3D v);

    // Computes the axis of rotation between the facet plane and the horizontal plane.
    Vector3D RotationAxis();

    // Computes angle of rotation between the facet plane and the horizontal plane.
    double RotationAngle();
    
public:

	// Initializes the transform class with the three facet vertices, setting facet_plane_
	// and horizontal_plane_.
    Transform(glob::VertexVector& facet_vertices);

    // Computes the rotation_matrix_.
    void RotationMatrix();

    // Perform the rotation, by multiplying every vertex in vertices with the rotation matrix.
    // Assumes vertices has a parameter (in this case, info()->rotated_point) to store the
    // restult in. If the coordinates of the vertex itself are rotated, the triangulation will
    // be altered, which is undesirable.
    void Rotate(glob::VertexVector& vertices);
};

#endif