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

#ifndef _particles_h_included_
#define _particles_h_included_

// -------------------------------------------------------------------------------- INCLUDES
// -----------------------------------------------------------------------------------------
#include <string>
#include <fstream>

#include "globals.h"

// --------------------------------------------------------------------------------- CLASSES
// -----------------------------------------------------------------------------------------

// Contains the assembly data, including particle positions and radii, limits of the
// assembly and simulation parameters, as well as functions related to the particles.
class Particles
{
private:
    // Stores data from the discrete-element method (DEM) simulation that generated the 
    // assembly. Only num_particles, particle_size_mean and particle_size_range are used
    // in the rest of the code, the rest is just stored in case the user wants to write this
    // data out again (along with the pore network or contact network, for example).
    struct DEMData
    {
        int         num_particles;

        // Mask used in Yade to determine the stress and strain settings. For more
        // information on Yade, see https://yade-dem.org/.
        int         stress_strain_mask;

        // Mean radius of the spherical particles
        double      particle_size_mean;

        // Range of the particle sizes, defined as the difference between the mean radius
        // and the minimum particle radius (assuming this is equal to the difference
        // between the mean radius and the maximum particle radius).
        double      particle_size_range;

        // Amount of stress or strain applied to the assembly in the DEM simulation.
        double      stress_strain_rate;

        // Young modulus of the particles in the DEM simulation.
        double      young_modulus;

        // Density of the particles in the DEM simulation.
        double      density;

        // Poisson ratio of the particles in the DEM simulation.
        double      poisson_ratio;

        // Friction angle of the particles in the DEM simulation.
        double      friction_angle;

        // Damping coefficient in the DEM simulation.
        double      damping_coef;
    } dem_data_;
    
    // It is assumed that the assembly is subsampled from a larger sample, to get rid of
    // wall effects etc. The limits of that subsample are stored here.
    struct Limits
    {
        double      x_min;
        double      x_max;
        double      y_min;
        double      y_max;
        double      z_min;
        double      z_max;
    } limits_;
    
    // This vector stores the particles as they are read in (coordinates + radii), in a
    // format that can be directly used in the triangulation initialization.
    glob::ParticleVector particles_;
    
    // Calculates the distance between point p and point c
    double Distance(glob::Point3D& p, glob::Point3D& c);
    
public:
    glob::ParticleVectorIter ParticlesBegin();
    glob::ParticleVectorIter ParticlesEnd();

    int NumParticles();
    double ParticleSizeMean();
    double ParticleSizeRange();
    double MaxParticleRadius();

    // Returns the mean particle radius divided by resolution. This function allows the
    // voxel size to be based on the particle radius, i.e. to have N voxels / particle.
    double ParticleVoxelSize(double resolution);

    // Tests whether the point p is inside or outside the limits given in limits_. An
    // additional padding can be added for OutsideLimits().
    bool InsideLimits(glob::Point3D p);
    bool OutsideLimits(glob::Point3D p, double padding);

    // Returns a string indicating the location of the particle, based on it's y coordinate. 
    // Can be "top" (above limits_), "bottom" (below limits_) or "middle" (anything else).
    // An additional padding can be added.
    std::string VerticalLocation(double y, double padding);
    
    // Returns a string indicating the location of the particle, based on it's x and
    // z coordinates. Can be "outside" (outside x and z limits_) or "inside" (anything
    // else). An additional padding can be added.
    std::string HorizontalLocation(double x, double z, double padding);
    
    // Approximates the volume of the intersection between a particle (defined by the 
    // VertexHandle) and the sample box, as given by the limits_, using voxels.
    double BoxParticleIntersectionVolume(glob::VertexHandle v);

    // Tests whether or not the point p is inside the sphere with centroid c and radius r.
    bool WithinSphere(glob::Point3D& p, glob::Point3D& c, double r);
    
    // Calculates the width of the gap between the particles at vertices u and v.
    double ParticleGapWidth(glob::VertexHandle& u, glob::VertexHandle& v);
    
    // Computes the radius of the cirulcar contact of two particles in contact, identified
    // by their respective VertexHandle's.
    double SphereContactRadius(glob::VertexHandle& u, glob::VertexHandle& v);
    
    // Computes the middle point of intersection of two overlapping particles. See also
    // http://stackoverflow.com/questions/5048701/finding-points-of-intersection-when-two-spheres-intersect
    double SphereCapHeight(glob::VertexHandle& u, glob::VertexHandle& v);
    
    // Computes the area of the spherical cap (with height h and radius r) created by the
    // intersection of two spherical particles.
    double SphereCapArea(double r, double h);
    
    // Computes the volume of the spherical cap (with height h and radius r) created by the
    // intersection of two spherical particles. The volume of the cap is equal to 1/2 of the
    // volume of the "lens" created by the intersection.
    double SphereCapVolume(double r, double h);
    
    // Computes the volume and area of the sphere at vertex v.
    double SphereVolume(glob::VertexHandle v);
    double SphereArea(glob::VertexHandle v);
    
    // Reads or writes the data, as stored in the private structs above.
    void ReadParticles(std::ifstream& file);
    void ReadData(std::ifstream& file);
    void ReadLimits(std::ifstream& file);
    void WriteCoordinates(std::ofstream& file, std::string c);
    void WriteRadii(std::ofstream& file);
};

#endif