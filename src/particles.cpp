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

#include "particles.h"

double Particles::Distance(glob::Point3D& p, glob::Point3D& c)
{
    return std::sqrt(std::pow(p.x() - c.x(), 2) +
                     std::pow(p.y() - c.y(), 2) +
                     std::pow(p.z() - c.z(), 2));
}

int Particles::NumParticles()
{
    return dem_data_.num_particles;
}

double Particles::ParticleSizeMean()
{
    return dem_data_.particle_size_mean;
}

double Particles::ParticleSizeRange()
{
    return dem_data_.particle_size_range;
}

double Particles::MaxParticleRadius()
{
    return dem_data_.particle_size_mean + dem_data_.particle_size_range;
}

double Particles::ParticleVoxelSize(double resolution)
{
    return dem_data_.particle_size_mean / resolution;
}

bool Particles::InsideLimits(glob::Point3D p)
{
    return (p.x() > limits_.x_min && p.x() < limits_.x_max &&
            p.y() > limits_.y_min && p.y() < limits_.y_max &&
            p.z() > limits_.z_min && p.z() < limits_.z_max );
}

bool Particles::OutsideLimits(glob::Point3D p, double padding)
{
    return (p.x() - padding < limits_.x_min || p.x() + padding > limits_.x_max ||
            p.y() - padding < limits_.y_min || p.y() + padding > limits_.y_max ||
            p.z() - padding < limits_.z_min || p.z() + padding > limits_.z_max );
}

std::string Particles::VerticalLocation(double y, double padding)
{
    if (y + padding > limits_.y_max)
    {
        return "top";
    }
    else if (y - padding < limits_.y_min)
    {
        return "bottom";
    }
    else
    {
        return "middle";
    }
}

std::string Particles::HorizontalLocation(double x, double z, double padding)
{
    if (x + padding > limits_.x_max || z + padding > limits_.z_max ||
        x - padding < limits_.x_min || z - padding < limits_.z_min)
    {
        return "outside";
    }
    else
    {
        return "inside";
    }
}


glob::ParticleVectorIter Particles::ParticlesBegin()
{
    return particles_.begin();
}

glob::ParticleVectorIter Particles::ParticlesEnd()
{
    return particles_.end();
}


double Particles::SphereVolume(glob::VertexHandle v)
{
    // Check if sphere exceeds limits, and if so, compute volume with voxels.
    if (OutsideLimits(v->point(), v->info().r))
    {
        return BoxParticleIntersectionVolume(v);
    }
    else
    {
        return (4/3.) * CGAL_PI * pow(v->info().r,3);
    }
}

double Particles::SphereArea(glob::VertexHandle v)
{
    // Limitation: does not account for spheres exceeding the subsample limits
    // (this is accounted for in sphere volume but not in sphere area)
    return 4 * CGAL_PI * pow(v->info().r,2);
}

double Particles::BoxParticleIntersectionVolume(glob::VertexHandle v)
{
    // Assuming 20 voxels per particle diameter
    double Dx = ParticleVoxelSize(10);
    double Dy = Dx;
    double Dz = Dx;
    double voxel_volume = Dx*Dy*Dz;
    
    double r = v->info().r;
    int Nx = round( 2 * r / Dx);
    int Ny = round( 2 * r / Dy);
    int Nz = round( 2 * r / Dz);
    
    // Using variable length arrays, even though it's not recommended
    // (speed gain over using std::vectors is significant)
    double* x = new double[Nx];
    for(int i=0; i<Nx; ++i)
    {
        x[i] = v->point().x() - r + Dx*i;
    }
    
    double* y = new double[Ny];
    for(int i=0; i<Ny; ++i)
    {
        y[i] = v->point().y() - r + Dy*i;
    }
    
    double* z = new double[Nz];
    for(int i=0; i<Nz; ++i)
    {
        z[i] = v->point().z() - r + Dz*i;
    }
    
    double solid_volume = 0.0;
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k=0; k<Nz; ++k)
            {
                // query_point is the center of the voxel
                glob::Point3D query_point = glob::Point3D(x[i] + Dx/2, y[j] + Dy/2, z[k] + Dz/2);
                
                // Test if query point is inside sphere and
                // test if query point is inside sample
                if (WithinSphere(query_point, v->point(), r) &&
                    InsideLimits(query_point))
                {
                    solid_volume += voxel_volume;
                }
            }
        }
    }
    
    delete [] x;
    delete [] y;
    delete [] z;
    
    return solid_volume;
}


bool Particles::WithinSphere(glob::Point3D& p, glob::Point3D& c, double r)
{
    return Distance(p,c) <= r;
}

double Particles::ParticleGapWidth(glob::VertexHandle& u,
                                 glob::VertexHandle& v)
{
    return Distance(u->point(), v->point()) - u->info().r - v->info().r;
}

double Particles::SphereContactRadius(glob::VertexHandle& u,
                                      glob::VertexHandle& v)
{
    // Compute distance from u to middlepoint of intersection
    double x = u->info().r - SphereCapHeight(u, v);

    // Compute area of circle (pi * radius^2 of intersecting circle)
    return sqrt(pow(u->info().r,2) - pow(x,2));
}

double Particles::SphereCapHeight(glob::VertexHandle& u,
                                   glob::VertexHandle& v)
{
    // Compute distance between u and v
    double distance = CGAL::sqrt(CGAL::squared_distance(u->point(), v->point()));
    
    // Compute distance to middle point (from u) with Pythagoras
    double x        = 0.5 * (distance - (std::pow(v->info().r,2) -
                                         std::pow(u->info().r,2)) / distance);
    
    return u->info().r - x;
}


double Particles::SphereCapArea(double r, double h)
{
    return 2 * CGAL_PI * r * h;
}


double Particles::SphereCapVolume(double r, double h)
{
    return CGAL_PI * std::pow(h,2) * (r - h / 3.0);
}


void Particles::ReadParticles(std::ifstream& file)
{
    particles_.resize(dem_data_.num_particles);
    
    int idx;
    double x,y,z,r;
    for(int i = 0; i < dem_data_.num_particles; ++i)
    {
        file >> idx >> x >> y >> z >> r;
        glob::VertexInfo info;
        info.r           = r;
        info.original_id = idx;
        particles_[i]    = std::make_pair(glob::Point3D(x,y,z), info);
    }
}


void Particles::ReadData(std::ifstream& file)
{
    file >> dem_data_.num_particles
         >> dem_data_.particle_size_mean
         >> dem_data_.particle_size_range
         >> dem_data_.stress_strain_mask
         >> dem_data_.stress_strain_rate
         >> dem_data_.young_modulus
         >> dem_data_.density
         >> dem_data_.poisson_ratio
         >> dem_data_.friction_angle
         >> dem_data_.damping_coef;
}

void Particles::ReadLimits(std::ifstream& file)
{
    file >> limits_.x_min >> limits_.x_max
         >> limits_.y_min >> limits_.y_max
         >> limits_.z_min >> limits_.z_max;
}

void Particles::WriteCoordinates(std::ofstream& file, std::string c)
{
    if (c == "x")
    {
        for(int i = 0; i < dem_data_.num_particles; ++i)
        {
            file << particles_[i].first.x() << std::endl;
        }
    }
    else if (c == "y")
    {
        for(int i = 0; i < dem_data_.num_particles; ++i)
        {
            file << particles_[i].first.y() << std::endl;
        }
    }
    else if (c == "z")
    {
        for(int i = 0; i < dem_data_.num_particles; ++i)
        {
            file << particles_[i].first.z() << std::endl;
        }
    }
    else
    {
        for(int i = 0; i < dem_data_.num_particles; ++i)
        {
            file << particles_[i].first.x() << " "
                 << particles_[i].first.y() << " "
                 << particles_[i].first.z() << std::endl;
        }
    }
}

void Particles::WriteRadii(std::ofstream& file)
{
    for(int i = 0; i < dem_data_.num_particles; ++i)
    {
        file << particles_[i].second.r << std::endl;
    }
}