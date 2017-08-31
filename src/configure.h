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

#ifndef _configure_h_included_
#define _configure_h_included_

// -------------------------------------------------------------------------------- INCLUDES
// -----------------------------------------------------------------------------------------
#include <fstream>
#include <iostream>
#include <vector>
#include <getopt.h>

// --------------------------------------------------------------------------------- CLASSES
// -----------------------------------------------------------------------------------------

// This class is
class Configure
{
private:
    // Path for list of directories to be used in network construction
    std::string path_;
    
    struct PoreParameters
    {
        // Are particles assumed to overlap? If so, Pack::CellParameters() will be faster.
        bool no_particle_overlap = false;
        
        // Void facet area/total facet area, i.e. threshold above which two cells are
        // merged in PoreNet (see step 7 in main.cpp). 0.4 is generally recommended.
        double porosity_threshold = 0.4;
        
    } pore_parameters_;
    
    // Converts a string to double while checking for invalid input
    double StringToDouble(char s[]);
    
    // Converts a string to bool while checking for invalid input
    bool StringToBool(std::string s);
    
    // Checks if valid path is provided
    std::string CheckPath(char s[]);
    
    // Print out the help information
    void ShowHelp(char *s);
    
public:
    
    // Reads in the command-line arguments using getopt
    void ReadArguments(int argc, char **argv);
    
    // Returns the path list. If path_ is not set, returns single entry ("/data")
    std::vector<std::string> GetPathList();
    
    // Getters for the parameters
    double GetPorosityThreshold();
    bool   GetNoParticleOverlap();
};


#endif