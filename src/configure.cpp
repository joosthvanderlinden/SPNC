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

double Configure::StringToDouble(char s[])
{
    char* end;
    double value = strtod(s, &end);
    if (end == s || *end != '\0')
    {
        std::printf("\nInvalid double argument: %s\n", s);
        exit(EXIT_FAILURE);
    }
    else
    {
        return value;
    }
}

bool Configure::StringToBool(std::string s)
{
    
    if (s.size() != 1 || s[0] < '0' || s[0] > '1')
    {
        std::printf("\nInvalid boolean argument: %s\n", s.c_str());
        std::printf("For booleans, use 1 (true) or 0 (false)\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        return s[0] == '1';
    }
}

std::string Configure::CheckPath(char s[])
{
    std::ifstream infile(s);
    if (infile.good())
    {
        std::printf("\nReading in file list: %s\n", s);
        return s;
    }
    else
    {
        std::printf("\nInvalid path: %s\n", s);
        exit(EXIT_FAILURE);
    }
}

void Configure::ShowHelp(char *s)
{
    std::cout << "\n\nUsage: " << s << " [-option] [argument]" << std::endl << std::endl;
    std::cout << "Options: " << std::endl;
    std::cout << "-h [--help]                  Show help information" << std::endl;
    std::cout << "-f [--file_path]             Path to file with a directory list" << std::endl;
    std::cout << "-p [--porosity_threshold]    Areal porosity merging threshold" << std::endl;
    std::cout << "-n [--no_particle_overlap]   Are particles assumed to overlap? (0 = yes, 1 = no)" << std::endl;

    std::cout<<"Example: " << s << " -f ../docs/example/dir_list.txt -p 0.4 -n 0" << std::endl << std::endl;
    
    exit(EXIT_SUCCESS);
}

void Configure::ReadArguments(int argc, char **argv)
{
    const struct option longopts[] =
    {
        {"help",                no_argument,        0, 'h'},
        {"file_path",           required_argument,  0, 'f'},
        {"porosity_threshold",  required_argument,  0, 'p'},
        {"no_particle_overlap", required_argument,  0, 'n'},
        {0,0,0,0},
    };
    
    opterr = 1; // turns off getopt error message
    int index;
    int iarg = 0;
    double v;
    while(iarg != -1)
    {
        iarg = getopt_long(argc, argv, "hf:p:n:", longopts, &index);

        switch (iarg)
        {
            case 'h':
                ShowHelp(argv[0]);
                break;
            case 'f':
                path_ = CheckPath(optarg);
                break;
            case 'p':
                v = StringToDouble(optarg);
                if ((v < 0) || (v > 1))
                {
                    std::fprintf(stderr, "porosity_threshold should be between 0.0 and 1.0.\n");
                    exit(EXIT_FAILURE);
                }
                pore_parameters_.porosity_threshold = v;
                break;
            case 'n':
                // Value for this parameter is checked in StringToBool() for validity
                pore_parameters_.no_particle_overlap = StringToBool(optarg);
                break;
            case ':':
                // missing option argument
                std::fprintf(stderr, "%s: option '-%c' requires an argument.\n",
                             argv[0], optopt);
                exit(EXIT_FAILURE);
            case '?':
                // invalid option
                std::fprintf(stderr, "%s: option '-%c' is invalid.\n",
                             argv[0], optopt);
                exit(EXIT_FAILURE);
        }
    }
}

std::vector<std::string> Configure::GetPathList()
{
    std::vector<std::string> path_list;
    if (path_.empty())
    {
        // By default, use example directory
        std::cout << "Using the assembly in docs/example/assembly_1." << std::endl;
        path_list.push_back("../docs/example/assembly_1");
    }
    else
    {
        std::ifstream file(path_);
        std::string   loc;
        
        while (std::getline(file, loc))
        {
            path_list.push_back(loc);
        }
    }
    return path_list;
}

double Configure::GetPorosityThreshold()
{
    return pore_parameters_.porosity_threshold;
}

bool Configure::GetNoParticleOverlap()
{
    return pore_parameters_.no_particle_overlap;
}

