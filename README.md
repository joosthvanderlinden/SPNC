# SPNC - Sphere Packing Network Construction
## Introduction
SPNC is a program to construct network from an assembly of spheres. The pore network represents pores (as nodes) and throats (as edges). For the contact network, particle centroids form the nodes and an edge is assigned if two particles are in contact. As such, constructing the contact network is relatively straightforward. For the construction of the pore network, SPNC uses the Delaunay triangulation. Because individual tetrahedra do not necessarily encapsulate pores, a merging criterion is adopted.

For more details, please read the MANUAL and the source code. The code is commented throughout and made deliberately verbose in some sections (e.g. main.cpp) to walk the reader through the steps of the network construction. If you use this program, or parts of it, please cite our publication [1].

[1] J.H. van der Linden, A. Sufian, G.A. Narsilio, A.R. Russell, A. Tordesillas, "A computational geometry approach to pore network construction for granular packings," Computers & Geosciences, vol. 112, pp. 133-143, 2018

## Requirements & Compilation
SPNC is implemented in C++(11) and relies heavily on the computational geometry library `CGAL`. In turn, `CGAL` requires `Boost`. To compile the code with `CMake`, an example `CMakeLists.txt` file is provided. Assuming the right compiler is available and `CMake`, `CGAL` and `Boost` are installed, compile SPNC on OSX/Unix as follows:

1. Create a directory `build` in the same directory as `src`
2. `cd build`
3. `cmake ../src`
4. `make`

For more information on building SPNC with CGAL, consult section 14 of the [CGAL installation manual](http://doc.cgal.org/latest/Manual/installation.html). For more information on how to use CGAL with Xcode, consult [this link](https://3d.bk.tudelft.nl/ken/en/2016/03/16/using-cgal-and-xcode.html). 

## Usage
Once compiled, SPNC can be run from the command line as
```shell
Usage: ./spnc [-option] [argument]

Options:
-h [--help]                  Show help information
-f [--file_path]             Path to file with a directory list
-p [--porosity_threshold]    Areal porosity merging threshold
-n [--no_particle_overlap]   Are particles assumed to overlap? (0 = yes, 1 = no)
Example: ./spnc -f ../docs/example/dir_list.txt -p 0.4 -n 0

```
An example sphere assembly is provided in /docs/example. See MANUAL for further details.

## Potential improvements
This software was written as part of a PhD thesis. The implementation can potentially be improved in many ways, such as:
- Use inheritance and implement a base class for general networks, and put the pore- and contact-network in child classes.
- Replace the window query with a search of neighbouring (and neighbour's neighbouring) tetrahedra, for further speedup.
- Integrate the disjoint-set in the network class, to avoid the (expensive) `Pack::AddPoreNetworkNodes()` function.
- More extensive exception handling and unit tests.

## Authors
Joost van der Linden <joosthvanderlinden@gmail.com>

The author gratefully acknowledges G.A. Narsilio, A. Tordesillas, A. Sufian, A.R. Russell and S.K. Matthai, for their feedback and support in the development of SPNC.

## License
This project is licensed under the GPL License - see the LICENSE.txt file for details.