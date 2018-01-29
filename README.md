# CUDA-crystal_surface
Outputs the surface atoms of a regular crystal lattice

INTRODUCTION

This simple program showcases the use of CUDA in extracting the surface of a finite crystal lattice. In atomic simulations, the surface of a crystalline lattice is defined as the number of elements whose coordination number (the number of nearest neighbors to a given site) is less than bulk coordination number, which is constant for near-equilibrium conditions, and for particular atoms of the crystal.

For example, Face-Centered Cubic lattices have a bulk coordination number of 12 across the whole bulk, meaning there are 12 atoms nearest to any one particular lattice atom. Any atom with less than 12 nearest sites is, therefore, a break in translation symmetry, corresponding to a surface region, be it an open surface or a closed surface.

The lattice input must consist of a .xyz file (VMD - http://www.ks.uiuc.edu/Research/vmd/ and related formats). The file looks like:

N

a

X1 x.xxxxx y.yyyyy z.zzzzz

X2 x.xxxxx y.yyyyy z.zzzzz

X3 x.xxxxx y.yyyyy z.zzzzz

X4 x.xxxxx y.yyyyy z.zzzzz

X5 x.xxxxx y.yyyyy z.zzzzz

X6 x.xxxxx y.yyyyy z.zzzzz

...

Where N is the total number of (x,y,z) triples, 'a' is a line containing the lattice parameter, and each of the following N lines contain the (x,y,z) triple of coordinates for each crystalline site, as well as each Xn chemical element label. Overlapping sites with the very same (x,y,z) triples have undefined behavior (they are not expected by the program).

COMPILING THE SOURCE

The code requires a CUDA capable machine of compute capability 2.x or higher, and a working version of the CUDA developer toolkit(https://developer.nvidia.com/cuda-toolkit). Recommended to compile with one of the following in terminal:

-> General optimization flags

nvcc -O3 -arch=sm_xx --ptxas-options=-v -o CUDA_surface2.0.x CUDA_Surface_vs2.0.cu

-> In case compiler can't find your libraries

nvcc -O3 -arch=sm_xx --ptxas-options=-v -o CUDA_surface2.0.x CUDA_Surface_vs2.0.cu -lm

-> In case of memcpy_inline errors (usually with Ubuntu 16.04 and most versions of CUDA)

nvcc -O3 -arch=sm_XX --ptxas-options=-v -o CUDA_surface2.0.x CUDA_Surface_vs2.0.cu -lm -D_FORCE_INLINES

In all of the above, change the XX in -arch=sm_XX to the compute capability of the GPU you're using, e.g for a Tesla K20x (compute capability 3.5) we would use -arch=sm_35.

USING THE PROGRAM

You need a valid .xyz file as input, as described in the introduction. For information on how to build these, either check up on Chemical Structure Databases, such as http://www.chemspider.com/, http://support.cas.org/content/cas-databases, http://cccbdb.nist.gov/ , all of which can be used to acquire structures and convert into .xyz using, for example, OpenBabel (http://openbabel.org/wiki/Main_Page). The developer of this software uses shell scripts and other means to build his own files however. Refer to DOI:10.5151/phypro-sic100-046 or to "Kittel, Introduction to Solid State Physics" for how to build your own arbitrary structures.

A sample coordinate file with a very large array of the FCC Cu crystal (N = 10 160 790, compressed with bzip2) is presented as test. Simply decompress the input (e.g using bunzip2 coord_z.xyz.bz2) and run the program with ./ on the folder where the coord_z.xyz file is located and it will output two results, the first being a surface.xyz file containing the list of surface atoms (in the same format as the input), and a second file called num_edge.dat, which is a histogram ordered as:

grp 0

grp 1

grp 2

grp 3

...

Where "grp i" corresponds to the number of atoms in the structure with coordination number i. For instance, if grp 7 is equal to 240, that means 240 atoms have coordination number of 7 (and if they were in a FCC structure, they are surface atoms). It should be observed that, in order to prevent branch divergences in code, this histogram counts the atom itself as being its own neighbor. Adjustments to fit other kinds of lattice as well as relaxing of parameters (segregating atoms that are within two lattice parameters from a surface, for example) can be made very simply directly in code for the time being.

This code may be used in combination with https://github.com/mgmonteiro/CUDA-Drill_Struct to create different kinds of structures by using different kinds of order parameters (other than coordination number), as well as changing the criteria for structure building in that other code.
