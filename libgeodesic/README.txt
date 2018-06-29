libgeodesic v0.1
Keenan Crane
January 31, 2013

INTRODUCTION
------------

This library is an implementation of the heat method as described
in the paper

   Crane, Weischedel, Wardetzky
   "Geodesics in Heat: A New Approach to Computing Distance Using Heat Flow"
   ACM Transactions on Graphics (to appear)

Currently only triangulated surfaces are supported.  For more
information, see the Doxygen documentation in

   doc/html/index.html
   doc/html/refman.pdf

The libgeodesic library is written in 100% pure ANSI C, but
adopts a simple object-oriented programming style.  All symbols
from the library are preceded by "hm" (which stands for "heat
method"), followed by the object name, and finally the method
or member name.  For instance, the method used to read a
triangle mesh in Wavefront OBJ format is named

   hmTriMeshReadOBJ()

where TriMesh is the object name and ReadOBJ is the method name.


BUILDING
--------

The libgeodesic library requires a linear algebra library
capable of solving sparse positive definite systems.
Currently, either of the following libraries may be used:

   SuiteSparse
   http://www.cise.ufl.edu/research/sparse/SuiteSparse/
   
   HSL_MA87
   http://www.hsl.rl.ac.uk/

Note that the latter library is more appropriate for highly
parallel machines (the CHOLMOD library found in SuiteSparse
is currently not parallelized).  See the links above for
installation instructions.

IMPORTANT: once the appropriate libraries have been installed,
           the chosen library must be specified by
           
           1. setting appropriate flags in the Makefile
           2. setting definitions in include/hmConstants.h

Finally, libgeodesic can be built by simply typing

   make

at the command line.  The current build process uses a standard
GNU Makefile; there is no support whatsoever for VisualStudio
or any other IDE.  However, if you are familiar with one of
these environments you should have no problem buildling a project
yourself!  Just add all the relevant source and include files,
and link to the appropriate libraries as necessary.


RUNNING
-------

The build system should produce a command-line executable called
"geodesic".  This executable can be used to compute the distance
on a triangulated surface, storing the result in a variety of
ways.  For information on usage, simply type

   ./geodesic -help

(Note that at present the help may not be 100% compatible with
the actual implementation -- use at your own risk!)



LIBRARY USAGE
-------------

The listing below gives an example of the most basic possible usage at the code
level -- for more advanced usage, see the Doxygen documentation in doc/

   hmContext context;
   hmTriMesh surface;
   hmTriDistance distance;

   /* initialize data */
   hmContextInitialize( &context );
   hmTriMeshInitialize( &surface );
   hmTriDistanceInitialize( &distance );

   /* read surface */
   hmTriMeshReadOBJ( surface, filename );
   distance->surface = mesh;

   /* set time for heat flow */
   hmTriDistanceEstimateTime( &distance );

   /* compute distance */
   hmTriDistanceBuild( &distance );

   /* set the nth vertex in the mesh as a source */
   size_t n = 0;
   hmClearArrayDouble( distance->isSource.values, nVertices, 0. );
   distance->isSource.values[0] = 1.;

   /* udpate distance function */
   hmTriDistanceUpdate( &distance );

   /* print distances */
   for( i = 0; i < distance->surface->nVertices; i++ )
   {
      printf( "%.20e\n", distance->distance.values[i] );
   }

   /* deallocate data */
   hmTriMeshDestroy( &surface );
   hmTriDistanceDestroy( &distance );
   hmContextDestroy( &context );


