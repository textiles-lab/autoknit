/**
 * @file
 * @author Keenan Crane <keenan@cs.caltech.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright 2012 Keenan Crane. All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE FREEBSD PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the author and should not be interpreted as representing official policies,
 * either expressed or implied, of any other person or institution.
 *
 */

#ifndef LIBGEODESIC_HMTRIMESH_H
#define LIBGEODESIC_HMTRIMESH_H

#include <stddef.h>
#include "hmVectorPairSizeTDouble.h"

/** \brief Triangle mesh data structure.
 *
 * Represents a collection of triangles and their vertex locations via a simple
 * vertex-face adjacency structure.
 *
 */
typedef struct hmTriMesh {
   /** \brief Number of vertices. */
   size_t nVertices;

   /** \brief Vertex coordinates as consecutive (x,y,z) triples. [3 x nVertices] */
   double* vertices;

   /** \brief Texture coordinates as consecutive (u,v) pairs. [2 x nVertices] */
   double* texCoords;

   /** \brief Number of triangles. */
   size_t nFaces;

   /** \brief Triangles as consecutive (i,j,k) triples of 0-based vertex indices. [3 x nFaces] */
   size_t* faces;

   /** \brief Indicates whether mesh data is owned by another object (nonzero if true; zero otherwise). */
   char referenced;

   /** \private
    * \brief Flags whether each vertex is on the surface boundary. [1 x surface->nVertices] */
   char* onBoundary;

   /** \private
    * \brief Vertex 1-ring neighbors and associated weights. [1 x surface->nVertices] */
   hmVectorPairSizeTDouble* vertexNeighbors;

   /** \private
    * \brief Dual areas associated with vertices. [1 x surface->nVertices] */
   double* vertexAreas;

   /** \private
    * \brief Cotangents of angles at each triangle corner. [3 x surface->nFaces] */
   double* weights;

   /** \brief Tangent vectors orthogonal to the edges in each triangle (used to accelerate gradient computation). [9 x nFaces] */
   double* edgeNormals;

   /** \brief Edge vectors weighted by opposing angle cotangents in each triangle (used to accelerate divergence computation). [9 x nFaces] */
   double* weightedEdges;

} hmTriMesh;

/** \brief Constructor.
 *
 * @param mesh Target object.
 * \memberof hmTriMesh
 *
 */
void hmTriMeshInitialize( hmTriMesh* mesh );

/** \brief Destructor.
 *
 * @param mesh Target object.
 * \memberof hmTriMesh
 *
 */
void hmTriMeshDestroy( hmTriMesh* mesh );

/** \brief Replaces mesh1 with mesh2.
 *
 * @param mesh1 Target object.
 * @param mesh2 Source object.
 * \memberof hmTriMesh
 *
 */
void hmTriMeshCopy( hmTriMesh* mesh1, const hmTriMesh* mesh2 );

/** \brief Copies mesh data.
 *
 * @param mesh Target object.
 * @param nVertices Number of vertices in source data.
 * @param vertices Vertex coordinates as consecutive (x,y,z) triples. [3 x nVertices]
 * @param texcoords Texture coordinates as consecutive (u,v) triples. [2 x nVertices]
 * @param nFaces Number of triangles in source data.
 * @param faces Triangles as consecutive (i,j,k) triples of 0-based vertex indices. [3 x nFaces]
 * \memberof hmTriMesh
 *
 */
void hmTriMeshCopyData( hmTriMesh* mesh,
                        size_t nVertices, const double* vertices,
                        size_t nFaces,    const size_t* faces );

/** \brief References external mesh data.
 *
 * Allows external mesh data to be referenced instead of making a copy.
 * This construction can be useful when working with large meshes or
 * a large collection of meshes because it avoids the cost of the copy
 * and allows computation to be done "in place."  Note however that
 * the external representation must be compatible with the internal
 * representation of hmTriMesh, i.e., a list of interleaved double-
 * valued vertex coordinates and a list of interleaved size_t-valued
 * vertex indices for each face.
 *
 * @param mesh Target object.
 * @param nVertices Number of vertices in source data.
 * @param vertices Vertex coordinates as consecutive (x,y,z) triples. [3 x nVertices]
 * @param nFaces Number of triangles in source data.
 * @param faces Triangles as consecutive (i,j,k) triples of 0-based vertex indices. [3 x nFaces]
 * \memberof hmTriMesh
 *
 */
void hmTriMeshReferenceData( hmTriMesh* mesh,
                             size_t nVertices, double* vertices,
                             size_t nFaces,    size_t* faces );

/** \brief Reads mesh from a Wavefront OBJ file.
 *
 * Input should be an ASCII text file containing lines of the form
 *
 *    v [x] [y] [z]
 *
 * specifying vertex locations (x,y,z) and
 *
 *    f [I0] [I1] [I2]
 *
 * specifying triangles (I0,I1,I2) as indices into the vertex
 * list, where indices start at 1.  Vertex normals and texture
 * coordinates are not currently supported.
 *
 * @param mesh Target object.
 * @param filename Path to input file.
 * \memberof hmTriMesh
 *
 */
void hmTriMeshReadOBJ( hmTriMesh* mesh, const char* filename );

/** \brief Writes mesh to a Wavefront OBJ file.
 *
 * Output will be an ASCII text file containing lines of the form
 *
 *    v [x] [y] [z]
 *
 * specifying vertex locations (x,y,z) and
 *
 *    f [I0] [I1] [I2]
 *
 * specifying triangles (I0,I1,I2) as indices into the vertex
 * list, where indices start at 1.  Vertex normals and texture
 * coordinates are not currently supported.
 *
 * @param mesh Source object.
 * @param filename Path to output file.
 * \memberof hmTriMesh
 *
 */
void hmTriMeshWriteOBJ( const hmTriMesh* mesh, const char* filename );

/** \private
 * \brief Updates mesh attributes (vertex areas, 1-ring neighbors, etc.).
 *
 * @param mesh Target object.
 * \memberof hmTriMesh
 *
 */
void hmTriMeshBuild( hmTriMesh* mesh );

/** \private
 * \brief Computes the 1-ring neighbors of each vertex.
 *
 * @param distance Target object.
 * \memberof hmTriMesh
 *
 */
void hmTriMeshFindOneRings( hmTriMesh* distance );

/** \private
 * \brief Computes dual areas associated with vertices.
 *
 * Uses barycentric dual areas, i.e., one-third the sum of the total
 * area of the triangles incident on each vertex.  These areas are
 * always well-defined and positive, and result in operators with the
 * same order of accuracy (\f$O(h)\f$) as those obtained by using
 * circumcentric dual cells.
 *
 * @param mesh Target object.
 * \memberof hmTriMesh
 *
 */
void hmTriMeshComputeVertexAreas( hmTriMesh* mesh );

/** \private
 * \brief Computes partial edge weights associated with each triangle corner.
 *
 * For piecewise linear finite elements, the Laplace operator is
 * expressed as the graph Laplacian with edge weights equal to
 * \f$\frac{1}{2}(\cot\alpha + \cot\beta)\f$, where \f$\alpha\f$
 * and \f$\beta\f$ are the two opposing angles.  This discretization
 * yields an \f$O(h)\f$ approximation where \f$h\f$ is the mean
 * edge length.
 *
 * @param distance Target object.
 * \memberof hmTriMesh
 *
 */
void hmTriMeshComputeWeights( hmTriMesh* mesh );

/** \brief Precomputes triangle edge data used to compute gradient and divergence.
 *
 * @param mesh Target object.
 * \memberof hmTriMesh
 *
 */
void hmTriMeshCacheEdgeData( hmTriMesh* mesh );

/** \brief \f$L^2\f$ distance between two scalar functions at vertices.
 *
 * For two functions \f$\phi_1\f$ and \f$\phi_2\f$, returns the quantity
 * \f$\left(\int_M (\phi_1-\phi_2)^2 \right)^{1/2}\f$.  For piecewise
 * linear functions this quantity can be computed exactly as
 * \f$\left(\sum_{i=1}^V \mathcal{A}_i (\phi_1^i-\phi_2^i)^2 \right)^{1/2}\f$
 * where \f$V\f$ is the number of vertices, \f$\phi_j^i\f$ is the function
 * value of the \f$j\f$th function at the \f$i\f$th vertex, and
 * \f$\mathcal{A}_i\f$ is the barycentric dual area associated with the
 * \f$i\f$th vertex.
 *
 * @param mesh Target object.
 * @param phi1 First function. [1 x nVertices]
 * @param phi2 Second function. [1 x nVertices]
 * @return Distance between functions.
 * \memberof hmTriMesh
 *
 */
double hmTriMeshL2Distance( hmTriMesh* mesh, double* phi1, double* phi2 );

/** \brief \f$L^\infty\f$ distance between two scalar functions at vertices.
 *
 * Returns the maximum (in magnitude) pointwise difference between two
 * functions \f$\phi_1\f$ and \f$\phi_2\f$.
 *
 * @param mesh Target object.
 * @param phi1 First function. [1 x nVertices]
 * @param phi2 Second function. [1 x nVertices]
 * @return Distance between functions.
 * \memberof hmTriMesh
 *
 */
double hmTriMeshLInfinityDistance( hmTriMesh* mesh, double* phi1, double* phi2 );

/** \brief Mean of pointwise relative error.
 *
 * Given two functions \f$\phi_1\f$ and \f$\phi_2\f$, returns the quantity
 * \f$\frac{1}{N}\sum_{i=1}^N |\phi_2/\phi_1|\f$.
 *
 * @param mesh Target object.
 * @param phi1 First function. [1 x nVertices]
 * @param phi2 Second function. [1 x nVertices]
 * @return Mean relative error.
 * \memberof hmTriMesh
 *
 */
double hmTriMeshMeanRelativeError( hmTriMesh* mesh, double* phi1, double* phi2 );

/** \brief Diameter of bounding sphere.
 *
 * Returns twice the maximum distance of any vertex from the center of mass.
 *
 * @param mesh Target object.
 * @return Mesh diameter.
 * \memberof hmTriMesh
 *
 */
double hmTriMeshDiameter( hmTriMesh* mesh );

#endif /* LIBGEODESIC_HMTRIMESH_H */
