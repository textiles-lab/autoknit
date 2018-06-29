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

#ifndef LIBGEODESIC_HMTRIDISTANCE_H
#define LIBGEODESIC_HMTRIDISTANCE_H

#include "hmTriMesh.h"
#include "hmVectorDouble.h"
#include "hmVectorPairSizeTDouble.h"
#include "hmDenseMatrix.h"
#include "hmSparseMatrix.h"
#include "hmCholeskyFactor.h"
#include <stddef.h>

/** \brief Distance function on a triangle mesh.
 *
 * Represents the geodesic distance function \f$\phi: M \rightarrow
 * \mathbb{R}\f$ to a given set \f$\gamma\f$ on a surface \f$M\f$.  The domain
 * \f$M\f$ is represented by a hmTriMesh -- note that no guarantees can be made
 * about the solution unless this mesh is _manifold_ (i.e., every edge
 * contained in at most two triangles; every vertex contained in a single ring
 * of triangles) and _oriented_ (i.e., all triangles specified in a consistent
 * winding order -- either clockwise or counter-clockwise).  The source set
 * is specified by setting values in hmTriDistance::isSource.  Once the
 * distance function has been initialized, it can be rapidly (re-)computed for
 * a new source set \f$\gamma\f$ by calling hmTriDistanceUpdate(); the
 * solution is stored in hmTriDistance::distance. Distance is computed using
 * the method described in Crane, Weischedel, and Wardetzky, _Geodesics in
 * Heat_, [arXiv:1204.6216v1](http://arxiv.org/abs/1204.6216) (2012).
 *
 */
typedef struct hmTriDistance {

   /* PUBLIC MEMBERS =================================================== */

   /** \brief Mesh of the domain \f$M\f$.  (__Note:__ since this member may reference an external object, it is not automatically deallocated upon destruction.) */
   hmTriMesh* surface;

   /** \brief Integration time \f$t\f$ for heat flow.
    *
    * Larger values of \f$t\f$ result in smoother approximations of
    * geodesic distance.  The value that best approximates the true
    * geodesic distance can be estimated by calling hmTriDistanceEstimateTime().
    * 
    * */
   double time;

   /** \brief Interpolates between domain boundary conditions; 0 means pure Neumann, 1 means pure Dirichlet. */
   double boundaryConditions;

   /** \brief Source set, equal to 1 at source vertices; 0 otherwise. [surface->nVertices x 1] */
   hmDenseMatrix isSource;

   /** \brief Distance values \f$\phi\f$. [surface->nVertices x 1] */
   hmDenseMatrix distance;

   /** \brief Flags whether to use verbose output (timing information, etc.). */
   char verbose;


   /* PRIVATE MEMBERS ================================================== */

   /** \private
    * \brief Solution \f$u\f$ to heat equation with Neumann boundary conditions. [surface->nVertices x 1] */
   hmDenseMatrix heatNeumann;

   /** \private
    * \brief Solution \f$u\f$ to heat equation with Dirichlet boundary conditions. [surface->nVertices x 1] */
   hmDenseMatrix heatDirichlet;

   /** \private
    * \brief Final solution to heat equation (points to entries of either hmTriDistance::heatNeumann or hmTriDistance::heatDirichlet). */
   double* heat;

   /** \private
    * \brief Divergence \f$\nabla \cdot X\f$ of normalized vector field \f$X\f$. [1 x surface->nVertices] */
   hmDenseMatrix potential;

   /** \private
    * \brief Conformal (i.e., weak) Laplacian \f$L_C\f$. */
   hmSparseMatrix laplacian;

   /** \private
    * \brief One-step heat flow operator \f$A-tL_C\f$ with zero Neumann boundary conditions. */
   hmSparseMatrix heatFlowNeumann;

   /** \private
    * \brief One-step heat flow operator \f$A-tL_C\f$ with zero Dirichlet boundary conditions. */
   hmSparseMatrix heatFlowDirichlet;

   /** \private
    * \brief Cholesky factor for hmTriDistance::Laplacian. */
   hmCholeskyFactor laplacianFactor;

   /** \private
    * \brief Cholesky factor for hmTriDistance::heatFlowNeumann. */
   hmCholeskyFactor heatFlowNeumannFactor;

   /** \private
    * \brief Cholesky factor for hmTriDistance::heatFlowDirichlet. */
   hmCholeskyFactor heatFlowDirichletFactor;

   /** \private
    * \brief Stack of start times in seconds with respect to processor time (for profiling). */
   hmVectorDouble startTimesProcessor;

   /** \private
    * \brief Stack of start times in seconds with respect to wall clock time (for profiling). */
   hmVectorDouble startTimesWallClock;

} hmTriDistance;


/* PUBLIC METHODS ====================================================== */

/** \brief Constructor.
 *
 * @param distance Target object.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceInitialize( hmTriDistance* distance );

/** \brief Destructor.
 *
 * @param distance Target object.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceDestroy( hmTriDistance* distance );

/** \brief Estimates a good value for the duration of heat flow.
 *
 * Uses the simple estimate \f$t = A/F\f$ where \f$A\f$ is the total
 * surface area and \f$F\f$ is the number of triangles in the mesh.  The
 * resulting value is stored in hmTriDistance::time.
 *
 * @param distance Target object.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceEstimateTime( hmTriDistance* distance );

/** \brief Specifies boundary conditions.
 *
 * Boundary conditions are specified via a parameter in the range \f$[0,1]\f$.
 * A value of zero yields pure Neumann conditions; a value of one yields
 * pure Dirichlet conditions.  Values in between interpolate between the two
 * boundary conditions.  For most purposes (namely: when the surface has no
 * boundary or when an accurate approximation of geodesic distance is
 * desired) this parameter has little effect and optimal performance will be
 * achieved by simply setting it to zero.  For surfaces with boundary and
 * large heat flow time (i.e., smoothed geodesic distance) a value of 1/2
 * tends to produce natural-looking behavior near the boundary.  __Note:__
 * in order for this parameter to take effect, its value must be set
 * _before_ calling hmTriDistanceBuild().
 *
 * @param distance Target object.
 * @param boundaryConditions Interpolation parameter.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceSetBoundaryConditions( hmTriDistance* distance,
                                         double boundaryConditions );

/** \brief Performs precomputation required for distance computation.
 *
 * The mesh hmTriDistance::surface and parameter hmTriDistance::time
 * must be defined.  This routine must be run before calling
 * hmTriDistanceUpdate().
 *
 * @param distance Target object.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceBuild( hmTriDistance* distance );

/** \brief Rebuilds data structures for a new flow time parameter \f$t\f$.
 *
 * This routine is useful primarily when computing solutions for a variety of
 * flow times \f$t\f$ -- it does not need to be called in standard usage.
 *
 * @param distance Target object.
 * @param time Heat flow integration time.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceUpdateTime( hmTriDistance* distance,
                              double time );

/** \brief Updates the distance function using the current source set \f$\gamma\f$.
 *
 * @param distance Target object.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceUpdate( hmTriDistance* distance );


/* PRIVATE METHODS ===================================================== */

/** \private
 * \brief Solves the heat equation.
 *
 * Uses one step of backward Euler for time \f$t\f$ with initial
 * conditions \f$u_0\f$ determined by the source set \f$\gamma\f$ in order
 * to obtain a function \f$u\f$ suitable as input for the heat method.
 * In particular, \f$u\f$ solves the linear elliptic problem
 * \f$(I-t\Delta)u = u_0\f$ where \f$\Delta\f$ is the negative-definite
 * Laplace-Beltrami operator on the domain \f$M\f$.  The Laplacian is
 * discretized using standard first-order piecewise linear finite elements
 * (i.e., the cotan-Laplacian), and the initial data is a Kronecker delta
 * over the source set.  The result is stored in hmTriDistance::heat as one
 * value per vertex.
 *
 * @param distance Target object.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceSolveHeatEquation( hmTriDistance* distance );

/** \private
 * \brief Computes the scalar potential of the distance function.
 *
 * Given an input function \f$u\f$, this routine computes the scalar
 * potential \f$\nabla \cdot X\f$, where \f$X = -\nabla u/|\nabla u|\f$
 * is the normalized gradient field of \f$u\f$.  The gradient and
 * divergence operators are discretized using standard piecewise linear
 * finite elements -- see _Crane et al_ for the relevant expressions.
 * The result is stored in hmTriDistance::potential as one value per vertex.
 *
 * @param distance Target object.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceComputePotential( hmTriDistance* distance );

/** \private
 * \brief Solves a Poisson equation.
 *
 * Solves the linear elliptic Poisson equation \f$\Delta\phi = f\f$
 * for the distance function \f$\phi\f$, where \f$\Delta\f$ is the
 * Laplace-Beltrami operator on the domain \f$M\f$ and \f$f\f$ is the scalar
 * potential computed by hmTriDistanceComputePotential().  The Laplacian is
 * discretized using standard first-order piecewise linear finite elements
 * (i.e., the cotan-Laplacian).  The result is stored in hmTriDistance::distance
 * as one value per vertex.
 *
 * @param distance Target object.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceSolvePoissonEquation( hmTriDistance* distance );

/** \private
 * \brief Builds matrices used for distance computation.
 *
 * @param distance Target object.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceBuildMatrices( hmTriDistance* distance );

/** \private
 * \brief Prefactors matrices used for distance computation.
 *
 * @param distance Target object.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceFactorMatrices( hmTriDistance* distance );

/** \private
 * \brief Starts timing a subroutine.
 *
 * To print elapsed time, call hmTriDistanceStopTiming().  Timings
 * are nested: the next call to hmTriDistanceStopTiming() will print
 * the time elapsed since the most recent call to
 * hmTriDistanceStartTiming().
 *
 * @param distance Target object.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceStartTiming( hmTriDistance* distance );

/** \private
 * \brief Stops timing a subroutine.
 *
 * Prints time elapsed in seconds together with the specified
 * label string.  Timings are nested: the next call to
 * hmTriDistanceStopTiming() will print the time elapsed since
 * the most recent call to hmTriDistanceStartTiming().
 *
 * @param distance Target object.
 * @param label Label for time interval -- pass NULL for no label.
 * \memberof hmTriDistance
 *
 */
void hmTriDistanceStopTiming( hmTriDistance* distance,
                              const char* label );

#endif /* LIBGEODESIC_HMTRIDISTANCE_H */

