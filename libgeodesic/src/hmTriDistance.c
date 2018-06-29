#include "hmTriDistance.h"
#include "hmVec3.h"
#include "hmUtility.h"
#include "hmConstants.h"
#include "hmVectorSizeT.h"
#include <time.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/time.h>

void hmTriDistanceInitialize( hmTriDistance* distance )
{
   distance->surface         = NULL;

   hmDenseMatrixInitialize( &distance->isSource,      0, 0 );
   hmDenseMatrixInitialize( &distance->distance,      0, 0 );
   hmDenseMatrixInitialize( &distance->heatNeumann,   0, 0 );
   hmDenseMatrixInitialize( &distance->heatDirichlet, 0, 0 );
   hmDenseMatrixInitialize( &distance->potential,     0, 0 );

   hmSparseMatrixInitialize( &distance->laplacian,         0, 0, 0 );
   hmSparseMatrixInitialize( &distance->heatFlowNeumann,   0, 0, 0 );
   hmSparseMatrixInitialize( &distance->heatFlowDirichlet, 0, 0, 0 );

   hmCholeskyFactorInitialize( &distance->laplacianFactor );
   hmCholeskyFactorInitialize( &distance->heatFlowNeumannFactor );
   hmCholeskyFactorInitialize( &distance->heatFlowDirichletFactor );

   /* by default, use pure Neumann boundary conditions */
   distance->boundaryConditions = 0.;

   hmVectorDoubleInitialize( &distance->startTimesProcessor );
   hmVectorDoubleInitialize( &distance->startTimesWallClock );
}

void hmTriDistanceDestroy( hmTriDistance* distance )
{
   hmDenseMatrixDestroy( &distance->isSource );
   hmDenseMatrixDestroy( &distance->distance );
   hmDenseMatrixDestroy( &distance->heatNeumann );
   hmDenseMatrixDestroy( &distance->heatDirichlet );
   hmDenseMatrixDestroy( &distance->potential );

   hmSparseMatrixDestroy( &distance->laplacian );
   hmSparseMatrixDestroy( &distance->heatFlowNeumann );
   hmSparseMatrixDestroy( &distance->heatFlowDirichlet );

   hmCholeskyFactorDestroy( &distance->laplacianFactor );
   hmCholeskyFactorDestroy( &distance->heatFlowNeumannFactor );
   hmCholeskyFactorDestroy( &distance->heatFlowDirichletFactor );

   hmVectorDoubleDestroy( &distance->startTimesProcessor );
   hmVectorDoubleDestroy( &distance->startTimesWallClock );
}

void hmTriDistanceEstimateTime( hmTriDistance* distance )
{
   size_t nFaces = distance->surface->nFaces;
   size_t* facesBegin = distance->surface->faces;
   size_t* facesEnd = facesBegin + 3*nFaces;
   size_t* f;
   double* vertices = distance->surface->vertices;
   double *p0, *p1, *p2;
   hmVec3 e01, e12, e20;
   double meanEdgeLength = 0.;
   double nEdges = 0.;

   /* iterate over faces */
   for( f = facesBegin; f != facesEnd; f += 3 )
   {
      /* get vertex coordinates p0, p1, p2 */
      p0 = &vertices[ f[0]*3 ];
      p1 = &vertices[ f[1]*3 ];
      p2 = &vertices[ f[2]*3 ];

      /* add edge lengths to mean */
      hmVec3Sub( e01, p1, p0 );
      hmVec3Sub( e12, p2, p1 );
      hmVec3Sub( e20, p0, p2 );
      meanEdgeLength += hmVec3Norm( e01 );
      meanEdgeLength += hmVec3Norm( e12 );
      meanEdgeLength += hmVec3Norm( e20 );

      nEdges += 3.;
   }
   meanEdgeLength /= nEdges;

   /* set t to square of mean edge length */
   distance->time = hmSquare( meanEdgeLength );
}

void hmTriDistanceSetBoundaryConditions( hmTriDistance* distance,
                                         double boundaryConditions )
{
   distance->boundaryConditions = boundaryConditions;
}

void hmTriDistanceBuild( hmTriDistance* distance )
{
   size_t nVertices = distance->surface->nVertices;

   if( distance->surface == NULL )
   {
      fprintf( stderr, "Error: hmTriDistanceBuild -- must specify a surface!\n" );
      exit( 1 );
   }

   hmTriDistanceDestroy( distance );

   hmVectorDoubleInitialize( &distance->startTimesProcessor );
   hmVectorDoubleInitialize( &distance->startTimesWallClock );
   hmTriDistanceStartTiming( distance );

   hmDenseMatrixInitialize( &distance->isSource,  nVertices, 1 );
   hmDenseMatrixInitialize( &distance->distance,  nVertices, 1 );
   hmDenseMatrixInitialize( &distance->potential, nVertices, 1 );

   /* only allocate space for both solutions if necessary */
   if( distance->boundaryConditions < 1. ) /* partial Neumann */
   {
      hmDenseMatrixInitialize( &distance->heatNeumann,   nVertices, 1 );
   }
   if( distance->boundaryConditions > 0. ) /* partial Dirichlet */
   {
      hmDenseMatrixInitialize( &distance->heatDirichlet, nVertices, 1 );
   }

   hmTriMeshBuild( distance->surface );
   hmTriDistanceBuildMatrices( distance );
   hmTriDistanceFactorMatrices( distance );

   hmTriDistanceStopTiming( distance, "Total build time" );
}

void hmTriDistanceUpdate( hmTriDistance* distance )
{
   hmTriDistanceStartTiming( distance );

   hmTriDistanceStartTiming( distance );
   hmTriDistanceSolveHeatEquation( distance );
   hmTriDistanceStopTiming( distance, "Solve heat equation" );

   hmTriDistanceStartTiming( distance );
   hmTriDistanceComputePotential( distance );
   hmTriDistanceStopTiming( distance, "Compute potential" );

   hmTriDistanceStartTiming( distance );
   hmTriDistanceSolvePoissonEquation( distance );
   hmTriDistanceStopTiming( distance, "Solve Poisson equation" );

   hmTriDistanceStopTiming( distance, "Total update time" );
}

void hmTriDistanceSolveHeatEquation( hmTriDistance* distance )
{
   size_t i;
   size_t nVertices = distance->surface->nVertices;
   const double BC = distance->boundaryConditions;
   const double rBC = 1.-BC;
   double *uN, *uD;

   /* only compute both solutions if necessary */
   if( BC < 1. ) /* partial Neumann */
   {
      hmCholeskyFactorBacksolve( &distance->heatFlowNeumannFactor,
                                 &distance->heatNeumann,
                                 &distance->isSource );
   }
   if( BC > 0. ) /* partial Dirichlet */
   {
      hmCholeskyFactorBacksolve( &distance->heatFlowDirichletFactor,
                                 &distance->heatDirichlet,
                                 &distance->isSource );
   }

   /* store the final solution in hmTriDistance::heat,
    * combining the two solutions if necessary */
   if( BC > 0. && BC < 1. )
   {
      /* get pointers to entries of Neumann and Dirichlet solutions */
      uN = distance->heatNeumann.values;
      uD = distance->heatDirichlet.values;

      /* write interpolated solution over Neumann solution */
      distance->heat = distance->heatNeumann.values;

      /* compute interpolated solution */
      for( i = 0; i < nVertices; i++ )
      {
         distance->heat[i] = rBC*uN[i] + BC*uD[i];
      }
   }
   else if( BC == 0. ) /* pure Neumann */
   {
      distance->heat = distance->heatNeumann.values;
   }
   else /* pure Dirichlet */
   {
      distance->heat = distance->heatDirichlet.values;
   }
}

void hmTriDistanceComputePotential( hmTriDistance* distance )
{
   /* array counters */
   size_t i;

   /* local data handles */
   int nFaces    = distance->surface->nFaces;
   int nVertices = distance->surface->nVertices;
   const size_t*         f = distance->surface->faces;
   const double*         w = distance->surface->weights;
   const double*      heat = distance->heat;
         double* potential = distance->potential.values;

   /* current triangle data */
   double u0, u1, u2; /* heat values */
   double rMag; /* reciprocal of magnitude */
   double *t0, *t1, *t2; /* edge normals */
   double *e0, *e1, *e2; /* cotan-weighted edge vectors */
   hmVec3 X; /* normalized gradient */
   double e0DotX, e1DotX, e2DotX;

   /* initialize potential to zero */
   hmClearArrayDouble( potential, distance->surface->nVertices, 0. );

   /* get pointers to first three edge normals */
   t0 = &distance->surface->edgeNormals[0];
   t1 = &distance->surface->edgeNormals[3];
   t2 = &distance->surface->edgeNormals[6];

   /* get pointers to first three weighted edges */
   e0 = &distance->surface->weightedEdges[0];
   e1 = &distance->surface->weightedEdges[3];
   e2 = &distance->surface->weightedEdges[6];

   /* add contribution from each face */
   for( i = 0; i < nFaces; i++ )
   {
      /* get heat values at three vertices */
      u0 = fabs( heat[ f[0] ] );
      u1 = fabs( heat[ f[1] ] );
      u2 = fabs( heat[ f[2] ] );

      /* normalize heat values so that they have roughly unit magnitude */
      rMag = 1./hmMaxDouble( hmMaxDouble( u0, u1 ), u2 );
      if( !isinf(rMag) )
      {
         u0 *= rMag;
         u1 *= rMag;
         u2 *= rMag;

         /* compute normalized gradient */
         X[0] = u0*t0[0] + u1*t1[0] + u2*t2[0];
         X[1] = u0*t0[1] + u1*t1[1] + u2*t2[1];
         X[2] = u0*t0[2] + u1*t1[2] + u2*t2[2];
         hmVec3Normalize( X );

         /* add contribution to divergence */
         e0DotX = hmVec3Dot( e0, X );
         e1DotX = hmVec3Dot( e1, X );
         e2DotX = hmVec3Dot( e2, X );
         potential[ f[0] ] -= e1DotX - e2DotX;
         potential[ f[1] ] -= e2DotX - e0DotX;
         potential[ f[2] ] -= e0DotX - e1DotX;

         if( isnan( potential[f[0]] ) ||
             isnan( potential[f[1]] ) ||
             isnan( potential[f[2]] ) )
         {
            fprintf( stderr, "NaN\n============\n" );
            fprintf( stderr, "heat: %e %e %e\n", heat[f[0]], heat[f[1]], heat[f[2]] );
            fprintf( stderr, " mag: %e\n",  hmMaxDouble( hmMaxDouble( u0, u1 ), u2 ));
            fprintf( stderr, "rMag: %e\n", rMag );
            fprintf( stderr, "   u: %e %e %e\n", u0, u1, u2 );
            fprintf( stderr, "   X: %e %e %e\n", X[0], X[1], X[2] );
            fprintf( stderr, "ei*X: %e %e %e\n", e0DotX, e1DotX, e2DotX );
            exit( 1 );
         }
      }
      
      /* move to next face */
      f += 3;
      w += 3;
      t0 += 9; t1 += 9; t2 += 9;
      e0 += 9; e1 += 9; e2 += 9;
   }

   /* remove mean value so that the potential is
    * in the range of the Laplace operator */
   hmRemoveMean( potential, nVertices );
}

void hmTriDistanceSolvePoissonEquation( hmTriDistance* distance )
{
   size_t i;
   const size_t nVertices = distance->surface->nVertices;
   double minDistance = DBL_MAX;
   double* phi;

   hmCholeskyFactorBacksolve( &distance->laplacianFactor,
                              &distance->distance,
                              &distance->potential );

   /* subtract the minimum value */
   phi = distance->distance.values;
   for( i = 0; i < nVertices; i++ )
   {
      minDistance = hmMinDouble( minDistance, phi[i] );
   }
   for( i = 0; i < nVertices; i++ )
   {
      phi[i] -= minDistance;
   }
}

void hmTriDistanceBuildMatrices( hmTriDistance* distance )
{
   size_t i;
   size_t nVertices = distance->surface->nVertices;
   size_t nNonZeros = 0;
   size_t nz; /* current nonzero */
   hmVectorSizeT columnStart;
   const hmVectorPairSizeTDouble* neighbors = distance->surface->vertexNeighbors;
   char* onBoundary = distance->surface->onBoundary;
   double columnSum;
   char diagonalSet;
   const hmPairSizeTDouble *neighborsBegin, *neighborsEnd, *currentNeighbor;
   const double boundaryConditions = distance->boundaryConditions;
   const double time = distance->time;
   double A; /* vertex area */

   hmSparseMatrix* laplacian         = &distance->laplacian;
   hmSparseMatrix* heatFlowNeumann   = &distance->heatFlowNeumann;
   hmSparseMatrix* heatFlowDirichlet = &distance->heatFlowDirichlet;

   /* determine the starting entry of nonzeros in
    * each column, keeping the lower triangle only */
   hmVectorSizeTInitialize( &columnStart );
   hmVectorSizeTPushBack( &columnStart, 0 );
   for( i = 0; i < nVertices; i++ )
   {
      /* add one for the diagonal entry */
      nNonZeros++;

      /* add one for each off-diagonal entry below the diagonal */
      neighborsBegin = neighbors[i].entries;
      neighborsEnd   = neighborsBegin + neighbors[i].size;
      for( currentNeighbor  = neighborsBegin;
           currentNeighbor != neighborsEnd;
           currentNeighbor ++ )
      {
         if( currentNeighbor->n > i )
         {
            nNonZeros++;
         }
      }

      /* keep track of the nonzeros so far */
      hmVectorSizeTPushBack( &columnStart, nNonZeros );
   }

   /* initialize matrices and copy column start pointers */
   hmSparseMatrixDestroy( laplacian );
   hmSparseMatrixInitialize( laplacian, nVertices, nVertices, nNonZeros );
   for( i = 0; i < nVertices+1; i++ )
   {
      distance->laplacian.columnStart[i] = columnStart.entries[i];
   }

   /* only build both heat flow operators if necessary */
   if( boundaryConditions < 1. ) /* partial Neumann */
   {
      hmSparseMatrixDestroy( heatFlowNeumann );
      hmSparseMatrixInitialize( heatFlowNeumann,  nVertices, nVertices, nNonZeros );
      for( i = 0; i < nVertices+1; i++ )
      {
         distance->heatFlowNeumann.columnStart[i] = columnStart.entries[i];
      }
   }
   if( boundaryConditions > 0. ) /* partial Dirichlet */
   {
      hmSparseMatrixDestroy( heatFlowDirichlet );
      hmSparseMatrixInitialize( heatFlowDirichlet,  nVertices, nVertices, nNonZeros );
      for( i = 0; i < nVertices+1; i++ )
      {
         distance->heatFlowDirichlet.columnStart[i] = columnStart.entries[i];
      }
   }

   /* fill nonzero entries */
   nz = 0;
   for( i = 0; i < nVertices; i++ )
   {
      /* get beginning and end of neighbor list for this vertex */
      neighborsBegin = neighbors[i].entries;
      neighborsEnd   = neighborsBegin + neighbors[i].size;

      /* get sum of neighbors' weights */
      columnSum = 0.;
      for( currentNeighbor  = neighborsBegin;
           currentNeighbor != neighborsEnd;
           currentNeighbor ++ )
      {
         columnSum += currentNeighbor->x;
      }

      /* set entries */
      diagonalSet = 0; /* the diagonal entry has not yet been set */
      for( currentNeighbor  = neighborsBegin;
           currentNeighbor != neighborsEnd+1;
           currentNeighbor ++ )
      {
         /* if we haven't yet set the diagonal entry and we've
          * already set all the off-diagonal entries OR we encounter
          * a neighbor with index smaller than the index of the
          * current vertex, then insert the diagonal entry here */
         if( !diagonalSet &&
             ( currentNeighbor == neighborsEnd || currentNeighbor->n > i ))
         {
            diagonalSet = 1;

            /* set diagonal entry of Laplacian, adding a small
               regularization term in order to get strict
               positive-definiteness (needed for CHOLMOD) */
            laplacian->values[ nz ] = columnSum + hmRegularization;
            laplacian->rowIndices[ nz ] = i;

            A = distance->surface->vertexAreas[ i ];

            if( boundaryConditions < 1. ) /* partial Neumann */
            {
               heatFlowNeumann->values[ nz ] = A + time*columnSum;
               heatFlowNeumann->rowIndices[ nz ] = i;
            }
            if( boundaryConditions > 0. ) /* partial Dirichlet */
            {
               if( onBoundary[ i ] )
               {
                  /* use the identity (times the mass matrix) for boundary
                   * rows/columns to enforce zero-Dirichlet conditions */
                  heatFlowDirichlet->values[ nz ] = A;
               }
               else
               {
                  heatFlowDirichlet->values[ nz ] = A + time*columnSum;
               }
               heatFlowDirichlet->rowIndices[ nz ] = i;
            }
            nz++;
         }

         /* set off-diagonal entries below the diagonal */
         if( currentNeighbor < neighborsEnd &&
             currentNeighbor->n > i )
         {
            laplacian->values[ nz ] = -currentNeighbor->x;
            laplacian->rowIndices[ nz ] = currentNeighbor->n;

            if( boundaryConditions < 1. ) /* partial Neumann */
            {
               heatFlowNeumann->values[ nz ] = -time*currentNeighbor->x;
               heatFlowNeumann->rowIndices[ nz ] = currentNeighbor->n;
            }
            if( boundaryConditions > 0. ) /* partial Dirichlet */
            {
               if( onBoundary[i] || onBoundary[ currentNeighbor->n ] )
               {
                  /* set off-diagonals to zero so that we retain
                   * the same sparsity pattern as other matrices */
                  heatFlowDirichlet->values[ nz ] = 0.;
               }
               else
               {
                  heatFlowDirichlet->values[ nz ] = -time*currentNeighbor->x;
               }
               heatFlowDirichlet->rowIndices[ nz ] = currentNeighbor->n;
            }
            nz++;
         }
      }
   }
}

void hmTriDistanceFactorMatrices( hmTriDistance* distance )
{
   /* Laplacian */
   hmCholeskyFactorDestroy    ( &distance->laplacianFactor );
   hmCholeskyFactorInitialize ( &distance->laplacianFactor );
   hmCholeskyFactorReorder    ( &distance->laplacianFactor, &distance->laplacian );
   hmCholeskyFactorSymbolic   ( &distance->laplacianFactor, &distance->laplacian );
   hmCholeskyFactorNumerical  ( &distance->laplacianFactor, &distance->laplacian );

   /* only factor both heat flow operators if necessary */
   /* (note that the symbolic factorization for Laplace can be reused in both
    * cases since all three matrices have the same sparsity pattern) */
   if( distance->boundaryConditions < 1. ) /* partial Neumann */
   {
      hmCholeskyFactorDestroy    ( &distance->heatFlowNeumannFactor );
      hmCholeskyFactorInitialize ( &distance->heatFlowNeumannFactor );
      hmCholeskyFactorCopy       ( &distance->heatFlowNeumannFactor, &distance->laplacianFactor );
#ifdef HM_USE_HSLMA87 /* currently no way to copy symbolic factorization in HSL_MA87... */
      hmCholeskyFactorSymbolic   ( &distance->heatFlowNeumannFactor, &distance->heatFlowNeumann );
#endif
      hmCholeskyFactorNumerical  ( &distance->heatFlowNeumannFactor, &distance->heatFlowNeumann );
   }
   if( distance->boundaryConditions > 0. ) /* partial Dirichlet */
   {
      hmCholeskyFactorDestroy    ( &distance->heatFlowDirichletFactor );
      hmCholeskyFactorInitialize ( &distance->heatFlowDirichletFactor );
      hmCholeskyFactorCopy       ( &distance->heatFlowDirichletFactor, &distance->laplacianFactor );
#ifdef HM_USE_HSLMA87
#else
      hmCholeskyFactorSymbolic   ( &distance->heatFlowDirichletFactor, &distance->heatFlowDirichlet );
#endif
      hmCholeskyFactorNumerical  ( &distance->heatFlowDirichletFactor, &distance->heatFlowDirichlet );
   }
}

void hmTriDistanceUpdateTime( hmTriDistance* distance, double time )
{
   hmTriDistanceStartTiming( distance );

   distance->time = time;

   hmTriDistanceBuildMatrices( distance );

   hmCholeskyFactorNumerical( &distance->laplacianFactor, &distance->laplacian );

   if( distance->boundaryConditions < 1. ) /* partial Neumann */
   {
      hmCholeskyFactorNumerical( &distance->heatFlowNeumannFactor, &distance->heatFlowNeumann );
   }

   if( distance->boundaryConditions > 0. ) /* partial Dirichlet */
   {
      hmCholeskyFactorNumerical( &distance->heatFlowDirichletFactor, &distance->heatFlowDirichlet );
   }

   hmTriDistanceStopTiming( distance, "Update t parameter" );
}

void hmTriDistanceStartTiming( hmTriDistance* distance )
{
   struct timeval t;
   struct timezone tz;
   double startTimeProcessor;
   double startTimeWallClock;

   if( !distance->verbose )
   {
      return;
   }

   startTimeProcessor = (double) clock() / (double) CLOCKS_PER_SEC;
   hmVectorDoublePushBack( &distance->startTimesProcessor, startTimeProcessor );
   
   gettimeofday( &t, &tz );
   startTimeWallClock = (double) t.tv_sec + 1e-6*(double) t.tv_usec;
   hmVectorDoublePushBack( &distance->startTimesWallClock, startTimeWallClock );
}

void hmTriDistanceStopTiming( hmTriDistance* distance,
                              const char* label )
{
   int i;
   struct timeval t;
   struct timezone tz;
   double startTimeProcessor, stopTimeProcessor;
   double startTimeWallClock, stopTimeWallClock;
   int nTabs = distance->startTimesProcessor.size-1;

   if( !distance->verbose )
   {
      return;
   }

   startTimeProcessor = hmVectorDoublePopBack( &distance->startTimesProcessor );
    stopTimeProcessor = (double) clock() / (double) CLOCKS_PER_SEC;
   
   gettimeofday( &t, &tz );
   startTimeWallClock = hmVectorDoublePopBack( &distance->startTimesWallClock );
    stopTimeWallClock = (double) t.tv_sec + 1e-6*(double) t.tv_usec;

   for( i = 0; i < nTabs; i++ ) printf( "\t" );
   printf( "%s\n", label );
   for( i = 0; i < nTabs; i++ ) printf( "\t" );
   printf( "--------------------------------------------\n" );
   for( i = 0; i < nTabs; i++ ) printf( "\t" );
   printf( " processor time: %f seconds\n", stopTimeProcessor-startTimeProcessor );
   for( i = 0; i < nTabs; i++ ) printf( "\t" );
   printf( "wall-clock time: %f seconds\n", stopTimeWallClock-startTimeWallClock );
   printf( "\n" );
}

