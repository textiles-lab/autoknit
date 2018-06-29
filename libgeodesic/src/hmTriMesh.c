#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "hmTriMesh.h"
#include "hmVec3.h"
#include "hmUtility.h"

void hmTriMeshInitialize( hmTriMesh* mesh )
{
   mesh->nVertices = 0;
   mesh->vertices  = NULL;
   mesh->texCoords = NULL;

   mesh->nFaces = 0;
   mesh->faces  = NULL;

   mesh->referenced = 0;

   mesh->onBoundary      = NULL;
   mesh->vertexNeighbors = NULL;
   mesh->vertexAreas     = NULL;
   mesh->weights         = NULL;
   mesh->edgeNormals     = NULL;
   mesh->weightedEdges   = NULL;
}

void hmTriMeshDestroy( hmTriMesh* mesh )
{
   if( !mesh->referenced )
   {
      if( mesh->vertices != NULL )
      {
         free( mesh->vertices );
      }

      if( mesh->faces != NULL )
      {
         free( mesh->faces );
      }
   }

   hmDestroy( mesh->onBoundary );
   hmDestroy( mesh->vertexNeighbors );
   hmDestroy( mesh->vertexAreas );
   hmDestroy( mesh->weights );
   hmDestroy( mesh->edgeNormals );
   hmDestroy( mesh->weightedEdges );
   hmDestroy( mesh->texCoords );

   mesh->vertices = NULL;
   mesh->texCoords = NULL;
   mesh->faces = NULL;
   mesh->referenced = 0;
}

void hmTriMeshCopy( hmTriMesh* mesh1, const hmTriMesh* mesh2 )
{
   hmTriMeshDestroy( mesh1 );
   hmTriMeshCopyData( mesh1, mesh2->nVertices, mesh2->vertices,
                             mesh2->nFaces,    mesh2->faces );

   mesh1->texCoords = malloc( mesh2->nVertices*2 * sizeof(double) );
   memcpy( mesh1->texCoords, mesh2->texCoords, mesh2->nVertices*2 * sizeof(double) );
}

void hmTriMeshCopyData( hmTriMesh* mesh,
                        size_t nVertices, const double* vertices,
                        size_t nFaces,    const size_t* faces )
{
   hmTriMeshDestroy( mesh );

   mesh->nVertices = nVertices;
   mesh->vertices = malloc( nVertices*3 * sizeof(double) );
   memcpy( mesh->vertices, vertices, nVertices*3 * sizeof(double) );

   mesh->nFaces = nFaces;
   mesh->faces = malloc( nFaces*3 * sizeof(size_t) );
   memcpy( mesh->faces, faces, nFaces*3 * sizeof(size_t) );

   mesh->referenced = 0;

   /* initialize texture coords to zero */
   mesh->texCoords = malloc( nVertices*2 * sizeof(double) );
   hmClearArrayDouble( mesh->texCoords, nVertices*2, 0. );
}

void hmTriMeshReferenceData( hmTriMesh* mesh,
                             size_t nVertices, double* vertices,
                             size_t nFaces,    size_t* faces )
{
   hmTriMeshDestroy( mesh );

   mesh->nVertices = nVertices;
   mesh->vertices = vertices;

   mesh->nFaces = nFaces;
   mesh->faces = faces;

   mesh->referenced = 1;

   /* initialize texture coords to zero */
   mesh->texCoords = malloc( nVertices*2 * sizeof(double) );
   hmClearArrayDouble( mesh->texCoords, nVertices*2, 0. );
}

#ifdef __linux__
char *fgetln(FILE *stream, size_t *len) {
	static char buffer[2000];
	char *line = fgets(buffer, sizeof(buffer), stream);
	if (line == NULL) {
		*len = 0;
		return NULL;
	}
	*len = strlen(line);
	return line;
}
#endif //__linux__

void hmTriMeshReadOBJ( hmTriMesh* mesh, const char* filename )
{
   FILE* in;
   char* line;
   char token[32];
   size_t size;
   double* v;
   size_t* f;
   unsigned long I0, I1, I2;

   hmTriMeshDestroy( mesh );
   hmTriMeshInitialize( mesh );

   if( !( in = fopen( filename, "r" )))
   {
      fprintf( stderr, "Error: could not read from file %s\n", filename );
      exit( 1 );
   }

   /* count the number of vertices, faces */
   while( !feof( in ))
   {
      line = fgetln( in, &size );
      if( line == NULL ) continue;

      sscanf( line, "%s", token );

      if( !strcmp( token, "v" ))
      {
         mesh->nVertices++;
      }
      else if( !strcmp( token, "f" ))
      {
         mesh->nFaces++;
      }
   }

   /* allocate storage */
   mesh->vertices  = malloc( mesh->nVertices*3 * sizeof(double) );
   mesh->texCoords = malloc( mesh->nVertices*2 * sizeof(double) );
   mesh->faces     = malloc(    mesh->nFaces*3 * sizeof(size_t) );

   /* read the mesh data */
   rewind( in );
   v = mesh->vertices;
   f = mesh->faces;
   while( !feof( in ))
   {
      line = fgetln( in, &size );
      if( line == NULL ) continue;

      sscanf( line, "%s", token );

      if( !strcmp( token, "v" ))
      {
         sscanf( line, "%*s %lf %lf %lf", &v[0], &v[1], &v[2] );
         v += 3;
      }
      else if( !strcmp( token, "f" ))
      {
         /* try reading triangle vertex indices, in several possible formats */
         if( 3 == sscanf( line, "%*s %lu %lu %lu", &I0, &I1, &I2 ) ||
             3 == sscanf( line, "%*s %lu/%*d %lu/%*d %lu/%*d", &I0, &I1, &I2 ) ||
             3 == sscanf( line, "%*s %lu/%*d/%*d %lu/%*d/%*d %lu/%*d/%*d", &I0, &I1, &I2 ) ||
             3 == sscanf( line, "%*s %lu//%*d %lu//%*d %lu//%*d", &I0, &I1, &I2 ) ||
             3 == sscanf( line, "%*s %lu// %lu// %lu//", &I0, &I1, &I2 ))
         {
            /* change 1-based indices to 0-based indices */
            f[0] = I0-1;
            f[1] = I1-1;
            f[2] = I2-1;
         }
         else
         {
            fprintf( stderr, "Error: could not parse face line %s\n", line );
            fprintf( stderr, "(in file %s)\n", filename );
            exit( 1 );
         }

         f += 3;
      }
   }

   fclose( in );
}

void hmTriMeshWriteOBJ( const hmTriMesh* mesh, const char* filename )
{
   FILE *out;
   double* v = mesh->vertices;
   double* vt = mesh->texCoords;
   size_t* f = mesh->faces;
   size_t i;
   unsigned long I0, I1, I2;

   if( !( out = fopen( filename, "w" )))
   {
      fprintf( stderr, "Warning: could not write to file %s\n", filename );
      return;
   }

   for( i = 0; i < mesh->nVertices; i++ )
   {
      fprintf( out, "v %.9f %.9f %.9f\n", v[0], v[1], v[2] );
      v += 3;
   }

   for( i = 0; i < mesh->nVertices; i++ )
   {
      fprintf( out, "vt %.9f %.9f\n", vt[0], vt[1] );
      vt += 2;
   }

   for( i = 0; i < mesh->nFaces; i++ )
   {
      I0 = 1+f[0];
      I1 = 1+f[1];
      I2 = 1+f[2];
      fprintf( out, "f %lu/%lu %lu/%lu %lu/%lu\n", I0, I0, I1, I1, I2, I2 );
      f += 3;
   }
}

void hmTriMeshBuild( hmTriMesh* mesh )
{
   hmTriMeshComputeWeights( mesh );
   hmTriMeshCacheEdgeData( mesh );
   hmTriMeshComputeVertexAreas( mesh );
   hmTriMeshFindOneRings( mesh );
}

void hmTriMeshFindOneRings( hmTriMesh* mesh )
{
   /* The 1-ring neighbors of each vertex are found by iterating over all triangles and keeping a
      list of which vertices appear next to which other vertices in each triangle.  Since vertices
      may share more than one face with each of their neighbors, this list will be redundant --
      the latter part of the routine extracts a list of unique neighbors, and also detects which
      vertices are on the surface boundary by checking whether each neighbor appears twice or just
      once.  Meanwhile, Laplacian edge weights are associated with each neighboring vertex in order
      to ease construction of the Laplace and heat flow matrices in hmTriDistanceBuildMatrices(). */

   size_t i, j;

   size_t nVertices = mesh->nVertices;
   size_t nFaces = mesh->nFaces;

   size_t* facesBegin = mesh->faces;
   size_t* facesEnd = facesBegin + 3*nFaces;
   size_t* f; /* current face */

   double* w = mesh->weights;

   hmVectorPairSizeTDouble *neighbors; /* temporary, redundant list of vertex neighbors */
   hmVectorPairSizeTDouble *uniqueNeighbors; /* final list of unique vertex neighbors */
   hmPairSizeTDouble neighbor; /* used to construct a record of the current neighbor */
   size_t lastNeighborIndex;
   size_t count; /* number of times a given neighbor appears */

   /* allocate a list of redundant neighbors for each vertex */
   neighbors = malloc( nVertices * sizeof( hmVectorPairSizeTDouble ));
   for( i = 0; i < nVertices; i++ )
   {
      hmVectorPairSizeTDoubleInitialize( &neighbors[i] );
   }

   /* allocate a list of unique neighbors for each vertex */
   hmDestroy( mesh->vertexNeighbors );
   mesh->vertexNeighbors = malloc( nVertices * sizeof(hmVectorPairSizeTDouble) );
   uniqueNeighbors = mesh->vertexNeighbors; /* short name */
   for( i = 0; i < nVertices; i++ )
   {
      hmVectorPairSizeTDoubleInitialize( &uniqueNeighbors[i] );
   }

   /* allocate an array of flags for boundary vertices */
   hmDestroy( mesh->onBoundary );
   mesh->onBoundary = malloc( nVertices * sizeof(char) );

   /* iterate over triangles */
   for( f = facesBegin; f != facesEnd; f+=3 )
   {
      /* for each triangle corner, append the two other corners and the associated weights to neighbor list */
      neighbor.n = f[1]; neighbor.x = w[2]; hmVectorPairSizeTDoublePushBack( &neighbors[ f[0] ], neighbor );
      neighbor.n = f[2]; neighbor.x = w[1]; hmVectorPairSizeTDoublePushBack( &neighbors[ f[0] ], neighbor );

      neighbor.n = f[2]; neighbor.x = w[0]; hmVectorPairSizeTDoublePushBack( &neighbors[ f[1] ], neighbor );
      neighbor.n = f[0]; neighbor.x = w[2]; hmVectorPairSizeTDoublePushBack( &neighbors[ f[1] ], neighbor );

      neighbor.n = f[0]; neighbor.x = w[1]; hmVectorPairSizeTDoublePushBack( &neighbors[ f[2] ], neighbor );
      neighbor.n = f[1]; neighbor.x = w[0]; hmVectorPairSizeTDoublePushBack( &neighbors[ f[2] ], neighbor );

      w += 3;
   }

   /* iterate over vertices */
   for( i = 0; i < nVertices; i++ )
   {
      /* sort neighbor list by index */
      hmVectorPairSizeTDoubleSort( &neighbors[i] );

      /* initially flag as an interior vertex */
      mesh->onBoundary[i] = 0;

      /* extract unique elements from neighbor list, summing weights */
      hmVectorPairSizeTDoubleResize( &uniqueNeighbors[i], 0 );
      lastNeighborIndex = -1;
      count = 0;
      for( j = 0; j < neighbors[i].size; j++ )
      {
         /* if we come across a new neighbor, add it to the list of unique neighbors */
         if( neighbors[i].entries[j].n != lastNeighborIndex )
         {
            /* if we encountered the previous neighbor only
             * once, this vertex must be on the surface boundary */
            if( count == 1 )
            {
               mesh->onBoundary[i] = 1;
            }
            count = 1;

            hmVectorPairSizeTDoublePushBack( &uniqueNeighbors[i], neighbors[i].entries[j] );
            lastNeighborIndex = neighbors[i].entries[j].n;
         }
         else
         {
            /* since we've seen this neighbor before, just accumulate its weight */
            uniqueNeighbors[i].entries[ uniqueNeighbors[i].size-1 ].x += neighbors[i].entries[j].x;
            count++;
         }
      }

      /* if the final neighbor was encountered only once, this is a boundary vertex */
      if( count == 1 )
      {
         mesh->onBoundary[i] = 1;
      }
   }

   /* deallocate storage for redundant neighbor list */
   for( i = 0; i < nVertices; i++ )
   {
      hmVectorPairSizeTDoubleDestroy( &neighbors[i] );
   }
   free( neighbors );
}


void hmTriMeshComputeVertexAreas( hmTriMesh* mesh )
{
   size_t nVertices = mesh->nVertices;
   size_t nFaces = mesh->nFaces;

   const size_t* facesBegin = mesh->faces;
   const size_t* facesEnd = facesBegin + 3*nFaces;
   const size_t* f;

   double* vertexAreas;
   double* vertices = mesh->vertices;
   double *p0, *p1, *p2;
   hmVec3 u, v, w;

   double A;

   /* initialize vertex areas to zero */
   hmDestroy( mesh->vertexAreas );
   mesh->vertexAreas = malloc( nVertices * sizeof( double ));
   hmClearArrayDouble( mesh->vertexAreas, nVertices, 0. );
   vertexAreas = mesh->vertexAreas;

   /* iterate over faces */
   for( f = facesBegin; f != facesEnd; f += 3 )
   {
      /* get vertex coordinates p0, p1, p2 */
      p0 = &vertices[ f[0]*3 ];
      p1 = &vertices[ f[1]*3 ];
      p2 = &vertices[ f[2]*3 ];

      /* compute (one-third of) the triangle area A = |u x v| */
      hmVec3Sub( u, p1, p0 );
      hmVec3Sub( v, p2, p0 );
      hmVec3Cross( w, u, v );
      A = hmVec3Norm( w ) / 6.;

      /* add contribution to each of the three corner vertices */
      vertexAreas[ f[0] ] += A;
      vertexAreas[ f[1] ] += A;
      vertexAreas[ f[2] ] += A;
   }
}

void hmTriMeshComputeWeights( hmTriMesh* mesh )
{
   /* array counters */
   int k;
   int j0, j1, j2;

   /* local data handles */
   size_t nFaces = mesh->nFaces;
   const size_t* facesBegin = mesh->faces;
   const size_t* facesEnd = facesBegin + 3*nFaces;
   const size_t* f;
   double* w; /* current weights */
   double* vertices = mesh->vertices;

   /* current triangle data */
   double* p[3]; /* vertex positions */
   hmVec3 u, v; /* edge vectors */
   hmVec3 N; /* triangle normal */
   double uvSinTheta, uvCosTheta;

   hmDestroy( mesh->weights );
   mesh->weights = malloc( 3*nFaces * sizeof(double) );
   w = mesh->weights;

   /* iterate over triangles */
   for( f = facesBegin; f != facesEnd; f += 3 )
   {
      /* get vertex coordinates */
      p[0] = &vertices[ f[0]*3 ];
      p[1] = &vertices[ f[1]*3 ];
      p[2] = &vertices[ f[2]*3 ];

      /* iterate over triangle corners */
      for( k = 0; k < 3; k++ )
      {
         /* get outgoing edge vectors u, v at current corner */
         j0 = (0+k) % 3;
         j1 = (1+k) % 3;
         j2 = (2+k) % 3;
         hmVec3Sub( u, p[j1], p[j0] );
         hmVec3Sub( v, p[j2], p[j0] );

         /* compute (one-half of) the cotangent weight */
         hmVec3Cross( N, u, v );
         uvSinTheta = hmVec3Norm( N );
         uvCosTheta = hmVec3Dot( u, v );
         w[k] = .5 * uvCosTheta / uvSinTheta;
      }

      /* move to the weights in the next triangle */
      w += 3;
   }
}

void hmTriMeshCacheEdgeData( hmTriMesh* mesh )
{
   /* array counters */
   size_t i;

   /* local data handles */
   size_t nFaces = mesh->nFaces;
   const size_t*         f = mesh->faces;
   const double*         w = mesh->weights;
         double*  vertices = mesh->vertices;

   /* current triangle data */
   double *p0, *p1, *p2; /* vertex coordinates */
   double *e0, *e1, *e2; /* edge vectors */
   double *t0, *t1, *t2; /* rotated edge vectors */
   hmVec3 N; /* triangle normal */

   /* allocate storage for edge data */
   hmDestroy( mesh->edgeNormals );
   hmDestroy( mesh->weightedEdges );
   mesh->edgeNormals   = malloc( 9*nFaces * sizeof( double ));
   mesh->weightedEdges = malloc( 9*nFaces * sizeof( double ));

   /* get pointers to first three edge normals */
   t0 = &mesh->edgeNormals[ 0 ];
   t1 = &mesh->edgeNormals[ 3 ];
   t2 = &mesh->edgeNormals[ 6 ];

   /* get pointers to first three weighted edges */
   e0 = &mesh->weightedEdges[ 0 ];
   e1 = &mesh->weightedEdges[ 3 ];
   e2 = &mesh->weightedEdges[ 6 ];

   /* iterate over faces */
   for( i = 0; i < mesh->nFaces; i++ )
   {
      /* get vertex coordinates */
      p0 = &vertices[ f[0]*3 ];
      p1 = &vertices[ f[1]*3 ];
      p2 = &vertices[ f[2]*3 ];

      /* compute edge vectors */
      hmVec3Sub( e0, p2, p1 );
      hmVec3Sub( e1, p0, p2 );
      hmVec3Sub( e2, p1, p0 );

      /* compute triangle normal */
      hmVec3Cross( N, e0, e1 );

      /* compute rotated edge vectors */
      hmVec3Cross( t0, N, e0 );
      hmVec3Cross( t1, N, e1 );
      hmVec3Cross( t2, N, e2 );

      /* scale edge vectors by cotangent weights at opposing corners */
      hmVec3Scale( e0, w[0] );
      hmVec3Scale( e1, w[1] );
      hmVec3Scale( e2, w[2] );

      /* move to the next face */
      f += 3;
      w += 3;
      t0 += 9; t1 += 9; t2 += 9;
      e0 += 9; e1 += 9; e2 += 9;
   }
}

double hmTriMeshL2Distance( hmTriMesh* mesh, double* phi1, double* phi2 )
{
   size_t i;
   double sum = 0.;

   hmTriMeshComputeVertexAreas( mesh );

   for( i = 0; i < mesh->nVertices; i++ )
   {
      sum += mesh->vertexAreas[i] * hmSquare( phi1[i] - phi2[i] );
   }

   return sqrt( sum );
}

double hmTriMeshLInfinityDistance( hmTriMesh* mesh, double* phi1, double* phi2 )
{
   size_t i;
   double maxDifference = 0.;

   for( i = 0; i < mesh->nVertices; i++ )
   {
      maxDifference = hmMaxDouble( maxDifference, fabs( phi1[i]-phi2[i] ));
   }

   return maxDifference;
}

double hmTriMeshMeanRelativeError( hmTriMesh* mesh, double* phi1, double* phi2 )
{
   size_t i;
   double mean = 0.;

   for( i = 0; i < mesh->nVertices; i++ )
   {
      if( phi1[i] != 0. )
      {
         mean += fabs( (phi2[i]-phi1[i]) / phi1[i] );
      }
   }

   mean /= (double) mesh->nVertices;

   return mean;
}

double hmTriMeshDiameter( hmTriMesh* mesh )
{
   size_t nVertices = mesh->nVertices;
   size_t nFaces = mesh->nFaces;

   const size_t* facesBegin = mesh->faces;
   const size_t* facesEnd = facesBegin + 3*nFaces;
   const size_t* f;

   const double* verticesBegin = mesh->vertices;
   const double* verticesEnd = verticesBegin + 3*nVertices;
   const double* v;

   double* vertices = mesh->vertices;
   double *p0, *p1, *p2;
   hmVec3 u1, u2, w;
   hmVec3 barycenter;

   double radius;
   double triangleArea, surfaceArea;

   hmVec3 centerOfMass;
   hmVec3Set( centerOfMass, 0., 0., 0. );

   /* iterate over faces */
   surfaceArea = 0.;
   for( f = facesBegin; f != facesEnd; f += 3 )
   {
      /* get vertex coordinates p0, p1, p2 */
      p0 = &vertices[ f[0]*3 ];
      p1 = &vertices[ f[1]*3 ];
      p2 = &vertices[ f[2]*3 ];

      /* compute barycenter (p0+p1+p2)/3 */
      hmVec3Set( barycenter, 0., 0., 0. );
      hmVec3Add( barycenter, barycenter, p0 );
      hmVec3Add( barycenter, barycenter, p1 );
      hmVec3Add( barycenter, barycenter, p2 );
      hmVec3Scale( barycenter, 1./3. );

      /* compute the triangle area A = |u x v| */
      hmVec3Sub( u1, p1, p0 );
      hmVec3Sub( u2, p2, p0 );
      hmVec3Cross( w, u1, u2 );
      triangleArea = hmVec3Norm( w ) / 2.;

      /* add contribution to total area */
      surfaceArea += triangleArea;

      /* add contribution to center of mass */
      hmVec3Scale( barycenter, triangleArea );
      hmVec3Inc( centerOfMass, barycenter );
   }

   hmVec3Scale( centerOfMass, 1./surfaceArea );

   /* iterate over vertices */
   radius = 0.;
   for( v = verticesBegin; v != verticesEnd; v+=3 )
   {
      hmVec3Sub( w, v, centerOfMass );
      radius = hmMaxDouble( radius, hmVec3Norm( w ));
   }

   return 2.*radius;
}
