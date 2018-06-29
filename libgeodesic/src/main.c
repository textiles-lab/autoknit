#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "hmTriDistance.h"
#include "hmContext.h"
#include "hmUtility.h"
#include "hmVectorSizeT.h"

int parseCommandLineArguments( int argc, char** argv, char** inputMesh, char** inputSources, char** outputMeshPrefix, char** outputDistancePrefix, char** inputReference, double* smoothness, double* boundaryConditions, char* verbose );
void printHelp( int argc, char** argv );
void printVersion( int argc, char** argv );
void readMesh( hmTriDistance* distance, hmTriMesh* mesh, const char* filename );
void readSources( int* nSourceSets, hmVectorSizeT** sourceSets, const char* filename );
void readReferenceValues( hmTriDistance* distance, int* nReferenceColumns, double*** referenceValues, char*** referenceNames, const char* filename );
void destroySourceSets( int nSourceSets, hmVectorSizeT** sourceSets );
void destroyReferenceValues( int nReferenceColumns, double*** referenceValues, char*** referenceNames );
void writeMesh( hmTriDistance* distance, const char* filename, int index );
void writeDistances( hmTriDistance* distance, const char* prefix, int index );
void setSources( hmTriDistance* distance, hmVectorSizeT* sourceSets, int index );
void compareDistance( hmTriDistance* distance, int nReferenceColumns, double** referenceValues, char** referenceNames, const char* prefix );

int main( int argc, char **argv )
{
   /* counters */
   size_t i;

   /* main data */
   hmContext context;
   hmTriMesh surface;
   hmTriDistance distance;

   int nSourceSets = 0;
   hmVectorSizeT* sourceSets = NULL;

   int nReferenceColumns = 0;
   double** referenceValues = NULL;
   char** referenceNames = NULL;

   /* parameters parsed from command-line arguments */
   char* inputMesh = NULL;
   char* inputSources = NULL;
   char* inputReference = NULL;
   char* outputMeshPrefix = NULL;
   char* outputDistancePrefix = NULL;
   double smoothness = -1.;
   double boundaryConditions = -1.;
   char verbose = 0;

   /* parse command-line arguments */
   if( !parseCommandLineArguments( argc, argv,
                                   &inputMesh,
                                   &inputSources,
                                   &outputMeshPrefix,
                                   &outputDistancePrefix,
                                   &inputReference,
                                   &smoothness,
                                   &boundaryConditions,
                                   &verbose ))
   {
      /* if something went wrong, bail out! */
      return 1;
   }

   /* initialize data */
   hmContextInitialize( &context );
   hmTriMeshInitialize( &surface );
   hmTriDistanceInitialize( &distance );

   /* read surface */
   readMesh( &distance, &surface, inputMesh );
   printf( "geodesic: %s\n", inputMesh );

   /* read source set(s) */
   readSources( &nSourceSets, &sourceSets, inputSources );
   if( inputReference != NULL )
   {
      if( nSourceSets != 1 )
      {
         fprintf( stderr, "Warning: comparison with reference values can be performed only\n" );
         fprintf( stderr, "for a single source set.  Reference values will be ignored.\n" );
      }
      else
      {
         /* read reference values */
         readReferenceValues( &distance, &nReferenceColumns, &referenceValues, &referenceNames, inputReference );
      }
   }

   /* set time for heat flow */
   hmTriDistanceEstimateTime( &distance );
   if( smoothness > 0. )
   {
      distance.time *= smoothness;
   }

   /* specify boundary conditions */
   if( boundaryConditions > 0. )
   {
      hmTriDistanceSetBoundaryConditions( &distance, boundaryConditions );
   }

   /* specify verbosity */
   distance.verbose = verbose;

   /* compute distance */
   hmTriDistanceBuild( &distance );
   for( i = 0; i < nSourceSets; i++ )
   {
      /* set current source set */
      setSources( &distance, sourceSets, i );

      /* udpate distance function */
      hmTriDistanceUpdate( &distance );

      /* write solution as plain text file */
      if( outputDistancePrefix != NULL )
      {
         writeDistances( &distance, outputDistancePrefix, i );
      }

      /* write solution as vertex coordinates */
      if( outputMeshPrefix != NULL )
      {
         writeMesh( &distance, outputMeshPrefix, i );
      }
   }

   /* compare with reference distances */
   if( inputReference != NULL && nSourceSets == 1 )
   {
      compareDistance( &distance,
                       nReferenceColumns,
                       referenceValues,
                       referenceNames,
                       outputMeshPrefix );
   }

   /* deallocate data */
   hmTriMeshDestroy( &surface );
   hmTriDistanceDestroy( &distance );
   hmContextDestroy( &context );
   destroySourceSets( nSourceSets, &sourceSets );
   destroyReferenceValues( nReferenceColumns, &referenceValues, &referenceNames );

   return 0;
}

int parseCommandLineArguments( int argc, char** argv,
                               char** inputMesh,
                               char** inputSources,
                               char** outputMeshPrefix,
                               char** outputDistancePrefix,
                               char** inputReference,
                               double* smoothness,
                               double* boundaryConditions,
                               char* verbose )
{
   int i;

   for( i = 1; i < argc; i++ )
   {
      if( !strcmp( argv[i], "-h" ) ||
          !strcmp( argv[i], "-help" )) /* help */
      {
         printHelp( argc, argv );
         return 0;
      }
      else if( !strcmp( argv[i], "-v" ) ||
               !strcmp( argv[i], "-version" )) /* version */
      {
         printVersion( argc, argv );
         return 0;
      }
      else if( !strcmp( argv[i], "-w" ) ||
               !strcmp( argv[i], "-verbose" )) /* verbose */
      {
         *verbose = 1;
      }
      else if( !strcmp( argv[i], "-i" ) ||
               !strcmp( argv[i], "-input" )) /* input mesh */
      {
         if( i == argc-1 || argv[i+1][0] == '-' )
         {
            fprintf( stderr, "-i: must specify input mesh!\n" );
            return 0;
         }
         *inputMesh = argv[i+1];
         i++;
      }
      else if( !strcmp( argv[i], "-s" ) ||
               !strcmp( argv[i], "-source" )) /* input sources */
      {
         if( i == argc-1 || argv[i+1][0] == '-' )
         {
            fprintf( stderr, "-s: must specify sources!\n" );
            return 0;
         }
         *inputSources = argv[i+1];
         i++;
      }
      else if( !strcmp( argv[i], "-o" ) ||
               !strcmp( argv[i], "-output" )) /* output mesh prefix */
      {
         if( i == argc-1 || argv[i+1][0] == '-' )
         {
            fprintf( stderr, "-o: must specify prefix for output meshes!\n" );
            return 0;
         }
         *outputMeshPrefix = argv[i+1];
         i++;
      }
      else if( !strcmp( argv[i], "-d" ) ||
               !strcmp( argv[i], "-distance" )) /* output distances prefix */
      {
         if( i == argc-1 || argv[i+1][0] == '-' )
         {
            fprintf( stderr, "-d: must specify prefix for output distance files!\n" );
            return 0;
         }
         *outputDistancePrefix = argv[i+1];
         i++;
      }
      else if( !strcmp( argv[i], "-r" ) ||
               !strcmp( argv[i], "-reference" )) /* reference distances */
      {
         if( i == argc-1 || argv[i+1][0] == '-' )
         {
            fprintf( stderr, "-r: must specify reference distances!\n" );
            return 0;
         }
         *inputReference = argv[i+1];
         i++;
      }
      else if( !strcmp( argv[i], "-m" ) ||
               !strcmp( argv[i], "-smoothness" )) /* smoothness */
      {
         if( i == argc-1 || argv[i+1][0] == '-' )
         {
            fprintf( stderr, "-m: must specify positive smoothness value!\n" );
            return 0;
         }
         *smoothness = strtod( argv[i+1], NULL );
         i++;
         if( *smoothness < 0. || isnan(*smoothness) || isinf(*smoothness) )
         {
            fprintf( stderr, "-m: must specify positive smoothness value!\n" );
            return 0;
         }
         if( *smoothness < 1. )
         {
            fprintf( stderr, "Warning: smoothness parameter below suggested value (i.e.,\n" );
            fprintf( stderr, "         less than 1.0) may result in numerical inaccuracy.\n" );
         }
      }
      else if( !strcmp( argv[i], "-b" ) ||
               !strcmp( argv[i], "-boundary" )) /* boundary conditions */
      {
         if( i == argc-1 || argv[i+1][0] == '-' )
         {
            fprintf( stderr, "-b: must specify a value between 0 and 1!\n" );
            return 0;
         }
         *boundaryConditions = strtod( argv[i+1], NULL );
         i++;
         if( *boundaryConditions < 0. || *boundaryConditions > 1. )
         {
            fprintf( stderr, "-b: must specify a value between 0 and 1!\n" );
            return 0;
         }
      }
      else if( argv[i][0] == '-' )
      {
         fprintf( stderr, "Unrecognized option \"%s\"!\n", argv[i] );
         return 0;
      }
      else
      {
         fprintf( stderr, "usage: %s [options] \n", argv[0] );
         fprintf( stderr, "(-help for more information)\n\n" );
         return 0;
      }
   }

   if( *inputMesh == NULL )
   {
      fprintf( stderr, "Error: must specify an input mesh!\n" );
      fprintf( stderr, "(-help for more information)\n\n" );
      return 0;
   }
   if( *inputSources == NULL )
   {
      fprintf( stderr, "Error: must specify a source set!\n" );
      fprintf( stderr, "(-help for more information)\n\n" );
      return 0;
   }
   if( *outputMeshPrefix == NULL && *outputDistancePrefix == NULL )
   {
      fprintf( stderr, "Warning: no output specified -- distances will not be saved!\n" );
   }

   return 1;
}

void printHelp( int argc, char** argv )
{
   fprintf( stdout, "geodesic %s -- Fast, robust geodesic distance computation on a variety of domains.\n", hmVersionString );
   fprintf( stdout, "\n" );
   fprintf( stdout, "usage: %s [options]\n", argv[0] );
   fprintf( stdout, "\n" );
   fprintf( stdout, "options:\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "   -v\n" );
   fprintf( stdout, "   -version\n" );
   fprintf( stdout, "         Print version information.\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "   -h\n" );
   fprintf( stdout, "   -help\n" );
   fprintf( stdout, "         Print this help message.\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "   -w\n" );
   fprintf( stdout, "   -verbose\n" );
   fprintf( stdout, "         Verbose output (timing information, etc.).\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "   -i     [inputMesh]\n" );
   fprintf( stdout, "   -input [inputMesh]\n" );
   fprintf( stdout, "         Specify the input mesh.  Guarantees on the output can be made only for\n" );
   fprintf( stdout, "         meshes that are triangulated 2-manifolds, i.e., every edge is contained\n" );
   fprintf( stdout, "         in at most two triangles and every vertex is contained in a sequence\n" );
   fprintf( stdout, "         of triangles where any two consecutive triangles share exactly one\n" );
   fprintf( stdout, "         edge.  Triangles in the input must be consistently oriented, i.e., the\n" );
   fprintf( stdout, "         three vertices should be specified in counter-clockwise order relative\n" );
   fprintf( stdout, "         to the outward-facing normal.  The file type will be automatically\n" );
   fprintf( stdout, "         determined from the file suffix -- currently supported formats include:\n" );
   fprintf( stdout, "         Wavefront OBJ.\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "   -s      [sources]\n" );
   fprintf( stdout, "   -source [sources]\n" );
   fprintf( stdout, "         Specify the source set as a plain ASCII text file -- multiple source\n" );
   fprintf( stdout, "         sets may be specified.  The first line should be an integer specifying\n" );
   fprintf( stdout, "         the number of sets; each set is specified by the number of points in\n" );
   fprintf( stdout, "         the set followed by a list of 1-based indices into the vertex set.  For\n" );
   fprintf( stdout, "         instance, to compute the distance functions relative to the two sets\n" );
   fprintf( stdout, "         {1,2,3,4} and {100}, the input would look like\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "            2\n" );
   fprintf( stdout, "            4 1 2 3 4\n" );
   fprintf( stdout, "            1 100\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "   -d        [distancePrefix]\n" );
   fprintf( stdout, "   -distance [distancePrefix]\n" );
   fprintf( stdout, "         Specify the output distance values as a plain ASCII text file -- the\n" );
   fprintf( stdout, "         command-line argument specifies the file prefix.  If multiple source\n" );
   fprintf( stdout, "         sets are specified, distances will be written as prefix.1.dist,\n" );
   fprintf( stdout, "         prefix.2.dist, etc.  Each file contains a list of distance values\n" );
   fprintf( stdout, "         (one per line) in the same order as the input mesh vertices.\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "   -o      [outputMeshPrefix]\n" );
   fprintf( stdout, "   -output [outputMeshPrefix]\n" );
   fprintf( stdout, "         Write an output mesh in Wavefront OBJ format, where files are named\n" );
   fprintf( stdout, "         using the specified prefix and numbered according to the source set\n" );
   fprintf( stdout, "         (prefix.1.obj, prefix.2.obj, etc.).  Distance values are stored in\n" );
   fprintf( stdout, "         vertex coordinates, i.e., each \"vt\" line contains two copies of the\n" );
   fprintf( stdout, "         distance value associated with that vertex:\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "            vt [d] [d]\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "         This format is convenient for visualizing isolines of the distance\n" );
   fprintf( stdout, "         function by mapping a repeating pattern onto the surface (several\n" );
   fprintf( stdout, "         such patterns have been provided in libgeodesic/data/textures).\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "   -r         [referenceValues]\n" );
   fprintf( stdout, "   -reference [referenceValues]\n" );
   fprintf( stdout, "         Compare with specified reference values (for instance, one might be\n" );
   fprintf( stdout, "         interested in comparing against the exact geodesic distances from the\n" );
   fprintf( stdout, "         corresponding smooth surface, or the exact polyhedral distance as\n" );
   fprintf( stdout, "         computed by the algorithm of Mitchell et al.)  This option can be\n" );
   fprintf( stdout, "         useful for testing the accuracy of the heat method.  If the file\n" );
   fprintf( stdout, "         contains multiple columns, values in the remaining columns will also be\n" );
   fprintf( stdout, "         compared against reference values.  Input should be a plain ASCII text\n" );
   fprintf( stdout, "         file.  The first line specifies the number of columns followed by the\n" );
   fprintf( stdout, "         number of distance values (which must equal the number of mesh\n" );
   fprintf( stdout, "         vertices), the second line specifies a name for each column (white\n" );
   fprintf( stdout, "         space delimited; each name should be no longer than 128 characters) and\n" );
   fprintf( stdout, "         the remainder of the file specifies distance values as floating point\n" );
   fprintf( stdout, "         numbers in the same order as the input mesh vertices.  For example, a\n" );
   fprintf( stdout, "         reference file with two columns might look like\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "            2 [nVertices]\n" );
   fprintf( stdout, "            reference Dijkstra\n" );
   fprintf( stdout, "            [referenceDistance0] [dijstraDistance0]\n" );
   fprintf( stdout, "            [referenceDistance1] [dijstraDistance1]\n" );
   fprintf( stdout, "            ...\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "         (Note that these values should correspond to the same source set as\n" );
   fprintf( stdout, "         the one specified using the -s option -- no automatic consistency\n" );
   fprintf( stdout, "         checking is performed!)  Once the distance has been computed, a\n" );
   fprintf( stdout, "         numerical comparison will be printed to standard output.  E.g.,\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "            Error relative to reference\n" );
   fprintf( stdout, "            ===========================\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "            libgeodesic\n" );
   fprintf( stdout, "            -----------\n" );
   fprintf( stdout, "               MEAN: [mean error]\n" );
   fprintf( stdout, "                MAX: [maximum error]\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "            Dijkstra\n" );
   fprintf( stdout, "            --------\n" );
   fprintf( stdout, "               MEAN: [mean error]\n" );
   fprintf( stdout, "                MAX: [maximum error]\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "         The MEAN and MAX values correspond to the mean relative and maximum\n" );
   fprintf( stdout, "         absolute error at each point, using the reference values as a baseline.\n" );
   fprintf( stdout, "         for comparison.  MAX error is normalized by mesh diameter.  The first\n" );
   fprintf( stdout, "         set of error values correspond to those computed by this executable\n" );
   fprintf( stdout, "         via libgeodesic; remaining values correspond to additional columns in\n" );
   fprintf( stdout, "         the input.  Note that error values are meaningful only with respect to\n" );
   fprintf( stdout, "         the reference values!  Please use this feature judiciously when making\n" );
   fprintf( stdout, "         comparisons.\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "         If the -o option is specified, distance values stored with the output\n" );
   fprintf( stdout, "         meshes will include the heat method solution (.0), the reference values\n" );
   fprintf( stdout, "         (.1), and any remaining columns from the reference input (.2-.n) in\n" );
   fprintf( stdout, "         that order.\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "   -m          [smoothnessValue]\n" );
   fprintf( stdout, "   -smoothness [smoothnessValue]\n" );
   fprintf( stdout, "         Specify the amount of regularization as a positive floating-point\n" );
   fprintf( stdout, "         value.  This value influences the smoothness of the solution.  In\n" );
   fprintf( stdout, "         particular, if m is the value of the command line argument then the\n" );
   fprintf( stdout, "         integration time used for heat flow is t = h^2 where h is the mean\n" );
   fprintf( stdout, "         length of all edges in the input mesh.  The default value m = 1\n" );
   fprintf( stdout, "         typically yields a good approximation of the exact geodesic distance;\n" );
   fprintf( stdout, "         larger values will produce smoother approximations.  Note that values\n" );
   fprintf( stdout, "         of m below 1 may increase the accuracy of the solution but may also\n" );
   fprintf( stdout, "         introduce additional numerical artifacts.\n" );
   fprintf( stdout, "    \n" );
   fprintf( stdout, "   -b        [boundaryConditions]\n" );
   fprintf( stdout, "   -boundary [boundaryConditions]\n" );
   fprintf( stdout, "         Specify boundary conditions as a floating-point value between 0 and\n" );
   fprintf( stdout, "         1, inclusive.  A value of 0 yields pure-Neumann conditions; a value\n" );
   fprintf( stdout, "         of 1 yields pure Dirichlet conditions.  For typical usage (small\n" );
   fprintf( stdout, "         flow time or surfaces without boundary) this parameter has little\n" );
   fprintf( stdout, "         effect and optimal performance will be achieved by setting it to the\n" );
   fprintf( stdout, "         default value of zero.  For surfaces with boundary where a large flow\n" );
   fprintf( stdout, "         time is desired (i.e., smoothed geodesic distance) a value of 0.5\n" );
   fprintf( stdout, "         tends to give natural-looking behavior near the boundary.\n" );
   fprintf( stdout, "\n" );
}

void printVersion( int argc, char** argv )
{
   fprintf( stdout, "geodesic -- %s\n", hmVersionString );
   fprintf( stdout, "\n" );
   fprintf( stdout, "Copyright 2012 Keenan Crane. All rights reserved.\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "Redistribution and use in source and binary forms, with or without modification,\n" );
   fprintf( stdout, "are permitted provided that the following conditions are met:\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "1. Redistributions of source code must retain the above copyright notice, this\n" );
   fprintf( stdout, "   list of conditions and the following disclaimer.\n" );
   fprintf( stdout, "2. Redistributions in binary form must reproduce the above copyright notice,\n" );
   fprintf( stdout, "   this list of conditions and the following disclaimer in the documentation\n" );
   fprintf( stdout, "   and/or other materials provided with the distribution.\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED\n" );
   fprintf( stdout, "WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF\n" );
   fprintf( stdout, "MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT\n" );
   fprintf( stdout, "SHALL THE FREEBSD PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,\n" );
   fprintf( stdout, "INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n" );
   fprintf( stdout, "LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR\n" );
   fprintf( stdout, "PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF\n" );
   fprintf( stdout, "LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE\n" );
   fprintf( stdout, "OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF\n" );
   fprintf( stdout, "ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n" );
   fprintf( stdout, "\n" );
   fprintf( stdout, "The views and conclusions contained in the software and documentation are those\n" );
   fprintf( stdout, "of the author and should not be interpreted as representing official policies,\n" );
   fprintf( stdout, "either expressed or implied, of any other person or institution.\n" );
}

void readMesh( hmTriDistance* distance,
               hmTriMesh* mesh,
               const char* filename )
{
   hmTriMeshReadOBJ( mesh, filename );
   distance->surface = mesh;
}

void readSources( int* nSourceSets,
                  hmVectorSizeT** sourceSets,
                  const char* filename )
{
   int i;
   int setSize;
   int rval;
   FILE* in;
   hmVectorSizeT* sourceSetsBegin;
   hmVectorSizeT* sourceSetsEnd;
   hmVectorSizeT* s;

   if( !( in = fopen( filename, "r" )))
   {
      fprintf( stderr, "Error: could not read from source set file %s\n", filename );
      exit( 1 );
   }

   /* get the number of sets */
   rval = fscanf( in, "%d", nSourceSets );
   if( rval == 0 || rval == EOF )
   {
      fprintf( stderr, "Error: malformatted source set file (%s)\n", filename );
      fclose( in );
      exit( 1 );
   }

   /* allocate an array of sets */
   if( !( *sourceSets = malloc( *nSourceSets * sizeof(hmVectorDouble))))
   {
      fprintf( stderr, "Error: could not allocate memory for source sets (check input format)!\n" );
      fclose( in );
      exit( 1 );
   }

   /* read each of the sets */
   sourceSetsBegin = *sourceSets;
   sourceSetsEnd = sourceSetsBegin + *nSourceSets;
   for( s = sourceSetsBegin; s != sourceSetsEnd; s++ )
   {
      /* get the set size */
      rval = fscanf( in, "%d", &setSize );
      if( rval == 0 || rval == EOF )
      {
         fprintf( stderr, "Error: malformatted source set file (%s)\n", filename );
         fclose( in );
         exit( 1 );
      }

      /* allocate space for the set */
      hmVectorSizeTInitialize( s );
      hmVectorSizeTResize( s, setSize );

      /* read the set indices */
      for( i = 0; i < setSize; i++ )
      {
         rval = fscanf( in, "%lu", &s->entries[i] );
         if( rval == 0 || rval == EOF )
         {
            fprintf( stderr, "Error: malformatted source set file (%s)\n", filename );
            fclose( in );
            exit( 1 );
         }

         /* convert 1-based indices to 0-based indices */
         s->entries[i]--;
      }
   }

   fclose( in );
}

void setSources( hmTriDistance* distance,
                 hmVectorSizeT* sourceSets,
                 int index )
{
   size_t i, n;
   size_t nVertices = distance->surface->nVertices;
   hmVectorSizeT* S = &sourceSets[index];

   /* initialize all vertices to zero, meaning "not a source" */
   hmClearArrayDouble( distance->isSource.values, nVertices, 0. );

   /* set the specified source vertices in the current set */
   for( i = 0; i < S->size; i++ )
   {

      /* make sure the source vertex index n is valid */
      n = S->entries[i];
      if( n >= nVertices )
      {
         /* print an error message, remembering that source
          * vertices were 1-based in the input */
         fprintf( stderr, "Error: source vertices must be in the range 1-nVertices!\n" );
         exit( 1 );
      }

      /* flag the current vertex as a source */
      distance->isSource.values[n] = 1.;
   }
}

void destroySourceSets( int nSourceSets,
                        hmVectorSizeT** sourceSets )
{
   hmVectorSizeT* sourceSetsBegin;
   hmVectorSizeT* sourceSetsEnd;
   hmVectorSizeT* s;

   /* deallocate individual vectors */
   sourceSetsBegin = *sourceSets;
   sourceSetsEnd = sourceSetsBegin + nSourceSets;
   for( s = sourceSetsBegin; s != sourceSetsEnd; s++ )
   {
      hmVectorSizeTDestroy( s );
   }

   /* deallocate array of vectors */
   hmDestroy( *sourceSets );
}

void writeMesh( hmTriDistance* distance,
                const char* prefix,
                int index )
{
   int i, j;
   char filename[2048];

   /* construct filename */
   sprintf( filename, "%s.%d.obj", prefix, index );

   /* copy distance values to both components of each vertex coordinate */
   for( i = 0; i < distance->surface->nVertices; i++ )
   {
      for( j = 0; j < 2; j++ )
      {
         distance->surface->texCoords[i*2+j] = distance->distance.values[i];
      }
   }

   /* write mesh */
   hmTriMeshWriteOBJ( distance->surface, filename );
}

void writeDistances( hmTriDistance* distance, const char* prefix, int index )
{
   int i;
   char filename[2048];
   FILE* out;

   /* construct filename */
   sprintf( filename, "%s.%d.dist", prefix, index );

   /* get file handle */
   if( !( out = fopen( filename, "w" )))
   {
      fprintf( stderr, "Error: could not write to file %s\n", filename );
      return;
   }

   /* write distances */
   for( i = 0; i < distance->surface->nVertices; i++ )
   {
      fprintf( out, "%.20e\n", distance->distance.values[i] );
   }

   /* close file */
   fclose( out );
}

void readReferenceValues( hmTriDistance* distance,
                          int* nReferenceColumns,
                          double*** referenceValues,
                          char*** referenceNames,
                          const char* filename )
{
   int i, j;
   const int maxNameSize = 128;
   size_t nVertices = distance->surface->nVertices;
   int nRows;
   int rval;
   FILE* in;

   /* try opening input */
   if( !( in = fopen( filename, "r" )))
   {
      fprintf( stderr, "Error: could not read from reference value file %s\n", filename );
      exit( 1 );
   }

   /* make sure the first entry specifies the number of columns */
   rval = fscanf( in, "%d", nReferenceColumns );
   if( rval == 0 || rval == EOF )
   {
      fprintf( stderr, "Error: malformatted reference value file (%s)\n", filename );
      fprintf( stderr, "       (First entry must specify number of columns!)\n" );
      fclose( in );
      exit( 1 );
   }

   /* make sure the second entry specifies the (correct) number of distance values */
   rval = fscanf( in, "%d", &nRows );
   if( rval == 0 || rval == EOF )
   {
      fprintf( stderr, "Error: malformatted reference value file (%s)\n", filename );
      fprintf( stderr, "       (Second entry must specify number of distance values!)\n" );
      fclose( in );
      exit( 1 );
   }
   if( nRows != nVertices )
   {
      fprintf( stderr, "Error: malformatted reference value file (%s)\n", filename );
      fprintf( stderr, "       (Number of distance values must equal number of vertices!)\n" );
      fclose( in );
      exit( 1 );
   }

   /* allocate storage for the distance values and associated names */
   if( !( *referenceValues = malloc( *nReferenceColumns * sizeof(double*))) ||
       !( *referenceNames  = malloc( *nReferenceColumns * sizeof(char*))))
   {
      fprintf( stderr, "Error: could not allocate space for %d columns of\n", *nReferenceColumns );
      fprintf( stderr, "       reference values. (Malformatted input?)\n" );
      fclose( in );
      exit( 1 );
   }
   for( i = 0; i < *nReferenceColumns; i++ )
   {
      if( !((*referenceValues)[i] = malloc( nVertices * sizeof(double))) ||
          !((*referenceNames)[i] = malloc( maxNameSize * sizeof(char))))
      {
         fprintf( stderr, "Error: could not allocate space for %lu reference values.\n", nVertices );
         fprintf( stderr, "       (Malformatted input?)\n" );
         fclose( in );
         exit( 1 );
      }
   }

   /* read names */
   for( i = 0; i < *nReferenceColumns; i++ )
   {
      rval = fscanf( in, "%s", (*referenceNames)[i] );
      if( rval == 0 || rval == EOF )
      {
         fprintf( stderr, "Error: unable to read column names in reference file.\n" );
         fprintf( stderr, "       (Malformatted input?)\n" );
         fclose( in );
         exit( 1 );
      }
   }

   /* read values */
   for( j = 0; j < nVertices; j++ )
   {
      for( i = 0; i < *nReferenceColumns; i++ )
      {
         rval = fscanf( in, "%lf", &(*referenceValues)[i][j] );
         if( rval == 0 || rval == EOF )
         {
            fprintf( stderr, "Error: unable to read reference values. (Malformatted input?)\n" );
            fclose( in );
            exit( 1 );
         }
      }
   }

   /* close input */
   fclose( in );
}

void destroyReferenceValues( int nReferenceColumns, double*** referenceValues, char*** referenceNames )
{
   int i;

   /* deallocate each column */
   for( i = 0; i < nReferenceColumns; i++ )
   {
      hmDestroy( (*referenceValues)[i] );
      hmDestroy( (*referenceNames)[i] );
   }

   /* deallocate array of columns */
   hmDestroy( *referenceValues );
   hmDestroy( *referenceNames );
}

void compareDistance( hmTriDistance* distance,
                      int nReferenceColumns,
                      double** referenceValues,
                      char** referenceNames,
                      const char* prefix )
{
   int i;
   size_t j;
   double lInfError, meanRelError;
   double diameter;

   diameter = hmTriMeshDiameter( distance->surface );

   printf( "Error relative to %s\n", referenceNames[0] );
   printf( "=============================================================\n\n" );

   /* compute error in solution computed by libgeodesic */
      lInfError = hmTriMeshLInfinityDistance( distance->surface, referenceValues[0], distance->distance.values );
   meanRelError = hmTriMeshMeanRelativeError( distance->surface, referenceValues[0], distance->distance.values );
   printf( "libgeodesic\n" );
   printf( "--------------------------------\n" );
   printf( " MAX: %.20f%%\n",   100.*lInfError/diameter );
   printf( "MEAN: %.20f%%\n\n", 100.*meanRelError );

   /* if mesh output is specified, write the reference solution in .1 */
   if( prefix != NULL )
   {
      for( j = 0; j < distance->surface->nVertices; j++ )
      {
         distance->distance.values[j] = referenceValues[0][j];
      }
      writeMesh( distance, prefix, 1 );
   }

   /* compute error in user-provided solutions */
   for( i = 1; i < nReferenceColumns; i++ )
   {
         lInfError = hmTriMeshLInfinityDistance( distance->surface, referenceValues[0], referenceValues[i] );
      meanRelError = hmTriMeshMeanRelativeError( distance->surface, referenceValues[0], referenceValues[i] );
      printf( "%s\n", referenceNames[i] );
      printf( "--------------------------------\n" );
      printf( " MAX: %.20f%%\n",   100.*lInfError/diameter );
      printf( "MEAN: %.20f%%\n\n", 100.*meanRelError );

      if( prefix != NULL )
      {
         for( j = 0; j < distance->surface->nVertices; j++ )
         {
            distance->distance.values[j] = referenceValues[i][j];
         }
         writeMesh( distance, prefix, 1+i );
      }
   }
}

