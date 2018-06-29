#include "hmVectorPairSizeTDouble.h"
#include "hmUtility.h"
#include "hmConstants.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void hmVectorPairSizeTDoubleInitialize( hmVectorPairSizeTDouble* vector )
{
   vector->size = 0;
   vector->storage = hmVectorDefaultStorage;
   vector->entries = malloc( hmVectorDefaultStorage * sizeof(hmPairSizeTDouble) );
}

void hmVectorPairSizeTDoubleDestroy( hmVectorPairSizeTDouble* vector )
{
   hmDestroy( vector->entries );
}

void hmVectorPairSizeTDoubleResize( hmVectorPairSizeTDouble* vector, size_t size )
{
   vector->size = size;
   vector->storage = hmMaxSizeT( hmVectorDefaultStorage, hmNextPowerOfTwo( size ));
   free( vector->entries );
   vector->entries = malloc( vector->storage * sizeof(hmPairSizeTDouble) );
}

void hmVectorPairSizeTDoublePushBack( hmVectorPairSizeTDouble* vector, hmPairSizeTDouble value )
{
   hmPairSizeTDouble* newEntries;

   if( vector->size == vector->storage )
   {
      vector->storage *= 2;
      newEntries = malloc( vector->storage*sizeof(hmPairSizeTDouble) );
      memcpy( newEntries, vector->entries, vector->size*sizeof(hmPairSizeTDouble) );
      free( vector->entries );
      vector->entries = newEntries;
   }

   vector->entries[vector->size] = value;
   vector->size++;
}

hmPairSizeTDouble hmVectorPairSizeTDoublePopBack( hmVectorPairSizeTDouble* vector )
{
   assert( vector->size > 0 );

   vector->size--;

   return vector->entries[ vector->size ];
}

void hmVectorPairSizeTDoubleSort( hmVectorPairSizeTDouble* vector )
{
   qsort( vector->entries,
          vector->size,
          sizeof(hmPairSizeTDouble),
          hmPairSizeTDoubleComparator );
}

