#include "hmVectorSizeT.h"
#include "hmUtility.h"
#include "hmConstants.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void hmVectorSizeTInitialize( hmVectorSizeT* vector )
{
   vector->size = 0;
   vector->storage = hmVectorDefaultStorage;
   vector->entries = malloc( hmVectorDefaultStorage * sizeof(size_t) );
}

void hmVectorSizeTDestroy( hmVectorSizeT* vector )
{
   hmDestroy( vector->entries );
}

void hmVectorSizeTResize( hmVectorSizeT* vector, size_t size )
{
   vector->size = size;
   vector->storage = hmMaxSizeT( hmVectorDefaultStorage, hmNextPowerOfTwo( size ));
   free( vector->entries );
   vector->entries = malloc( vector->storage * sizeof(size_t) );
}

void hmVectorSizeTPushBack( hmVectorSizeT* vector, size_t value )
{
   size_t* newEntries;

   if( vector->size == vector->storage )
   {
      vector->storage *= 2;
      newEntries = malloc( vector->storage*sizeof(size_t) );
      memcpy( newEntries, vector->entries, vector->size*sizeof(size_t) );
      free( vector->entries );
      vector->entries = newEntries;
   }

   vector->entries[vector->size] = value;
   vector->size++;
}

size_t hmVectorSizeTPopBack( hmVectorSizeT* vector )
{
   assert( vector->size > 0 );

   vector->size--;

   return vector->entries[ vector->size ];
}

void hmVectorSizeTSort( hmVectorSizeT* vector )
{
   qsort( vector->entries,
          vector->size,
          sizeof(size_t),
          hmSizeTComparator );
}

