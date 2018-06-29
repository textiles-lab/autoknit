#include "hmVectorDouble.h"
#include "hmUtility.h"
#include "hmConstants.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void hmVectorDoubleInitialize( hmVectorDouble* vector )
{
   vector->size = 0;
   vector->storage = hmVectorDefaultStorage;
   vector->entries = malloc( hmVectorDefaultStorage * sizeof(double) );
}

void hmVectorDoubleDestroy( hmVectorDouble* vector )
{
   hmDestroy( vector->entries );
}

void hmVectorDoubleResize( hmVectorDouble* vector, size_t size )
{
   vector->size = size;
   vector->storage = hmMaxSizeT( hmVectorDefaultStorage, hmNextPowerOfTwo( size ));
   free( vector->entries );
   vector->entries = malloc( vector->storage * sizeof(double) );
}

void hmVectorDoublePushBack( hmVectorDouble* vector, double value )
{
   double* newEntries;

   if( vector->size == vector->storage )
   {
      vector->storage *= 2;
      newEntries = malloc( vector->storage*sizeof(double) );
      memcpy( newEntries, vector->entries, vector->size*sizeof(double) );
      free( vector->entries );
      vector->entries = newEntries;
   }

   vector->entries[vector->size] = value;
   vector->size++;
}

double hmVectorDoublePopBack( hmVectorDouble* vector )
{
   assert( vector->size > 0 );

   vector->size--;

   return vector->entries[ vector->size ];
}

void hmVectorDoubleSort( hmVectorDouble* vector )
{
   qsort( vector->entries,
          vector->size,
          sizeof(double),
          hmDoubleComparator );
}

