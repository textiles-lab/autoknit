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
 * @section DESCRIPTION
 *
 * This file contains common utility functions used throughout.
 *
 */

#ifndef LIBGEODESIC_HMUTILITY_H
#define LIBGEODESIC_HMUTILITY_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/* Safe(r) array deallocation. */
#define hmDestroy( x ) if( x != NULL ) { free( x ); x = NULL; }

/** \brief Returns the square of the argument.
 *
 * Useful for squaring long expressions without having to first
 * name an intermediate variable.
 *
 * @param x Value.
 * @return Squared value.
 *
 */
static __inline__ double hmSquare( double x )
{
   return x*x;
}

/** \brief Returns the smaller value.
 *
 * @param x First value.
 * @param y Second value.
 * @return Smaller value.
 *
 */
static __inline__ size_t hmMinSizeT( size_t x, size_t y )
{
   if( x <= y )
   {
      return x;
   }

   return y;
}

/** \brief Returns the smaller value.
 *
 * @param x First value.
 * @param y Second value.
 * @return Smaller value.
 *
 */
static __inline__ double hmMinDouble( double x, double y )
{
   if( x <= y )
   {
      return x;
   }

   return y;
}

/** \brief Returns the larger value.
 *
 * @param x First value.
 * @param y Second value.
 * @return Larger value.
 *
 */
static __inline__ size_t hmMaxSizeT( size_t x, size_t y )
{
   if( x >= y )
   {
      return x;
   }

   return y;
}

/** \brief Returns the larger value.
 *
 * @param x First value.
 * @param y Second value.
 * @return Larger value.
 *
 */
static __inline__ double hmMaxDouble( double x, double y )
{
   if( x >= y )
   {
      return x;
   }

   return y;
}

/** \brief Removes the mean value.
 *
 * @param array Target object; must be no shorter than the specified length.
 * @param length Array length.
 *
 */
static __inline__ void hmRemoveMean( double* array, size_t length )
{
   size_t i;
   double mean = 0.;

   /* compute mean */
   for( i = 0; i < length; i++ )
   {
      mean += array[i];
   }
   mean /= (double) length;

   /* subtract mean from each entry */
   for( i = 0; i < length; i++ )
   {
      array[i] -= mean;
   }
}

/** \brief Sets all array elements to the specified default value.
 *
 * @param array Target object; must be no shorter than the specified length.
 * @param length Array length.
 * @param value Default value.
 *
 */
static __inline__ void hmClearArrayDouble( double* array, size_t length, double value )
{
   size_t i;

   assert( array != NULL );

   for( i = 0; i < length; i++ )
   {
      array[i] = value;
   }
}

/** \brief Returns the next largest power of two.
 *
 * @param x Argument.
 * @return Next largest power of two.
 *
 */
static __inline__ int hmNextPowerOfTwo( size_t x )
{
   if( x == 0 ) return 1;

   /* Muahahahaha!! */
   x--;
   x |= ( x >> 1  );
   x |= ( x >> 2  );
   x |= ( x >> 4  );
   x |= ( x >> 8  );
   x |= ( x >> 16 );
#ifdef __LP64__ /* only if we're on a 64-bit architecture */
   x |= ( x >> 32 );
#endif
   x++;

   return x;
}

/** \brief Comparator for size_t values.
 *
 * Returns a negative, zero, or positive value depending on
 * whether the first element is less than, equal to, or
 * greater than the second element.  Compatible with qsort
 * (from the standard library).
 *
 * @param e1 First element.
 * @param e2 Second element.
 * @return Ordering.
 *
 */
static __inline__ int hmSizeTComparator( const void* e1,
                                         const void* e2 )
{
   double n1 = *((size_t*) e1 );
   double n2 = *((size_t*) e2 );

   if( n1 < n2 )
   {
      return -1;
   }
   
   if( n1 > n2 )
   {
      return 1;
   }

   return 0;
}

/** \brief Comparator for doubles.
 *
 * Returns a negative, zero, or positive value depending on
 * whether the first element is less than, equal to, or
 * greater than the second element.  Compatible with qsort
 * (from the standard library).
 *
 * @param e1 First element.
 * @param e2 Second element.
 * @return Ordering.
 *
 */
static __inline__ int hmDoubleComparator( const void* e1,
                                          const void* e2 )
{
   double x1 = *((double*) e1 );
   double x2 = *((double*) e2 );

   if( x1 < x2 )
   {
      return -1;
   }
   
   if( x1 > x2 )
   {
      return 1;
   }

   return 0;
}

/** \brief Samples a value from the specified interval uniformly at random.
 *
 * @param minValue Smallest value.
 * @param maxValue Largest value.
 * @return Sample.
 *
 */
static __inline__ double hmRandomReal( double minValue,
                                       double maxValue )
{
   const double rRandMax = 1. / (double) RAND_MAX;
   double u = (double) rand() * (double) rRandMax;

   return u*(maxValue-minValue) + minValue;
}

#endif /* LIBGEODESIC_HMUTILITY_H */

