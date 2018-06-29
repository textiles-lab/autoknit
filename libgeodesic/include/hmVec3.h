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

#ifndef LIBGEODESIC_HMVEC3_H
#define LIBGEODESIC_HMVEC3_H

#include <stdio.h>
#include <stddef.h>
#include <math.h>

/** \class hmVec3
 *  \brief Vector in Euclidean 3-space.
 */
typedef double hmVec3[3];

/** \brief Initializes vector with specified components.
 *
 * @param u Target object.
 * @param x x-coordinate.
 * @param y y-coordinate.
 * @param z z-coordinate.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Set( hmVec3 u, double x, double y, double z )
{
   u[0] = x;
   u[1] = y;
   u[2] = z;
}

/** \brief Initializes vector from the first three components of an array.
 *
 * @param u Target object.
 * @param x Array of at least three doubles.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3FromArray( hmVec3 u, const double* x )
{
   u[0] = x[0];
   u[1] = x[1];
   u[2] = x[2];
}

/** \brief Computes sum \f$w=u+v\f$ of two vectors \f$u\f$ and \f$v\f$.
 *
 * @param w Sum of u and v.
 * @param u First summand.
 * @param v Second summand.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Add( hmVec3 w, const hmVec3 u, const hmVec3 v )
{
   w[0] = u[0] + v[0];
   w[1] = u[1] + v[1];
   w[2] = u[2] + v[2];
}

/** \brief Replaces the vector \f$u\f$ with the vector \f$u+v\f$.
 *
 * @param u Initial vector.
 * @param v Increment.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Inc( hmVec3 u, const hmVec3 v )
{
   u[0] += v[0];
   u[1] += v[1];
   u[2] += v[2];
}

/** \brief Computes difference \f$w=u-v\f$ of two vectors \f$u\f$ and \f$v\f$.
 *
 * @param w Difference.
 * @param u Minuend.
 * @param v Subtrahend.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Sub( hmVec3 w, const hmVec3 u, const hmVec3 v )
{
   w[0] = u[0] - v[0];
   w[1] = u[1] - v[1];
   w[2] = u[2] - v[2];
}

/** \brief Replaces the vector \f$u\f$ with the vector \f$u-v\f$.
 *
 * @param u Initial vector.
 * @param v Decrement.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Dec( hmVec3 u, const hmVec3 v )
{
   u[0] -= v[0];
   u[1] -= v[1];
   u[2] -= v[2];
}

/** \brief Replaces the vector \f$u\f$ with the vector \f$-u\f$.
 *
 * @param u Argument.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Negate( hmVec3 u )
{
   u[0] = -u[0];
   u[1] = -u[1];
   u[2] = -u[2];
}

/** \brief Multiples a vector \f$u\f$ by a constant factor \f$a\f$ to get the new vector \f$w=au\f$.
 *
 * @param w Scaled vector.
 * @param a Scale.
 * @param u Initial vector.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Mul( hmVec3 w, double a, const hmVec3 u )
{
   w[0] = a*u[0];
   w[1] = a*u[1];
   w[2] = a*u[2];
}

/** \brief Replaces a vector \f$u\f$ with the rescaled vector \f$au\f$.
 *
 * @param u Vector.
 * @param a Scale.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Scale( hmVec3 u, double a )
{
   u[0] *= a;
   u[1] *= a;
   u[2] *= a;
}

/** \brief Computes Euclidean inner product \f$a = u \cdot v \f$ of two vectors \f$u\f$ and \f$v\f$.
 *
 * @param u First argument.
 * @param v Second argument.
 * @return Inner product.
 * \memberof hmVec3
 *
 */
static __inline__ double hmVec3Dot( const hmVec3 u, const hmVec3 v )
{
   return u[0]*v[0] +
          u[1]*v[1] +
          u[2]*v[2] ;
}

/** \brief Computes cross product \f$w = u \times v\f$ of two vectors \f$u\f$ and \f$v\f$.
 *
 * @param w Cross product.
 * @param u First argument.
 * @param v Second argument.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Cross( hmVec3 w, const hmVec3 u, const hmVec3 v )
{
   w[0] = u[1]*v[2] - u[2]*v[1];
   w[1] = u[2]*v[0] - u[0]*v[2];
   w[2] = u[0]*v[1] - u[1]*v[0];
}

/** \brief Computes Euclidean length \f$a = ||u|| = \sqrt{u \cdot u}\f$ of a vector \f$u\f$.
 *
 * @param u Argument.
 * @return Norm.
 * \memberof hmVec3
 *
 */
static __inline__ double hmVec3Norm( const hmVec3 u )
{
   return sqrt( u[0]*u[0] +
                u[1]*u[1] +
                u[2]*u[2] );
}

/** \brief Computes squared Euclidean length \f$a = ||u||^2 = u \cdot u\f$ of a vector \f$u\f$.
 *
 * @param u Argument.
 * @return Norm.
 * \memberof hmVec3
 *
 */
static __inline__ double hmVec3Norm2( const hmVec3 u )
{
   return u[0]*u[0] +
          u[1]*u[1] +
          u[2]*u[2] ;
}

/** \brief Divides a vector \f$u\f$ by its length.
 *
 * @param u Argument.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Normalize( hmVec3 u )
{
   double rNorm = 1. / sqrt( u[0]*u[0] + u[1]*u[1] + u[2]*u[2] );

   u[0] *= rNorm;
   u[1] *= rNorm;
   u[2] *= rNorm;
}

/** \brief Computes unit length vector \f$w = u/|u|\f$ parallel to \f$u\f$.
 *
 * @param w Normalized vector.
 * @param u Argument.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Unit( hmVec3 w, const hmVec3 u )
{
   double rNorm = 1. / sqrt( u[0]*u[0] + u[1]*u[1] + u[2]*u[2] );

   w[0] = u[0]*rNorm;
   w[1] = u[1]*rNorm;
   w[2] = u[2]*rNorm;
}

/** \brief Prints the vector \f$u\f$ to stdout.
 *
 * @param u Argument.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3Print( const hmVec3 u )
{
   printf( "[ %e %e %e ]\n", u[0], u[1], u[2] );
}

/** \brief Prints the vector \f$u\f$ to stderr.
 *
 * @param u Argument.
 * \memberof hmVec3
 *
 */
static __inline__ void hmVec3PrintErr( const hmVec3 u )
{
   fprintf( stderr, "[ %e %e %e ]\n", u[0], u[1], u[2] );
}

#endif /* LIBGEODESIC_HMVEC3_H */

