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

#ifndef LIBGEODESIC_HMCHOLESKYFACTOR_H
#define LIBGEODESIC_HMCHOLESKYFACTOR_H

#include <stddef.h>
#include "hmDenseMatrix.h"
#include "hmSparseMatrix.h"
#include "hmConstants.h"

/* ======================= CHOLMOD IMPLEMENTATION ======================= */
#ifdef HM_USE_CHOLMOD
#include <cholmod.h>

/** \brief Cholesky factor \f$L\f$ for a real symmetric positive-definite sparse matrix \f$A = LL^T\f$.
 *
 */
typedef struct hmCholeskyFactor {

   /** \brief Cholesky factor in CHOLMOD format. */
   cholmod_factor *data;

} hmCholeskyFactor;

#endif /* HM_USE_CHOLMOD */
/* ===================== END CHOLMOD IMPLEMENTATION ===================== */

/* ======================= HSLMA87 IMPLEMENTATION ======================= */
#ifdef HM_USE_HSLMA87
#include <hsl_mc68i.h>
#include <hsl_ma87d.h>

/** \brief Cholesky factor \f$L\f$ for a real symmetric positive-definite sparse matrix \f$A = LL^T\f$.
 *
 */
typedef struct hmCholeskyFactor {

   /** \brief Dimension of factored matrix. */
   int degree;

   /** \brief Cholesky factor in HSLMA87 format. */
   void *keep;
   
   /** \brief Matrix reordering. [1 x nColumns] */
   int *order;

   /** \brief Factorization parameters. */
   struct ma87_control control;

   /** \brief Diagnostic information. */
   struct ma87_info info;

} hmCholeskyFactor;

#endif /* HM_USE_HSLMA87 */
/* ===================== END HSLMA87 IMPLEMENTATION ===================== */

/** \brief Constructor.
 *
 * @param factor Target object.
 * \memberof hmCholeskyFactor
 *
 */
void hmCholeskyFactorInitialize( hmCholeskyFactor* factor );

/** \brief Destructor.
 *
 * @param factor Target object.
 * \memberof hmCholeskyFactor
 *
 */
void hmCholeskyFactorDestroy( hmCholeskyFactor* factor );

/** \brief Replaces factor1 with a copy of factor2.
 *
 *  Useful for reusing a symbolic factorization from a matrix with
 *  the same sparsity pattern.
 *
 * @param factor1 Target object.
 * @param factor2 Source object.
 * \memberof hmCholeskyFactor
 *
 */
void hmCholeskyFactorCopy( hmCholeskyFactor* factor1,
                           hmCholeskyFactor* factor2 );

/** \brief Computes a fill-reducing matrix reordering.
 *
 * When using HSL_MA87, the reordering scheme can be specified
 * in the file hmConstants.h.  When using CHOLMOD, this routine
 * is ignored since the reordering is computed at the same time
 * as the symbolic factorization (in hmCholeskyFactorSymbolic()).
 *
 * @param factor Target object.
 * @param matrix Matrix to be reordering.
 * \memberof hmCholeskyFactor
 *
 */
void hmCholeskyFactorReorder( hmCholeskyFactor* factor,
                              hmSparseMatrix* matrix );

/** \brief Performs symbolic Cholesky factorization.
 *
 * @param factor Target object.
 * @param matrix Matrix to be factored.
 * \memberof hmCholeskyFactor
 *
 */
void hmCholeskyFactorSymbolic( hmCholeskyFactor* factor,
                               hmSparseMatrix* matrix );

/** \brief Computes numerical Cholesky factorization.
 *
 * Should be called only after symbolic factorization.
 *
 * @param factor Target object.
 * @param matrix Matrix to be factored.
 * \memberof hmCholeskyFactor
 *
 */
void hmCholeskyFactorNumerical( hmCholeskyFactor* factor,
                                hmSparseMatrix* matrix );

/** \brief Solves a system \f$Ax=b\f$ by applying backsubstitution on the Cholesky factor \f$L\f$ of \f$A\f$.
 *
 * @param factor Target object.
 * @param x Solution vector.  Must be initialized, but not necessarily with the dimensions of the solution.
 * @param b Data vector.  Must have dimensions compatible with the Cholesky factor.
 * \memberof hmCholeskyFactor
 *
 */
void hmCholeskyFactorBacksolve( hmCholeskyFactor* factor,
                                      hmDenseMatrix* x,
                                const hmDenseMatrix* b );

#endif /* LIBGEODESIC_HMCHOLESKYFACTOR_H */

