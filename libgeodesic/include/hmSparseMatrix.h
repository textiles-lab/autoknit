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

#ifndef LIBGEODESIC_HMSPARSEMATRIX_H
#define LIBGEODESIC_HMSPARSEMATRIX_H

#include <stddef.h>
#include "hmConstants.h"

/* ======================= CHOLMOD IMPLEMENTATION ======================= */
#ifdef HM_USE_CHOLMOD
#include <cholmod.h>

/** \brief Real sparse matrix. */
typedef struct hmSparseMatrix {

   /** \brief Number of rows (equal to data->nrow). */
   size_t nRows;

   /** \brief Number of columns (equal to data->ncol). */
   size_t nColumns;

   /** \brief Number of nonzero entries (equal to data->nzmax). */
   size_t nNonZeros;

   /** \brief Values of nonzero entries (points to data->x). */
   double* values;

   /** \brief Row indices of nonzero entries (points to data->i). */
   SuiteSparse_long* rowIndices;

   /** \brief Starting indices for nonzeros in each column (points to data->p). */
   SuiteSparse_long* columnStart;

   /** \private
    * \brief CHOLMOD sparse matrix. */
   cholmod_sparse* data;

} hmSparseMatrix;

#endif /* HM_USE_CHOLMOD */
/* ===================== END CHOLMOD IMPLEMENTATION ===================== */


/* ======================= HSLMA87 IMPLEMENTATION ======================= */
#ifdef HM_USE_HSLMA87
#include <hsl_ma87d.h>

/** \brief Real sparse matrix. */
typedef struct hmSparseMatrix {

   /** \brief Number of rows. */
   size_t nRows;

   /** \brief Number of columns. */
   size_t nColumns;

   /** \brief Number of nonzero entries. */
   size_t nNonZeros;

   /** \brief Values of nonzero entries. */
   double* values;

   /** \brief Row indices of nonzero entries. */
   int* rowIndices;

   /** \brief Starting indices for nonzeros in each column. */
   int* columnStart;

} hmSparseMatrix;


#endif /* HM_USE_HSLMA87 */
/* ===================== END HSLMA87 IMPLEMENTATION ===================== */

/** \brief Constructor.
 *
 * Builds an empty matrix.
 *
 * @param matrix    Target object.
 * @param nRows     Number of rows.
 * @param nColumns  Number of columns.
 * @param nNonZeros Number of nonzero entries.
 * \memberof hmSparseMatrix
 *
 */
void hmSparseMatrixInitialize( hmSparseMatrix* matrix,
                               size_t nRows,
                               size_t nColumns,
                               size_t nNonZeros );

/** \brief Destructor.
 *
 * @param matrix Target object.
 * \memberof hmSparseMatrix
 *
 */
void hmSparseMatrixDestroy( hmSparseMatrix* matrix );

#endif /* LIBGEODESIC_HMSPARSEMATRIX_H */

