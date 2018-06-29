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

#ifndef LIBGEODESIC_HMDENSEMATRIX_H
#define LIBGEODESIC_HMDENSEMATRIX_H

#include <stddef.h>
#include "hmConstants.h"

/* ======================= CHOLMOD IMPLEMENTATION ======================= */
#ifdef HM_USE_CHOLMOD
#include <cholmod.h>

/** \brief Real dense matrix. */
typedef struct hmDenseMatrix {

   /** \brief Number of rows (equal to data->nrow). */
   size_t nRows;

   /** \brief Number of columns (equal to data->ncol). */
   size_t nColumns;

   /** \brief Values of entries (points to data->x). */
   double* values;

   /** \private
    * \brief CHOLMOD dense matrix. */
   cholmod_dense* data;

} hmDenseMatrix;

#endif /* HM_USE_CHOLMOD */
/* ===================== END CHOLMOD IMPLEMENTATION ===================== */

/* ======================= HSLMA87 IMPLEMENTATION ======================= */
#ifdef HM_USE_HSLMA87
#include <hsl_ma87d.h>

/** \brief Real dense matrix. */
typedef struct hmDenseMatrix {

   /** \brief Number of rows (equal to data->nrow). */
   size_t nRows;

   /** \brief Number of columns (equal to data->ncol). */
   size_t nColumns;

   /** \brief Values of entries (points to data->x). */
   double* values;

} hmDenseMatrix;

#endif /* HM_USE_HSLMA87 */
/* ===================== END HSLMA87 IMPLEMENTATION ===================== */

/** \brief Constructor.
 *
 * @param matrix Target object.
 * @param nRows Number of rows.
 * @param nColumns Number of columns.
 * \memberof hmDenseMatrix
 *
 */
void hmDenseMatrixInitialize( hmDenseMatrix* matrix,
                              size_t nRows,
                              size_t nColumns );

/** \brief Destructor.
 *
 * @param matrix Target object.
 * \memberof hmDenseMatrix
 *
 */
void hmDenseMatrixDestroy( hmDenseMatrix* matrix );

/** \brief Replaces matrix1 with a copy of matrix2.
 *
 * @param matrix1 Target object.
 * @param matrix2 Source object.
 * \memberof hmDenseMatrix
 *
 */
void hmDenseMatrixCopy(       hmDenseMatrix* matrix1,
                        const hmDenseMatrix* matrix2 );

#endif /* LIBGEODESIC_HMDENSEMATRIX_H */

