#include "hmDenseMatrix.h"
#include "hmUtility.h"
#include <stdlib.h>
#include <string.h>

/* ======================= CHOLMOD IMPLEMENTATION ======================= */
#ifdef HM_USE_CHOLMOD
extern cholmod_common hmCHOLMODCommon;

void hmDenseMatrixInitialize( hmDenseMatrix* matrix,
                              size_t nRows,
                              size_t nColumns )
{
   size_t leadingDimension = nRows;

   matrix->data =
   cholmod_l_allocate_dense(
         nRows,
         nColumns,
         leadingDimension,
         CHOLMOD_REAL,
         &hmCHOLMODCommon
      );

   matrix->nRows    = nRows;
   matrix->nColumns = nColumns;
   matrix->values   = matrix->data->x;
}

void hmDenseMatrixDestroy( hmDenseMatrix* matrix )
{
   if( matrix->data != NULL )
   {
      cholmod_l_free_dense( &matrix->data,
                            &hmCHOLMODCommon );

      matrix->data   = NULL;
      matrix->values = NULL;

      matrix->nRows    = 0;
      matrix->nColumns = 0;
   }
}

void hmDenseMatrixCopy(       hmDenseMatrix* matrix1,
                        const hmDenseMatrix* matrix2 )
{
   hmDenseMatrixDestroy( matrix1 );

   matrix1->data = cholmod_l_copy_dense( matrix2->data, &hmCHOLMODCommon );

   matrix1->nRows    = matrix1->data->nrow;
   matrix1->nColumns = matrix1->data->ncol;
   matrix1->values   = matrix1->data->x;
}

#endif /* HM_USE_CHOLMOD */
/* ===================== END CHOLMOD IMPLEMENTATION ===================== */

/* ======================= HSLMA87 IMPLEMENTATION ======================= */
#ifdef HM_USE_HSLMA87
extern struct ma87_control hmHSLMA87control;
extern struct ma87_info    hmHSLMA87info;

void hmDenseMatrixInitialize( hmDenseMatrix* matrix,
                              size_t nRows,
                              size_t nColumns )
{
   matrix->nRows    = nRows;
   matrix->nColumns = nColumns;

   matrix->values = malloc( nRows*nColumns * sizeof(double));
}

void hmDenseMatrixDestroy( hmDenseMatrix* matrix )
{
   matrix->nRows    = 0;
   matrix->nColumns = 0;

   hmDestroy( matrix->values );
}

void hmDenseMatrixCopy(       hmDenseMatrix* matrix1,
                        const hmDenseMatrix* matrix2 )
{
   int nEntries;

   hmDenseMatrixDestroy( matrix1 );

   matrix1->nRows    = matrix2->nRows;
   matrix1->nColumns = matrix2->nColumns;
   nEntries = matrix2->nRows * matrix2->nColumns;

   matrix1->values = malloc( nEntries * sizeof( double ));
   memcpy( matrix1->values, matrix2->values, nEntries * sizeof(double));
}

#endif /* HM_USE_HSLMA87 */
/* ===================== END HSLMA87 IMPLEMENTATION ===================== */

