#include <stdlib.h>
#include "hmSparseMatrix.h"
#include "hmUtility.h"

#ifdef HM_USE_CHOLMOD
/* ======================= CHOLMOD IMPLEMENTATION ======================= */
extern cholmod_common hmCHOLMODCommon;

void hmSparseMatrixInitialize( hmSparseMatrix* matrix,
                               size_t nRows,
                               size_t nColumns,
                               size_t nNonZeros )
{
   int sorted = 1;
   int packed = 1;
   int stype = -1; /* assume matrix is symmetric and use entries only in the lower triangle */

   matrix->data =
   cholmod_l_allocate_sparse(
         nRows,
         nColumns,
         nNonZeros,
         sorted,
         packed,
         stype,
         CHOLMOD_REAL,
         &hmCHOLMODCommon
      );

   matrix->nRows       = nRows;
   matrix->nColumns    = nColumns;
   matrix->nNonZeros   = nNonZeros;
   matrix->values      = matrix->data->x;
   matrix->rowIndices  = matrix->data->i;
   matrix->columnStart = matrix->data->p;
}

void hmSparseMatrixDestroy( hmSparseMatrix* matrix )
{
   if( matrix->data != NULL )
   {
      cholmod_l_free_sparse( &matrix->data,
                             &hmCHOLMODCommon );

      matrix->data        = NULL;
      matrix->values      = NULL;
      matrix->rowIndices  = NULL;
      matrix->columnStart = NULL;

      matrix->nRows     = 0;
      matrix->nColumns  = 0;
      matrix->nNonZeros = 0;
   }
}

#endif /* HM_USE_CHOLMOD */
/* ===================== END CHOLMOD IMPLEMENTATION ===================== */

/* ======================= HSLMA87 IMPLEMENTATION ======================= */
#ifdef HM_USE_HSLMA87
extern struct ma87_control hmHSLMA87control;
extern struct ma87_info    hmHSLMA87info;

void hmSparseMatrixInitialize( hmSparseMatrix* matrix,
                               size_t nRows,
                               size_t nColumns,
                               size_t nNonZeros )
{
   matrix->nRows       = nRows;
   matrix->nColumns    = nColumns;
   matrix->nNonZeros   = nNonZeros;
   matrix->values      = malloc( nNonZeros    * sizeof(double));
   matrix->rowIndices  = malloc( nNonZeros    * sizeof(int));
   matrix->columnStart = malloc( (nColumns+1) * sizeof(int));
}

void hmSparseMatrixDestroy( hmSparseMatrix* matrix )
{
   hmDestroy( matrix->values      );
   hmDestroy( matrix->rowIndices  );
   hmDestroy( matrix->columnStart );
}

#endif /* HM_USE_HSLMA87 */
/* ===================== END HSLMA87 IMPLEMENTATION ===================== */

