#include <stdio.h>
#include <string.h>
#include "hmCholeskyFactor.h"
#include "hmUtility.h"
#include "hmConstants.h"

/* ======================= CHOLMOD IMPLEMENTATION ======================= */
#ifdef HM_USE_CHOLMOD
extern cholmod_common hmCHOLMODCommon;

void hmCholeskyFactorInitialize( hmCholeskyFactor* factor )
{
   factor->data = NULL;
}

void hmCholeskyFactorDestroy( hmCholeskyFactor* factor )
{
   if( factor->data != NULL )
   {
      cholmod_l_free_factor( &factor->data, &hmCHOLMODCommon );
      factor->data = NULL;
   }
}

void hmCholeskyFactorCopy( hmCholeskyFactor* factor1,
                           hmCholeskyFactor* factor2 )
{
   hmCholeskyFactorDestroy( factor1 );

   factor1->data = cholmod_l_copy_factor( factor2->data, &hmCHOLMODCommon );
}

void hmCholeskyFactorReorder( hmCholeskyFactor* factor,
                              hmSparseMatrix* matrix )
{}

void hmCholeskyFactorSymbolic( hmCholeskyFactor* factor,
                              hmSparseMatrix* matrix )
{
   /* clear the old factor if already built */
   hmCholeskyFactorDestroy( factor );

   /* build symbolic factorization (this step is likely the
    * most expensive part of the entire algorithm!) */
   factor->data = cholmod_l_analyze( matrix->data, &hmCHOLMODCommon );
}

void hmCholeskyFactorNumerical( hmCholeskyFactor* factor,
                            hmSparseMatrix* matrix )
{
   /* build numerical factorization */
   cholmod_l_factorize( matrix->data, factor->data, &hmCHOLMODCommon );
}

void hmCholeskyFactorBacksolve( hmCholeskyFactor* factor,
                                      hmDenseMatrix* x,
                                const hmDenseMatrix* b )
{
   hmDenseMatrixDestroy( x );

   x->data = cholmod_l_solve(
         CHOLMOD_A,
         factor->data,
         b->data,
         &hmCHOLMODCommon
      );

   x->nRows    = x->data->nrow;
   x->nColumns = x->data->ncol;
   x->values   = x->data->x;
}

#endif /* HM_USE_CHOLMOD */
/* ===================== END CHOLMOD IMPLEMENTATION ===================== */

/* ======================= HSLMA87 IMPLEMENTATION ======================= */
#ifdef HM_USE_HSLMA87
extern struct ma87_control hmHSLMA87control;
extern struct ma87_info    hmHSLMA87info;
extern struct mc68_control hmHSLMC68control;
extern struct mc68_info    hmHSLMC68info;

void hmCholeskyFactorInitialize( hmCholeskyFactor* factor )
{
   factor->degree = 0;
   factor->keep  = NULL;
   factor->order = NULL;
}

void hmCholeskyFactorDestroy( hmCholeskyFactor* factor )
{
   if( factor->keep != NULL )
   {
      ma87_finalise( &factor->keep, &hmHSLMA87control );
      factor->keep = NULL;
   }

   hmDestroy( factor->order );
}

void hmCholeskyFactorCopy( hmCholeskyFactor* factor1,
                           hmCholeskyFactor* factor2 )
{
   hmCholeskyFactorDestroy( factor1 );

   factor1->degree = factor2->degree;
   factor1->order = malloc( factor2->degree * sizeof(int));
   memcpy( factor1->order, factor2->order, factor2->degree * sizeof(int));
}

void hmCholeskyFactorReorder( hmCholeskyFactor* factor,
                              hmSparseMatrix* matrix )
{
   /* clear the old factor if already built */
   hmCholeskyFactorDestroy( factor );
   factor->degree = matrix->nColumns;

   /* compute a matrix reordering to reduce fill-in */
   factor->order = malloc( factor->degree * sizeof(int));
   mc68_order(
         hmReorderingScheme,
         factor->degree,
         matrix->columnStart,
         matrix->rowIndices,
         factor->order,
         &hmHSLMC68control,
         &hmHSLMC68info
      );
}

void hmCholeskyFactorSymbolic( hmCholeskyFactor* factor,
                               hmSparseMatrix* matrix )
{
   /* build symbolic factorization (this step is likely the
    * most expensive part of the entire algorithm!) */
   ma87_analyse(
         factor->degree,
         matrix->columnStart,
         matrix->rowIndices,
         factor->order,
         &factor->keep,
         &hmHSLMA87control,
         &hmHSLMA87info
      );
}

void hmCholeskyFactorNumerical( hmCholeskyFactor* factor,
                                hmSparseMatrix* matrix )
{
   /* build numerical factorization */
   ma87_factor(
         factor->degree,
         matrix->columnStart,
         matrix->rowIndices,
         matrix->values,
         factor->order,
         &factor->keep,
         &hmHSLMA87control,
         &hmHSLMA87info
      );
}

void hmCholeskyFactorBacksolve( hmCholeskyFactor* factor,
                                      hmDenseMatrix* x,
                                const hmDenseMatrix* b )
{
   int job  = 0; /* just solve Ax=b */
   int nrhs = 1; /* only one right-hand side */

   /* ma87_solve() overwrites the data vector with the
    * solution vector, so we need to make a copy first */
   hmDenseMatrixDestroy( x );
   hmDenseMatrixCopy( x, b );

   ma87_solve(
         job,
         nrhs,
         factor->degree,
         x->values,
         factor->order,
         &factor->keep,
         &hmHSLMA87control,
         &hmHSLMA87info
      );
}

#endif /* HM_USE_HSLMA87 */
/* ===================== END HSLMA87 IMPLEMENTATION ===================== */

