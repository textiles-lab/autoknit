#include "hmContext.h"
#include "hmConstants.h"

/* ======================= CHOLMOD IMPLEMENTATION ======================= */
#ifdef HM_USE_CHOLMOD
#include <cholmod.h>

cholmod_common hmCHOLMODCommon;

void hmContextInitialize( hmContext* context )
{
   cholmod_l_start( &hmCHOLMODCommon );
   context->initialized = 1;
}

void hmContextDestroy( hmContext* context )
{
   cholmod_l_finish( &hmCHOLMODCommon );
   context->initialized = 0;
}

#endif /* HM_USE_CHOLMOD */
/* ===================== END CHOLMOD IMPLEMENTATION ===================== */

/* ======================= HSLMA87 IMPLEMENTATION ======================= */
#ifdef HM_USE_HSLMA87
#include <hsl_mc68i.h>
#include <hsl_ma87d.h>

struct ma87_control hmHSLMA87control;
struct ma87_info    hmHSLMA87info;

struct mc68_control hmHSLMC68control;
struct mc68_info    hmHSLMC68info;

void hmContextInitialize( hmContext* context )
{
   ma87_default_control( &hmHSLMA87control );
   mc68_default_control( &hmHSLMC68control );
   context->initialized = 1;
}

void hmContextDestroy( hmContext* context )
{
   context->initialized = 0;
}

#endif /* HM_USE_CHOLMOD */
/* ===================== END CHOLMOD IMPLEMENTATION ===================== */

