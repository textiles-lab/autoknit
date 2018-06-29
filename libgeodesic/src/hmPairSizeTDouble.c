#include "hmPairSizeTDouble.h"

int hmPairSizeTDoubleComparator( const void* e1,
                                 const void* e2 )
{
   size_t n1 = ((hmPairSizeTDouble*) e1 )->n;
   size_t n2 = ((hmPairSizeTDouble*) e2 )->n;

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

