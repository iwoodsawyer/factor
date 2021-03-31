/* Starting from version 7.8, MATLAB BLAS expects ptrdiff_t arguments for integers */
#if MATLAB_VERSION >= 0x0708
#include <stddef.h>
#include <stdlib.h>
#endif 
#include <string.h>
        
/* Define MX_HAS_INTERLEAVED_COMPLEX for version <9.4 */
#ifndef MX_HAS_INTERLEAVED_COMPLEX
#define MX_HAS_INTERLEAVED_COMPLEX 0
#endif

/* Starting from version 7.6, MATLAB BLAS is seperated */
#if MATLAB_VERSION >= 0x0705
#include <blas.h>
#endif
#include <lapack.h>

#ifndef min
#define min(a,b) ((a) <= (b) ? (a) : (b))
#endif
