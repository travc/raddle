

#ifndef __JN_MEM_ALLOC_FLAG__
#define __JN_MEM_ALLOC_FLAG__

#include <stdlib.h>

#ifdef UN_EXTERN
	#define JN_GLOB
#else 
	#define JN_GLOB extern
#endif


double gJN_total_bytes_requested;

#define JN_MALLOC(X)   	jn_malloc(X,#X)
#define JN_CALLOC(X,Y)	jn_calloc(X,Y,#X,#Y)

extern void *jn_malloc(size_t bytes,char *VarName);
extern void *jn_calloc(size_t Num, size_t bytes, char *NumName,char *VarName);


#endif
