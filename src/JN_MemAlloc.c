#include <stdio.h>
#include <stdlib.h>
#include "JN_MemAlloc.h"



void * jn_malloc(size_t bytes, char *VarName){
  
  void *tmp;
  if((tmp=malloc(bytes))==NULL){
    
    fprintf(stderr,"\n\nMemory Allocation Failure.\n");
    fprintf(stderr,"malloc(%s)",VarName);
    fprintf(stderr,"%s = %d bytes",VarName,(int)bytes);
    fprintf(stderr,"Total bytes requested via jn_malloc and jn_calloc:\n");
    fprintf(stderr,"%0.f bytes = %.3f Kb = %3.f Mb\n",gJN_total_bytes_requested,gJN_total_bytes_requested/1000.0,gJN_total_bytes_requested/1000000.0);
    fprintf(stderr,"The program needs more memory.\n");
    
    exit(1);
  }else{
    gJN_total_bytes_requested+=(double)bytes;
    
  };
  return(tmp);
}
     
extern void *jn_calloc(size_t Num, size_t bytes, char *NumName,char *VarName){

  void *tmp;
  if((tmp=calloc(Num,bytes))==NULL){
    
    fprintf(stderr,"\n\nMemory Allocation Failure.\n");
    fprintf(stderr,"calloc(%s, %s)\n",NumName,VarName);
    fprintf(stderr,"%s = %d items\n",NumName,(int)Num);
    fprintf(stderr,"%s = %d bytes\n",VarName,(int)bytes);
    
    fprintf(stderr,"Total bytes requested via jn_malloc and jn_calloc:\n");
    fprintf(stderr,"%0.f bytes = %.3f Kb = %3.f Mb\n",gJN_total_bytes_requested,gJN_total_bytes_requested/1000.0,gJN_total_bytes_requested/1000000.0);
    fprintf(stderr,"The program needs more memory.\n");
    
    exit(1);
  }else{
    gJN_total_bytes_requested+=(double)bytes;
    
  };
  return(tmp);
  
}



