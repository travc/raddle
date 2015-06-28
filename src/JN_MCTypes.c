#include <stdio.h>
#include <stdlib.h>
#include "JN_MCTypes.h"
#include "JN_MemAlloc.h"
struct intval *AllocIntval(int Min,int Max,int Step){
  
  struct intval * Ptr;
  Ptr=(struct intval *) JN_MALLOC(sizeof(intval));
  Ptr->H=AllocDiscreteHisto(Min,Max,Step);
  InitIntvalToZero(Ptr);
  
  return(Ptr);
};

struct doubval *AllocDoubval(double Min,double Max,double Step){
  struct doubval * Ptr;
  Ptr=(doubval* ) JN_MALLOC(sizeof(doubval));
  Ptr->H=AllocContHisto(Min,Max,Step);
  InitDoubvalToZero(Ptr);
  return(Ptr);
};


struct DiscreteHisto *AllocDiscreteHisto(int Min,int Max,int Step){
  
  struct DiscreteHisto * Ptr;
  if(Step<=0) 
    Ptr=NULL;
  else{
    Ptr=(struct DiscreteHisto *) JN_MALLOC(sizeof(struct DiscreteHisto));
    Ptr->Min=Min;
    Ptr->Max=Max;
    Ptr->Step=Step;
    Ptr->Remainder=((Max-Min+1)%Step!=0);
    
    Ptr->nBins=(int)(Max-Min+1)/Step+2+Ptr->Remainder;  /* 2 extra for values < Min and >Max */
    Ptr->Counts=(int *)JN_MALLOC((size_t)Ptr->nBins*sizeof(int));
    InitDiscreteHisto(Ptr);
  }
  return(Ptr);

};

struct ContHisto *AllocContHisto(double Min,double Max,double Step){
  struct ContHisto * Ptr;

  if(Step<=0.0){
    Ptr=NULL;
  } else{
    Ptr=(struct ContHisto *) JN_MALLOC(sizeof(struct ContHisto));
    Ptr->Min=Min;
    Ptr->Max=Max;
    Ptr->Step=Step;
    Ptr->nBins=(int)(Max-Min)/Step;  /* Initial Calculation */
    if(fabs((Min+(double)Ptr->nBins*Step)-Max)>(1e-4)*Step)
      Ptr->Remainder=1; /* Check for remainder */
    else
      Ptr->Remainder=0;
    Ptr->nBins+=2+Ptr->Remainder; /* Update calculation and add bins for instances <Min or >Max */
    Ptr->Counts=(int *)JN_MALLOC((size_t)Ptr->nBins*sizeof(int));
    InitContHisto(Ptr);

  }
  return(Ptr);
}


void InitIntvalToZero(struct intval *t){
  t->v=0;
  t->Mean=0.0;
  t->SS=0.0;
  t->Var=0.0;
  t->nObs=0.0;
  if(t->H!=NULL)
    InitDiscreteHisto(t->H);
  

};

void InitIntvalSummaryToZero(struct intval *t){
  
  // Init the Summary but NOT the initial value
  t->Mean=0.0;
  t->SS=0.0;
  t->Var=0.0;
  t->nObs=0.0;
  if(t->H!=NULL)
    InitDiscreteHisto(t->H);
  

};


void InitDoubvalToZero(struct doubval *t){
  t->v=0.0;
  t->Mean=0.0;
  t->SS=0.0;
  t->Var=0.0;
  t->nObs=0.0;
  if(t->H!=NULL)
    InitContHisto(t->H);
  
};



void InitDoubvalSummaryToZero(struct doubval *t){

  // Init the Summary but NOT the initial value
  t->Mean=0.0;
  t->SS=0.0;
  t->Var=0.0;
  t->nObs=0.0;
  if(t->H!=NULL)
    InitContHisto(t->H);
  

};

void InitDiscreteHisto(struct DiscreteHisto  *D){
  int i;
  if(D!=NULL){
    for(i=0;i<D->nBins;i++)
    D->Counts[i]=0;
    D->nObs=0;
  }

};

void InitContHisto(struct ContHisto * C){
  int i;
  if(C!=NULL){
    for(i=0;i<C->nBins;i++)
      C->Counts[i]=0;
    C->nObs=0;
  }
};



/* INCREMENTING FUNCTIONS */
void IncrementIntval(struct intval *t){
  double n;

  t->nObs++;

  n=(double)t->nObs;

  t->Mean+=((double)t->v-t->Mean)/n;
  t->SS+=(double)(t->v*t->v);

  if(t->nObs>1){
    t->Var=(t->SS-n*(t->Mean*t->Mean)/(n*(n-1.0)));
  }

  IncrementDiscreteHisto(t->H,t->v);
};

void IncrementDoubval(struct doubval *t){
  double n;

  t->nObs++;

  n=(double)t->nObs;

  t->Mean+=(t->v-t->Mean)/n;
  t->SS+=(t->v*t->v);

  if(t->nObs>1){
    t->Var=(t->SS-n*(t->Mean*t->Mean)/(n*(n-1.0)));
  }

  IncrementContHisto(t->H,t->v);
};


void IncrementDiscreteHisto(struct DiscreteHisto * H,int v){

  if(H!=NULL){
    H->nObs++;
    
    if(v>H->Max){
      H->Counts[H->nBins-1]++;
    } else if (v<H->Min){
      H->Counts[0]++;
    } else{
      H->Counts[(int)(v-H->Min)/H->Step+1]++;
    }
  }

};

void IncrementContHisto(struct ContHisto * H, double v){

  if(H!=NULL){
    H->nObs++;
    
    if(v>=H->Max){
      H->Counts[H->nBins-1]++;
    } else if (v<H->Min){
      H->Counts[0]++;
    } else{
      H->Counts[(int)((v-H->Min)/H->Step+1)]++;
    }
  }

};


struct DoubvalCovar *AllocDoubvalCovar(int n){
  int i,j;
  struct DoubvalCovar *Ptr;
  Ptr=(DoubvalCovar *)JN_MALLOC(sizeof(DoubvalCovar));
  Ptr->Num=n;
  Ptr->Square=(doubval ***)JN_MALLOC((size_t)n*sizeof(doubval**));
  
  for(i=0;i<n;i++){
    Ptr->Square[i]=(doubval**)JN_MALLOC((size_t)n*sizeof(doubval*));
    for(j=0;j<n;j++){
      Ptr->Square[i][j]=AllocDoubval(0.0,0.0,0.1);
    }
  }
  InitDoubvalCovar(Ptr);
  return(Ptr);
  
};

void InitDoubvalCovar(struct DoubvalCovar *D){
  
  int i,j;
  for(i=0;i<D->Num;i++){
    for(j=0;j<D->Num;j++){
      InitDoubvalToZero(D->Square[i][j]);
    }
  }
};

void IncrementDoubvalCovar(struct DoubvalCovar *D, struct doubval **Vals){

  int i,j;
  int n;
  n=Vals[0]->nObs;

  for(i=0;i<D->Num;i++){

    
    for(j=0;j<D->Num;j++){
      
      if(n!=Vals[i]->nObs)
	fprintf(stderr,"\n\nnObs not identical for variables being compared in DoubvalCovar!\n\n");
      
      /* Increment sum of cross-products */ 
      D->Square[i][j]->SS+=Vals[i]->v*Vals[j]->v;
      
      /* After cross-products calculated store covariance in the Var field of the doubval structure */ 
      /* NOTE: CODE DIFFERS FROM ERIC'S BY A FACTOR OF N-1.0 */
      if(n>1){
	D->Square[i][j]->Var=(D->Square[i][j]->SS-n*Vals[i]->Mean*Vals[j]->Mean)/(n*(n-1.0));
      }
    }

  }
  
};


/* OUTPUT FUNCTIONS */
void PrintDiscreteHistKey(FILE *outfile,struct DiscreteHisto *H){
  int i;
  int t;
  
  /* First bin */
  fprintf(outfile,"<%d  ",H->Min);

  for(i=0;i<H->nBins-3;i++){
    if(H->Step==1)
      fprintf(outfile,"%d      ",(H->Min+i));
    else
      fprintf(outfile,"[%d,%d] ",(H->Min+i*H->Step),(H->Min+(i+1)*H->Step-1));
  }
  
  /* Last two bins */
  if(H->Remainder){
    t=H->Min+(H->nBins-3)*H->Step;
    if(t==H->Max)
      fprintf(outfile,">=%d    ",t);
    else{
      fprintf(outfile,"[%d,%d] ",t,H->Max-1);
      fprintf(outfile,">%d     ",H->Max);
    }
  }  else {
    fprintf(outfile,"[%d,%d] ",H->Max-H->Step+1,H->Max);
    fprintf(outfile,">%d",H->Max);
  }
  fprintf(outfile,"\n");
};

void PrintContHistKey(FILE *outfile,struct ContHisto * H){

  int i;  
  double t;

  
  /* First bin */
  fprintf(outfile,"<%.2f) ",H->Min);
  
  for(i=0;i<H->nBins-3;i++){
    fprintf(outfile,"[%.2f,%.2f) ",(H->Min+i*H->Step),(H->Min+(i+1)*H->Step));
  }
  
  /* Last two bins */
  if(H->Remainder){
    t=H->Min+(double)(H->nBins-3)*H->Step;
    fprintf(outfile,"[%.2f,%.2f)",t,H->Max);
    fprintf(outfile,">=%.2f     ",H->Max);
  }  else {
    fprintf(outfile,"[%.2f,%.2f) ",H->Max-H->Step,H->Max);
    fprintf(outfile,">=%.2f",H->Max);
  }
  fprintf(outfile,"\n");


}


void PrintDiscreteHistKeySimple(FILE *outfile,struct DiscreteHisto *H){
  int i;
  int t;
  
  /* First bin */
  fprintf(outfile,"<%d",H->Min);

  for(i=0;i<H->nBins-3;i++){
    fprintf(outfile,"%d ",(H->Min+i*H->Step));
  }
  
  /* Last two bins */
  if(H->Remainder){
    t=H->Min+(H->nBins-3)*H->Step;
    if(t==H->Max)
      fprintf(outfile,">=%d ",t);
    else{
      fprintf(outfile,"%d ",t);
      fprintf(outfile,">%d ",H->Max);
    }
  }  else {
    fprintf(outfile,"%d ",H->Max-H->Step+1);
    fprintf(outfile,">%d ",H->Max);
  }
  fprintf(outfile,"\n");
};

void PrintContHistKeySimple(FILE *outfile,struct ContHisto * H){

  int i;  
  double t;
  
  /* First bin */
  fprintf(outfile,"<%.2f ",H->Min);
  
  for(i=0;i<H->nBins-3;i++){
    fprintf(outfile,"%.2f ",(H->Min+i*H->Step));
  }
  
  /* Last two bins */
  if(H->Remainder){
    t=H->Min+(double)(H->nBins-3)*H->Step;
    fprintf(outfile,"%.2f ",t);
    fprintf(outfile,">=%.2f ",H->Max);
  }  else {
    fprintf(outfile,"%.2f ",H->Max-H->Step);
    fprintf(outfile,">=%.2f",H->Max);
  }
  fprintf(outfile,"\n");


}





void PrintDiscreteHistCount(FILE *outfile,struct DiscreteHisto *H){
  int i;
  for(i=0;i<H->nBins;i++)
    fprintf(outfile,"%d ",H->Counts[i]);
  fprintf(outfile,"\n");
};

void PrintDiscreteHistProp(FILE *outfile,struct DiscreteHisto *H){
  int i;
  for(i=0;i<H->nBins;i++)
    fprintf(outfile,"%.8f ",H->Counts[i]/(double)H->nObs);
  fprintf(outfile,"\n");
};


void PrintContHistCount(FILE *outfile,struct ContHisto *H){
  int i;
  for(i=0;i<H->nBins;i++)
    fprintf(outfile,"%d ",H->Counts[i]);
  fprintf(outfile,"\n");
};

void PrintContHistProp(FILE *outfile,struct ContHisto *H){
  int i;
  for(i=0;i<H->nBins;i++)
    fprintf(outfile,"%.8f ",H->Counts[i]/(double)H->nObs);
  fprintf(outfile,"\n");
};

