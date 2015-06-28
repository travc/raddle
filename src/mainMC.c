/*
 *  mainMC.c
 *  Raddle
 *
 *  Created by John Novembre on 1/24/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */


#include <stdio.h>


#include "mainMC.h"
struct SamplerInfo InitSettingsAndData(int argc,char **argv);
void MC_CalcPtsig2MarginSurf(struct LsurfPoint *P,struct SamplerInfo * S,int LocusNo,double * s);
void MC_CalcPtsig2MarginAcrossLociSurf(struct LsurfPoint *P,struct SamplerInfo * S,double * s);
void MC_CalcPtsig2MarginAcrossLociSurfFast(struct LsurfPoint *P,struct SamplerInfo * S,struct LsurfPoint ** PerLocusSurf,int sCtr);
void MC_LikelihoodSurface1d(struct SamplerInfo * S);
void MC_LikelihoodSurfaceHelper(struct SamplerInfo * S, int D, double * svals,int * Ctr);
void MC_LikelihoodSurface(struct SamplerInfo * S);
void MC_DrawReplicates(struct SamplerInfo * S);
void OutputLikelihoodSurfaces(FILE *outfile,struct SamplerInfo S);

void OutputMLEs(FILE * outfile,struct SamplerInfo S);
struct MC_gslwrap_Params { struct SamplerInfo *S; int LocusNo; double offset;int CalcCode; double sig2base;};
double MC_gslwrap_1dPerLocus(double sig2,void * p);
double MC_gslwrap_1dAcrossLoci(double sig2,void * params);
void MC_1dMaximize(struct SamplerInfo * S, double (*f)(double,void *), struct MC_gslwrap_Params p, int CalcCI,struct MLE *);

double MC_gslwrap_f_2dPerLocus(const gsl_vector * x, void * p);
void MC_gslwrap_df_2dPerLocus(const gsl_vector * x, void * p, gsl_vector *g);
void MC_gslwrap_fdf_2dAcrossLoci(const gsl_vector * x, void * p, double * f, gsl_vector *g);
void MC_gslwrap_fdf_2dPerLocus(const gsl_vector * x, void * p, double * f, gsl_vector *g);
void MC_2dMaximize(struct SamplerInfo * S, double(*f)(const gsl_vector *,void *),void (*df)(const gsl_vector *x,void *,gsl_vector *),
 void (*fdf)(const gsl_vector *, void *,double *, gsl_vector *), double(*fCI)(double, void *),
 struct MC_gslwrap_Params p, double sig2x,double sig2y,int CalcCI,struct MLE * MLE);
void MC_PointEstimates(struct SamplerInfo * S);
/*void CalcSurface(struct BaseInfo B,struct LocusData LD, struct LocusParams L,struct ModelParamsPtr M, struct LocusParamsPtr LP, struct LsurfPoint * Lsurf,int n);
void OutputLikelihoodSurface(FILE * outfile,struct BaseInfo B,struct LocusData LD,struct LsurfPoint * Lsurf);*/
int main(int argc,char **argv){

  struct SamplerInfo S;
  
  int d;
  FILE * outfile;
  char filename[100];
   int i,k;
  
  
  /* Read Settings and Data */
  S=InitSettingsAndData(argc,argv);
  
  k=0;

  /* Initialize driving values */
  for(d=0;d<S.B.Dim;d++){
    S.MP_IS.muPtr[d]->v=0.0;
    S.MP_IS.sig2Ptr[d]->v=S.B.sig2Drive;
  }
  
  /* Initialize topo and times */
  if(S.B.UserInputTopoTimes){
    ReadTreeTimes(S.B,S.LD,S.L);
    for(i=0;i<S.B.nLoci;i++)
      Calc_xbar_vprime(&S.L[i].NodeList[0],S.B,S.LD[i],&S.L[i]);
  }else{
    for(i=0;i<S.B.nLoci;i++){
      if(S.B.FixTimes)
        SampleTimesRandom(S.B,S.LD[i],&S.L[i]);
      if(S.B.FixTopo)
        SampleTopoRandom(S.B,S.LD[i],&S.L[i]);  
    }
  }
    
    
  //getchar();
  MC_DrawReplicates(&S);
  if(S.B.CalcMultiDimSurf)
    MC_LikelihoodSurface(&S);  
  else if(S.B.Calc1DimSurf)
    MC_LikelihoodSurface1d(&S);  

  if(S.B.Calc1DimSurf ||S.B.CalcMultiDimSurf){
    strcpy(filename,S.B.FileNameBase);
    strcat(filename,".MCout");
    if(!(outfile=fopen(filename,"w"))){
      fprintf(stderr,"Error opening output file.\n");
      exit(1);
    }

    
    OutputLikelihoodSurfaces(stdout, S);
    OutputLikelihoodSurfaces(outfile, S);
    fclose(outfile);
    
  }
  
  
  MC_PointEstimates(&S);
  strcpy(filename,S.B.FileNameBase);
  strcat(filename,".MLEout");
  if(!(outfile=fopen(filename,"w"))){
    fprintf(stderr,"Error opening output file.\n");
    exit(1);
  }

  OutputMLEs(stdout,S);
  OutputMLEs(outfile,S);
  fclose(outfile);
  //getchar();
 
    /* Normalize the surface *//*
    for(j=0;j<B.nLsurfPoints;j++){
      S.Lsurf[i][j].LogL-=log(S.B.MaxReps);
      for(k=2;k<=S.LD[i].n;k++){
        S.Lsurf[i][j].LogL+=log(k);
      }
      for(k=0;k<B.Dim;k++){
          Lsurf[i][j].LogL+=log(B.sig2Max[k]-B.sig2Min[k]);
        if(B.muMax[k]-B.muMin[k]>0)
          Lsurf[i][j].LogL+=log(B.muMax[k]-B.muMin[k]);
        if(!B.Marginalize_x0)
            if(B.x0Max[k]-B.x0Min[k]>0)
              Lsurf[i][j].LogL+=log(B.x0Max[k]-B.x0Min[k]);
      }
    }
    
  

  
  for(i=0;i<S.B.nLoci;i++)
    OutputLikelihoodSurface(stdout,S.B,S.LD[i],S.Lsurf[i]);
 
  if(!(outfile=fopen(strcat(S.B.FileNameBase,".MCout"),"w"))){
    fprintf(stderr,"Error opening output file.\n");
    exit(1);
  }
  OutputLikelihoodSurface(outfile,S.B,S.LD[0],S.Lsurf[0]);
  */

  getsd(&S.B.seed1,&S.B.seed2);
  strcpy(filename,S.B.FileNameBase);
  strcat(filename,".seeds");
  outfile=fopen(filename,"w");
  fprintf(outfile,"%ld %ld\n",S.B.seed1,S.B.seed2);
  fclose(outfile);
  return(0);
 
}



void MC_DrawReplicates(struct SamplerInfo * S){

  int i,d,n;
   
  double LogISweight;

  
  for(i=0;i<S->B.nLoci;i++){
    /* Loop over replicates */
    for(n=0;n<S->B.MaxReps;n++){
  
      if(n%(S->B.MaxReps/10+1)==0)
        printf(" --Locus %d Replicate %d--\n",i,n);

      LogISweight=0.0;

      /* Loop over loci */
    
            
      
      
      /* Draw topos, times */
      if(!S->B.FixTimes)
        SampleTimesRandom(S->B,S->LD[i],&S->L[i]);
      if(S->B.FixTopo){
        Calc_xbar_vprime(&S->L[i].NodeList[0],S->B,S->LD[i],&S->L[i]);
      }else{
        if(S->B.Debug>5)
          printf("Drawing topologies... %d \n",n);

        if(S->B.UseIS){
          SampleTopoIS(S->B.Dim,S->LD[i],S->LP[i],S->LP[i],S->MP_IS,S->B.TopoISHeat);
        }else{
          SampleTopoRandom(S->B,S->LD[i],&S->L[i]);
        }
      }
      LogISweight=S->L[i].NodeList[0].Data->Logwr; /* Set IS weight */
      S->meanISweight[i]+=exp(LogISweight);
      S->meanISweightSquare[i]+=exp(2.0*LogISweight);
      
    
      /* Calculate summary statistic for replicate */
      S->L[i].NodeList[0].Data->CalcFlag=0; /* So that x0 term is not calculated in CalcLikelhoodExpo*/
      for(d=0;d<S->B.Dim;d++){
        S->RepSumStats[i][n].sumSig[d]=CalcLikelihoodExpo(S->LP[i].NodeListPtr[0],d,S->MP,S->LP[i]);
        S->RepSumStats[i][n].xbar[d]=S->L[i].NodeList[0].lChild->Data->xbar[d];
      }
      S->RepSumStats[i][n].LogDet=CalcLogVarCovarDet(&S->L[i].NodeList[0],S->MP,S->LP[i]);
      S->L[i].NodeList[0].Data->CalcFlag=1; /* Reset CalcFlag */
      
      S->RepSumStats[i][n].vMRCA=S->LP[i].t1Ptr->v-S->L[i].NodeList[0].lChild->t;
      S->RepSumStats[i][n].S2_MRCA=S->L[i].NodeList[0].lChild->Data->S2;
      
      S->RepSumStats[i][n].t1=S->LP[i].t1Ptr->v;
      S->RepSumStats[i][n].LogISweight=LogISweight;
      
      /* Store sumSamSizes,LogDet,xbar,t1-t2,S2, ISweight */
        
      if(S->B.Debug>5){
        for(d=0;d<S->B.Dim;d++){
          printf("sumSig %g (%g), xbar %g ",
            S->RepSumStats[i][n].sumSig[d],S->RepSumStats[i][n].sumSig[d]/2.0,S->RepSumStats[i][n].xbar[d]);
        }
        printf("LogDet %g LogISweight %g %g %g %g\n",S->RepSumStats[i][n].LogDet,S->RepSumStats[i][n].LogISweight,S->RepSumStats[i][n].vMRCA,S->RepSumStats[i][n].S2_MRCA,S->RepSumStats[i][n].t1);
      }
   
    } /* End replicates*/
    S->meanISweight[i]/=S->B.MaxReps;
    S->meanISweightSquare[i]/=S->B.MaxReps;
    
  } /* End loop over loci */
}





void MC_CalcPtsig2MarginSurf(struct LsurfPoint * P,struct SamplerInfo * S,int LocusNo,double * s){

  double LogL,Ltot,LsumSquare,ISnormConstant;
  int n,d;
  
  /*for(d=0;d<S->B.Dim;d++)
    printf("%g ",s[d]);
  printf(" In CalcPtsig2mMarginSurf for LocusNo %d %d %d %d\n",LocusNo,P->sig2,P->mu,P->x0);*/
  
  for(d=0;d<S->B.Dim;d++){
    if(P->sig2)
      P->sig2[d]=s[d];
    if(P->mu)
      P->mu[d]=0.0/0.0;
    if(P->x0)
      P->x0[d]=0.0/0.0;
  }
  P->f=S->L[LocusNo].f.v;
  P->LocusNo=LocusNo;

  Ltot=0.0;
  LsumSquare=0.0;
  
  for(n=0;n<S->B.MaxReps;n++){
    LogL=0.0;
    for(d=0;d<S->B.Dim;d++){
      LogL+=-((S->LD[LocusNo].n-1)*0.5)*log(2.0*PI*s[d])-0.5*S->RepSumStats[LocusNo][n].LogDet-0.5*S->RepSumStats[LocusNo][n].sumSig[d]/s[d];
    }
    Ltot+=exp(LogL+S->RepSumStats[LocusNo][n].LogISweight);
    LsumSquare+=exp(2*(LogL+S->RepSumStats[LocusNo][n].LogISweight));
    //printf("%d %g\n",sCtr,exp(LogL+S->RepSumStats[i][n].LogISweight+LogFactorial(S->LD[i].n)),LsumSquare*exp(LogFactorial(S->LD[i].n)));
  }
  if(S->B.UnbiasedISestimate)
    ISnormConstant=S->B.MaxReps; /* the number of replicates */
  else
    ISnormConstant=(S->meanISweight[LocusNo])*S->B.MaxReps; /* the sum of the IS weights */
  
  P->LogL=log(Ltot)-log(ISnormConstant)+LogFactorial(S->LD[LocusNo].n);
  
  /* Compute sampling variance of the likelihood point estimate */
  /* the square of the factorial term is necessary to complete LsumSquare */
  /* c^2*E[X^2] */
  /* Following Eqn on p. 35 of Liu 2004 */
  P->LsumSquare=(exp(2.0*LogFactorial(S->LD[LocusNo].n))*LsumSquare)/ISnormConstant;
  P->LSE=sqrt((P->LsumSquare-exp(2.0*P->LogL))/S->B.MaxReps);
        
  return;
}


void MC_CalcPtsig2MarginAcrossLociSurf(struct LsurfPoint * P,struct SamplerInfo * S,double * s){
    
  struct LsurfPoint x;
  int d,i;
  double z;
  /* Set to null pointers we won't use */
  x.mu=0x0;
  x.sig2=0x0;
  x.x0=0x0;
  
  for(d=0;d<S->B.Dim;d++){
    if(P->sig2)
      P->sig2[d]=s[d];
    if(P->mu)
      P->mu[d]=0.0/0.0;
    if(P->x0)
      P->x0[d]=0.0/0.0;
  }


  P->f=0.0/0.0;
  P->LocusNo=-1;
  P->LsumSquare=1.0;
  P->LogL=0;
  z=0.0;

  for(i=0;i<S->B.nLoci;i++){
    MC_CalcPtsig2MarginSurf(&x,S,i,s);
    P->LogL+=x.LogL;
    z+=x.LSE*x.LSE/exp(2.0*x.LogL);

  }
  P->LSE=sqrt(P->LsumSquare-exp(2.0*P->LogL));
  /* SE based on based on propogation of error, taking advantage of fact that Monte Carlo error from locus to locus is indpt */
  P->LSE=sqrt(z)*exp(P->LogL);
  P->LsumSquare=0.0/0.0;
  return;
}

void MC_CalcPtsig2MarginAcrossLociSurfFast(struct LsurfPoint * P,struct SamplerInfo * S,struct LsurfPoint ** PerLocusSurf,int sCtr){
  int d,i;
  double z;
  for(d=0;d<S->B.Dim;d++){
    P->sig2[d]=PerLocusSurf[0][sCtr].sig2[d];
    P->mu[d]=0.0/0.0;
    P->x0[d]=0.0/0.0;
  }

  P->f=0.0/0.0;
  P->LocusNo=-1;
  P->LsumSquare=1.0;
  P->LogL=0;
  z=0.0;
  for(i=0;i<S->B.nLoci;i++){
      P->LogL+=PerLocusSurf[i][sCtr].LogL;
      z+=PerLocusSurf[i][sCtr].LSE*PerLocusSurf[i][sCtr].LSE/exp(2.0*PerLocusSurf[i][sCtr].LogL);
  }
  P->LsumSquare=0.0/0.0;
  /* SE based on based on propogation of error, taking advantage of fact that Monte Carlo error from locus to locus is indpt */
  P->LSE=sqrt(z)*exp(P->LogL);
  return;
}


  
void MC_LikelihoodSurfaceHelper(struct SamplerInfo * S, int D, double * svals,int * Ctr){
  int sCtr,i;
  double s;
  if(D>=0){
    for(sCtr=0,s=S->B.sig2Min[D];sCtr<S->B.sig2nvals[D];s+=S->B.sig2Delt[D],sCtr++){
      svals[D]=s;
      //printf("Making a recursive call...\t%d %d ",D,*Ctr);
      //for(d=0;d<S->B.Dim;d++){
      //  printf("%g ",svals[d]);
      //}
      //printf("\n");
 
      MC_LikelihoodSurfaceHelper(S,D-1,svals,Ctr);
    }    
  }else {
    //printf("In main section... \t\t%d %d ",D, *Ctr);
    //for(d=0;d<S->B.Dim;d++){
    //  printf("%g ",svals[d]);
    //}
    //printf("\n");
    for(i=0;i<S->B.nLoci;i++){
      MC_CalcPtsig2MarginSurf(&S->sig2MarginSurf[i][*Ctr],S,i,svals);
    }
    MC_CalcPtsig2MarginAcrossLociSurfFast(&S->sig2MarginAllLociSurf[*Ctr],S,S->sig2MarginSurf,*Ctr);
    *Ctr=*Ctr+1;
  }
}

void MC_LikelihoodSurface(struct SamplerInfo * S){

  double * svals;
  svals=(double *) calloc((size_t)S->B.Dim,sizeof(double));
  int Ctr=0;
  int D=S->B.Dim-1;
  printf("Calculating multi-dimensional likelihood surface...\n");
  MC_LikelihoodSurfaceHelper(S,D,svals,&Ctr);

  free(svals);
}



void MC_LikelihoodSurface1d(struct SamplerInfo * S){

  int d,sCtr,i;
  double s;
  double * svals;
  svals=(double *) calloc((size_t)S->B.Dim,sizeof(double));

  for(sCtr=0,s=S->B.sig2Min1d;sCtr<S->B.sig2nvals1d;s+=S->B.sig2Delt1d,sCtr++){
    for(d=0;d<S->B.Dim;d++){
      svals[d]=s;
    }

    for(i=0;i<S->B.nLoci;i++){
      MC_CalcPtsig2MarginSurf(&S->sig2Margin1dSurf[i][sCtr],S,i,svals);
    }
    MC_CalcPtsig2MarginAcrossLociSurfFast(&S->sig2MarginAllLoci1dSurf[sCtr],S,S->sig2Margin1dSurf,sCtr);
    
  }
  free(svals);
}




void MC_LikelihoodSurfaceOld(struct SamplerInfo * S){
  int i,d,n;

 
  double LogL,Ltot,LsumSquare;
  double ISnormConstant;
  int sCtr;
  double s;
  LogL=0.0;
  /* Post process stored values to build surfaces... */
  printf("Building sig2Margin likelihood surfaces\n");
  /* First create sig2MarginSurf */
  for(i=0;i<S->B.nLoci;i++){
   
    for(sCtr=0,s=S->B.sig2Min1d;sCtr<S->B.sig2nvals1d;s+=S->B.sig2Delt1d,sCtr++){
        for(d=0;d<S->B.Dim;d++){
          S->sig2Margin1dSurf[i][sCtr].sig2[d]=s;
          S->sig2Margin1dSurf[i][sCtr].mu[d]=0.0/0.0;
          S->sig2Margin1dSurf[i][sCtr].x0[d]=0.0/0.0;
        }
        S->sig2Margin1dSurf[i][sCtr].f=S->L[i].f.v;
        S->sig2Margin1dSurf[i][sCtr].LocusNo=i;
        Ltot=0.0;
        LsumSquare=0.0;
        for(n=0;n<S->B.MaxReps;n++){
          //LogL=-((S->LD[i].n-1)*0.5)*log(2.0*PI*s)-0.5*log(exp(S->RepSumStats[i][n].LogDet)*2.0)-0.5*S->RepSumStats[i][n].sumSig[d]/(s*2.0);
          for(d=0;d<S->B.Dim;d++){
            LogL=-((S->LD[i].n-1)*0.5)*log(2.0*PI*s)-0.5*S->RepSumStats[i][n].LogDet-0.5*S->RepSumStats[i][n].sumSig[d]/s;
          }
          Ltot+=exp(LogL+S->RepSumStats[i][n].LogISweight);
          LsumSquare+=exp(2*(LogL+S->RepSumStats[i][n].LogISweight));
          //printf("%d %g\n",sCtr,exp(LogL+S->RepSumStats[i][n].LogISweight+LogFactorial(S->LD[i].n)),LsumSquare*exp(LogFactorial(S->LD[i].n)));
        }
        if(S->B.UnbiasedISestimate)
          ISnormConstant=S->B.MaxReps; /* the number of replicates */
        else
          ISnormConstant=(S->meanISweight[i])*S->B.MaxReps; /* the sum of the IS weights */
        
        S->sig2Margin1dSurf[i][sCtr].LogL=log(Ltot)-log(ISnormConstant)+LogFactorial(S->LD[i].n);
        
        /* Compute sampling variance of the likelihood point estimate */
        /* the square of the factorial term is necessary to complete LsumSquare */
        /* c^2*E[X^2] */
        S->sig2Margin1dSurf[i][sCtr].LsumSquare=(exp(2.0*LogFactorial(S->LD[i].n))*LsumSquare)/ISnormConstant;
        //S->sig2Margin1dSurf[i][sCtr].LsVar=S->sig2Margin1dSurf[i][sCtr].LsumSquare-exp(2.0*S->sig2Margin1dSurf[i][sCtr].LogL);
        
        //printf("%g %g %g\n",S->sig2MarginSurf[i][d][sCtr].LogL,S->sig2MarginSurf[i][d][sCtr].LsumSquare,S->sig2MarginSurf[i][d][sCtr].LsVar);
      
    }
  }
  printf("Summing across loci\n"); 
  /* Then sig2MarginAllLoci */
  for(sCtr=0,s=S->B.sig2Min1d;sCtr<S->B.sig2nvals1d;s+=S->B.sig2Delt1d,sCtr++){
    for(d=0;d<S->B.Dim;d++){
      S->sig2MarginAllLoci1dSurf[sCtr].sig2[d]=s;
      S->sig2MarginAllLoci1dSurf[sCtr].mu[d]=0.0/0.0;
      S->sig2MarginAllLoci1dSurf[sCtr].x0[d]=0.0/0.0;
    }
    S->sig2MarginAllLoci1dSurf[sCtr].f=S->L[i].f.v;
    S->sig2MarginAllLoci1dSurf[sCtr].LocusNo=-1;
    S->sig2MarginAllLoci1dSurf[sCtr].LsumSquare=1.0;
    for(i=0;i<S->B.nLoci;i++){
      S->sig2MarginAllLoci1dSurf[sCtr].LogL+=S->sig2Margin1dSurf[i][sCtr].LogL;
      S->sig2MarginAllLoci1dSurf[sCtr].LsumSquare*=S->sig2Margin1dSurf[i][sCtr].LsumSquare;
    }
    //S->sig2MarginAllLoci1dSurf[sCtr].LsVar=S->sig2MarginAllLoci1dSurf[sCtr].LsumSquare-exp(2.0*S->sig2MarginAllLoci1dSurf[sCtr].LogL);
  }
 
 /*if(S->B.Outputx0sig2muSurf){
    // LEGACY CODE...
    printf("Building x0sig2mu likelihood surfaces\n");
    // Then x0sig2musurf 
    for(d=0;d<S->B.Dim;d++){
      for(i=0;i<S->B.nLoci;i++){
    
        for(q=0,sCtr=0,s=S->B.sig2Min[d];q<S->B.sig2nvals[d];s+=S->B.sig2Delt[d],q++){
          for(p=0,m=S->B.muMin[d];p<S->B.munvals[d];m+=S->B.muDelt[d],p++){
            for(r=0,x=S->B.x0Min[i][d];r<S->B.x0nvals[i][d];x+=S->B.x0Delt[i][d],sCtr++,r++){
              S->x0sig2muSurf[i][d][sCtr].sig2=s;
              S->x0sig2muSurf[i][d][sCtr].mu=m;
              S->x0sig2muSurf[i][d][sCtr].x0=x;
              S->x0sig2muSurf[i][d][sCtr].f=S->L[i].f.v;
              S->x0sig2muSurf[i][d][sCtr].LocusNo=i;
              S->x0sig2muSurf[i][d][sCtr].d=d;
              Ltot=0.0;
              for(n=0;n<S->B.MaxReps;n++){
              
                sumSigp=S->RepSumStats[i][n].sumSig[d]+pow(S->RepSumStats[i][n].xbar[d]-x-m*(S->RepSumStats[i][n].vMRCA+S->RepSumStats[i][n].S2_MRCA),2)/(S->RepSumStats[i][n].vMRCA+S->RepSumStats[i][n].S2_MRCA);
                sumSamSizesp=S->LD[i].n;
                LogDetp=S->RepSumStats[i][n].LogDet+log(S->RepSumStats[i][n].vMRCA+S->RepSumStats[i][n].S2_MRCA);
                LogL=-(0.5*sumSamSizesp)*log(2*PI*s)-0.5*LogDetp-0.5*sumSigp/s;
                Ltot+=exp(LogL+S->RepSumStats[i][n].LogISweight);
                if(S->B.Debug>6)  
                  printf("m: %g s: %g x:%g sumSigp: %g (%g) sumSizesp %g LogDetp %g NormFactor %g LogL %g LogISweight %g\n",m,s,x,sumSigp,sumSigp/2.0,sumSamSizesp,LogDetp,exp(-(0.5*sumSamSizesp)*log(2.0*PI)-0.5*LogDetp),LogL,S->RepSumStats[i][n].LogISweight);
              }
              S->x0sig2muSurf[i][d][sCtr].LogL=log(Ltot)-log(S->B.MaxReps)+LogFactorial(S->LD[i].n);
           
            }
            
          }
        }

      }
    }
  }*/
            
          
          /* Use Log Sum Exp trick: First find Bp */
  /*        Bp=S->sig2muSurf[d][sCtr][0].LogL;
          for(j=1;j<S->B.MaxReps;j++){
            if(S->sig2muSurf[d][sCtr][j].LogL>Bp)
              Bp=S->sig2muSurf[d][sCtr][j].LogL;
	      }*/
          /* Calculate sum of exponentials */
  /*         SumExp=0.0;
          for(j=0;j<S->B.MaxReps;j++){
            SumExp+=exp(S->sig2muSurf[d][sCtr][j].LogL-Bp);
	    }*/
          /* Take log */
  //        S->sig2muSurf[d][sCtr][0].LogL=log(SumExp)+Bp-log(S->B.MaxReps);
    
}


void OutputLikelihoodSurfaces(FILE * outfile,struct SamplerInfo S){
  
  int i,j,d;
  double SampleVarISweight;
  fprintf(outfile,"# MaxReps %d IS: %d ISheat: %g sig2Drive: %g UnbiasedISestimate: %d FixTopo: %d FixTimes: %d Seeds: %ld %ld\n",S.B.MaxReps,S.B.UseIS,S.B.TopoISHeat,S.B.sig2Drive,S.B.UnbiasedISestimate,S.B.FixTopo,S.B.FixTimes,S.B.seed1,S.B.seed2);
  for(i=0;i<S.B.nLoci;i++){
    SampleVarISweight=(S.meanISweightSquare[i]-pow(S.meanISweight[i],2))*(S.B.MaxReps/(S.B.MaxReps-1.0));
    fprintf(outfile,"# Locus %s meanISweight %g VarISweight: %g ESS: %g\n",S.LD[i].ID,S.meanISweight[i],SampleVarISweight,S.B.MaxReps/(1+SampleVarISweight));
  }
  
  if(S.B.CalcMultiDimSurf){
    if(S.B.Outputsig2MarginSurf){
      for(i=0;i<S.B.nLoci;i++){
        for(j=0;j<S.B.sig2nvalsTot;j++){
          fprintf(outfile,"%s ",S.LD[i].ID);
          for(d=0;d<S.B.Dim;d++){
            fprintf(outfile,"%g ",S.sig2MarginSurf[i][j].sig2[d]);
          }
          fprintf(outfile,"%g %g %g\n",S.sig2MarginSurf[i][j].f,S.sig2MarginSurf[i][j].LogL,
            log(S.sig2MarginSurf[i][j].LSE));
        }
      }
    }
     
       
    if(S.B.Outputsig2MarginAllLociSurf){
      for(j=0;j<S.B.sig2nvalsTot;j++){
          fprintf(outfile,"AllLoci ");
          for(d=0;d<S.B.Dim;d++){
            fprintf(outfile,"%g ",S.sig2MarginAllLociSurf[j].sig2[d]);
          }
          fprintf(outfile,"%g %g %g\n",S.sig2MarginAllLociSurf[j].f,S.sig2MarginAllLociSurf[j].LogL,
            log(S.sig2MarginAllLociSurf[j].LSE));
      }
    }
  }else if(S.B.Calc1DimSurf){

    
    if(S.B.Outputsig2MarginSurf){
      for(i=0;i<S.B.nLoci;i++){
        for(j=0;j<S.B.sig2nvals1d;j++){
          fprintf(outfile,"%s %g %g %g %g\n",S.LD[i].ID,S.sig2Margin1dSurf[i][j].sig2[0],
           S.sig2Margin1dSurf[i][j].f,S.sig2Margin1dSurf[i][j].LogL,log(S.sig2Margin1dSurf[i][j].LSE));
        }
      }
    }
      
    if(S.B.Outputsig2MarginAllLociSurf){
      for(j=0;j<S.B.sig2nvals1d;j++){
          fprintf(outfile,"AllLoci %g %g %g %g\n",S.sig2MarginAllLoci1dSurf[j].sig2[0],
            S.sig2MarginAllLoci1dSurf[j].f,S.sig2MarginAllLoci1dSurf[j].LogL,log(S.sig2MarginAllLoci1dSurf[j].LSE));
      }
    }
  }
      
          
/*    if(S.B.Outputx0sig2muSurf){

      for(i=0;i<S.B.nLoci;i++){
        for(j=0;j<S.B.nx0sig2muSurfPoints[i][d];j++){
          fprintf(outfile,"%s %d",S.LD[i].ID,d);
          fprintf(outfile," %g %g %g  ",S.x0sig2muSurf[i][d][j].sig2,S.x0sig2muSurf[i][d][j].mu,S.x0sig2muSurf[i][d][j].x0);
          fprintf(outfile,"%g %g\n",S.x0sig2muSurf[i][d][j].f,(S.x0sig2muSurf[i][d][j].LogL));
        }
      }
    }
    */
  
}


void OutputMLEs(FILE * outfile,struct SamplerInfo S){
  int i,d;
  double SampleVarISweight;
  fprintf(outfile,"# MaxReps %d IS: %d ISheat: %g sig2Drive: %g UnbiasedISestimate: %d FixTopo: %d FixTimes: %d Seeds: %ld %ld\n",S.B.MaxReps,S.B.UseIS,S.B.TopoISHeat,S.B.sig2Drive,S.B.UnbiasedISestimate,S.B.FixTopo,S.B.FixTimes,S.B.seed1,S.B.seed2);
  for(i=0;i<S.B.nLoci;i++){
    SampleVarISweight=(S.meanISweightSquare[i]-pow(S.meanISweight[i],2))*(S.B.MaxReps/(S.B.MaxReps-1.0));
    fprintf(outfile,"# Locus %s meanISweight %g VarISweight: %g ESS: %g\n",S.LD[i].ID,S.meanISweight[i],SampleVarISweight,S.B.MaxReps/(1+SampleVarISweight));
  }
  fprintf(outfile,"# GSL_Params: lower_bound %g upper_bound %g Relative_Error %g Max_Iterations %g Eps_for_finite_difference %g\n",S.B.gsla,S.B.gslb,S.B.gslRelErr,S.B.gslMaxIter,S.B.gsldfeps);
  fprintf(outfile,"# Elapsed_Time: %ld Start_Time: %s",(long)(time(NULL)-S.B.StartTime),ctime(&S.B.StartTime));
//  fprintf("LocusID");

  if(S.B.CalcMLEperLocus){
    for(i=0;i<S.B.nLoci;i++){
    
      fprintf(outfile,"%s ",S.LD[i].ID);
      fprintf(outfile,"%g %g %g ",S.sig2Margin1dMLE[i].v,S.sig2Margin1dMLE[i].lowerCI,S.sig2Margin1dMLE[i].upperCI);
      
      for(d=0;d<S.B.Dim;d++){
        fprintf(outfile,"%g %g %g ",S.sig2MarginMLE[i][d].v,S.sig2MarginMLE[i][d].lowerCI,S.sig2MarginMLE[i][d].upperCI);
      }
      fprintf(outfile,"%g %g %g %g",S.sig2Margin1dMLE[i].LogL,S.sig2Margin1dMLE[i].LogLSE,S.sig2MarginMLE[i][0].LogL,S.sig2MarginMLE[i][0].LogLSE);
      fprintf(outfile,"\n");
    }
  }
  fprintf(outfile,"AllLoci ");
  fprintf(outfile,"%g %g %g ",S.sig2MarginAllLoci1dMLE.v,S.sig2MarginAllLoci1dMLE.lowerCI,S.sig2MarginAllLoci1dMLE.upperCI);
    
  for(d=0;d<S.B.Dim;d++){
    fprintf(outfile,"%g %g %g ",S.sig2MarginAllLociMLE[d].v,S.sig2MarginAllLociMLE[d].lowerCI,S.sig2MarginAllLociMLE[d].upperCI);
  }
  fprintf(outfile,"%g %g ",S.sig2MarginAllLoci1dMLE.LogL,S.sig2MarginAllLoci1dMLE.LogLSE);
  fprintf(outfile,"%g %g ",S.sig2MarginAllLociMLE[0].LogL,S.sig2MarginAllLociMLE[0].LogLSE);

  fprintf(outfile,"\n");

  
};
  

struct SamplerInfo InitSettingsAndData(int argc,char **argv){ 
  
  
  
  struct SamplerInfo S;
  int i,j,k,d;

  S.B.FixTopo=0;
  S.B.FixTimes=0;
  S.B.Debug=1;
  S.B.MaxReps=5000;
  S.B.Dim=2;
  S.B.Prior_mu=0;
  S.B.Prior_sig2=0;
  S.B.Prior_x0=0;
  S.B.Outputsig2MarginSurf=1;
  S.B.Outputsig2MarginAllLociSurf=1;
  S.B.Outputx0sig2muSurf=0;
  S.B.LongLikelihoodCalc=0;
  S.B.UserInputTopoTimes=0;
  S.B.Calc1DimSurf=0;
  S.B.CalcMultiDimSurf=0;
  S.B.CalcMLEperLocus=0;
  S.B.sMin=0;
  S.B.sMax=300;
  S.B.snvals=11;
  
  S.B.gsla=1e-10;
  S.B.gslb=1e10;
  S.B.gslRelErr=0.01;
  S.B.gslMaxIter=5000;
  S.B.gslCIperLocus=1;

  S.B.UseIS=1;
  S.B.UnbiasedISestimate=1;
  S.B.TopoISHeat=1.0;
  S.B.sig2Drive=25;
  S.B.muDrive=0.0;
  
  
  S.B.StartTime=time(NULL);
  S.B.seed1=(long)S.B.StartTime / 834  + 7737;
  S.B.seed2=(S.B.seed1 + 438977)  / 264;

  
  printf("Reading commandline...\n");
  ReadCommandLine(argc,argv,&S.B);
  strcpy(S.B.FileNameBase,S.B.InFileName);

  S.B.gsldfeps=(1.0/S.B.sig2Drive)*1e-2;

  S.B.sig2Min=(double *)calloc((size_t)S.B.Dim,sizeof(double));
  S.B.sig2Max=(double *)calloc((size_t)S.B.Dim,sizeof(double));
  S.B.sig2nvals=(double *)calloc((size_t)S.B.Dim,sizeof(double));
  S.B.sig2Delt=(double *)calloc((size_t)S.B.Dim,sizeof(double));

  S.B.muMin=(double *)calloc((size_t)S.B.Dim,sizeof(double));
  S.B.muMax=(double *)calloc((size_t)S.B.Dim,sizeof(double));
  S.B.munvals=(double *)calloc((size_t)S.B.Dim,sizeof(double));
  S.B.muDelt=(double *)calloc((size_t)S.B.Dim,sizeof(double));


  S.B.x0Min=(double **)calloc((size_t)S.B.nLoci,sizeof(double));
  S.B.x0Max=(double **)calloc((size_t)S.B.nLoci,sizeof(double));
  S.B.x0nvals=(double **)calloc((size_t)S.B.nLoci,sizeof(double));
  S.B.x0Delt=(double **)calloc((size_t)S.B.nLoci,sizeof(double));
  for(i=0;i<S.B.nLoci;i++){
    S.B.x0Min[i]=(double *)calloc((size_t)S.B.Dim,sizeof(double));
    S.B.x0Max[i]=(double *)calloc((size_t)S.B.Dim,sizeof(double));
    S.B.x0nvals[i]=(double *)calloc((size_t)S.B.Dim,sizeof(double));
    S.B.x0Delt[i]=(double *)calloc((size_t)S.B.Dim,sizeof(double));
  }

  S.B.sig2Min1d=S.B.sMin;
  S.B.sig2Max1d=S.B.sMax;
  S.B.sig2nvals1d=S.B.snvals;
  if(S.B.sig2nvals1d>1)
    S.B.sig2Delt1d=(S.B.sig2Max1d-S.B.sig2Min1d)/((double)S.B.sig2nvals1d-1.0);
  else 
    S.B.sig2Delt1d=1.0;
      
  for(d=0;d<S.B.Dim;d++){


    
    for(i=0;i<S.B.nLoci;i++){
      S.B.x0Min[i][d]=0.0;
      S.B.x0Max[i][d]=0.0;
      //S.B.x0Min[i][d]=2861.0;
      //S.B.x0Max[i][d]=2861.0;
      
      S.B.x0nvals[i][d]=1;
      if(S.B.x0nvals[i][d]>1)
        S.B.x0Delt[i][d]=(S.B.x0Max[i][d]-S.B.x0Min[i][d])/((double)S.B.x0nvals[i][d]-1.0);
      else
        S.B.x0Delt[i][d]=1.0;
    }

    S.B.sig2Min[d]=S.B.sMin;
    S.B.sig2Max[d]=S.B.sMax;
    S.B.sig2nvals[d]=S.B.snvals;
    if(S.B.sig2nvals[d]>1)
      S.B.sig2Delt[d]=(S.B.sig2Max[d]-S.B.sig2Min[d])/((double)S.B.sig2nvals[d]-1.0);
    else
      S.B.sig2Delt[d]=1.0;

    
    S.B.muMin[d]=0.0;
    S.B.muMax[d]=0.0;
    S.B.munvals[d]=1;
    if(S.B.munvals[d]>1)
      S.B.muDelt[d]=(S.B.muMax[d]-S.B.muMin[d])/((double)S.B.munvals[d]-1.0);
    else
      S.B.muDelt[d]=1.0;
    

  }
  
  printf("Reading data...\n");

  // Declare Memory for Arrays
  S.LD=ReadData(S.B);

  printf("nLoci: %d\n ",S.B.nLoci);


  /* Allocate memory */
  S.L=AllocLocusParams(S.B,S.LD);
  S.LP=AllocLocusParamsPtr(S.B,S.LD);
  S.M=AllocModelParams(S.B);
  S.MP=AllocModelParamsPtr(S.B);
  S.M_IS=AllocModelParams(S.B);
  S.MP_IS=AllocModelParamsPtr(S.B);

  S.RepSumStats=(struct RepSummaryStats **) calloc((size_t)S.B.nLoci,sizeof(struct RepSummaryStats *));
  for(i=0;i<S.B.nLoci;i++){
    S.RepSumStats[i]=(struct RepSummaryStats *) calloc((size_t)S.B.MaxReps,sizeof(struct RepSummaryStats));
    if(!S.RepSumStats[i]){
      fprintf(stderr,"Malloc error\n");
      exit(-1);
    }
  }
  for(i=0;i<S.B.nLoci;i++){
    for(j=0;j<S.B.MaxReps;j++){
      S.RepSumStats[i][j].sumSig=(double *)calloc((size_t)S.B.Dim,sizeof(double));
      S.RepSumStats[i][j].xbar=(double *)calloc((size_t)S.B.Dim,sizeof(double));
    }
  }



  /* Declare mem sigs2Margin1dSurf */
  S.sig2Margin1dSurf=(struct LsurfPoint **) calloc((size_t)S.B.nLoci,sizeof(struct LsurfPoint *));
  for(i=0;i<S.B.nLoci;i++){
    S.sig2Margin1dSurf[i]=(struct LsurfPoint *) calloc((size_t)S.B.sig2nvals1d,sizeof(struct LsurfPoint));
    if(!S.sig2Margin1dSurf[i]){
      fprintf(stderr,"Malloc error\n");
      exit(-1);
    }
    for(j=0;j<S.B.sig2nvals1d;j++){
      S.sig2Margin1dSurf[i][j].sig2=(double *) calloc((size_t)S.B.Dim,sizeof(double));
      S.sig2Margin1dSurf[i][j].mu=(double *) calloc((size_t)S.B.Dim,sizeof(double));
      S.sig2Margin1dSurf[i][j].x0=(double *) calloc((size_t)S.B.Dim,sizeof(double));
    }
  }
    
  /* Declare mem sigs2MarginAllLoci */
  S.sig2MarginAllLoci1dSurf=(struct LsurfPoint *) calloc((size_t)S.B.sig2nvals1d,sizeof(struct LsurfPoint));
  for(j=0;j<S.B.sig2nvals1d;j++){
    S.sig2MarginAllLoci1dSurf[j].sig2=(double *) calloc((size_t)S.B.Dim,sizeof(double));
    S.sig2MarginAllLoci1dSurf[j].mu=(double *) calloc((size_t)S.B.Dim,sizeof(double));
    S.sig2MarginAllLoci1dSurf[j].x0=(double *) calloc((size_t)S.B.Dim,sizeof(double));
  }
  
  /* Now for multi-dimensional likelihood surfaces */

  S.B.sig2nvalsTot=1;
  for(d=0;d<S.B.Dim;d++){
    S.B.sig2nvalsTot*=S.B.sig2nvals[d];
  }
  printf("Total size of multi-dimensional MLE grid: %d\n",S.B.sig2nvalsTot);
  /* Declare mem sigs2MarginSurf */
  S.sig2MarginSurf=(struct LsurfPoint **) calloc((size_t)S.B.nLoci,sizeof(struct LsurfPoint *));
  for(i=0;i<S.B.nLoci;i++){
    S.sig2MarginSurf[i]=(struct LsurfPoint *) calloc((size_t)S.B.sig2nvalsTot,sizeof(struct LsurfPoint));
    if(!S.sig2MarginSurf[i]){
      fprintf(stderr,"Malloc error\n");
      exit(-1);
    }
    for(j=0;j<S.B.sig2nvalsTot;j++){
      S.sig2MarginSurf[i][j].sig2=(double *) calloc((size_t)S.B.Dim,sizeof(double));
      S.sig2MarginSurf[i][j].mu=(double *) calloc((size_t)S.B.Dim,sizeof(double));
      S.sig2MarginSurf[i][j].x0=(double *) calloc((size_t)S.B.Dim,sizeof(double));
    }
  }
    
  /* Declare mem sigs2MarginAllLoci */
  S.sig2MarginAllLociSurf=(struct LsurfPoint *) calloc((size_t)S.B.sig2nvalsTot,sizeof(struct LsurfPoint));
  for(j=0;j<S.B.sig2nvalsTot;j++){
    S.sig2MarginAllLociSurf[j].sig2=(double *) calloc((size_t)S.B.Dim,sizeof(double));
    S.sig2MarginAllLociSurf[j].mu=(double *) calloc((size_t)S.B.Dim,sizeof(double));
    S.sig2MarginAllLociSurf[j].x0=(double *) calloc((size_t)S.B.Dim,sizeof(double));
  }
  
  /* Declare mem for MLEs */
  S.sig2Margin1dMLE=(struct MLE *)calloc((size_t)S.B.nLoci,sizeof(struct MLE));
  S.sig2MarginMLE=(struct MLE**)calloc((size_t)S.B.nLoci,sizeof(struct MLE *));
  for(j=0;j<S.B.nLoci;j++){
    S.sig2MarginMLE[j]=(struct MLE*)calloc((size_t)S.B.Dim,sizeof(struct MLE));
  }
  S.sig2MarginAllLociMLE=(struct MLE *)calloc((size_t)S.B.Dim,sizeof(struct MLE));
  
  /* Declare mem x0sig2muSurf */
  /* LEGACY CODE 
  S.x0sig2muSurf=(struct LsurfPoint ***)calloc((size_t)S.B.nLoci,sizeof(struct LsurfPoint ** ));
  S.B.nx0sig2muSurfPoints=(int **)calloc((size_t)S.B.nLoci,sizeof(int *));

  for(i=0;i<S.B.nLoci;i++){
    S.x0sig2muSurf[i]=(struct LsurfPoint **)calloc((size_t)S.B.Dim,sizeof(struct LsurfPoint *));
    S.B.nx0sig2muSurfPoints[i]=(int *)calloc((size_t)S.B.Dim,sizeof(int));

    if(S.x0sig2muSurf[i]==0){
      fprintf(stderr,"Mem alloc error...\n");
      exit(-1);
    }
  }
  for(i=0;i<S.B.nLoci;i++){
    for(d=0;d<S.B.Dim;d++){
      S.B.nx0sig2muSurfPoints[i][d]=S.B.sig2nvals[d]*S.B.munvals[d]*S.B.x0nvals[i][d];
      S.x0sig2muSurf[i][d]=(struct LsurfPoint *)calloc((size_t)S.B.nx0sig2muSurfPoints[i][d],sizeof(struct LsurfPoint));
      if(S.x0sig2muSurf[i][d]==0){
        fprintf(stderr,"Mem alloc error...\n");
      exit(-1);
      }
    }
  }
  */

  /* Declare mem meanISweight meanISweightSquare */
  
  S.meanISweight=(double *)calloc((size_t)S.B.nLoci,sizeof(double));
  S.meanISweightSquare=(double *)calloc((size_t)S.B.nLoci,sizeof(double));
  for(i=0;i<S.B.nLoci;i++){
    S.meanISweight[i]=0.0;
    S.meanISweightSquare[i]=0.0;
  }

  /* Set-up pointers */
  for(i=0;i<S.B.Dim;i++){
    S.MP.sig2Ptr[i]=&(S.M.sig2[i]);
    S.MP.muPtr[i]=&(S.M.mu[i]);
    S.MP_IS.sig2Ptr[i]=&(S.M_IS.sig2[i]);
    S.MP_IS.muPtr[i]=&(S.M_IS.mu[i]);

  }
  for(i=0;i<S.B.nLoci;i++){
    for(j=0;j<S.B.Dim;j++){
      S.LP[i].x0Ptr[j]=&S.L[i].x0[j];
    }
    S.LP[i].t1Ptr=&S.L[i].t1;
    S.LP[i].fPtr=&S.L[i].f;
    S.L[i].f.v=S.LD[i].f;
    for(j=0;j<S.LD[i].n;j++){
      S.LP[i].NodeTimesPtr[j]=&S.L[i].NodeTimes[j];
      S.L[i].NodeList[j+1].key=j+1;
    }
    for(j=0;j<2*S.LD[i].n;j++){
      S.LP[i].NodeListPtr[j]=&S.L[i].NodeList[j];
    }
  }



  // Output Data
  if(S.B.Debug>0){
    for(i=0;i<S.B.nLoci;i++){
      printf("===Locus %d=== f:%g\n",i,S.LD[i].f);
      for(j=0;j<S.B.Dim;j++){
        printf("====Dimension: %d====\n",j+1);
        for(k=0;k<S.LD[i].n;k++){
          printf("%f ",S.LD[i].x[j][k]);
        }
        printf("\n");
      }
      printf("\n");
    }
  }

  return(S);
}



double MC_gslwrap_1dPerLocus(double sig2,void * p){
  struct MC_gslwrap_Params * params=(struct MC_gslwrap_Params *) p;
  int i;
  struct LsurfPoint x;
  double * s;
  
  /* Set Lsurf pointers that we won't use to null */
  //printf(" In MC_gslwrap_1dPerLocus for LocusNo %d %d %d %d\n",params->LocusNo,x.sig2,x.mu,x.x0);
  x.mu=0x0;
  x.sig2=0x0;
  x.x0=0x0;
  
  //printf(" In MC_gslwrap_1dPerLocus for LocusNo %d %d %d %d\n",params->LocusNo,x.sig2,x.mu,x.x0);
  
  s=(double *) calloc((size_t)params->S->B.Dim,sizeof(double));
  if(params->CalcCode==0){
    s[0]=sig2;
    s[1]=params->sig2base;
  } else if (params->CalcCode==1){
    s[0]=params->sig2base;
    s[1]=sig2;
  }else{
    /* Set all components of s to sig2 */
    for(i=0;i<params->S->B.Dim;i++){
      s[i]=sig2;
    }
  };
   
  MC_CalcPtsig2MarginSurf(&x,params->S,params->LocusNo,s);
  
  free(s);

  if(-x.LogL+params->offset==-x.LogL+params->offset && -x.LogL+params->offset!=1.0/0.0)
    return(-x.LogL+params->offset);
  else
    return(DBL_MAX);
}

double MC_gslwrap_1dAcrossLoci(double sig2,void * p){
  struct MC_gslwrap_Params * params=(struct MC_gslwrap_Params *) p;
  int i;
  struct LsurfPoint x;
  double * s;

  /* Set Lsurf pointers that we won't use to null */
  x.sig2=0x0;
  x.mu=0x0;
  x.x0=0x0;

  s=(double *) calloc((size_t)params->S->B.Dim,sizeof(double));
  if(params->CalcCode==0){
    s[0]=sig2;
    s[1]=params->sig2base;
  } else if (params->CalcCode==1){
    s[0]=params->sig2base;
    s[1]=sig2;
  }else{
    /* Set all components of s to sig2 */
    for(i=0;i<params->S->B.Dim;i++){
      s[i]=sig2;
    }
  };
  //printf("%g %g %g\n",s[0],s[1],-x.LogL+params->offset);
  MC_CalcPtsig2MarginAcrossLociSurf(&x,params->S,s);
  free(s);
  
  if(-x.LogL+params->offset==-x.LogL+params->offset && -x.LogL+params->offset!=1.0/0.0)
    return(-x.LogL+params->offset);
  else
    return(DBL_MAX);

 
}

void MC_1dMaximize(struct SamplerInfo * S, double (*f)(double, void*), struct MC_gslwrap_Params p, int CalcCI,struct MLE * MLE){
  int status;
  int iter = 0, max_iter = S->B.gslMaxIter;
  int iterMLE=0,iterCIl=0,iterCIu=0;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  const gsl_root_fsolver_type *Tr;
  gsl_root_fsolver *sr;
  double m = S->B.sig2Drive;
  double r1=0.0, r2=0.0;
  double a = S->B.gsla;
  double b = S->B.gslb;
  gsl_function F;
     
     
  /* Find minimum of -LogL */  
  p.offset=0.0; 
  p.CalcCode=2;
  F.function = f;
  F.params = (void *) &p;
     
  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m,a,b);

  if(S->B.Debug>8){
    printf ("Minimizing using %s method\n", gsl_min_fminimizer_name (s));
     
    //printf ("%5s [%9s, %9s] %9s %10s %9s\n",
    //        "iter", "lower", "upper", "min", "err(est)");
     
    printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
           iter, a, b, m, b - a);
  }
  do {
    iter++;
    status = gsl_min_fminimizer_iterate (s);
     
    m = gsl_min_fminimizer_x_minimum (s);
    a = gsl_min_fminimizer_x_lower (s);
    b = gsl_min_fminimizer_x_upper (s);
     
    status = gsl_min_test_interval (a, b, 0.0, S->B.gslRelErr);
    if(S->B.Debug>8){ 
      if (status == GSL_SUCCESS)
        printf ("Converged:\n");
     
      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
          iter, a, b, m, b - a);
    }
  } while (status == GSL_CONTINUE && iter < max_iter);
  
  if(status==GSL_CONTINUE){
    fprintf(stderr,"Failed to converge while searching for MLE in MC_1dMaximize.\n");
   
    m=0.0/0.0;
    r1=0.0/0.0;
    r2=0.0/0.0;
    MLE->v=m;
    MLE->LogL=0.0/0.0;
    
    CalcCI=0;
  }else{
    MLE->v=m;
    MLE->LogL=-(*f)(m,(void *) &p);
    iterMLE=iter;
  }
  
  gsl_min_fminimizer_free(s);
  
  if(CalcCI){
    /* Find left confidence interval */
    p.offset=0.0;
    p.offset=-(*f)(m,&p)-2.0; 
    a=S->B.gsla;
    b=m;
    Tr = gsl_root_fsolver_brent;
    sr = gsl_root_fsolver_alloc (Tr);
    gsl_root_fsolver_set (sr, &F,a,b);

    if(S->B.Debug>8){
      printf ("Finding CI using %s method\n", gsl_root_fsolver_name (sr));
       
      //printf ("%5s [%9s, %9s] %9s %10s %9s\n",
      //        "iter", "lower", "upper", "min", "err(est)");
       
      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
             iter, a, b, r1, b - a);
    }
    do {
      iter++;
      status = gsl_root_fsolver_iterate (sr);
       
      r1 = gsl_root_fsolver_root (sr);
      a = gsl_root_fsolver_x_lower (sr);
      b = gsl_root_fsolver_x_upper (sr);
       
      status = gsl_root_test_interval (a, b, 0.0,S->B.gslRelErr);
      if(S->B.Debug>8){ 
        if (status == GSL_SUCCESS)
          printf ("Converged:\n");
       
        printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
            iter, a, b, r1, b - a);
      }
    } while (status == GSL_CONTINUE && iter < max_iter);

    if(status==GSL_CONTINUE){
      fprintf(stderr,"Failed to converge while searching for left confidence interval boundary in MC_1dMaximize.\n");
      r1=0.0/0.0;
    }
    iterCIl=iter;
    
    /* Find right confidence interval */
    p.offset=0.0;
    p.offset=-(*f)(m,&p)-2.0; 
    a=m;
    b=S->B.gslb;
    Tr = gsl_root_fsolver_brent;
    sr = gsl_root_fsolver_alloc (Tr);
    gsl_root_fsolver_set (sr, &F,a,b);

    if(S->B.Debug>8){
      printf ("Finding CI using %s method\n", gsl_root_fsolver_name (sr));
       
      //printf ("%5s [%9s, %9s] %9s %10s %9s\n",
      //        "iter", "lower", "upper", "min", "err(est)");
       
      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
             iter, a, b, r2, b - a);
    }
    do {
      iter++;
      status = gsl_root_fsolver_iterate (sr);
       
      r2 = gsl_root_fsolver_root (sr);
      a = gsl_root_fsolver_x_lower (sr);
      b = gsl_root_fsolver_x_upper (sr);
       
      status = gsl_root_test_interval (a, b, 0.0,S->B.gslRelErr);
      if(S->B.Debug>8){ 
        if (status == GSL_SUCCESS)
          printf ("Converged:\n");
        printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
            iter, a, b, r2, b - a);
      }
    } while (status == GSL_CONTINUE && iter < max_iter);
    if(status==GSL_CONTINUE){
      fprintf(stderr,"Failed to converge while searching for right confidence interval boundary in MC_1dMaximize.\n");
      r2=0.0/0.0;
    }
    iterCIu=iter;
    gsl_root_fsolver_free(sr);
    p.offset=0.0;
  }

  MLE->lowerCI=r1;
  MLE->upperCI=r2;
  printf("MLEwithCI: [%g, %g, %g] L: %g NumIterations: %d %d %d\n",MLE->lowerCI,MLE->v,MLE->upperCI,MLE->LogL,iterMLE,iterCIl,iterCIu);

}

double MC_gslwrap_f_2dAcrossLoci(const gsl_vector * x, void * p){

  int i;
  double s[2];
  struct LsurfPoint f;
  f.sig2=0;f.mu=0;f.x0=0;
  struct MC_gslwrap_Params * params=(struct MC_gslwrap_Params *) p;

  for(i=0;i<2;i++){
    s[i]=gsl_vector_get(x,i);
     if(s[i]<0.1)
      s[i]=0.1;
  }
  MC_CalcPtsig2MarginAcrossLociSurf(&f,params->S,s);
  if(-f.LogL==-f.LogL) // test for nan;
    return(-f.LogL);
  else
    return(FLT_MAX);
};

void MC_gslwrap_df_2dAcrossLoci(const gsl_vector * x, void * p, gsl_vector *g){
  /* Calcs the function at x and approximates the gradient at x using a central difference equation */
  struct MC_gslwrap_Params * params=(struct MC_gslwrap_Params *) p;

  int i;
  struct LsurfPoint f0a, f0b, f1a, f1b;
  double s[2];
  double eps=params->S->B.gsldfeps;

  /* Set Lsurf pointers that we won't use to null */
  f0a.sig2=0; f0a.mu=0; f0a.x0=0;
  f0b.sig2=0; f0b.mu=0; f0b.x0=0;
  f1a.sig2=0; f1a.mu=0; f1a.x0=0;
  f1b.sig2=0; f1b.mu=0; f1b.x0=0;



  for(i=0;i<2;i++){
    s[i]=gsl_vector_get(x,i);
     if(s[i]<0.1)
      s[i]=0.1;
  }
  s[0]+=eps;
  MC_CalcPtsig2MarginAcrossLociSurf(&f0a,params->S,s);
  s[0]-=2*eps;
  MC_CalcPtsig2MarginAcrossLociSurf(&f0b,params->S,s);
  s[0]+=eps;
  
  s[1]+=eps;
  MC_CalcPtsig2MarginAcrossLociSurf(&f1a,params->S,s);
  s[1]-=2*eps;
  MC_CalcPtsig2MarginAcrossLociSurf(&f1b,params->S,s);
  s[1]+=eps;
  
  gsl_vector_set(g,0,(-f0a.LogL+f0b.LogL)/(2*eps));
  gsl_vector_set(g,1,(-f1a.LogL+f1b.LogL)/(2*eps));

};

void MC_gslwrap_fdf_2dAcrossLoci(const gsl_vector * x, void * p, double * f, gsl_vector *g){
  /* Calcs the function at x and approximates the gradient at x using a central difference equation */
  struct MC_gslwrap_Params * params=(struct MC_gslwrap_Params *) p;

  int i;
  struct LsurfPoint f0a, f0b, f1a, f1b;
  double s[2];
  double eps=params->S->B.gsldfeps;

  /* Set Lsurf pointers that we won't use to null */
  f0a.sig2=0; f0a.mu=0; f0a.x0=0;
  f0b.sig2=0; f0b.mu=0; f0b.x0=0;
  f1a.sig2=0; f1a.mu=0; f1a.x0=0;
  f1b.sig2=0; f1b.mu=0; f1b.x0=0;



  for(i=0;i<2;i++){
    s[i]=gsl_vector_get(x,i);
     if(s[i]<.1)
      s[i]=0.1;
  }
  s[0]+=eps;
  MC_CalcPtsig2MarginAcrossLociSurf(&f0a,params->S,s);
  s[0]-=2*eps;
  MC_CalcPtsig2MarginAcrossLociSurf(&f0b,params->S,s);
  s[0]+=eps;
  
  s[1]+=eps;
  MC_CalcPtsig2MarginAcrossLociSurf(&f1a,params->S,s);
  s[1]-=2*eps;
  MC_CalcPtsig2MarginAcrossLociSurf(&f1b,params->S,s);
  s[1]+=eps;
  
  gsl_vector_set(g,0,(-f0a.LogL+f0b.LogL)/(2*eps));
  gsl_vector_set(g,1,(-f1a.LogL+f1b.LogL)/(2*eps));
  *f=-(f0a.LogL+f0b.LogL)/2;
}



double MC_gslwrap_f_2dPerLocus(const gsl_vector * x, void * p){

  int i;
  double s[2];
  struct LsurfPoint f;
  f.sig2=0;f.mu=0;f.x0=0;
  struct MC_gslwrap_Params * params=(struct MC_gslwrap_Params *) p;

  for(i=0;i<2;i++){
    s[i]=gsl_vector_get(x,i);
     if(s[i]<0.1)
      s[i]=0.1;
  }
  MC_CalcPtsig2MarginSurf(&f,params->S,params->LocusNo,s);
  return(-f.LogL);
};

void MC_gslwrap_df_2dPerLocus(const gsl_vector * x, void * p, gsl_vector *g){
  /* Calcs the function at x and approximates the gradient at x using a central difference equation */
  struct MC_gslwrap_Params * params=(struct MC_gslwrap_Params *) p;

  int i;
  struct LsurfPoint f0a, f0b, f1a, f1b;
  double s[2];
  double eps=params->S->B.gsldfeps;

  /* Set Lsurf pointers that we won't use to null */
  f0a.sig2=0; f0a.mu=0; f0a.x0=0;
  f0b.sig2=0; f0b.mu=0; f0b.x0=0;
  f1a.sig2=0; f1a.mu=0; f1a.x0=0;
  f1b.sig2=0; f1b.mu=0; f1b.x0=0;



  for(i=0;i<2;i++){
    s[i]=gsl_vector_get(x,i);
     if(s[i]<0.1)
      s[i]=0.1;
  }
  s[0]+=eps;
  MC_CalcPtsig2MarginSurf(&f0a,params->S,params->LocusNo,s);
  s[0]-=2*eps;
  MC_CalcPtsig2MarginSurf(&f0b,params->S,params->LocusNo,s);
  s[0]+=eps;
  
  s[1]+=eps;
  MC_CalcPtsig2MarginSurf(&f1a,params->S,params->LocusNo,s);
  s[1]-=2*eps;
  MC_CalcPtsig2MarginSurf(&f1b,params->S,params->LocusNo,s);
  s[1]+=eps;
  
  gsl_vector_set(g,0,(-f0a.LogL+f0b.LogL)/(2*eps));
  gsl_vector_set(g,1,(-f1a.LogL+f1b.LogL)/(2*eps));

};

void MC_gslwrap_fdf_2dPerLocus(const gsl_vector * x, void * p, double * f, gsl_vector *g){
  /* Calcs the function at x and approximates the gradient at x using a central difference equation */
  struct MC_gslwrap_Params * params=(struct MC_gslwrap_Params *) p;

  int i;
  struct LsurfPoint f0a, f0b, f1a, f1b;
  double s[2];
  double eps=params->S->B.gsldfeps;

  /* Set Lsurf pointers that we won't use to null */
  f0a.sig2=0; f0a.mu=0; f0a.x0=0;
  f0b.sig2=0; f0b.mu=0; f0b.x0=0;
  f1a.sig2=0; f1a.mu=0; f1a.x0=0;
  f1b.sig2=0; f1b.mu=0; f1b.x0=0;



  for(i=0;i<2;i++){
    s[i]=gsl_vector_get(x,i);
     if(s[i]<.1)
      s[i]=0.1;
  }
  s[0]+=eps;
  MC_CalcPtsig2MarginSurf(&f0a,params->S,params->LocusNo,s);
  s[0]-=2*eps;
  MC_CalcPtsig2MarginSurf(&f0b,params->S,params->LocusNo,s);
  s[0]+=eps;
  
  s[1]+=eps;
  MC_CalcPtsig2MarginSurf(&f1a,params->S,params->LocusNo,s);
  s[1]-=2*eps;
  MC_CalcPtsig2MarginSurf(&f1b,params->S,params->LocusNo,s);
  s[1]+=eps;
  
  gsl_vector_set(g,0,(-f0a.LogL+f0b.LogL)/(2*eps));
  gsl_vector_set(g,1,(-f1a.LogL+f1b.LogL)/(2*eps));
  *f=-(f0a.LogL+f0b.LogL)/2;
}


void MC_2dMaximize(struct SamplerInfo * S, double(*f)(const gsl_vector *,void *),void (*df)(const gsl_vector *x,void *,gsl_vector *),
 void (*fdf)(const gsl_vector *, void *,double *, gsl_vector *), double(*fCI)(double, void *),
 struct MC_gslwrap_Params p, double sig2x,double sig2y,int CalcCI,struct MLE * MLE){

   size_t iter = 0;
   int iterMLE=0,iterCIl0=0,iterCIu0=0,iterCIl1=0,iterCIu1=0;
   int status;
   double size;
   int error=0;
   const gsl_multimin_fminimizer_type *T;
   gsl_multimin_fminimizer *s;

   const gsl_root_fsolver_type *Tr;
   gsl_root_fsolver *sr;
   double a,b,r1,r2;
   
   gsl_vector *ss, *x;
   gsl_multimin_function my_func;
   gsl_function F;
   
   /* Initialize functions */
   my_func.f = f;
   my_func.n = 2;
   my_func.params = &p;
   F.function = fCI;
   F.params = (void *) &p;
   r1=0.0;
   r2=0.0;
   /* Set initial step size for NM simplex algorithm */
   ss=gsl_vector_alloc(2);
   gsl_vector_set_all(ss,sig2x*1e-2);

   /* Set Starting point, x = (sig2x,sig2y) */   
   x = gsl_vector_alloc (2);
   gsl_vector_set (x, 0, sig2x);
   gsl_vector_set (x, 1, sig2y);
 
   T = gsl_multimin_fminimizer_nmsimplex;
   s = gsl_multimin_fminimizer_alloc (T, 2);
 
   gsl_multimin_fminimizer_set (s, &my_func, x, ss);
   
   if(p.S->B.Debug>8)
    printf("Finding minimum using Nelder Mead simplex algorithm...\n");
     
   do
     {
       iter++;
       status = gsl_multimin_fminimizer_iterate (s);
  
       if (status){
         fprintf(stderr,"***ERROR*** Error in Multidimensional maximization: %d.\n\tConsider finding MLE by calculating a grid-based approximation \n\tto the likelihood surface using -s option on command-line.\n",status);
         error=1;
         break;
        }
       //gsl_vector_fprintf (stdout, s->gradient,"%g");
       size=gsl_multimin_fminimizer_size(s);
       status = gsl_multimin_test_size (size,sig2x*1e-3);
      if(p.S->B.Debug>8){
         if (status == GSL_SUCCESS)
           printf ("Minimum found at:\n");
   
         printf ("%5d %.5f %.5f %10.5f %.3f %d %d\n", (int)iter,
                 gsl_vector_get (s->x, 0),
                 gsl_vector_get (s->x, 1),
                 s->fval,size,status,iter<p.S->B.gslMaxIter);
      }
     }
   while (status == GSL_CONTINUE && iter < p.S->B.gslMaxIter);
  
   if(status==GSL_CONTINUE||error){
      if(status==GSL_CONTINUE)
        fprintf(stderr,"Failed to converge while searching for MLE in MC_2dMaximize.\n");
      gsl_vector_set(s->x,0,0.0/0.0);
      gsl_vector_set(s->x,1,0.0/0.0);
       
      MLE[0].v=0.0/0.0;
      MLE[0].LogL=0.0/0.0;       
      MLE[0].upperCI=0.0/0.0;
      MLE[0].lowerCI=0.0/0.0;
      MLE[1].v=0.0/0.0;
      MLE[1].LogL=0.0/0.0;       
      MLE[1].upperCI=0.0/0.0;
      MLE[1].lowerCI=0.0/0.0;

      CalcCI=0;
      error=0;
   }else{
  
    MLE[0].v= gsl_vector_get (s->x, 0);
    MLE[0].LogL=-s->fval;
    MLE[1].v= gsl_vector_get (s->x, 1);
    MLE[1].LogL=-s->fval;
    iterMLE=iter;
   }
   
   if(CalcCI){
    /* Find left confidence interval for dimension 0 */
    p.CalcCode=0;
    p.offset=-s->fval-2.0;
    p.sig2base=MLE[1].v; 
    a=S->B.gsla;
    b=MLE[0].v;
    //printf("[%g, %g] %g %d\n",MLE[0].v,MLE[1].v,MLE[0].LogL,iter);

    //printf("%g %g %g %g\n",a,b,(*fCI)(a,(void *) &p),(*fCI)(b,(void *) &p));
    Tr = gsl_root_fsolver_brent;
    sr = gsl_root_fsolver_alloc (Tr);
    gsl_root_fsolver_set (sr, &F,a,b);

    if(S->B.Debug>8){
      printf ("Finding CI using %s method\n", gsl_root_fsolver_name (sr));
       
      //printf ("%5s [%9s, %9s] %9s %10s %9s\n",
      //        "iter", "lower", "upper", "min", "err(est)");
       
      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
             (int)iter, a, b, r1, b - a);
    }
    do {
      iter++;
      status = gsl_root_fsolver_iterate (sr);
       
      r1 = gsl_root_fsolver_root (sr);
      a = gsl_root_fsolver_x_lower (sr);
      b = gsl_root_fsolver_x_upper (sr);
       
      status = gsl_root_test_interval (a, b, 0.0,S->B.gslRelErr);
      if(S->B.Debug>8){ 
        if (status == GSL_SUCCESS)
          printf ("Converged:\n");
       
        printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
            (int)iter, a, b, r1, b - a);
      }
    } while (status == GSL_CONTINUE && iter < S->B.gslMaxIter);

    if(status==GSL_CONTINUE){
      fprintf(stderr,"Failed to converge while searching for left confidence interval boundary in MC_2dMaximize.\n");
      r1=0.0/0.0;
    }
    
    iterCIl0=iter;
   
    
    /* Find right confidence interval for dimension 0 */
    p.offset=-s->fval-2.0;
    a=MLE[0].v;
    b=S->B.gslb;
    Tr = gsl_root_fsolver_brent;
    sr = gsl_root_fsolver_alloc (Tr);
    gsl_root_fsolver_set (sr, &F,a,b);

    if(S->B.Debug>8){
      printf ("Finding CI using %s method\n", gsl_root_fsolver_name (sr));
       
      //printf ("%5s [%9s, %9s] %9s %10s %9s\n",
      //        "iter", "lower", "upper", "min", "err(est)");
       
      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
             (int)iter, a, b, r2, b - a);
    }
    do {
      iter++;
      status = gsl_root_fsolver_iterate (sr);
       
      r2 = gsl_root_fsolver_root (sr);
      a = gsl_root_fsolver_x_lower (sr);
      b = gsl_root_fsolver_x_upper (sr);
       
      status = gsl_root_test_interval (a, b, 0.0,S->B.gslRelErr);
      if(S->B.Debug>8){ 
        if (status == GSL_SUCCESS)
          printf ("Converged:\n");
        printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
            (int)iter, a, b, r2, b - a);
      }
    } while (status == GSL_CONTINUE && iter < S->B.gslMaxIter);
    if(status==GSL_CONTINUE){
      fprintf(stderr,"Failed to converge while searching for right confidence interval boundary in MC_2dMaximize.\n");
      r2=0.0/0.0;
    }

    iterCIu0=iter;
    
    MLE[0].upperCI=r2;
    MLE[0].lowerCI=r1;
    

    /* Find left confidence interval for dimension 1 */
    p.CalcCode=1;
    p.sig2base=MLE[0].v; 
    p.offset=-s->fval-2.0; 
    a=S->B.gsla;
    b=MLE[1].v;
    Tr = gsl_root_fsolver_brent;
    sr = gsl_root_fsolver_alloc (Tr);
    gsl_root_fsolver_set (sr, &F,a,b);

    if(S->B.Debug>8){
      printf ("Finding CI using %s method\n", gsl_root_fsolver_name (sr));
       
      //printf ("%5s [%9s, %9s] %9s %10s %9s\n",
      //        "iter", "lower", "upper", "min", "err(est)");
       
      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
             (int)iter, a, b, r1, b - a);
    }
    do {
      iter++;
      status = gsl_root_fsolver_iterate (sr);
       
      r1 = gsl_root_fsolver_root (sr);
      a = gsl_root_fsolver_x_lower (sr);
      b = gsl_root_fsolver_x_upper (sr);
       
      status = gsl_root_test_interval (a, b, 0.0,S->B.gslRelErr);
      if(S->B.Debug>8){ 
        if (status == GSL_SUCCESS)
          printf ("Converged:\n");
       
        printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
            (int)iter, a, b, r1, b - a);
      }
    } while (status == GSL_CONTINUE && iter < S->B.gslMaxIter);

    if(status==GSL_CONTINUE){
      fprintf(stderr,"Failed to converge while searching for left confidence interval boundary in MC_2dMaximize.\n");
      r1=0.0/0.0;
    }
    iterCIu1=iter;

    
    /* Find right confidence interval for dimension 1 */
    p.offset=-s->fval-2.0; 
    a=MLE[1].v;
    b=S->B.gslb;
    Tr = gsl_root_fsolver_brent;
    sr = gsl_root_fsolver_alloc (Tr);
    gsl_root_fsolver_set (sr, &F,a,b);

    if(S->B.Debug>8){
      printf ("Finding CI using %s method\n", gsl_root_fsolver_name (sr));
       
      //printf ("%5s [%9s, %9s] %9s %10s %9s\n",
      //        "iter", "lower", "upper", "min", "err(est)");
       
      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
             (int)iter, a, b, r2, b - a);
    }
    do {
      iter++;
      status = gsl_root_fsolver_iterate (sr);
       
      r2 = gsl_root_fsolver_root (sr);
      a = gsl_root_fsolver_x_lower (sr);
      b = gsl_root_fsolver_x_upper (sr);
       
      status = gsl_root_test_interval (a, b, 0.0,S->B.gslRelErr);
      if(S->B.Debug>8){ 
        if (status == GSL_SUCCESS)
          printf ("Converged:\n");
        printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
            (int)iter, a, b, r2, b - a);
      }
    } while (status == GSL_CONTINUE && iter < S->B.gslMaxIter);
    if(status==GSL_CONTINUE){
      fprintf(stderr,"Failed to converge while searching for right confidence interval boundary in MC_2dMaximize.\n");
      r2=0.0/0.0;
    }

    iterCIl1=iter;

    
    MLE[1].upperCI=r2;
    MLE[1].lowerCI=r1;
    



    
    gsl_root_fsolver_free(sr);
    p.offset=0.0;
  }
  
  printf("MLEwithCIs: [%g, %g, %g] [%g, %g, %g] L: %g NumIterations: %d %d %d %d %d\n",MLE[0].lowerCI,MLE[0].v,MLE[0].upperCI,
    MLE[1].lowerCI,MLE[1].v,MLE[1].upperCI,MLE[0].LogL,iterMLE,iterCIl0,iterCIu0,iterCIl1,iterCIu1);
      
   gsl_multimin_fminimizer_free (s);
   gsl_vector_free (x);
   gsl_vector_free(ss);
      
}



void MC_PointEstimates(struct SamplerInfo * S){
  int i,d;
  struct MC_gslwrap_Params p;
  struct LsurfPoint x;
  double * dummy;
  
  double * s;
  
  s=(double *) calloc((size_t)S->B.Dim,sizeof(double));
  dummy=(double *) calloc((size_t)S->B.Dim,sizeof(double));

  p.S=S;
  x.sig2=dummy;
  x.mu=dummy;
  x.x0=dummy;
  
  //gsl_set_error_handler_off();
  if(S->B.CalcMLEperLocus){
    for(i=0;i<S->B.nLoci;i++){
      printf("Finding MLE for locus %s\n",S->LD[i].ID);
      p.LocusNo=i;
      /* Find MLEs */
      MC_1dMaximize(S,&MC_gslwrap_1dPerLocus,p,S->B.gslCIperLocus,&S->sig2Margin1dMLE[i]);
      if(S->B.Dim==2) 
        MC_2dMaximize(S,&MC_gslwrap_f_2dPerLocus,&MC_gslwrap_df_2dPerLocus,&MC_gslwrap_fdf_2dPerLocus,&MC_gslwrap_1dPerLocus,p,S->sig2Margin1dMLE[i].v,S->sig2Margin1dMLE[i].v,S->B.gslCIperLocus,S->sig2MarginMLE[i]);
      /* Calculate standard errors on MLEs */
      for(d=0;d<S->B.Dim;d++) s[d]=S->sig2Margin1dMLE[i].v;
      MC_CalcPtsig2MarginSurf(&x,S,i,s);
      S->sig2Margin1dMLE[i].LogLSE=log(x.LSE);
      for(d=0;d<S->B.Dim;d++) s[d]=S->sig2MarginMLE[i][d].v; 
      MC_CalcPtsig2MarginSurf(&x,S,i,s);
      for(d=0;d<S->B.Dim;d++) S->sig2MarginMLE[i][d].LogLSE=log(x.LSE);
    }
  }
  printf("Finding MLE for full dataset.\n");
  MC_1dMaximize(S,&MC_gslwrap_1dAcrossLoci,p,1,&(S->sig2MarginAllLoci1dMLE));
  if(S->B.Dim==2)
    MC_2dMaximize(S,&MC_gslwrap_f_2dAcrossLoci,&MC_gslwrap_df_2dAcrossLoci,&MC_gslwrap_fdf_2dAcrossLoci,&MC_gslwrap_1dAcrossLoci,p,
      S->sig2MarginAllLoci1dMLE.v,S->sig2MarginAllLoci1dMLE.v,1,S->sig2MarginAllLociMLE);
  /* Calculate standard errors on MLEs */
  for(d=0;d<S->B.Dim;d++) s[d]=S->sig2MarginAllLoci1dMLE.v;
  MC_CalcPtsig2MarginAcrossLociSurf(&x,S,s);
  S->sig2MarginAllLoci1dMLE.LogLSE=log(x.LSE);
  for(d=0;d<S->B.Dim;d++) s[d]=S->sig2MarginAllLociMLE[d].v; 
  MC_CalcPtsig2MarginAcrossLociSurf(&x,S,s);
  
  for(d=0;d<S->B.Dim;d++) S->sig2MarginAllLociMLE[d].LogLSE=log(x.LSE);
  free(s);
  
  
}
