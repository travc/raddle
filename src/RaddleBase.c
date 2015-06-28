/*
 *  RaddleBase.c
 *  Raddle
 *
 *  Created by John Novembre on 12/6/05.
 *
 */

#include <unistd.h> /* For Malloc Debug */

#include "RaddleBase.h"
#include "mainMC.h"

#include<math.h>


struct LocusData * ReadData(struct BaseInfo B){
  
  struct LocusData *LocusArray;
  int i,j,k;
  FILE *in;
  char filename[80];
  // Declare Memory for Arrays
  LocusArray=(struct LocusData *) calloc((size_t)B.nLoci,sizeof(struct LocusData));

  strcpy(filename,B.InFileName);
  strcat(filename,".xy");



  if(!(in=fopen(filename,"r"))){
    fprintf(stderr,"Error opening input file %s.\n",filename);
    exit(-1);
  }


  fscanf(in,"%*s%*s%*s");
  for(i=0;i<B.Dim;i++)
    fscanf(in,"%*s");


  for(i=0;i<B.nLoci;i++){

    // Read LocusArray[i].n
    fscanf(in,"%s%lf%d",LocusArray[i].ID,&LocusArray[i].f,&LocusArray[i].n);
    
    
    // Declare Mem for Positions
    LocusArray[i].x=(double **) calloc((size_t)B.Dim,sizeof(double *));
    for(j=0;j<B.Dim;j++){
      LocusArray[i].x[j]=(double *) calloc((size_t)LocusArray[i].n,sizeof(double));
      
    }

    /* Read first row data points */
    if(B.Dim==2)
      fscanf(in,"%lf%lf\n",&LocusArray[i].x[0][0],&LocusArray[i].x[1][0]);
    else
      for(k=0;k<B.Dim;k++)
        fscanf(in,"%lf ",&LocusArray[i].x[k][0]);


    // Read in Remaining Postions
    for(j=1;j<LocusArray[i].n;j++){
      if(B.Dim==2)
        fscanf(in,"%*s%*f%*d%lf%lf\n",&LocusArray[i].x[0][j],&LocusArray[i].x[1][j]);
      else{
        fscanf(in,"%*s%*f%*d");
        for(k=0;k<B.Dim;k++)
          fscanf(in,"%lf ",&LocusArray[i].x[k][j]);
    
      }
    }
  }

  fclose(in);
  return(LocusArray);

}



struct LocusData * ReadDataOldFormat(struct BaseInfo B){
  
  struct LocusData *LocusArray;
  int i,j,k;
  FILE *in;
  char filename[80];
  // Declare Memory for Arrays
  LocusArray=(struct LocusData *) calloc((size_t)B.nLoci,sizeof(struct LocusData));

  strcpy(filename,B.InFileName);
  strcat(filename,".xy");



  if(!(in=fopen(filename,"r"))){
    fprintf(stderr,"Error opening input file %s.\n",filename);
    exit(-1);
  }

  for(i=0;i<B.nLoci;i++){

    // Read LocusArray[i].n
    fscanf(in,"%s%lf%d\n",LocusArray[i].ID,&LocusArray[i].f,&LocusArray[i].n);
    
    
    // Declare Mem for Positions
    LocusArray[i].x=(double **) calloc((size_t)B.Dim,sizeof(double *));
    for(j=0;j<B.Dim;j++){
      LocusArray[i].x[j]=(double *) calloc((size_t)LocusArray[i].n,sizeof(double));
      
    }

    // Read in Postions
    for(j=0;j<LocusArray[i].n;j++){
      if(B.Dim==2)
        fscanf(in,"%lf,%lf\n",&LocusArray[i].x[0][j],&LocusArray[i].x[1][j]);
      else
        for(k=0;k<B.Dim;k++)
          fscanf(in,"%lf ",&LocusArray[i].x[k][j]);
    }
  }

  fclose(in);
  return(LocusArray);

}






// Code from Monty to generate times from reconstructed 
// bd phylogeny



void bdtimes(double *times, double f, double r, double s, int samsize) {
	double t1choice(double f, double r, double s, int j);
	double trandom(double u, double f, double xi, double t1);
	//double unif(void);
	int j, compare(const void *e1, const void *e2);

	times[0] = times[1] = BIGTIME;
	times[samsize+1] = 0.0;  // so qsort can be used, they should not change
	times[1] = t1choice(f, r, s, samsize);
	for (j=2; j<=samsize; j++)
	  times[j] = trandom((double)ranf(), f, r+s, times[1]);
	qsort((void *) times, (size_t) samsize+1, sizeof(double), compare);

}  // end, bdtimes

double t1choice(double f, double r, double s, int j) {
	  double x, t=0.0, xi;
	  //double unif();
	  if (fabs(r) < 1e-7) { // r = 0
	    x = exp(log((double)ranf()) / j);
	    if (fabs(s) < 1e-7) // s = 0
	      return 2 * x / (f * (1 - x));
	    else if (s > 0.0)   // s > 0
	      return -log(f*(1-x)/(f-(f-2*s)*x))/s;
	    else if (s < 0.0)   // s < 0
	      return -log((f*(1-x)-2*s)/((f-2*s)*(1-x)))/s;
	  } // end, r = 0
	  else if (r > 0.0) {
	    xi = r + s;
	    if (fabs(xi) < 1e-7)  // xi = 0
	      do {
		x = exp(log((double)ranf()) / j);
		t = 2 * x / (f * (1 - x));
	      } while ((double)ranf() > exp(-r * t));
	    else if (xi > 0.0)      // xi = r + s > 0
	      do {
		x = exp(log((double)ranf()) / j);
		t = -log(f*(1-x)/(f-(f-2*xi)*x))/xi;
	      } while ((double)ranf() > exp(-r * t));
	    else if (xi < 0.0)  //  xi < 0
	      do {
		x = exp(log((double)ranf()) / j);
		t = -log((f*(1-x)-2*xi)/((f-2*xi)*(1-x)))/xi;
	      } while((double)ranf() > exp(-r * t));
	    return t;
	  } // end, r > 0
    return -1;
	} // end, t1choice

double trandom(double u, double f, double xi, double t1) {
	double capC;
	
	if (xi != 0.0) {
		capC = (f - (f-2*xi) * exp(-xi*t1))/(1-exp(-xi*t1));
		return log(1 + 2*u*xi/(capC-f*u))/xi;
		}
	else {
		capC = (2 + f * t1) / t1;
		return 2 * u / (capC - f * u);
		}
	}


/* needed for qsort, modified to sort in reverse order */
int compare(const void *e1, const void *e2)  {
	const double *x1, *x2;
	x1 = e1;
	x2 = e2;
	return (*x1 > *x2) ? -1 : (*x1 < *x2) ? 1: 0;
	}





struct ModelParams AllocModelParams(struct BaseInfo B){
  
  struct ModelParams MP;
  MP.sig2=(doubval *)calloc((size_t)B.Dim,sizeof(doubval));
  MP.mu=(doubval *)calloc((size_t)B.Dim,sizeof(doubval));
  MP.LogPriorPrsig2=(doubval *)calloc((size_t)B.Dim,sizeof(doubval));
  MP.LogPriorPrmu=(doubval *)calloc((size_t)B.Dim,sizeof(doubval));
  
  return(MP);
}

struct ModelParamsPtr AllocModelParamsPtr(struct BaseInfo B){
  struct ModelParamsPtr MP;
  MP.sig2Ptr=(doubval **)calloc((size_t)B.Dim,sizeof(doubval *));
  MP.muPtr=(doubval **)calloc((size_t)B.Dim,sizeof(doubval *));
  MP.LogPriorPrsig2Ptr=(doubval **)calloc((size_t)B.Dim,sizeof(doubval *));
  MP.LogPriorPrmuPtr=(doubval **)calloc((size_t)B.Dim,sizeof(doubval *));

  return(MP);
  
};


struct LocusParams * AllocLocusParams(struct BaseInfo B,struct LocusData * LD){

  int i,j,k;
  struct LocusParams * LP;

  LP=(struct LocusParams *)calloc((size_t)B.nLoci,sizeof(struct LocusParams));

  // Declare Mem for ancestral location
  for(i=0;i<B.nLoci;i++){
    LP[i].x0=(doubval *) calloc((size_t)B.Dim,sizeof(doubval));
    LP[i].LogPriorPrx0=(doubval *)calloc((size_t)B.Dim,sizeof(doubval));
    
    // Declare Mem For NodeList
    LP[i].NodeTimes=(doubval *)calloc((size_t)(LD[i].n-1),sizeof(doubval));
    LP[i].NodeList=(struct bNode *)calloc((size_t)(2*LD[i].n),sizeof(struct bNode));  
    for(j=0;j<2*LD[i].n;j++){
      LP[i].NodeList[j].t=0.0;
      LP[i].NodeList[j].v=0.0;
      LP[i].NodeList[j].Data=(struct bNodeData *)calloc((size_t)1,sizeof(struct bNodeData));
      LP[i].NodeList[j].Data->CalcFlag=1;
      LP[i].NodeList[j].Data->LogL=0.0;
      LP[i].NodeList[j].Data->xbar=(double *)calloc((size_t)B.Dim,sizeof(double));
      LP[i].NodeList[j].Data->S2=0.0;
      LP[i].NodeList[j].Data->vPrime=0.0;
      LP[i].NodeList[j].Data->c=(double *)calloc((size_t)LD[i].n,sizeof(double));
    }
    
    for(j=1;j<=LD[i].n;j++){
      for(k=0;k<B.Dim;k++){
        LP[i].NodeList[j].Data->xbar[k]=LD[i].x[k][j-1];
      }    
    }
    
    
    LP[i].LevelsList=(struct bNode ***)calloc((size_t)(LD[i].n),sizeof(struct bNode **));
    for(j=0;j<LD[i].n;j++){
      LP[i].LevelsList[j]=(struct bNode **)calloc((size_t)(j+1),sizeof(struct bNode *));
    }
    
  }
  return(LP);
}


struct LocusParamsPtr * AllocLocusParamsPtr(struct BaseInfo B,struct LocusData * LD){


  int i,j;
  struct LocusParamsPtr * LP;

  LP=(struct LocusParamsPtr *)calloc((size_t)B.nLoci,sizeof(struct LocusParamsPtr));

  // Declare Mem for ancestral location
  for(i=0;i<B.nLoci;i++){
    LP[i].x0Ptr=(doubval **) calloc((size_t)B.Dim,sizeof(doubval *));
    LP[i].LogPriorPrx0Ptr=(doubval **)calloc((size_t)B.Dim,sizeof(doubval *));
    
    // Declare Mem For NodeList
    LP[i].NodeTimesPtr=(doubval **)calloc((size_t)(LD[i].n-1),sizeof(doubval *));
    
    LP[i].NodeListPtr=(struct bNode **)calloc((size_t)(2*LD[i].n),sizeof(struct bNode *));  
    
    LP[i].LevelsListPtr=(struct bNode ****)calloc((size_t)(LD[i].n),sizeof(struct bNode ***));
    for(j=0;j<LD[i].n;j++){
      LP[i].LevelsListPtr[j]=(struct bNode ***)calloc((size_t)(j+1),sizeof(struct bNode **));
    }
    
  }
  return(LP);


};



void CopyLocusParams(struct BaseInfo B, struct LocusData LD,struct LocusParams LO,struct LocusParams *LC){

  int i,j;
  for(i=0;i<B.Dim;i++){
    LC->x0[i]=LO.x0[i];
    LC->LogPriorPrx0[i]=LO.LogPriorPrx0[i];
  }
  LC->t1=LO.t1;
  LC->f=LO.f;

  LC->LogPrt1Givenjf=LO.LogPrt1Givenjf;
  LC->LogPrNodeTimesGivenjf=LO.LogPrNodeTimesGivenjf;

  for(i=0;i<LD.n-1;i++){
    LC->NodeTimes[i]=LO.NodeTimes[i];
  }
  for(i=0;i<2*LD.n;i++){
    LC->NodeList[i].key=LO.NodeList[i].key;
    LC->NodeList[i].Tip=LO.NodeList[i].Tip;
    LC->NodeList[i].t=LO.NodeList[i].t;
    LC->NodeList[i].v=LO.NodeList[i].v;
    LC->NodeList[i].Data->LogL=LO.NodeList[i].Data->LogL;
    LC->NodeList[i].Data->CalcFlag=LO.NodeList[i].Data->CalcFlag;
    LC->NodeList[i].Data->vPrime=LO.NodeList[i].Data->vPrime;    
    LC->NodeList[i].Data->S2=LO.NodeList[i].Data->S2;
    for(j=0;j<B.Dim;j++){
      LC->NodeList[i].Data->xbar[j]=LO.NodeList[i].Data->xbar[j];
    }
    for(j=0;j<LD.n;j++){
      LC->NodeList[i].Data->c[j]=LO.NodeList[i].Data->c[j];
    }
  }


  LC->NodeList[0].lChild=&LC->NodeList[2*LD.n-1];
  
  for(i=LD.n+1;i<2*LD.n;i++){
    for(j=0;j<2*LD.n;j++){
      if(LC->NodeList[j].key==LO.NodeList[i].lChild->key)
	LC->NodeList[i].lChild=&LC->NodeList[j];
      if(LC->NodeList[j].key==LO.NodeList[i].rChild->key)
	LC->NodeList[i].rChild=&LC->NodeList[j];
      if(LC->NodeList[j].key==LO.NodeList[i].Anc->key)
	LC->NodeList[i].Anc=&LC->NodeList[j];
    }
  }


  for(i=0;i<LD.n;i++){
    for(j=0;j<i+1;j++)
      LC->LevelsList[i][j]=LO.LevelsList[i][j];
  }

}


/* The goal of this function is to sample a topology using the importance sampling algorithm */
void SampleTopoIS(int Dim, struct LocusData LD, struct LocusParamsPtr LP, struct LocusParamsPtr LPprime,struct ModelParamsPtr MP, double TopoISHeat){
  
  int i,j;
  
  /* w: The weight for the whole replicate */
  double Logwr;
  
  /* wk: The partial weight for the node */
  double Logwrk;
  
  /* LogLk: The term the node will contribute to the log likelihood */
  double LogLk;
  
  /* L: The log likelihood function for the replicate */
  double LogL;
  
  /* Active List of nodes */
  struct bNode ** ActiveList;
  
  /* Number on Active List of nodes */
  int nActive;
  
  int Node1;
  int Node2;
  int Ctr;

  double vPrime1, vPrime2;
  nActive=LD.n;
  Ctr=LD.n;
    
  /* Initialize wr */
  Logwr=0.0;
  LogL=0.0;
  /* Initialize ActiveList */
  ActiveList=(struct bNode **)calloc((size_t)(LD.n),sizeof(struct Node *) );

  for(j=1;j<=LD.n;j++){
    LPprime.NodeListPtr[j]->t=0.0;
    LPprime.NodeListPtr[j]->v=0.0;
    LPprime.NodeListPtr[j]->Data->CalcFlag=1;
    LPprime.NodeListPtr[j]->Data->LogL=0.0;
    LPprime.NodeListPtr[j]->Data->vPrime=0.0;
  }


  for(i=1;i<=LD.n;i++){
    ActiveList[i-1]=LPprime.NodeListPtr[i];
  }

  
  /* Build tree starting from current data and working back in time*/
  while(nActive>1){
    
    //printf("nActive: %d CoalTime: %g\n",nActive,LP.NodeTimesPtr[nActive-2]->v);
    SampleNodeISwLogs(Dim,nActive,LP.NodeTimesPtr[nActive-2]->v,ActiveList,MP,
            &LogLk,&Logwrk,&Node1,&Node2,TopoISHeat);
    Logwr+=Logwrk;
    LogL+=LogLk;

    Ctr++;
    //printf("%d %d %d %g\n",Ctr,Node1,Node2,LogL);
    
    
    /* Connect Pointers to new node */
    LPprime.NodeListPtr[Ctr]->lChild=ActiveList[Node1];
    ActiveList[Node1]->Anc=LP.NodeListPtr[Ctr];
    LPprime.NodeListPtr[Ctr]->rChild=ActiveList[Node2];
    ActiveList[Node2]->Anc=LP.NodeListPtr[Ctr];
    LPprime.NodeListPtr[Ctr]->key=Ctr;
    
    /* Initialize raw branch length information */
    LPprime.NodeListPtr[Ctr]->t=LP.NodeTimesPtr[nActive-2]->v;
    LPprime.NodeListPtr[Ctr]->lChild->v=LPprime.NodeListPtr[Ctr]->t-LPprime.NodeListPtr[Ctr]->lChild->t;
    LPprime.NodeListPtr[Ctr]->rChild->v=LPprime.NodeListPtr[Ctr]->t-LPprime.NodeListPtr[Ctr]->rChild->t;

    /* Initialize vPrime branch length information */
    vPrime1=LPprime.NodeListPtr[Ctr]->lChild->v+LPprime.NodeListPtr[Ctr]->lChild->Data->S2;
    vPrime2=LPprime.NodeListPtr[Ctr]->rChild->v+LPprime.NodeListPtr[Ctr]->rChild->Data->S2;
    LPprime.NodeListPtr[Ctr]->lChild->Data->vPrime=vPrime1;
    LPprime.NodeListPtr[Ctr]->rChild->Data->vPrime=vPrime2;

    /* Calculate S2 */
    LPprime.NodeListPtr[Ctr]->Data->S2=vPrime1*vPrime2/(vPrime1+vPrime2);

    /* Calculate xbar */
    for(i=0;i<Dim;i++){
      LPprime.NodeListPtr[Ctr]->Data->xbar[i]=(vPrime2*ActiveList[Node1]->Data->xbar[i]+vPrime1*ActiveList[Node2]->Data->xbar[i])/(vPrime1+vPrime2);
    }
    
    /* Assign L */
    LPprime.NodeListPtr[Ctr]->Data->LogL=LogL; 
    LPprime.NodeListPtr[Ctr]->Data->LogLk=LogLk; 
    LPprime.NodeListPtr[Ctr]->Data->Logwr=Logwr; 
    LPprime.NodeListPtr[Ctr]->Data->Logwrk=Logwrk; 
      
    //printf("SampleTopoIS: Joining %d %d to make node %d\n",LPprime.NodeListPtr[Ctr]->lChild->key,LPprime.NodeListPtr[Ctr]->rChild->key,LPprime.NodeListPtr[Ctr]->key);
    //printf("Logwrk: %g Logwr: %g\n",LPprime.NodeListPtr[Ctr]->Data->Logwrk,LPprime.NodeListPtr[Ctr]->Data->Logwr);
  

    // Add a pointer to it to the ActiveList; and overwrite Node1
    ActiveList[MIN(Node1,Node2)]=LPprime.NodeListPtr[Ctr];

  
    // Shift list down to overwrite Node2
    for(i=MAX(Node1,Node2);i<nActive-1;i++){
      ActiveList[i]=ActiveList[i+1];
    }
    nActive--;
    
      
  }
  
  /* Handle the root lineage */
  LogLk=0.0;
  for(i=0;i<Dim;i++){
    LogLk+=LogfN(LPprime.NodeListPtr[Ctr]->Data->xbar[i],LP.x0Ptr[i]->v+MP.muPtr[i]->v*LP.t1Ptr->v,(LP.t1Ptr->v-LPprime.NodeListPtr[Ctr]->t+LPprime.NodeListPtr[Ctr]->Data->S2)*MP.sig2Ptr[i]->v);
    
    LPprime.NodeListPtr[0]->Data->xbar[i]=LPprime.NodeListPtr[Ctr]->Data->xbar[i];    
  }
  
  LPprime.NodeListPtr[0]->lChild=ActiveList[0];
  ActiveList[0]->Anc=LPprime.NodeListPtr[0];

  LPprime.NodeListPtr[0]->t=LP.t1Ptr->v;
  
  LPprime.NodeListPtr[0]->lChild->v=LPprime.NodeListPtr[0]->t-LPprime.NodeListPtr[0]->lChild->t;   
  /* Initialize vPrime branch length information */
  LPprime.NodeListPtr[0]->lChild->Data->vPrime=(LP.t1Ptr->v-LPprime.NodeListPtr[Ctr]->t+LPprime.NodeListPtr[Ctr]->Data->S2);
   
  
  LogL+=LogLk;
  
  LPprime.NodeListPtr[0]->Data->LogL=LogL;
  LPprime.NodeListPtr[0]->Data->LogLk=LogLk;
  LPprime.NodeListPtr[0]->Data->Logwr=LPprime.NodeListPtr[Ctr]->Data->Logwr;
  //printf("%g\n",LPprime.NodeListPtr[0]->Data->LogL);
  //printf("SampleTopoIS Logwr: %g\n",LPprime.NodeListPtr[0]->Data->Logwr);
  
  free(ActiveList);
}




void SampleNodeISwLogs(int Dim, int k,double CoalTime, struct bNode** ActiveList, struct ModelParamsPtr MP,
      double * LogLk,double *wrk,int * Node1, int * Node2,double heat){

  double * bvec;
  double * probvec;
  double * probID1vec;
  double * probID2vec;
  double * LogLkvec;
  double psum=0.0;
  double z;
  int i,j;
  int l1, l2;

  bvec=(double *)calloc((size_t)k*(k-1)/2+1,sizeof(double));
  probvec=(double *)calloc((size_t)k*(k-1)/2+1,sizeof(double));
  probID1vec=(double *)calloc((size_t)k*(k-1)/2+1,sizeof(double));
  probID2vec=(double *)calloc((size_t)k*(k-1)/2+1,sizeof(double));
  LogLkvec=(double *)calloc((size_t)k*(k-1)/2+1,sizeof(double));

  /* Init cumulative distribution */
  probvec[0]=-1e300;
 
  for (i=1, l2=1; l2<k; l2++) { //  lp = l' in notes
    for (l1=0; l1<l2; l1++, i++)  {
  
      bvec[i]=0.0;
      for(j=0;j<Dim;j++){
        /* Calculate probability of joining nodes l and lp */ 
        bvec[i] += LogfN(ActiveList[l1]->Data->xbar[j]-ActiveList[l2]->Data->xbar[j],0.0,((CoalTime-ActiveList[l1]->t)+ActiveList[l1]->Data->S2+(CoalTime-ActiveList[l2]->t)+ActiveList[l2]->Data->S2)*MP.sig2Ptr[j]->v);
        //printf("x1-x2: %g sd:%g Sig2Val: %g\n",ActiveList[l1]->Data->xbar[j]-ActiveList[l2]->Data->xbar[j],((CoalTime-ActiveList[l1]->t)+ActiveList[l1]->Data->S2+(CoalTime-ActiveList[l2]->t)+ActiveList[l2]->Data->S2)*MP.sig2Ptr[j]->v,MP.sig2Ptr[j]->v);
      
      }
      bvec[i]/=heat; /* To flatten distribution */

      /* Calculate Log of sum of exponentials using trick */
      //printf("i:%d bvec:%g ",i,bvec[i]);

      if(i==1){
        psum=bvec[i];
        //printf(" %g\n",psum);
     
      }else{
        if(bvec[i]-psum>300){
          psum=bvec[i];
        }else{
          psum=psum+log(1.0+exp(bvec[i]-psum));
        }

        //printf(" %g\n",psum);
     
      }
      
      probvec[i] = psum;
      probID1vec[i]=l1;
      probID2vec[i]=l2;
      LogLkvec[i]=bvec[i]*heat;
    }
  }
	
  /*  printf("\nLogLikelihood of possible nodes: ");
    for(i=1;i<=k*(k-1)/2;i++){
      printf("%g ",LogLkvec[i]);
    }     
    printf("\n");
    printf("Relative log probability of choosing possible nodes: ");
    for(i=1;i<=k*(k-1)/2;i++){
      printf("%g ",probvec[i]);
    }     
    printf("\n");*/
  /*  Choose the coalescent event */
  if (k>2){	  
    
    z = log(ranf())+probvec[k*(k-1)/2];
   
     j = 0;
    while(z > probvec[j])
      j++;
    
    if(j==0){
    
      fprintf(stderr,"Underflow error in SampleNodeIS().\n");
      exit(-1);
    }
    
    l1=probID1vec[j];
    l2=probID2vec[j];
    
    *wrk = log(2.0)-log((double)k*(k-1));
    *wrk -= bvec[j]- probvec[k*(k-1)/2];
    *LogLk = LogLkvec[j];
  } else {  // k == 2
    l1 = 0;
    l2 = 1;
    *wrk =0.0;
    *LogLk=LogLkvec[1];
  }
		  
  *Node1=l1;
  *Node2=l2;
  
  free(bvec);
  free(probvec);
  free(probID1vec);
  free(probID2vec);
  free(LogLkvec);
  return;
}


void CalcISweight(struct BaseInfo B, struct LocusData LD, struct LocusParamsPtr LP,struct ModelParamsPtr MP, double TopoISHeat){

  
  int i,j,l1,l2;
  
  /* w: The weight for the whole replicate */
  double Logwr;
  
 
  /* L: The log likelihood function for the replicate */
  double LogL;
  
  /* Active List of nodes */
  struct bNode ** ActiveList;
  /* Number on Active List of nodes */
  int nActive;  
  int Ctr;
  int Node1=0, Node2=0;
  double CoalTime;
  double * bvec;
  double psum=0.0, Bp; 

  nActive=LD.n;
  Ctr=LD.n+1;


    
  /* Initialize wr */
  Logwr=0.0;
  LogL=0.0;


  /* Initialize ActiveList */
  ActiveList=(struct bNode **)calloc((size_t)(LD.n),sizeof(struct Node *) );
  bvec=(double *)calloc((size_t)LD.n*(LD.n-1)/2.0+1,sizeof(double));

  for(i=1;i<=LD.n;i++){
    ActiveList[i-1]=LP.NodeListPtr[i];
  }

  while(nActive>1){
  
    /*for(i=0;i<nActive;i++){
      printf("%d ",ActiveList[i]->key);
    }
    printf("\n");
    printf("key: %d rChild: %d lChild: %d ",LP.NodeListPtr[Ctr]->key,LP.NodeListPtr[Ctr]->rChild->key,LP.NodeListPtr[Ctr]->lChild->key);*/
    
    CoalTime=LP.NodeListPtr[Ctr]->t;
    /* Compute Logwrk and LogLk */

    
  
    for (i=1, l2=1; l2<nActive; l2++) { //  lp = l' in notes
      for (l1=0; l1<l2; l1++, i++)  {
  
        bvec[i]=0.0;
        for(j=0;j<B.Dim;j++){
          /* Calculate probability of joining nodes l and lp */ 
          bvec[i] += LogfN(ActiveList[l1]->Data->xbar[j]-ActiveList[l2]->Data->xbar[j],0.0,((CoalTime-ActiveList[l1]->t)+ActiveList[l1]->Data->S2+(CoalTime-ActiveList[l2]->t)+ActiveList[l2]->Data->S2)*MP.sig2Ptr[j]->v);
        }
        bvec[i]/=TopoISHeat; /* To flatten distribution */

        /* Calculate Log of sum of exponentials using trick */
        if(i==1){
          psum=bvec[i];
          Bp=bvec[i];
        }else{
          if(bvec[i]>Bp){ /* Calculate B */
            Bp=bvec[i];      
          }
          if(bvec[i]-psum>300){
            psum=bvec[i];
          }else{
            psum=psum+log(1.0+exp(bvec[i]-psum));
          }
        }
      
        //printf("%d %d %g %g %g %g\n",ActiveList[l1]->key,ActiveList[l2]->key,bvec[i],LogL,Logwr, psum);
          
        if(((ActiveList[l1]==LP.NodeListPtr[Ctr]->lChild)||(ActiveList[l1]==LP.NodeListPtr[Ctr]->rChild))
        &&((ActiveList[l2]==LP.NodeListPtr[Ctr]->lChild)||(ActiveList[l2]==LP.NodeListPtr[Ctr]->rChild))){
          //printf("**");
          Node1=l1;
          Node2=l2;
          Logwr+=log(2.0)-bvec[i]-log((double)nActive*(nActive-1));
          LogL+=bvec[i]*TopoISHeat;
        
        }
        
        //printf("\n");
      }
    }
    Logwr+=psum;

    ActiveList[MIN(Node1,Node2)]=LP.NodeListPtr[Ctr];

   
    // Shift list down to overwrite Node2
    for(i=MAX(Node1,Node2);i<nActive-1;i++){
      ActiveList[i]=ActiveList[i+1];
    }

   
    nActive--;
    Ctr++;
	}
  //printf("CalcISweight Logwr: %g\n",Logwr);
    
    
   /* Handle the root lineage */

  Ctr--;
  for(i=0;i<B.Dim;i++){
    LogL+=LogfN(ActiveList[0]->lChild->Data->xbar[i],LP.x0Ptr[i]->v+MP.muPtr[i]->v*LP.t1Ptr->v,(LP.t1Ptr->v-ActiveList[0]->t+ActiveList[0]->Data->S2)*MP.sig2Ptr[i]->v);
  }   
    
  LP.NodeListPtr[0]->Data->LogL=LogL;
  LP.NodeListPtr[0]->Data->Logwr=Logwr;
  //printf("%g\n",LP.NodeListPtr[0]->Data->LogL);
  
  
  free(ActiveList);
  free(bvec);

}



double CalcLikelihoodGiven_xbar_vprime(struct bNode *p,struct BaseInfo B, struct ModelParamsPtr M,struct LocusParamsPtr LP){

  double x0,x1,x2,v1,v2,mu,sig2x;
  double L,v;
  int i;
  
  if(p->lChild==0){
    L=1.0;
  }else if (p->Anc==0){
    L=CalcLikelihoodGiven_xbar_vprime(p->lChild,B,M,LP);
    v=p->lChild->Data->vPrime;
    for(i=0;i<B.Dim;i++){
      x1=p->lChild->Data->xbar[i];
      x0=LP.x0Ptr[i]->v;
      mu=M.muPtr[i]->v;
      sig2x=M.sig2Ptr[i]->v;   
      L*=fN(x1,x0+mu*v,v*sig2x);
   
    }
    p->Data->LogL=log(L);
  } else{
    L=CalcLikelihoodGiven_xbar_vprime(p->lChild,B,M,LP);
    L*=CalcLikelihoodGiven_xbar_vprime(p->rChild,B,M,LP);
    v1=p->lChild->Data->vPrime;
    v2=p->rChild->Data->vPrime;
    for(i=0;i<B.Dim;i++){
      x1=p->lChild->Data->xbar[i];
      x2=p->rChild->Data->xbar[i];
      mu=M.muPtr[i]->v;
      sig2x=M.sig2Ptr[i]->v;
      L*=fN(x1-x2+mu*(v2-v1),0.0,(v1+v2)*sig2x);
    }
  }

  return(L);
}


double CalcLogLikelihood(struct bNode *p,struct BaseInfo B,struct LocusData LD, struct ModelParamsPtr M,struct LocusParamsPtr LP){

  double LogL;
  int i;
  
  if(p->lChild==0){
    LogL=0.0;
    p->t=0.0;
    for(i=0;i<B.Dim;i++){
      p->Data->xbar[i]=LD.x[i][p->key-1];
    }
    p->Data->S2=0.0;
  }else if (p->Anc==0){
    LogL=CalcLogLikelihood(p->lChild,B,LD,M,LP);

    p->lChild->v=LP.t1Ptr->v-p->lChild->t;
    p->lChild->Data->vPrime=p->lChild->v+p->lChild->Data->S2;
    //printf("%d %g %g %d\n",p->key,p->t,LogL,p->Data->CalcFlag);
    
    if(p->Data->CalcFlag){
      for(i=0;i<B.Dim;i++){
        LogL+=LogfN(p->lChild->Data->xbar[i],LP.x0Ptr[i]->v+M.muPtr[i]->v*LP.t1Ptr->v,p->lChild->Data->vPrime*M.sig2Ptr[i]->v);
      }
    }

    //printf("CalcLogLikelihood: %d %g %g\n",p->key,p->t,LP.x0Ptr[0]->v);
    
    p->Data->LogL=LogL;
  } else{
    p->t=LP.NodeTimesPtr[2*LD.n-1-p->key]->v;
    LogL=CalcLogLikelihood(p->lChild,B,LD,M,LP);
    LogL+=CalcLogLikelihood(p->rChild,B,LD,M,LP);
    
 
    p->lChild->Data->vPrime=p->lChild->v+p->lChild->Data->S2;
    p->rChild->Data->vPrime=p->rChild->v+p->rChild->Data->S2;
        
    p->Data->S2=p->lChild->Data->vPrime*p->rChild->Data->vPrime/(p->lChild->Data->vPrime+p->rChild->Data->vPrime);
    
    for(i=0;i<B.Dim;i++){
      
      LogL+=LogfN(p->lChild->Data->xbar[i]-p->rChild->Data->xbar[i],0.0,(p->lChild->Data->vPrime+p->rChild->Data->vPrime)*M.sig2Ptr[i]->v);
      p->Data->xbar[i]=(p->rChild->Data->vPrime*p->lChild->Data->xbar[i]+p->lChild->Data->vPrime*p->rChild->Data->xbar[i])/(p->rChild->Data->vPrime+p->lChild->Data->vPrime);      
      //printf("%d %d %d %g %g",p->key,p->lChild->key,p->rChild->key,p->Data->xbar[i],LogL);
    }
    //printf("\n");
  }
  return(LogL);
}


void Calc_xbar_vprime(struct bNode *p,struct BaseInfo B,struct LocusData LD, struct LocusParams * LP){

  int i;
  
  if(p->lChild==0){
    p->t=0.0;
    for(i=0;i<B.Dim;i++){
      p->Data->xbar[i]=LD.x[i][p->key-1];
    }
    p->Data->S2=0.0;
  }else if (p->Anc==0){
    Calc_xbar_vprime(p->lChild,B,LD,LP);
    p->lChild->v=LP->t1.v-p->lChild->t;
    p->lChild->Data->vPrime=p->lChild->v+p->lChild->Data->S2;
    //printf("%d %g %d\n",p->key,p->t,p->Data->CalcFlag);
    
    
    //printf("CalcLogLikelihood: %d %g %g\n",p->key,p->t,LP.x0Ptr[0]->v);
    
  } else{
  
    Calc_xbar_vprime(p->lChild,B,LD,LP);
    Calc_xbar_vprime(p->rChild,B,LD,LP);
    p->t=LP->NodeTimes[2*LD.n-1-p->key].v;
    p->lChild->v=p->t-p->lChild->t;
    p->rChild->v=p->t-p->rChild->t;
    p->lChild->Data->vPrime=p->lChild->v+p->lChild->Data->S2;
    p->rChild->Data->vPrime=p->rChild->v+p->rChild->Data->S2;
        
    p->Data->S2=p->lChild->Data->vPrime*p->rChild->Data->vPrime/(p->lChild->Data->vPrime+p->rChild->Data->vPrime);
    
    for(i=0;i<B.Dim;i++){
      p->Data->xbar[i]=(p->rChild->Data->vPrime*p->lChild->Data->xbar[i]+p->lChild->Data->vPrime*p->rChild->Data->xbar[i])/(p->rChild->Data->vPrime+p->lChild->Data->vPrime);      
      //printf("%d %d %d %g ",p->key,p->lChild->key,p->rChild->key,p->Data->xbar[i]);
    }
    //printf("\n");
  }
}


double CalcLogLikelihoodGiven_xbar_vprime(struct bNode *p,struct BaseInfo B,struct LocusData LD, struct ModelParamsPtr M,struct LocusParamsPtr LP){

  double LogL;
  int i;
  
  if(p->lChild==0){
    LogL=0.0;
  }else if (p->Anc==0){
    LogL=CalcLogLikelihoodGiven_xbar_vprime(p->lChild,B,LD,M,LP);    
    if(p->Data->CalcFlag){
      for(i=0;i<B.Dim;i++){
        LogL+=LogfN(p->lChild->Data->xbar[i],LP.x0Ptr[i]->v+M.muPtr[i]->v*LP.t1Ptr->v,p->lChild->Data->vPrime*M.sig2Ptr[i]->v);
      }
    }
    //printf("%d %g %g\n",p->key,p->t,LogL);
    
    p->Data->LogL=LogL;
  } else{
    LogL=CalcLogLikelihoodGiven_xbar_vprime(p->lChild,B,LD,M,LP);
    LogL+=CalcLogLikelihoodGiven_xbar_vprime(p->rChild,B,LD,M,LP);
    
    for(i=0;i<B.Dim;i++){
      
      LogL+=LogfN(p->lChild->Data->xbar[i]-p->rChild->Data->xbar[i],0.0,(p->lChild->Data->vPrime+p->rChild->Data->vPrime)*M.sig2Ptr[i]->v);
      //printf("%d %g %g",p->key,p->Data->xbar[i],LogL);
    }
    //printf("\n");
  }
  return(LogL);
}


inline double LogfN(double x, double y,double z){
    return -(x - y) * (x - y)*0.5 /z-0.5*log(2.0 * 3.14159 * z);
}


/* This function derived from Gavin Cawley's
"On a Fast, Compact Approximation of the exponential function"
Neural Computation 2000 12:2009-2012
*/
inline double fastexp(double y){

  union{ double d;
  #ifdef LITTLE_ENDIAN
    struct {int j,i;}n;
  #else
    struct {int i,j;}n;
  #endif
  }_eco;
  printf("I found conditions where the approximation implied in the fastexp function fails!  Miserably!\n");

  _eco.n.i=(int)(EXP_A*(y))+(1072693248-EXP_C);
  _eco.n.j=0;
  return _eco.d;
  
}


double fN(double x, double y, double z) {  //  Normal distribution
	double xx, safeexp(double x);

	if (z > 0.0) {
		xx = (x - y) * (x - y) / (2.0 * z);
		//return safeexp(-xx) / sqrt(2.0 * 3.14159 * z);
		return exp(-xx) / sqrt(2.0 * 3.14159 * z);
		}
	else  {
		printf("*** Error in fN, x = %g, y = %g, z = %g ***\n", x, y, z);
		exit(0);
		}
	}

double safeexp(double x) { 
	if (x > -300.0) return exp(x);
	else return 0.0;
}


void CalcCvector(struct bNode *p, int nTips){
  int i;
  double vNorm;
  // If a tip;
  if(p->lChild==0){
    for(i=0;i<nTips;i++){
      if(i+1==p->key)
        p->Data->c[i]=1.0;
      else
        p->Data->c[i]=0.0;

      //printf("%g ",p->Data->c[i]);
    }
    //printf("%d <-Tip \n",p->key);
  }else if(p->Anc==0){
    CalcCvector(p->lChild,nTips);
    for(i=0;i<nTips;i++){
      p->Data->c[i]=(p->lChild->Data->c[i])/p->lChild->Data->vPrime;
      //printf("%g ",p->Data->c[i]);

    }
    //printf("%d <-Root C vector \n",p->key);
  } else{
    CalcCvector(p->lChild,nTips);
    CalcCvector(p->rChild,nTips);

  
    vNorm=p->lChild->Data->vPrime+p->rChild->Data->vPrime;
    for(i=0;i<nTips;i++){
        if(p->lChild->Data->c[i]){
          p->Data->c[i]=(p->lChild->Data->c[i])*p->rChild->Data->vPrime/vNorm;          
        }
        if(p->rChild->Data->c[i]){
          p->Data->c[i]=(p->rChild->Data->c[i])*p->lChild->Data->vPrime/vNorm;
        }
    //printf("%g ",p->Data->c[i]);

    }
    //printf("%d %d %d<-Interior\n",p->key,p->lChild->key,p->rChild->key);
  }
}



double CalcLikelihoodExpo(struct bNode *p,int d,struct ModelParamsPtr M,struct LocusParamsPtr LP){

  double L1,L2,L3;
 
  
  if(p->lChild==0){
    return(0.0);
  }else if (p->Anc==0){
    L1=CalcLikelihoodExpo(p->lChild,d,M,LP);
    if(p->Data->CalcFlag)
      L2=pow(p->lChild->Data->xbar[d]-LP.x0Ptr[d]->v-M.muPtr[d]->v*LP.t1Ptr->v,2)/p->lChild->Data->vPrime;
    else
      L2=0;
    //printf("%d x1:%g x0:%g mu:%g v:%g L: %g \n",p->key,p->lChild->Data->xbar[d],LP.x0Ptr[d]->v,M.muPtr[d]->v,p->lChild->Data->vPrime,L1+L2);      
    return(L1+L2);
  
  } else{
    L1=CalcLikelihoodExpo(p->lChild,d,M,LP);
    L2=CalcLikelihoodExpo(p->rChild,d,M,LP);
    L3=pow(p->lChild->Data->xbar[d]-p->rChild->Data->xbar[d],2)/(p->lChild->Data->vPrime+p->rChild->Data->vPrime);
    //printf("%d x1:%g x0:%g mu:%g v1:%g v2:%g L: %g\n",p->key,p->lChild->Data->xbar[d],p->rChild->Data->xbar[d],M.muPtr[d]->v,p->lChild->Data->vPrime,p->rChild->Data->vPrime,L1+L2+L3);      


    return(L1+L2+L3);
  }
}

double CalcLogVarCovarDet(struct bNode *p, struct ModelParamsPtr M,struct LocusParamsPtr LP ){

  double D1,D2,D3;
  if(p->lChild==0){
    return(0.0);
  }else if (p->Anc==0){
    D1=CalcLogVarCovarDet(p->lChild,M,LP);
    if(p->Data->CalcFlag)
      D2=log(p->lChild->Data->vPrime);
    else
      D2=0;
    //printf("%d vPrime:%g D:%g\n",p->key,p->lChild->Data->vPrime,D1+D2);      
    return(D1+D2);
  
  } else{
    D1=CalcLogVarCovarDet(p->lChild,M,LP);
    D2=CalcLogVarCovarDet(p->rChild,M,LP);
    D3=log(p->lChild->Data->vPrime+p->rChild->Data->vPrime);
    //printf("%d vPrime:%g vPrime:%g D: %g\n",p->key,p->lChild->Data->vPrime,p->rChild->Data->vPrime,D1+D2+D3);      
    return(D1+D2+D3);

    
  }
 };

void CalcLogPrt1Givenjf(struct BaseInfo B, struct LocusData LD, struct LocusParamsPtr * L){

  L->LogPrt1GivenjfPtr->v=(LD.n-1)*log(L->fPtr->v*L->t1Ptr->v)-(LD.n+1)*log(2+L->fPtr->v*L->t1Ptr->v);
  /* Next line: normalizing constant (optional for many MCMC calculations) */    
  L->LogPrt1GivenjfPtr->v+=log(2*L->fPtr->v*LD.n);
}

void CalcLogPrNodeTimesGivenjf(struct BaseInfo B, struct LocusData LD, struct LocusParamsPtr * L){


  int j;
  L->LogPrNodeTimesGivenjfPtr->v=0.0;
  for (j=2;j<=LD.n;j++){
    L->LogPrNodeTimesGivenjfPtr->v-=2*log(2+L->fPtr->v*L->NodeTimesPtr[j-2]->v);
    /* Next line: normalizing constant (optional for many MCMC calculations) */
    L->LogPrNodeTimesGivenjfPtr->v+=log(j-1)+log(2*(2+L->fPtr->v*L->t1Ptr->v))-log(L->t1Ptr->v);

  }
}


void CalcLogPriorPrmu(struct BaseInfo B,struct ModelParamsPtr M){
  printf("Fix prior on mu to be true conjugate prior\n");
  int i;
  for(i=0;i<B.Dim;i++){
    M.LogPriorPrmuPtr[i]->v=-pow((M.muPtr[i]->v-B.mu0[i]),2.0)/(2.0*(B.tau2mu0[i]));
    /* Next line: normalizing constant (which is optional for many MCMC calculations) */
    M.LogPriorPrmuPtr[i]->v-=0.5*log(2.0*PI*B.tau2mu0[i]);
  }
}

void CalcLogPriorPrx0(struct BaseInfo B,struct LocusParamsPtr * L){

  int i;
  for(i=0;i<B.Dim;i++){
    L->LogPriorPrx0Ptr[i]->v=-pow(L->x0Ptr[i]->v-B.x00[i],2)/(2.0*B.tau2x00[i]);
    /* Next line: normalizing constant (which is optional for many MCMC calculations) */
    L->LogPriorPrx0Ptr[i]->v-=0.5*log(2.0*PI*B.tau2x00[i]);
  }
}


void CalcLogPriorPrsig2(struct BaseInfo B, struct ModelParamsPtr M){
  int i;
  for (i=0;i<B.Dim;i++){
    
    M.LogPriorPrsig2Ptr[i]->v=-(B.nu0[i]/2.0+1.0)*log(M.sig2Ptr[i]->v)-B.nu0[i]*B.sig20[i]/(2.0*M.sig2Ptr[i]->v);
    /* Next line: normalizing constant (which is optional for many MCMC calculations) */
    M.LogPriorPrsig2Ptr[i]->v+=(B.nu0[i]/2.0)*log(B.nu0[i]/2.0)+(B.nu0[i]/2.0)*log(B.sig20[i])-LogGamma(B.nu0[i]/2.0);
    //printf("sig2: %g nu0:%g sig20: %g LogPriorPrsig2: %g\n",M.sig2Ptr[i]->v,B.nu0[i],B.sig20[i],M.LogPriorPrsig2Ptr[i]->v);

  }
}

void SampleTimesRandom(struct BaseInfo B, struct LocusData LD, struct LocusParams *LP){

  int i;
  double times[MAX_RARE_ALLELES];
  // Initialize the times with a random draw from the 
  // prior
  bdtimes(times,LP->f.v,0.0,0.0,LD.n);
  // These shifts are in place to deal with the format
  // of times vector from Monty's bdtimes();
  LP->t1.v=times[1];
  for(i=0;i<LD.n-1;i++)
    LP->NodeTimes[i].v=times[i+2];
  if(B.Debug>5 ||B.FixTimes){
    printf("Times: %g ",LP->t1.v);
    for(i=0;i<LD.n-1;i++)
      printf("%g ",LP->NodeTimes[i].v);
    printf("\n");
  }

};


void SampleTopoRandom(struct BaseInfo B, struct LocusData LD, struct LocusParams *LP){

  int i,j;
  
  int Node1;
  int Node2;
  int Ctr;

  struct bNode** ActiveList;
  int nActive;

  double vPrime1, vPrime2,x1,x2;


  // Initialize the Topo with a random topology
  // To-do: Write code to initialize with UPGMA topology
  
  // Declare Mem for Active List
  ActiveList=(struct bNode **)calloc((size_t)(LD.n),sizeof(struct Node *) );


  // Init Active List

  for(i=1;i<=LD.n;i++){
    LP->NodeList[i].key=i;
    LP->NodeList[i].Data->LogL=0;
    LP->NodeList[i].t=0;
    for(j=0;j<B.Dim;j++){
      LP->NodeList[i].Data->xbar[j]=LD.x[j][i-1];
    }
    ActiveList[i-1]=&LP->NodeList[i];
  }

  nActive=LD.n;
  Ctr=LD.n;

  // Store LevelsList
  for(i=0;i<nActive;i++)
    LP->LevelsList[nActive-1][i]=ActiveList[i];


  while(nActive>1){

    // Join two from ActiveList to make a new node
    // Note: Node 0 is reserved for the root 
    Node1=UniformRV(0,nActive-1);
    Node2=DuwoxRV(0,nActive-1,Node1);
    Ctr++;

    // Initialize the node in the NodeList
    LP->NodeList[Ctr].lChild=ActiveList[Node1];
    ActiveList[Node1]->Anc=&LP->NodeList[Ctr];
    LP->NodeList[Ctr].rChild=ActiveList[Node2];
    ActiveList[Node2]->Anc=&LP->NodeList[Ctr];
    LP->NodeList[Ctr].key=Ctr;
    
    /* Initialize raw branch lenght information */
    LP->NodeList[Ctr].t=LP->NodeTimes[nActive-2].v;
    LP->NodeList[Ctr].lChild->v=LP->NodeList[Ctr].t-LP->NodeList[Ctr].lChild->t;
    LP->NodeList[Ctr].rChild->v=LP->NodeList[Ctr].t-LP->NodeList[Ctr].rChild->t;
    
    /* Initialize vPrime branch length information */
    vPrime1=LP->NodeList[Ctr].lChild->v+LP->NodeList[Ctr].lChild->Data->S2;
    vPrime2=LP->NodeList[Ctr].rChild->v+LP->NodeList[Ctr].rChild->Data->S2;
    LP->NodeList[Ctr].lChild->Data->vPrime=vPrime1;
    LP->NodeList[Ctr].rChild->Data->vPrime=vPrime2;

    /* Calculate S2 */
    LP->NodeList[Ctr].Data->S2=vPrime1*vPrime2/(vPrime1+vPrime2);
    
    /* Calculate xbar */
    for(i=0;i<B.Dim;i++){
      x1=ActiveList[Node1]->Data->xbar[i];
      x2=ActiveList[Node2]->Data->xbar[i];
      LP->NodeList[Ctr].Data->xbar[i]=(vPrime2*x1+vPrime1*x2)/(vPrime1+vPrime2);
    }

    

    if(B.Debug>6)
      printf("Joining %d %d to make node %d\n",ActiveList[Node1]->key,ActiveList[Node2]->key,Ctr);


    // Add a pointer to it to the ActiveList; and overwrite Node1
    ActiveList[MIN(Node1,Node2)]=&LP->NodeList[Ctr];

    // Shift list down to overwrite Node2
    for(i=MAX(Node1,Node2);i<nActive-1;i++){
      ActiveList[i]=ActiveList[i+1];
    }
    nActive--;
    // Store LevelsList
    for(i=0;i<nActive;i++)
      LP->LevelsList[nActive-1][i]=ActiveList[i];
    
  }

  
  LP->NodeList[0].lChild=ActiveList[0];
  ActiveList[0]->Anc=&LP->NodeList[0];
  LP->NodeList[0].key=0;
  
  LP->NodeList[0].t=LP->t1.v;
  LP->NodeList[0].lChild->v=LP->NodeList[0].t-LP->NodeList[0].lChild->t;
  LP->NodeList[0].lChild->Data->vPrime=LP->t1.v-LP->NodeList[Ctr].t+LP->NodeList[Ctr].Data->S2;
  LP->NodeList[0].Data->CalcFlag=1; 
  
  LP->NodeList[0].Data->Logwr=0;
 
  free(ActiveList);

  if(B.Debug>9||B.FixTopo==1){
    printf("Output of tree level-by-level:\n");
    for(i=LD.n-1;i>=0;i--){
      printf("Level %d: ",i);
      for(j=0;j<=i;j++){
	printf("%d ",LP->LevelsList[i][j]->key);
      }
      printf("\n");
    }
  }
  

};

void CalcSig2MarginParams(struct BaseInfo B,int d,struct LocusData *LD,struct LocusParamsPtr *L,struct ModelParamsPtr * M,double * sumSamSizes, double *sumSig,double * Det){

    int i;
    double s,r,t;
    
    s=0.0;
    r=0.0;
    t=0.0;
    for(i=0;i<B.nLoci;i++){
      s+=CalcLikelihoodExpo(L[i].NodeListPtr[0],d,*M,L[i]);
      r+=LD[i].n-1;
      if(L[i].NodeListPtr[0]->Data->CalcFlag)
        r++;
      if(Det)
        t+=CalcLogVarCovarDet(L[i].NodeListPtr[0],*M,L[i]);
    }
    
    *sumSig=s;
    *sumSamSizes=r;
    if(Det)
      *Det=t;

}





void ReadTreeTimes(struct BaseInfo B, struct LocusData * LD,struct LocusParams *LP){

  char *treeString;
  int i,l;
  FILE * infile;
  char filename[80];

  strcpy(filename,B.InFileName);
  strcat(filename,".trees_in");

  if(!(infile=fopen(filename,"r"))){
    fprintf(stderr,"Error opening input file %s.\n",filename);
  }


 
  // Declare memory for treeString
  treeString=(char *)calloc((size_t)2*1000*40,sizeof(char));
  if(!treeString){
    fprintf(stderr,"Error declaring memory for tree.  Too many nodes?\n");
    exit(1);
  }
  
  for(l=0;l<B.nLoci;l++){
    fscanf(infile,"%s",treeString);
    ReadbTree(treeString,LP[l].NodeList,&ReadbNodeDataSimple);

    // Set true root
    LP[l].NodeList[0].key=0;
   
    //LP->NodeList[0].lChild=MRCA;
    //LP->t1=MRCA->t+MRCA->v;
  
    if(B.Debug>10){
      PrintbTreeWithData(stdout,&LP[l].NodeList[0],&PrintbNodeDataSimple);
      printf("\n");
    }
    // Set Node Times
    for(i=0;i<LD[l].n-1;i++){
      LP[l].NodeTimes[i].v=LP[l].NodeList[2*LD[l].n-i-1].t;
    }
    LP[l].t1.v=LP[l].NodeList[0].t;

  }
     
  free(treeString);
  fclose(infile);
};

void OutputTreeTimes(struct BaseInfo B, struct LocusData * LD,struct LocusParams *LP){


  int l;
  FILE * outfile;
  char filename[80];

  strcpy(filename,B.InFileName);
  strcat(filename,".trees_out");

  if(!(outfile=fopen(filename,"w"))){
    fprintf(stderr,"Error opening input file %s.\n",filename);
  }


  for(l=0;l<B.nLoci;l++){
  
    PrintbTreeWithData(outfile,&LP[l].NodeList[0],&PrintbNodeDataSimple);
  
  }
  fclose(outfile);
};

