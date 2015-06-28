#include "mainMC.h"
#include "RaddleBase.h"


void ReadCommandLine(int argc,char **argv, struct BaseInfo* B){
  int i,tmp;
  void usage();

  i=1;
  if(argc<2){
    usage();
    exit(0);
  }
  while(i<argc){
    switch(argv[i][1]){
    case 'S':
      i++;
      B->seed1=(long)atoi(argv[i]);
      B->seed2=(long)atoi(argv[i+1]);
      i++;
      break;
    case 'd':
      i++;
      B->Dim=atoi(argv[i]);
      break;
    case 'n':
      i++;
      B->nLoci=atoi(argv[i]);
      break;
    case 'M':
      i++;
      B->MaxReps=atoi(argv[i]);
      break;    
    case 'f':
      i++;
      sprintf(B->InFileName,"%s",argv[i]);
      break;
    case 'o':
      i++;
      sprintf(B->FileNameBase,"%s",argv[i]);
      break;
    case 'D':
      i++;
      B->Debug=atoi(argv[i]);
      break;
    case 'H':
      i++;
      B->TopoISHeat=atof(argv[i]);
      break;
    case 'm':
      B->UseIS=0;
      break;
    case 's':
      i++;
      B->sMin=atof(argv[i]);
      i++;
      B->sMax=atof(argv[i]);
      i++;
      B->snvals=atoi(argv[i]);
      i++;
      tmp=atoi(argv[i]);
      if(tmp){
        B->CalcMultiDimSurf=1;
        B->Calc1DimSurf=0;
      }else{
        B->CalcMultiDimSurf=0;
        B->Calc1DimSurf=1;
      }
      break;
    case 'v':
      i++;
      B->sig2Drive=atof(argv[i]);
      break;
    case 'a':
      i++;
      B->gsla=atof(argv[i]);
      break;
    case 'b':
      i++;
      B->gslb=atof(argv[i]);
      break;
    case 'P':
      
      B->CalcMLEperLocus=1;
      break;
    default:
      fprintf(stderr,"Unknown option: %s\n",argv[i]);
      usage();
      break;
    }
    i++;
  }

  setall(B->seed1,B->seed2);  
  
}

void usage(){

  fprintf(stderr,"RaddleMC <options>\n");
  fprintf(stderr,"\t=====Basic parameters=====\n");
  fprintf(stderr,"\t-S <integer> <integer>\t\t\tRandom number seeds\n");
  fprintf(stderr,"\t-d <integer> \t\t\t\tDimensionality of habitat\n");
  fprintf(stderr,"\t-n <integer> \t\t\t\tNumber of loci\n");
  fprintf(stderr,"\t-M <integer> \t\t\t\tMax Number of iterations\n");
  fprintf(stderr,"\t-f <filename> \t\t\t\tInput filename\n");
  /*fprintf(stderr,"\t-o <filename> \t\t\t\tOutput filenamebase\n");*/
  fprintf(stderr,"\t=====Parameters affecting approximation of likelihood function=====\n");
  fprintf(stderr,"\t-v <double> \t\t\t\tDriving value for sig2\n");
  fprintf(stderr,"\t-H <double> \t\t\t\tHeat for IS sampler\n");
  fprintf(stderr,"\t-m \t\t\t\t\tUse vanilla Monte Carlo\n");
  fprintf(stderr,"\t=====Parameters affecting finding MLEs=====\n");
  fprintf(stderr,"\t-P \t\t\t\t\tCalculate MLEs for each locus individually\n");
  fprintf(stderr,"\t-a <double> \t\t\t\tLower bound for sig2 search\n");
  fprintf(stderr,"\t-b <double> \t\t\t\tUpper bound for sig2 search\n");   
  fprintf(stderr,"\t-s <double> <double> <int> <bool> \tMin Sig2, Max Sig2, nGridPoints, CalcMultiDimSurf\n");
  exit(0);
}


