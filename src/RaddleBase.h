/*
 *  RaddleBase.h
 *  Raddle
 *
 *  Created by John Novembre on 12/6/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __RaddleBase__
#define __RaddleBase__


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "ranlib.h"
#include "MathStatRand.h"
#include "ECA_utilities.h"
//#include "ECA_MemAlloc.h"


#include "JN_MemAlloc.h"
#include "JN_MCTypes.h"
#include "JN_RVutils.h"
#include "JN_TreeUtils.h"


//#include RADDLE_H 
/* RADDLE_H is either RaddleMCMC.h or mainMC.h */


#include<stdio.h>
#include<math.h>
#define min(x, y)  (((x) < (y)) ? x : y)
#define max(x, y)  (((x) > (y)) ? x : y)

#define BIGTIME 1e300
#define PI 3.1415926535897932384626433832795

struct bNodeData{

  //double L;
  double LogL;
  double LogLk;
  double Logwr;
  double Logwrk;
  int CalcFlag;
  double *xbar;
  double S2;  
  double vPrime;
  double *c;
};


struct LocusData{
  char ID[80];
  double f;
  int n; // the # of rare alleles for the locus
  double **x;  // multidim array: n * Dim.   

};

struct LocusParams{
  /* True parameters */
  struct doubval *x0;  // Dim-dimensional root position
  struct bNode *NodeList;  // 2n-1 array representing genealogy
  struct bNode ***LevelsList;
  struct doubval t1;
  struct doubval *NodeTimes;
  struct doubval f;
  
  /* Prior probs on parameters */
  struct doubval LogPrt1Givenjf;
  struct doubval LogPrNodeTimesGivenjf;
  struct doubval * LogPriorPrx0;
};

struct LocusParamsPtr{
  
  struct doubval **x0Ptr;  // Dim-dimensional root position
  struct bNode **NodeListPtr;  // 2n-1 array representing genealogy
  struct bNode ****LevelsListPtr;
  struct doubval *t1Ptr;
  struct doubval **NodeTimesPtr;
  struct doubval *fPtr;
  
  struct doubval * LogPrt1GivenjfPtr;
  struct doubval * LogPrNodeTimesGivenjfPtr;
  struct doubval ** LogPriorPrx0Ptr;
};


struct ModelParams{
  /* True parameters */
  struct doubval *mu;
  struct doubval *sig2;
  
  /* Prior probs of parameters */
  struct doubval *LogPriorPrmu;
  struct doubval *LogPriorPrsig2;
  };

struct ModelParamsPtr{

  struct doubval **muPtr;
  struct doubval **sig2Ptr;
  
  struct doubval **LogPriorPrmuPtr;
  struct doubval **LogPriorPrsig2Ptr;
  
};

struct BaseInfo{

  int Dim;
  int nLoci;
  
  char InFileName[100];
  char FileNameBase[100];
  FILE * outfile;
  long seed1;
  long seed2;
  int Debug;
  time_t StartTime;
  int MaxReps;
  double sig2Min1d, sig2Max1d, sig2Delt1d, sig2nvals1d;
  double gsla,gslb,gslRelErr,gslMaxIter,gslCIperLocus,gsldfeps;
  double * sig2Min, * sig2Max, * sig2Delt, * sig2nvals;
  int sig2nvalsTot;
  
  /* Next four lines for options not currently really being used */
  int * nsig2muSurfPoints;
  int ** nx0sig2muSurfPoints;
  double ** x0Min, ** x0Max, ** x0nvals, ** x0Delt;
  double * muMin, * muMax, * munvals, *muDelt;
  
  double sMin,sMax,snvals;
  double muDrive,sig2Drive;
  double TopoISHeat;
  int UseIS;
  int UnbiasedISestimate; /* 1=Eqn 2.4 Liu 2001 0=Eqn 2.3 Liu 2001 */
  int LongLikelihoodCalc;
  int Calc1DimSurf;
  int CalcMultiDimSurf;
  int CalcMLEperLocus;
  int OutputTimes;
  int OutputTopo;
  int Outputsig2MarginSurf;
  int Outputsig2MarginAllLociSurf;
  int Outputx0sig2muSurf;
  
  /* Variables to control whether the following values are updated during 
  the MC iteration */
  int FixTopo;
  int FixTimes;
  
  int UserInputTopoTimes;
  
  /* Variables to control the prior used for mu,sig2,x0 */
  int Prior_mu;
  int Prior_sig2;
  int Prior_x0;
  
   
  double * mu0;  /* Mean of the gaussian prior on mu */
  double * tau2mu0; /* Variance of the gaussian prior on mu */
  double * nu0; /* For the df of the inverse chi-square prior on sig2 */  
  double * sig20; /* For the sigma2 parameter of the inverse chi-square prior on sig2 */
  double * x00; /* Mean of the gaussian prior on x0. Note: constant across loci */
  double * tau2x00; /* Variance of the guassian prior on x0.  Note: constant across loci */ 


  
};




/* Probability functions */
double CalcJoint(struct BaseInfo B,struct ModelParams M, struct LocusData * LD,struct LocusParams ** LP, int N);
double CalcLikelihood(struct bNode *p,struct BaseInfo B,struct LocusData LD, struct ModelParamsPtr M,struct LocusParamsPtr LP);
double CalcLikelihoodGiven_xbar_vprime(struct bNode *p,struct BaseInfo B,struct ModelParamsPtr M,struct LocusParamsPtr LP);
double CalcLogLikelihood(struct bNode *p,struct BaseInfo B,struct LocusData LD, struct ModelParamsPtr M,struct LocusParamsPtr LP);
double CalcLogLikelihoodGiven_xbar_vprime(struct bNode *p,struct BaseInfo B,struct LocusData LD, struct ModelParamsPtr M,struct LocusParamsPtr LP);
void Calc_xbar_vprime(struct bNode *p,struct BaseInfo B,struct LocusData LD, struct LocusParams * LP);
void CalcISweight(struct BaseInfo B, struct LocusData LD, struct LocusParamsPtr LP,struct ModelParamsPtr MP, double TopoISHeat);

double CalcLikelihoodExpo(struct bNode *p,int d,struct ModelParamsPtr M,struct LocusParamsPtr LP);
double CalcLogVarCovarDet(struct bNode *p,struct ModelParamsPtr M,struct LocusParamsPtr LP);

void CalcLogPrt1Givenjf(struct BaseInfo B, struct LocusData LD, struct LocusParamsPtr * L);
void CalcLogPrNodeTimesGivenjf(struct BaseInfo B, struct LocusData LD, struct LocusParamsPtr * L);
void CalcLogPriorPrmu(struct BaseInfo B,struct ModelParamsPtr M);
void CalcLogPriorPrx0(struct BaseInfo B,struct LocusParamsPtr * L);
void CalcLogPriorPrsig2(struct BaseInfo B, struct ModelParamsPtr M);
void CalcSig2MarginParams(struct BaseInfo B,int d,struct LocusData *LD,struct LocusParamsPtr *L,struct ModelParamsPtr * M,double * nu, double *nus2,double *Det);

double Prsig2(struct BaseInfo B,struct ModelParams M);
double PrTopo(struct LocusParams L);
double Prf(struct BaseInfo B,struct ModelParams  M);
double Prx0(struct BaseInfo B,struct LocusParams L);
double Prt1(struct BaseInfo B,struct ModelParams M,struct LocusData LD,struct LocusParams LP);
double PrDataGivenAllElse(struct BaseInfo B,struct ModelParams M,struct LocusData LD,struct LocusParams LP);
double PrNodeTimes(struct BaseInfo B,struct ModelParams M,struct LocusData LD,struct LocusParams LP);
double fN(double x, double y, double z);
inline double LogfN(double x, double y,double z);

/* For fastexp */
# define EXP_A (1048576/M_LN2)
# define EXP_C 60801

inline double fastexp(double y);

/* Sampling functions */
void SampleTopoIS(int Dim, struct LocusData LD, struct LocusParamsPtr LP, struct LocusParamsPtr LPprime,struct ModelParamsPtr MP, double TopoISHeat);
void SampleNodeISwLogs(int Dim, int k,double CoalTime, struct bNode** ActiveList, struct ModelParamsPtr MP,
      double * Lk,double *wrk,int * Node1, int * Node2,double heat);
void SampleTimesRandom(struct BaseInfo B,struct LocusData LD, struct LocusParams *LP);
void SampleTopoRandom(struct BaseInfo B, struct LocusData LD, struct LocusParams* LP);
void bdtimes(double *times, double f, double r, double s, int samsize);

      
      
/* Utility functions */
void ReadCommandLine(int argc,char **argv,struct BaseInfo *B);
struct LocusData * ReadData(struct BaseInfo B);
void ReadTreeTimes(struct BaseInfo B, struct LocusData * LD,struct LocusParams *LP);
void OutputTreeTimes(struct BaseInfo B, struct LocusData * LD,struct LocusParams *LP);

struct ModelParams AllocModelParams(struct BaseInfo B);
struct LocusParams * AllocLocusParams(struct BaseInfo B,struct LocusData * LD);
struct ModelParamsPtr AllocModelParamsPtr(struct BaseInfo B);
struct LocusParamsPtr * AllocLocusParamsPtr(struct BaseInfo B,struct LocusData * LD);

void CopyLocusParams(struct BaseInfo B, struct LocusData LD,struct LocusParams LO,struct LocusParams *LC);
void CalcCvector(struct bNode *p, int nTips);



#endif
