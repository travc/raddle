#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#ifndef __RaddleMCMC__
#define __RaddleMCMC__

#include "RaddleBase.h"
#include "JN_MCTypes.h"

#define MAX_RARE_ALLELES 200
int gECA_mall_call_SetFlag;
double gECA_mall_call_bytes;


struct BaseInfo{

  int Dim;
  int nLoci;
  
  
  double * mu0;  /* Mean of the gaussian prior on mu */
  double * tau2mu0; /* Variance of the gaussian prior on mu */
  double * nu0; /* For the df of the inverse chi-square prior on sig2 */  
  double * sig20; /* For the sigma2 parameter of the inverse chi-square prior on sig2 */
  double * x00; /* Mean of the gaussian prior on x0. Note: constant across loci */
  double * tau2x00; /* Variance of the guassian prior on x0.  Note: constant across loci */ 

 
  char InFileName[100];
  char FileNameBase[100];
  FILE * outfile;
  long seed1;
  long seed2;
  int Debug;

  int BurnIn;
  int ChainSampleInterval;
  int MaxReps;
  
  int OutputTimes;
  int OutputTopo;
  int OutputScaleAccept;
  int OutputTopoTimesAccept;
  /* Variables to control whether the following values are updated during 
  the MCMC iteration */
  int FixTopo;
  int FixTimes;
  int Fix_mu;
  int Fix_sig2;
  int Fix_x0;
  
  int ScaleTimesMuSig2;
  int ScaleMux0;
  int UserInputTopoTimes;
  
  /* Variables to control the prior used for mu,sig2,x0 */
  int Prior_mu;
  int Prior_sig2;
  int Prior_x0;

  double TimesProposeSD;
  double fProposeVar;
  double t1ProposeSD;
  double TopoISHeat;
  double alphaProposeSD;
  double betaProposeSD;
  double Initx0SD;
  
  
};

struct ChainInfo{

  struct BaseInfo B;
  
  struct LocusData * LD;

  int NumSweepsAfterBurnIn;  
  int NumSweeps;



  /* Two sets of locus params: for current and proposed (but not necc in that order)  */
  struct LocusParams *LP;
  struct LocusParams * LP2; 

  /* Two sets of model params: for current and proposed (but not necc in that order)  */
  struct ModelParams MP;
  struct ModelParams MP2;


  /* Pointers for current versus proposed parameters */
  struct LocusParamsPtr *LPcurrent;  // Array of pointers to the current locus params
  struct LocusParamsPtr *LPprime;  // Array of pointers to proposed locus params
  struct ModelParamsPtr MPcurrent;
  struct ModelParamsPtr MPprime;

  /* Workspace for UpdateTopoAndTimes */
  double WorkSpace[MAX_RARE_ALLELES*(MAX_RARE_ALLELES-1)/2];

  struct doubval *MHTopoAccept;
  struct doubval *MHTimesAccept;
  struct doubval MHScaleTimesMuSig2Accept;
  struct doubval MHScaleMux0;

};


/* Initialization and Mem Alloc Functions */
struct ChainInfo InitSettingsAndData(int argc, char **argv);

/* MCMC Initialization Functions */
void InitChain(struct ChainInfo * C);
void Initx0(struct BaseInfo B, struct LocusData LD,struct LocusParams* LP);
void InitTimes(struct BaseInfo B,struct LocusData LD, struct LocusParams *LP,struct ModelParams M);
void InitTopo(struct BaseInfo B, struct LocusData LD, struct LocusParams* LP,struct ModelParams MP);
void Initsig2(struct BaseInfo B,struct ModelParamsPtr MP, struct LocusData *LD,struct LocusParamsPtr *LP);

/* MCMC functions */
void DoASingleSweep(struct ChainInfo *C);

void IncrementValues(struct ChainInfo *C);
void OutputHistograms(struct ChainInfo *C);
void InitOutputChainState(struct ChainInfo * C);
void OutputChainState(int IterNo,struct ChainInfo C);

/* Proposal Functions */
void UpdateTopoAndTimes(struct BaseInfo B,struct LocusData LD, struct ModelParamsPtr M, struct LocusParamsPtr LPcurrent,struct LocusParamsPtr LPprime, double *Workspace);
void Updatex0(struct BaseInfo B,struct LocusData LD, struct ModelParamsPtr M, struct LocusParamsPtr * LPcurrent, struct LocusParamsPtr * LPprime);
void Updatemu(struct BaseInfo B, struct LocusData * LD,struct LocusParamsPtr * LP, struct ModelParamsPtr * M,struct ModelParamsPtr * MP);

void ProposeTopoIS(struct BaseInfo B,struct LocusData LD, struct ModelParamsPtr M, struct LocusParamsPtr LPcurrent,struct LocusParamsPtr LPprime);
void DoMHStepTopoIS(struct BaseInfo B,struct doubval * Accept,struct LocusData LD, struct ModelParamsPtr M, struct LocusParamsPtr LPcurrent,struct LocusParamsPtr LPprime);

void ProposeScalingTimesMuSig2(struct BaseInfo B,struct LocusData * LD,struct LocusParamsPtr * L, struct LocusParamsPtr * LP,struct ModelParamsPtr M,struct ModelParamsPtr MP);
void DoMHStepScalingTimesMuSig2(struct BaseInfo B,struct doubval * Accept, struct LocusData * LD,struct LocusParamsPtr * L, struct LocusParamsPtr * LP,struct ModelParamsPtr * M,struct ModelParamsPtr * MP);

void ProposeScalingMux0(struct BaseInfo B,int d,struct LocusData * LD,struct LocusParamsPtr * L, struct LocusParamsPtr * LP,struct ModelParamsPtr M,struct ModelParamsPtr MP);
void DoMHStepScalingMux0(struct BaseInfo B,struct doubval * Accept,int d,struct LocusData * LD,struct LocusParamsPtr * L, struct LocusParamsPtr * LP,struct ModelParamsPtr * M,struct ModelParamsPtr * MP);


void ProposeTimes(struct BaseInfo B, struct LocusData LD,struct LocusParamsPtr L,struct LocusParamsPtr LP);
void DoMHStepTimes(struct BaseInfo B, struct doubval * Accept,struct LocusData LD,struct LocusParamsPtr * L,struct LocusParamsPtr * LP, struct ModelParamsPtr M);


void Updatesig2(struct BaseInfo B, struct LocusData * LD,struct LocusParamsPtr * LP,struct ModelParamsPtr * M,struct ModelParamsPtr *MP);

void PrintbNodeData(FILE * outfile,struct bNode p);

#endif




