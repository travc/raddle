/*
 *  mainMC.h
 *  Raddle
 *
 *  Created by John Novembre on 1/24/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#define DEBUG 1

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_multimin.h>
#include<gsl/gsl_min.h>
#include<gsl/gsl_roots.h>
#ifndef __RaddleMC__
#define __RaddleMC__

#include "RaddleBase.h"
#include "JN_MCTypes.h"

#define MAX_RARE_ALLELES 200
int gECA_mall_call_SetFlag;
double gECA_mall_call_bytes;



struct MLE{
  double v;
  double lowerCI;
  double upperCI;
  double LogL;
  double LogLSE;
};


struct SamplerInfo{
  
  struct BaseInfo B;
  struct LocusData * LD;
  int NumSweeps;
  
  struct LocusParams *L;
  struct ModelParams M;
  struct ModelParams M_IS; /* Driving values of IS sampler */

  struct LocusParamsPtr *LP;  // Array of pointers to the current locus params
  struct ModelParamsPtr MP;
  struct ModelParamsPtr MP_IS;
  struct LsurfPoint ** sig2Margin1dSurf;    // nLoci * sig2vals1d 
  struct LsurfPoint * sig2MarginAllLoci1dSurf;    // sig2vals1d 
  struct MLE * sig2Margin1dMLE; // nLoci
  struct MLE sig2MarginAllLoci1dMLE; 
  
  struct LsurfPoint ** sig2MarginSurf;    // nLoci * sig2valsTot 
  struct LsurfPoint * sig2MarginAllLociSurf;    // sig2valsTot
  struct MLE ** sig2MarginMLE;  // nLoci * Dim
  struct MLE * sig2MarginAllLociMLE;  // Dim
  
  struct PostMeanEst ** sig2PostMean; // nLoci*nDim;
  double sig2PostMeanAllLoci;
  
  struct LsurfPoint *** x0sig2muSurf; // nLoci * nDim * (sig2nvals * munvals * x0nvals)
  struct RepSummaryStats ** RepSumStats;
  
  double * meanISweight;
  double * meanISweightSquare;
};

struct RepSummaryStats{

  double * sumSig;
  double * xbar;  
  double LogDet;
  double vMRCA;
  double S2_MRCA;

  double t1;
  double LogISweight;

};


struct LsurfPoint{
  int LocusNo;
  double * sig2;
  double * mu;
  double * x0;
  double f;
  double LogL;
  double LsumSquare; /* Sum of squares of L */
  double LSE; /* MC sampling error of L */
  double L;
};




#endif
