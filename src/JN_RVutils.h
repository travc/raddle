
#ifndef __JN_RV_UTILS_FLAG__
#define __JN_RV_UTILS_FLAG__

#include <math.h>
#include "ranlib.h"
#include "MathStatRand.h"
#include "ECA_utilities.h"
//#include "ECA_MemAlloc.h"

#define MIN(x, y)  (((x) < (y)) ? x : y)
#define MAX(x, y)  (((x) > (y)) ? x : y)


void Calc2Moments(double * x, int l, double * m);
double LogOfNormalPDF(double mean, double var, double x);
double LogOfBetaPDF(double alpha,double beta, double x);
double LogOfInvGammaPDF(double theta,double alpha,double beta);
double TruncatedNormalRV(double mean, double var, double a, double b);
double ReflectTruncatedNormalRV(double mean, double var, double a, double b);
double ltqnorm(double p);

#endif
