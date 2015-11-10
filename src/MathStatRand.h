/***********    MathStatRand
*
*	This is a c header file for Eric Anderson's 
*	"MathStatRand" library which includes routines that he has
*	commonly used in his programs for the last couple of years
*	and which he wanted to put together in one place
*
*	The types of functions herein are:
*		Math:	Mathematical Special functions like the gamma function
*		Stat:   Probability densities and mass functions and cdf's
*		Rand:	Functions for generating random variates that are not widely available.
*				These functions rely on the functions available in the ranlib.c
*				random number generator library
*
************/





/* MATHEMATICAL FUNCTIONS */
extern double LogGamma(double x);
extern double IncGammaCF(double a, double x);
extern double IncGammaSer(double a, double x);
extern double GammaP(double a, double x);
extern double GammaQ(double a, double x) ;
extern double ErrFunct(double x);
extern double ErrFunctComp(double x);
extern double LogFactorial(int x);
double RecursiveLogStirling3(double ***M, int n, int i, int MaxN);

/* STATISTICAL FUNCTIONS */
extern double LogMultinomialPMF(int n, double *p, int *x, int inNum);
extern double LogCMDPMF(int *N, double *a, int K);
extern double LogPEUrnModelPMF(int *N, int *X, double s, int K);
extern double LogMultinomUrnPMF(int *N, int *X, int K);
extern double LogHGUrnPMF(int *N, int *X, int K);
extern double NormalPDF(double mean, double var, double x);
extern double LogDirichletPDF(double *pars, double *rvs, int numComp);
extern double NormalCDF(double mean, double var, double x);
extern double LogisticCDF(double mu, double var, double x);
extern double NormalInvCDF(double mean, double var, double x);


/* RANDOM NUMBER GENERATION FUNCTIONS */
extern int UniformRV(int a, int n);
extern int DuwoxRV(int a, int n, int x);
extern void OrderedPairRV(int Lo, int Hi, int *i, int *j);
int IntFromProbsRV(double *Probs, int lo, int n);
extern void DirichletRV(double *par, int n, double *out);
