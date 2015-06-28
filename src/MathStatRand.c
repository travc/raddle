/***********    MathStatRand
*
*	This is a c source code file to make Eric Anderson's 
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ranlib.h"
#include "ECA_utilities.h"
#include "ECA_MemAlloc.h"
#include "MathStatRand.h"



/******************************************************************************************/
/*                                                                                        */
/*						          MATHEMATICAL FUNCTIONS                                  */
/*                                                                                        */  
/******************************************************************************************/




/*************
*  FUNCTION:  double ErrFunct
*  PURPOSE:   Compute the ErrorFunction
			  
			
*  INPUT:     double x
   
   	Requires the function GammaP.  Based on Numerical Recipes in C, page 220.
   	
   AUTHOR:  Eric C. Anderson
   DATE:  18 AUGUST 2000
*************/
double 
ErrFunct(double x) 
{
	if(x < 0.0) {
		return( -GammaP(.5, x*x));
	}
	else {
		return( GammaP(.5, x*x));
	}
}

/*************
*  FUNCTION:  double ErrFunctComp
*  PURPOSE:   Compute the ErrorFunction
			  
			
*  INPUT:     double x
   
   	Requires the function GammaP and GammaQ  Based on Numerical Recipes in C, page 220.
   	
   AUTHOR:  Eric C. Anderson
   DATE:  18 AUGUST 2000
*************/
double 
ErrFunctComp(double x) 
{
	if(x < 0.0) {
		return( 1.0 + GammaP(.5, x*x));
	}
	else {
		return( GammaQ(.5, x*x));
	}
}

/*************
*  FUNCTION:  double GammaP
*  PURPOSE:   Compute the incomplete Gamma function.  P(a,x)
			  
			
*  INPUT:     double a;
			  double x;
   
   	Requires the functions IncGammaCF and IncGammaSer
   	
   AUTHOR:  Eric C. Anderson
   DATE:  18 AUGUST 2000
*************/
double 
GammaP(double a, double x) 
{
	if(x<0.0 || a <= 0.0)  {
		printf("\n\nx<0.0 or a<=0.0 in GammaP\n\nExiting to system...");
		exit(1);
	}
	if(x <= a + 1.0) {
		return(IncGammaSer(a,x));
	}
	else {
		return(1.0 - IncGammaCF(a,x));
	}
}


/*************
*  FUNCTION:  double GammaQ
*  PURPOSE:   Compute the incomplete Gamma function.  Q(a,x).
			  This is the part integrated from x to infinity.
			  
			
*  INPUT:     double a;
			  double x;
   
   	Requires the functions IncGammaCF and IncGammaSer
   	
   AUTHOR:  Eric C. Anderson
   DATE:  18 AUGUST 2000
*************/
double 
GammaQ(double a, double x) 
{
	if(x<0.0 || a <= 0.0)  {
		printf("\n\nx<0.0 or a<=0.0 in GammaP\n\nExiting to system...");
		exit(1);
	}
	if(x <= a + 1.0) {
		return(1.0 - IncGammaSer(a,x));
	}
	else {
		return(IncGammaCF(a,x));
	}
}

/*************
*  FUNCTION:  double IncGammaSer
*  PURPOSE:   Compute the incomplete Gamma function by its series representation.  
			  Note that this computes the part from 0 to x (called P(a,x) )
			  
			  This is the preferred method when x < (a+1)
			
*  INPUT:     double x

	This function is based on the series representation described
	in "Numerical Recipes in C" (Cambridge University Press), 
   
   	Requires the function GammaLog
   	
   AUTHOR:  Eric C. Anderson
   DATE:  18 AUGUST 2000
*************/
#define MAXIT 100
#define PRECIS 1.0e-7

double 
IncGammaSer(double a, double x) 
{
	double old_funct, new_funct;  /*  this is what we ultimately want to return; */
	double term;  /*  the next term to be added to the series */
	double i;
	double numer;  /*  the normalizing constant---Gamma(a).  Also turns out to be the numerator of the terms that */
				   /*  get summed. */
	double denom;  /*  the denominator of the terms that get summed */
	
	/*  set value of the log of the normalizing constant */
	numer = exp(LogGamma(a));
	denom = a * numer;
	
	term = numer/denom;
	old_funct = term;
	
	for(i=1;i<MAXIT;i++)  {
		term *= 1.0/(a + i) * x;  /*  the terms may be computed recursively */
		new_funct = old_funct + term;
		
		/*  check for convergence */
		if( fabs(new_funct - old_funct) < PRECIS)  {
			break;
		}
		else {
			old_funct = new_funct;
		}
	}
	
	if(i >= MAXIT - .001) {
		/*  done here we warn if it failed to converge.  */
		printf("\n\nFailed to converge in IncGammaSer\n\n");
		exit(1);
	}
	
	/*  otherwise, we return the result */
	/*  first we have to add the coefficients and normalize it: */
	new_funct *= exp(-x) * pow(x,a) * 1.0/numer;
	return(new_funct);
}
#undef MAXIT
#undef PRECIS

/*************
*  FUNCTION:  double IncGammaCF
*  PURPOSE:   Compute the incomplete Gamma function by its continued
			  fraction representation.  Note that this computes the 
			  part from x to infinity (called Q(a,x) )
			  
			  This is the preferred method for x > (a + 1)
			
*  INPUT:     double x

	This function is based on the continued fraction representation described
	in "Numerical Recipes in C" (Cambridge University Press), and evaluated using
	Lentz's method for continued fractions.  Basically following the pseudocode on
	page 171.
   
   	Requires the function GammaLog
   	
   AUTHOR:  Eric C. Anderson
   DATE:  18 AUGUST 2000
*************/
#define SMIDGEN 1.0e-30
#define MAXIT 100
#define PRECIS 2.0e-7

double 
IncGammaCF(double a, double x) 
{
	double funct;  /*  this is what we ultimately want to return; */
	double A,B,C,D;
	double Delta;
	double j;
	double log_normo;  /*  the log of the normalizing constant---Gamma(a) */
	
	/*  set value of the log of the normalizing constant */
	log_normo = LogGamma(a);
	
	/*  here for the zero subscript part: */
	B = 0.0;  /*  b_0 is really zero, so we just make it tiny */
	funct = SMIDGEN;
	C = funct;
	D = 0.0;
	
	/*  then for the 1 subscript part things are sort of different still */
	A = 1.0;
	B = x + 1.0 - a;
	D = B + A * D;
	if(fabs(D) < SMIDGEN) D = SMIDGEN;
	C = B +  A/C;
	if(fabs(C) < SMIDGEN) C = SMIDGEN;
	D = 1.0/D;
	Delta = C * D;
	funct = Delta * funct;
	/*  don't even bother checking for convergence after this first "iteration" */
	
	
	/*  now we iterate through the 2,3... and so forth subscript parts until converged */
	/*  so each loop corresponds to the j+1-th subscript of the continued fraction */
	for(j=1.0;j<MAXIT;j++)  {
		/*  now, we define the new values for A and B on the j-th level of the continued fraction */
		A = -j * (j - a);
		B = x + 1.0 + (2 * j) - a;	
		
		D = B + A * D;
		if(fabs(D) < SMIDGEN) D = SMIDGEN;
		C = B + A/C;
		if(fabs(C) < SMIDGEN) C = SMIDGEN;
		D = 1.0/D;
		Delta = C * D;
		funct = Delta * funct;

		/*  here we check for convergence.  If it has, we return the appropriate result */
		if(fabs(Delta - 1.0) < PRECIS) {
			/*  add the coefficients and the normalizing constant */
			funct *= exp(-x - log_normo) * pow(x,a);
			return(funct);
		}
	}
	
	/*  done here we warn if it failed to converge.  Should never get here */
	printf("\n\nFailed to converge in IncGammaCF\n\n");
	exit(1);
	
	return(-55.55);  /*  just put this in so that it doesn't give a compile warning */
					 /*  about having no return value. */

}
#undef SMIDGEN
#undef MAXIT
#undef PRECIS

/*************
*  FUNCTION:  double LogGamma
*  PURPOSE:   Compute log of the Gamma Function
			
*  INPUT:     double x

	This function is based on the six term series of Lanczos as described
	in "Numerical Recipes in C" (Cambridge University Press).  I use the 
	choice of gamma = 5 and N = 6.  
	
*************/
double 
LogGamma(double x) 
{
	/*  declare the coefficients in the series as a static double array */
	static double Coeff[6] = {76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	double Prelim;  /*  for the part before the series */
	double Series;
	int i;
	double denom;	
	
	/*  We start off by computing Gamma(x+1) by the series: */
	
	/*  compute the preliminary part */
	Prelim = (x+.5) * log(x+5.5) - (x+5.5);

	/*  add the log of sqrt(2\pi) to that: */
	Prelim += 0.91893853320467274;
	
	/*  now compute the series.  Start with the constant term c_0 */
	Series = 1.000000000190015;

	denom = x;
	for(i=0;i<6;i++)  {
		Series += Coeff[i]/(++denom);
	}
	
	/*  now we just have to return the right thing.  Since we just computed */
	/*  Recall Gamma(x+1) = x * Gamma(x), so we have to subtract log(x) from this. */
	/*  We can do this with only one call to log by dividing Series by x. */
	return(Prelim + log(Series/x));


}

/*************
*  FUNCTION:  double LogFactorial
*  PURPOSE:   Compute the log of x!
			  
			
*  INPUT:     double x
   
   	Requires the function LogGamma
   	
   	Based on the implementation in Numerical Recipes in C.
   	
   AUTHOR:  Eric C. Anderson
   DATE:  18 AUGUST 2000
*************/
double 
LogFactorial(int x) 
{
	static double PreStored[201];  /*  we will precompute and store value for up to 200 */
	int i;
	
	if(x<0) {
		printf("\n\nNegative argument in call of LogFactorial\n\nExiting...");
		exit(1);
	}
	if(x<=1) {
		return(0.0);
	}
	if(x <= 200)  {
		if(!PreStored[x])  {  /*  if PreStored[x] is zero, then we must fill in the table */
			for(i=2;i<201;i++)  {
				PreStored[i] = LogGamma((double)i + 1.0);
			}
		}
		return(PreStored[x]);
	}
	else    /*  otherwise we have to compute it. */
		return(LogGamma((double)x + 1.0));
	
}



/*
	This recursive function returns the Stirling number of the third kind 
	as described on page 6 of Johnson, Kotz, and Balakrishnan's "Discrete
	Multivariate Distributions."  The Stirling number of the third kind,
	s(n,i) is defined as those coefficients which satisfy:
	
	x(x+1)...(x+n-1) = \sum_{i=0}^n s(n,i)x^i
	
	Which is pretty cool. 
	
	s(n,i) also satisfies the recursions:
	
		s(n,1) = (n-1)!
		
		s(n,n) = 1
		
		s(n,i) = s(n-1,i-1) + (n-1)s(n-1,i)     for i = 1,2,...,n-1
		
	These recursions are used in this function.
	
	The way you use this function:  The variable you pass in as M should be a double ***
	should be the address of a double ** variable. it
	might be simplest to give it global scope.  It is within M that this function will
	store previous results to be used in later recursions.  If M==NULL, then this function
	will allocate memory to M as needed.  If M != NULL but it doesn't have sufficient
	memory allocated to it, you could be screwed.  MaxN is the largest n that will ever
	get passed to this function. 
	
	M is subscripted by [n][i]
*/
double RecursiveLogStirling3(double ***M, int n, int i, int MaxN)
{
	int j,k;
	double norm;  /* a factor to be able to do the calculation without overflow */
	
	if(n > MaxN) {
		fprintf(stderr, "\n\nn = %d > MaxN = %d in RecursiveStirling3.",n,MaxN);
		fprintf(stderr, "\n...Exiting to system...\n\n");
		exit(1);
	}
	if(n<1 || i < 1 || i>n) {
		fprintf(stderr, "\n\nError in RecursiveStirling3. n = %d, i = %d",n,i);
		fprintf(stderr, "\n...Exiting to system...\n\n");
		exit(1);
	}
	if(*M==NULL)  {
		(*M) = (double **)ECA_CALLOC( (size_t)(MaxN+1), sizeof(double *));
		for(j=0;j<MaxN+1;j++)  {
			(*M)[j] = (double *)ECA_CALLOC( (size_t)j+1, sizeof(double));
		}
		/* initialize to -1.0's */
		for(j=0;j<=MaxN;j++)  {
			for(k=0;k<=j;k++)  {
				(*M)[j][k] = -1.0;
			}
		}
		/* do a few simple initializations */
		for(j=1;j<=MaxN;j++)  {
			(*M)[j][j] = 0.0;
			(*M)[j][1] = LogFactorial(j-1);
		}
	}
	
	/* take care of the easy cases first */
	if(n==i) {
		return(0.0);
	}
	if(i==1) {
		return((*M)[n][i]);
		
	}
	/*  now implement the recursion. */
	if((*M)[n][i]==-1.0) {
		if((*M)[n-1][i]==-1.0) {
			RecursiveLogStirling3(M,n-1,i,MaxN);
		}
		if((*M)[n-1][i-1]==-1.0) {
			RecursiveLogStirling3(M,n-1,i-1,MaxN);
		}
		
		/* now that these have been taken care of, we can return the result, but only AFTER
		   recording the result in M  */
		norm = (*M)[n-1][i-1] < (*M)[n-1][i] ? (*M)[n-1][i-1]  : (*M)[n-1][i];
		
		(*M)[n][i] = norm + log( exp((*M)[n-1][i-1] - norm) + (double)(n-1) * 
							exp((*M)[n-1][i] - norm) );
	}
	return((*M)[n][i]);
}


/******************************************************************************************/
/*                                                                                        */
/*						          STATISTICAL FUNCTIONS                                   */
/*                                                                                        */  
/******************************************************************************************/


/*************
*  FUNCTION:  double LogMultinomialPMF(int n, double *p, int *x, int inNum)
*  PURPOSE:   Return the log of the multinomial probability function
*
*  INPUT: 	  int n  -->  number of trials
			  double *p  --> the probability of success
*  			  int *x  -->  the number of successes
			  int inNum --> the number of categories
			  
	Originally function log_multi_prob
*************/
double 
LogMultinomialPMF(int n, double *p, int *x, int inNum)
{
	int k;  /*  subscript over categories */
	double log_prob; 
	int check_n;
	double check_p;
	
	
	log_prob = LogFactorial(n);  /*  initialize to log n! to accumulate the sum */
	check_n = 0;
	check_p = 0.0;
	
	for(k=0;k<inNum;k++) {
		if(p[k] != 0 ) {
			log_prob += x[k] * log(p[k]) - LogFactorial(x[k]);
		}
		else if(x[k] != 0) {
			printf("\n\nYou're bumming my friend.  You have a zero probability event \nin LogMultinomialPMF!\n\n");	
			printf("\nCrashing and burning\n");
			printf("\n...Exiting to system...\n");
			exit(1);
		}
		check_n += x[k];
		check_p += p[k];
	}
	
	if(check_n != n) {
		printf("\n\nSum of x's is not n in LogMultinomialPMF\n\n");
		printf("\nThis cannot be tolerated!!  Get me out of here!");
		printf("\n...Exiting to system...\n");
		exit(1);
	}
	if(check_p < .99999999999 || check_p > 1.00000000001) {
		printf("\n\nSum of p's is not 1.0 in LogMultinomialPMF\n\n");
		printf("\nThis is a grave problem!!  Going away now!");
		printf("\n...Exiting to system...\n");
		exit(1);
	}
	return(log_prob);
}




double 
LogCMDPMF(int *N, double *a, int K)
{
/*************
*  FUNCTION:  double LogCMDPMF(int *N, double *a, int K)
*  PURPOSE:   To compute the log of the probability of a Dirichlet-compound multinomial random
			  vector N, given the parameters (a_1,...,a_K) of the Dirichlet dsn.
			  
			  This function follows the pmf given in "Discrete Multivariate Distributions"
			  Equation 35.152 on page 80, though we denote n_i by N_i and n by SumN
			  and alpha_i by a_i and a_\dot by Sum_a, and k by K (the number of categories)
   
   			  A big note here is that sometimes it will be convenient to allow for categories
   			  that "don't exist."  In other words, if one of the a_i's is zero, then the corresponding
   			  N_i had better be zero, and if it is, then we will want to go ahead and 
   			  calculate the probability of a random vector with only the K-1 components.
   
   AUTHOR:  Eric C. Anderson
   DATE:  11/02/00 Created
***************/
	int i;
	int SumN=0;
	double Sum_a=0.0;
	double Result = 0.0;  /*  this will accumulate the result (a sum of logs of  */
						  /*  the factors in Equation 35.152) */
	
	/*  First we compute SumN and Sum_a, and check to make sure that there aren't any negative */
	/*  numbers there  */
	for(i=0;i<K;i++)  {
		SumN += N[i];
		Sum_a += a[i];
		if(a[i] < 0.0 || N[i] < 0) {  /*  we do not allow negative parameters or variables */
			printf("\n\na[%d] or N[%d] < zero in LogCMDPMF()", i,i);
			printf("Exiting to system...");
			exit(1);
		}
	}
	
	/*  then we compute the log of the pmf */
	/*  first do the coefficient part: */
	Result = LogFactorial(SumN) + LogGamma(Sum_a) - LogGamma(SumN + Sum_a);
	for(i=0;i<K;i++)  {
		if(a[i] > 0.0)  {
			Result +=  LogGamma(N[i] + a[i]) - LogFactorial(N[i]) - LogGamma(a[i]);
		}
		else {  /*  in this case, a[i] == 0.0 */
			if(N[i] > 0) {
				printf("\n\nN[%d] > 0 while a[i] == 0 in LogCompMultProb()", i);
				printf("Exiting to system...");
				exit(1);
			}
		}
	}  /*  closes loop over i */
	
	return(Result);
}


double 
LogPEUrnModelPMF(int *N, int *X, double s, int K)
{
/*************
*  FUNCTION:  
*  PURPOSE:   To compute the log of the probability of a Dirichlet-compound multinomial random
			  vector N, given the parameters in the form of a Polya Urn Model.
			  
			  X is an array of ivals giving the initial number of balls of K different colors in the urn 
			  s is a double which is the number of balls of like color to the one drawn
			  that get stochastically replaced (along with the original drawn ball). 
			  N is an array of ivals of the number of balls drawn of different colors 
			  
			  Hence s == 0.0 is Multinomial, but we can't use this function to compute 
			  the multinomial probability function because the LogGamma function will croak
			  on the thing.
			  
			  We allow for colors to exist amongst the K, even if there aren't any
			  initial balls of that color.  This is done just by calling LogCompMultProb. 
   
*  AUTHOR:  Eric C. Anderson
*  DATE:  11/02/00 Created.  
***************/

	int i;
	double *alpha;  /*  for an array of alphas */
	double result;
	
	alpha = dvector(0,K-1);
	

	for(i=0;i<K;i++)  {
		alpha[i] = (double)X[i]/s;
	}
	
	result = LogCMDPMF(N,alpha,K);
	
	return(result);
}



double 
LogMultinomUrnPMF(int *N, int *X, int K)
{
/*************
*  FUNCTION:  
*  PURPOSE:   To compute the log of the probability of a multinomial random
			  vector N, given the parameters in the form of a Polya Urn Model.
			  (i.e. with the counts in the urn from which one is drawing, and
			  assuming that the stochastic replacement quantity is zero.---so
			  sampling is with replacement)
			  
			  X is an array of ints giving the initial number of balls of K different colors in the urn 
 			  N is an array of ints of the number of balls drawn of different colors 
			  K is the number of categories
			  
			  We allow for colors to exist amongst the K, even if there aren't any
			  initial balls of that color, so long as there aren't any balls of that
			  color drawn!   
   
*  AUTHOR:  Eric C. Anderson
*  DATE:    Created 4/16/01
***************/

	int i;
	double *p;  /*  for an array of proportions */
	double X_sum=0.0;
	int N_sum = 0;
	double result;
	
	p = dvector(0,K-1);
	
	for(i=0;i<K;i++)  {
		X_sum += (double)X[i];
		N_sum += N[i];
	}
	
	for(i=0;i<K;i++)  {
		p[i] = (double)X[i]/X_sum;
	}
	
	/*  start with the factorial of the number of draws: */
	result = LogFactorial(N_sum);
	
	/*  then cycle over the different categories */
	for(i=0;i<K;i++)  {
		if(p[i] > 0.0)  {
			result +=     (double)N[i] 
					   *  log(p[i])
					   -  LogFactorial(N[i]);
		}
		else {  /*  in this case, p[i] == 0.0 */
			if(N[i] > 0) {
				printf("\n\nN[%d] = %d (> 0) while p[%d] == 0.0 in LogMultinomUrnPMF()",
						i,N[i],i);
				printf("Exiting to system...");
				exit(1);
			}
		}
	}  /*  closes loop over i */
	
	
	free_dvector(p,0,K-1);
	
	return(result);
}



double 
LogHGUrnPMF(int *N, int *X, int K)
{
/*************
*  FUNCTION:  LogHGUrnPMF(ival *N, ival *X, int K)
*  PURPOSE:   To compute the log of the probability of a multivariate hypergeometric random
			  vector N, given the parameters in the form of a Polya Urn Model.
			  (i.e. with the counts in the urn from which one is drawing, and
			  assuming that the stochastic replacement quantity is -1.---so
			  sampling is without replacement)
			  
			  X is an array of ints giving the initial number of balls of K different colors in the urn 
 			  N is an array of ints of the number of balls drawn of different colors 
			  K is the number of categories
			    
			  We allow for colors to exist amongst the K, even if there aren't any
			  initial balls of that color, so long as there aren't any balls of that
			  color drawn!   
   
*  AUTHOR:  Eric C. Anderson
*  DATE:    Created 4/16/01
***************/

	int i;
	int X_sum=0;
	int N_sum = 0;
	double result;
	
	
	for(i=0;i<K;i++)  {
		X_sum += X[i];
		N_sum += N[i];
		if(N[i] > X[i])  {
			printf("\n\nN[%d] = %d > X[%d] = %d in LogHGUrnPMF()",
					i,N[i],i,X[i]);
			printf("Exiting to system...");
			exit(1);
		}
	}

	
	/*  start with the normalizing constant: */
	result = LogFactorial(X_sum-N_sum) + LogFactorial(N_sum) - LogFactorial(X_sum);
	
	/*  then cycle over the different categories */
	for(i=0;i<K;i++)  {
		result +=   LogFactorial(X[i])
				  - LogFactorial(X[i] - N[i])
				  - LogFactorial(N[i]);
	}  /*  closes loop over i */
		
	return(result);
}














/*************
*  FUNCTION:  NormalPDF(double mean, double var, double x)
*  PURPOSE:   Return the value of the density function for a given
			  normal random variable with mean and var.  
*
*  INPUT:     double mean -->  the mean of the normal distribution 
			  double var  -->  the variance of the normal distribution
*  			  double x  -->  the value whose density we shall return

	Originally "normal_density"


*************/



double
NormalPDF(double mean, double var, double x)
{
	return(  (1.0/(sqrt(2.0*3.14159*var) ) ) * 
			exp ( (-pow(x-mean,2))/ (2.0*var) ) );
}



/*

	Returns the log of the density of a dirichlet random vector

*/
double
LogDirichletPDF(double *pars, double *rvs, int numComp)
{
	int i;
	double sum=0.0;
	double output=0.0;
	
	for(i=0;i<numComp;i++)  {
		output += ((pars[i]-1.0) * log(rvs[i]) ) - LogGamma(pars[i]);
		sum += pars[i];
	}
	
	return(output + LogGamma(sum));

}

/*************
*  FUNCTION:  double NormalCDF(double mean, double var, double x)
*  PURPOSE:   NORM(al) LOW(er) TAIL (probability)
*			  Compute the probability that X < x when
			  X is a N(mean,var) random variable.
*
*  INPUT:     double mean -->  the mean of the normal distribution to draw from
			  double var  -->  the variance of the normal distribution to draw from
			  double x  -->  the value of little x
			  
	Originally norm_low_tail			  
			  
*************/



double 
NormalCDF(double mean, double var, double x)
{
	double q;
	
	q = (x-mean)/(sqrt(2.0*var));   /*  q is N(0,1/2). */
	
	
	if(q <= 0.0) {
		return(.5 * ErrFunctComp(-q) );
	}
	else {
		return(.5 + .5 *ErrFunct(q) );
	}

}


/*************
*  FUNCTION:  LogisticCDF(double mu, double beta, double x)
*  PURPOSE:   Return the value F(x;mu,beta) of the cdf for a logistic dsn which is 
			  parametrized in terms of its mean and variance.
			  
*	Originally "norm_up_tail"
*************/

double
LogisticCDF(double mu, double var, double x)
{
	double beta;
	
	beta = sqrt(3.0*var)/3.1415926535897932385;
	
	return( 1.0 /  
		( 1.0 + exp( -(x-mu)/beta))    );
}


/*************
*  FUNCTION:  double NormalInvCDF(double mean, double var, double x)
*  PURPOSE:   NORM(al) LOW(er) TAIL (probability)
*			  Compute the probability that X > x when
			  X is a N(mean,var) random variable.
*
*  INPUT:     double mean -->  the mean of the normal distribution to draw from
			  double var  -->  the variance of the normal distribution to draw from
			  double x  -->  the value of little x
			  
	Originally norm_low_tail
*************/



double 
NormalInvCDF(double mean, double var, double x)
{
	double q;
	
	q = (x-mean)/(sqrt(2.0*var));  /*  q is N(0,1/2). */
	
	
	if(q >= 0.0) {
		return(.5 * ErrFunctComp(q) );
	}
	else {
		return(.5 + .5 * ErrFunct(-q) );
	}

}



/******************************************************************************************/
/*                                                                                        */
/*								RANDOM NUMBER GENERATION FUNCTIONS                        */
/*                                                                                        */  
/******************************************************************************************/




/*************
*  FUNCTION:  uniform(int a, int n)
*  PURPOSE:   Generate a random integer uniformly distributed between a and n,
*			  inclusive.  This uses ranlib's ignuin function which does exactly
*			  the same thing, but I wrote this so it accepts ints instead of longs
*
*  INPUT:     int a  -->  low endpoint of the interval (inclusive)
*			  int n  -->  high endpoint of the interval (inclusive)

*	Originally "uniform"
*************/


int 
UniformRV(int a, int n)
{
	return( (int)ignuin( (long)a, (long)n ) );
}



/*************
				
				Draws a discrete uniform between a and n, excluding x

	*  INPUT:     int a  -->  low endpoint of the interval (inclusive)
	*			  int n  -->  high endpoint of the interval (inclusive)
	*  			  int x  -->  the value between a and n to be excluded
	
		Taken from earlier funtion "duwox"
*************/
int 
DuwoxRV(int a, int n, int x)
{

	int r;
	
	if(x > n || x < a) {
		printf("\nWARNING: x is not between a and n in duwox(a,n,x)\n");
	}
	
	r = UniformRV(a,n-1);
	
	if(r >= x) {
		r += 1;
	}
	
	return(r);
}



/*
	return the pair (i,j) where i < j and they are both between Lo and Hi
	
	taken from "UrnOverlap's" draw_ordered_pair
*/
void 
OrderedPairRV(int Lo, int Hi, int *i, int *j)
{


	int a,b;
	
	if(Lo >= Hi) {
		printf("Lo = %d >= Hi = %d, in function draw_ordered_pair()",Lo,Hi);
		printf("\nExiting to system...\n\n");
		exit(1);
	}
	else {
		a = UniformRV(Lo,Hi);
		b = DuwoxRV(Lo,Hi,a);
		
		if(a<b) {
			*i = a;
			*j = b;
		}
		else if(a>b) {
			*i = b;
			*j = a;
		}
		else {  /*  in this case a == b */
			printf("\n\ndraw_ordered_pair working incorrectly.\nExiting to System...\n\n");
			exit(1);
		}
	}

}



/*
	Given a vector of n probabilities summing to one, and each
	one corresponding to an index between 0 and n-1 and an integer
	between  L and L + n - 1, this 
	function returns the integer drawn randomly from the probabilities.
	
	DOES NOT CHECK TO VERIFY THE PROBS SUM TO ONE!!
*/
int IntFromProbsRV(double *Probs, int lo, int n)
{
	int i;
	double cum = 0.0;
	double rando;
	
	rando = (double)ranf();
	
	for(i=0;i<n;i++)  {
		cum += Probs[i];
		if(cum >= rando)
			break;
	}
	
	if(i==lo+n) i--;  /*  if for some reason it ran over the bounds, reign it in. */
	
	return(lo + i);
}





/*

	generates a dirichlet random vector by normalized gammas


*/
void
DirichletRV(double *par, int n, double *out)
{
	int i;
	double normo=0.0;
	
	for(i=0;i<n;i++)  {
	
		/*  have it print an error message if any of the parameters are less than  */
		/*  _or_ equal to zero. */
		if(par[i] <= 0.0)  {
			printf("\n\nDRAT!! par[%d] <= 0.0 in DirichletRV",i);
		}	
		out[i] = (double)gengam(.5f,(float)par[i]);   /*  draw a gamma with proper scale and shape pars */
										/*  use scale parameter .010, because it seems that I was occasionally */
										
		
		/*  check here---if any of them are less than or equal to zero, there is */
		/*  a problem.  gengam returns floats that are sometimes equal to zero.  Very stupid! */
		/*  so, the fix here is going to to be, convert those zeroes into 1.0e-44. That is very */
		/*  small (about the limit of floating point precision on many machines.  So, this is not */
		/*  very elegant, but it should work. */
		if(out[i] <= 0.0) {
			/* fprintf(stderr,"\nRealized a component of Dirichlet Random Vector <= 0.0\n"); */
			out[i] = 1.0e-44;
		}
		
		normo += out[i];
	}

	/*  one more loop to normalize those so they sum to one: */
	for(i=0;i<n;i++)  {
		out[i] /= normo;  
		
		
	}
}













