
#include "JN_RVutils.h"
void Calc2Moments(double * x, int l, double * m){

  int i;
  double sum=0.0;
  
  for(i=0;i<l;i++){
    sum+=x[i];    
  }
  // Mean
  m[0]=sum/(double)l;
  
  sum=0.0;
  for(i=0;i<l;i++){
    sum+=(x[i]-m[0])*(x[i]-m[0]);
  }
  // Varaince
  m[1]=sum/(double)(l-1);
  

}

double
LogOfNormalPDF(double mean, double var, double x)
{
  return(  log((1.0/(sqrt(2.0*3.14159*var) ) )) + 
	   ( (-pow(x-mean,2))/ (2.0*var) ) );
}

double LogOfBetaPDF(double alpha,double beta, double x){

  return ( (alpha-1.0)*log(x)+(beta-1.0)*log(1.0-x) );

}

double ReflectTruncatedNormalRV(double mean, double var, double a, double b){

  // Note: Uses reflecting boundaries
  double x;

  // if a and b 
  x=gennor((float)mean,(float)sqrt(var));
  while((x<a)||(x>b)){
    if(x<a)
      x=a+(a-x);
    if(x>b)
      x=b-(x-b);
  }
 
  return(x);
}
  
  
  
double LogOfInvGammaPDF(double theta,double beta,double alpha){
  double x;
  x=-alpha*log(beta)-(alpha+1.0)*log(theta)-beta/theta-LogGamma(alpha);
  return(x);
}
/* double TruncatedNormalRV(double mean, double var, double a, double b){ */


/*   double Fa,Fb; */
/*   double x; */
/*   double n=0; */

/*   // If rejection scheme will be efficient */

/*   if(MIN(mean-b,a-mean)/sqrt(var)>1){ */
/*     // rejection scheme */
/*     do{ */
/*       x=gennor((float)mean,(float)sqrt(var)); */
/*       n++; */
/*     }while(((x<a)||(x>b))&&(n<10)); */
/*     if(n<10) */
/*       return(x); */
/*   } */

/*   // Otherwise: inverse cdf method */
/*   Fa=NormalCDF(mean,var,a); */
/*   Fb=NormalCDF(mean,var,b); */
/*   return(sqrt(var)*ltqnorm((double)genunf((float)Fa,(float)Fb))+mean); */
/* } */

