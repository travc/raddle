
#ifndef __JN_MCTYPES_FLAG__
#define __JN_MCTYPES_FLAG__

#include<math.h>
#include<stdio.h>
/* Structure and associated functions that are convenient for use in Monte Carlo simulations  */
/* The basic principle is that any random variable is declared as a structure that carries    */
/* summary statistics and a histogram of it's distribution                                    */

/* Please NOTE: This code bears EXTREME debt to Eric Anderson's MCtypes code.  Much of this   */
/* is directly adapted from his clever and inspiring work!                                    */


typedef struct DiscreteHisto{
  int nBins;
  int Min;
  int Max;
  int Step;
  int Remainder;  /* Used if (Max-Min)/Step is not an integer value */ 
  int *Counts;
  int nObs;

} DiscreteHisto;

typedef struct ContHisto{
  int nBins;
  double Min;
  double Max;
  double Step;
  int Remainder;  /* Used if (Max-Min)/Step is not an integer value */ 

  int *Counts;
  int nObs;

} ContHisto;


typedef struct intval {

  int v;      /* Current value */ 
  double Mean;
  double SS;  /* Sum of squares */
  double Var;
  int nObs;
  struct DiscreteHisto * H;
  
} intval;



typedef struct doubval{

  double v;      /* Current value */ 
  double Mean;
  double SS;  /* Sum of squares */
  double Var;
  int nObs;
  struct ContHisto * H;


} doubval;

typedef struct DoubvalCovar{
  int Num;
  struct doubval ***Square;

} DoubvalCovar;


/* MEM ALLOC FUNCTIONS */
extern struct intval *AllocIntval(int Min, int Max,int Step);
extern struct doubval *AllocDoubval(double Min, double Max,double Step);
extern struct DiscreteHisto *AllocDiscreteHisto(int Min,int Max,int Step);
extern struct ContHisto *AllocContHisto(double Min,double Max,double Step);
extern struct DoubvalCovar *AllocDoubvalCovar(int n);

/* INITIALIZATION FUNCTIONS */
extern void InitIntvalToZero(struct intval *t);
extern void InitIntvalSummaryToZero(struct intval *t);
extern void InitDoubvalToZero(struct doubval *t);
extern void InitDoubvalSummaryToZero(struct doubval *t);
extern void InitDiscreteHisto(struct DiscreteHisto * H);
extern void InitContHisto(struct ContHisto * H);
extern void InitDoubvalCovar(struct DoubvalCovar *D);

/* INCREMENTING FUNCTIONS */
extern void IncrementIntval(struct intval *t);
extern void IncrementDoubval(struct doubval *t);
extern void IncrementDiscreteHisto(struct DiscreteHisto *D,int v);
extern void IncrementContHisto(struct ContHisto *C, double v);
extern void IncrementDoubvalCovar(struct DoubvalCovar *D, struct doubval **Values);
/* For IncrementDvalCovar the CODE DIFFERS FROM ERIC'S ON ONE LINE BY A FACTOR OF N-1.0 */

/* OUTPUT FUNCTIONS */
extern void PrintDiscreteHistKey(FILE *outfile,struct DiscreteHisto *H);
extern void PrintContHistKey(FILE *outfile,struct ContHisto * H);
extern void PrintDiscreteHistKeySimple(FILE *outfile,struct DiscreteHisto *H);
extern void PrintContHistKeySimple(FILE *outfile,struct ContHisto * H);

extern void PrintDiscreteHistCount(FILE *outfile,struct DiscreteHisto *H);
extern void PrintDiscreteHistProp(FILE *outfile,struct DiscreteHisto *H);
extern void PrintContHistCount(FILE *outfile,struct ContHisto *H);
extern void PrintContHistProp(FILE *outfile,struct ContHisto *H);


#endif 
