
//Now, my goal is to write a program that will read in a set of trees,
//drop mutations, take the subset of those that are in low-frequency,
//and produce the appropriate IS input file from the resultant mutant
//subtree.  I'll call this program MakeMutants.c.  The input will be a
//filename for the output of simCoalLatticeAll, the lower bound on the
//frequency and upperbound on the frequency for loci to output to the
//IS, and an output filename.


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "ranlib.h"
#include "MathStatRand.h"
#include "ECA_utilities.h"
#include "ECA_MemAlloc.h"


typedef struct node{
  struct node *left;  // pointer for tree structure
  struct node *right; // pointer for tree structure
  struct node *next;  // pointer for list structure
  int key;            // identifier
  int x0;             // init x position
  int y0;             // init y position
  int x;              // x position
  int y;              // y position
  int vecPos;         // x-y position as a single number
  int tB;             // time of node's birth (looking back in time)
  int tD;             // time of node's death 
  int nDescendants;   // number of descendants to the node
  
} node;

typedef struct Params{
  // To be read from commandline
  long seed1,seed2;
  double lowerBoundFreq;
  double upperBoundFreq;;
  double lowerGridLim;
  double upperGridLim;
  double delGrid;
  char infilename[80];
  char outfilebase[80];
  int verb;

  // To be read from input file
  int dim;
  double vx;
  double vy;
  double f;
  int nloci;
  int nlowFreqLoci;

  // To be used in program
  FILE* infile;
  FILE* outfile;
  struct node *NodeList;
  int nTips;
  int nNodes;

} Params;

int gECA_mall_call_SetFlag;
double gECA_mall_call_bytes;

int main(int argc, char **argv){

  int i;
  struct Params P;
  struct Params GetParams(int argc, char **argv);
  void InitOutFiles(struct Params *P);
  void TermOutFiles(struct Params *P);

  
  struct node* ReadTree(struct Params *P);
  struct node *SampleRoot;
  struct node *MutantTreeRoot;
  

  struct node* DropMutation(struct node *root,struct Params *P);
  int Age;
  void PrintTipData(struct node *p,struct Params *P);

  // Input filename, upper frequency, lower frequency, output filename
  // Open input file
  // Read parameters, dim, vx,vy,f,nloci

  P=GetParams(argc,argv);

  // Open output files and write headers
  InitOutFiles(&P);


  // Begin loop over nloci
  for(i=0;i<P.nloci;i++) {

    // Read first tree from input file
    SampleRoot=ReadTree(&P);


    // Drop mutation
    MutantTreeRoot=DropMutation(SampleRoot,&P);
    // if derived allele fits in frequency range, 
    if((double)MutantTreeRoot->nDescendants>=P.lowerBoundFreq*P.nTips && (double)MutantTreeRoot->nDescendants<=P.upperBoundFreq*P.nTips){

      //output tree data to output file
      fprintf(P.outfile,"%d ",MutantTreeRoot->nDescendants);
      PrintTipData(MutantTreeRoot,&P);
      fprintf(P.outfile,"\n");
      if(P.verb>0)
	printf("\n");
      P.nlowFreqLoci++;

      // Get an age of the mutation
      // Note: use assumption that if a single event occurs among
      // a series of binomial trials (discrete generations)
      // the exact trial on which it occurs is uniformly 
      // distributed

      Age=UniformRV(MutantTreeRoot->tB,MutantTreeRoot->tD);      
      // Output Age and Frequency to stdout so it can be greped out
      printf("Age: %d Freq: %f\n",Age,MutantTreeRoot->nDescendants/(double)P.nTips);

    }
    

  }

  // Close files
  TermOutFiles(&P);

  exit(0);

}

struct Params GetParams(int argc, char **argv){

  // Input filename, upper frequency, lower frequency, output filename
  // Open input file
  // Read parameters, dim, vx,vy,f,nloci
  void usage();
  int SeedsSpec, BoundsSpec,InFileSpec,OutFileSpec;
  int i;
  struct Params P;
  SeedsSpec=BoundsSpec=InFileSpec=OutFileSpec=0;
  i=1;
  P.verb=0;
  P.upperGridLim=-1.0;
  P.lowerGridLim=-1.0;
  while(i<argc) {
    switch(argv[i][1]){
    case 'S':
      i++;
      P.seed1=(long)atoi(argv[i]);
      P.seed2=(long)atoi(argv[i+1]);
      i++;
      SeedsSpec=1;
      break;
    case 'o':
      i++;
      sprintf(P.outfilebase,argv[i]);
      OutFileSpec=1;
      break;
    case 'i':
      i++;
      sprintf(P.infilename,argv[i]);
      InFileSpec=1;
      break;

    case 'A':
      i++;
      P.lowerGridLim=atof(argv[i]);
      i++;
      P.upperGridLim=atof(argv[i]);
      i++;
      P.delGrid=atof(argv[i]);
      break;

    case 'B':
      i++;
      P.lowerBoundFreq=atof(argv[i]);
      i++;
      P.upperBoundFreq=atof(argv[i]);
      BoundsSpec=1;
      break;
          
    case 'V':
      i++;
      P.verb=atoi(argv[i]);
      break;
    default:
      fprintf(stderr,"Unrecognized option: %s\n",argv[i]);
      usage();
      exit(1);
    }
  
    i++;

  }


  if(!SeedsSpec)
    fprintf(stderr,"Seeds not specified...\n");
  if(!InFileSpec)
    fprintf(stderr,"Input file not specified...\n");
  if(!OutFileSpec)
    fprintf(stderr,"Output file not specified...\n");
  if(!BoundsSpec)
    fprintf(stderr,"Mutant allele frequency bounds not specified...\n");
  if(!SeedsSpec || !InFileSpec || !OutFileSpec || !BoundsSpec){
    usage();
    exit(1);
  }
  
  
  setall(P.seed1,P.seed2);  
  
  // Open input file
  P.infile=fopen(P.infilename,"r");
  if(!P.infile){
    fprintf(stderr,"Problem opening input file %s\n",P.infilename);
    exit(1);
  }

  // Read parameters, dim, vx,vy,f,nloci
  fscanf(P.infile,"dim: %d\n",&P.dim);
  fscanf(P.infile,"Vx: %lf\n",&P.vx);
  if(P.dim==2)
    fscanf(P.infile,"Vy: %lf\n",&P.vy);
  else 
    P.vy=0;
  fscanf(P.infile,"f: %lf\n",&P.f);
  fscanf(P.infile,"nloci: %d\n",&P.nloci);

  if(P.f<0 || P.f>1){
    fprintf(stderr,"Error in input f value: %f",P.f);
    exit(1);
  }
  if(P.nloci<0) {
    fprintf(stderr,"Error in input nloci: %d",P.nloci);
    exit(1);
  }
  if(P.vx<0) {
    fprintf(stderr,"Error in input vx: %f",P.vx);
    exit(1);
  }
  if(P.vy<0) {
    fprintf(stderr,"Error in input vy: %f",P.vy);
    exit(1);
  }

  P.nlowFreqLoci=0;
  
  return(P);

};

void usage(){
  fprintf(stderr,"usage: MakeMutants <Options> \n");
  fprintf(stderr,"\tOptions: \n");
  fprintf(stderr,"\t\t-S <integer1> <integer2>: Random number seeds\n");
  fprintf(stderr,"\t\t-B <float1> <float2>: lower/upper bounds on mut. allele freq.\n");
  fprintf(stderr,"\t\t-A <float1> <float2> <float3>: lower/upper bounds and increment for MLE grid in IS input file.\n");
  fprintf(stderr,"\t\t-i <filename>: Input filename\n");
  fprintf(stderr,"\t\t-o <filneame>: Output base filename\n");
  fprintf(stderr,"\t\t-V <value>: Verbosity\n");
  
}


void InitOutFiles(struct Params *P){
  // Open output files and write headers
  // Open Data output file
  char filename[100];
  strcpy(filename,P->outfilebase);
  strcat(filename,".xy");
  P->outfile=fopen(filename,"w");
  if(!P->outfile){
    fprintf(stderr,"Problem opening output file %s.\n",filename);
    exit(1);
  }

  

};

void TermOutFiles(struct Params *P){

  char filename[100];



  // Get seeds for output
  getsd(&P->seed1,&P->seed2);
  strcpy(filename,P->outfilebase);
  strcat(filename,".seeds");
  P->outfile=fopen(filename,"w");
  fprintf(P->outfile,"%ld %ld\n",P->seed1,P->seed2);
  fclose(P->outfile);


  // Prepare input file for IS 
  strcpy(filename,P->outfilebase); 
  strcat(filename,".in"); 
  P->outfile=fopen(filename,"w");   
  fprintf(P->outfile,"initseed %ld\n",P->seed1); 
  fprintf(P->outfile,"f %f\n",P->f); 

  // if grid lims not specified
  if(P->lowerGridLim<0){
    if(P->dim==2){ 
      fprintf(P->outfile,"vxlow %lf\n",P->vx-.3*P->vx); 
      fprintf(P->outfile,"vxhigh %lf\n",P->vx+.3*P->vx); 
      fprintf(P->outfile,"delvx %lf\n",P->vx*.15);     
      fprintf(P->outfile,"vylow %lf\n",P->vy-.3*P->vy); 
      fprintf(P->outfile,"vyhigh %lf\n",P->vy+.3*P->vy); 
      fprintf(P->outfile,"delvy %lf\n",P->vy*.15); 
    }else{ 
      fprintf(P->outfile,"vxlow %lf\n",P->vx-.25*P->vx); 
      fprintf(P->outfile,"vxhigh %lf\n",P->vx+.25*P->vx); 
    fprintf(P->outfile,"delvx %lf\n",P->vx*.03125);     
    } 
  }else{
    if(P->dim==2){ 
      fprintf(P->outfile,"vxlow %lf\n",P->lowerGridLim); 
      fprintf(P->outfile,"vxhigh %lf\n",P->upperGridLim); 
      fprintf(P->outfile,"delvx %lf\n",P->delGrid);     
      fprintf(P->outfile,"vylow %lf\n",P->lowerGridLim); 
      fprintf(P->outfile,"vyhigh %lf\n",P->upperGridLim); 
      fprintf(P->outfile,"delvy %lf\n",P->delGrid); 
    }else{ 
      fprintf(P->outfile,"vxlow %lf\n",P->lowerGridLim); 
      fprintf(P->outfile,"vxhigh %lf\n",P->upperGridLim); 
    fprintf(P->outfile,"delvx %lf\n",P->delGrid);     
    } 



  }
  fprintf(P->outfile,"maxcase %d\n",P->nlowFreqLoci); 
  fprintf(P->outfile,"end\n"); 
  fclose(P->outfile); 
    
  
}


struct node * ReadTree(struct Params *P){

  struct node * ReadNodes(char *treeString,int *ctr,int *n,struct Params * P);
  void PrintData(struct node *p,int dim);

  char *treeString;
  struct node * root;
  int ctr1,ctr2;


  // Read P.nTips
  fscanf(P->infile,"%*d%*d%d",&P->nTips);
  P->nNodes=2*P->nTips-1;

  // Declare memory for NodeList and treeString
  P->NodeList=(struct node *)calloc((size_t)(P->nNodes),sizeof(struct node) );
  treeString=(char *)calloc((size_t)P->nNodes*40,sizeof(char));
  if(!P->NodeList || !treeString){
    fprintf(stderr,"Error declaring memory for tree.  Too many nodes?\n");
    exit(1);
  }
  
  fscanf(P->infile,"%s",treeString);
  ctr1=0;
  ctr2=0;
  root=ReadNodes(treeString,&ctr1,&ctr2,P);
  if(P->verb>0)
    PrintData(root,P->dim);
  return(root);

};

struct node* ReadNodes(char *treeString,int *ctr,int *n,struct Params *P){

  struct node *nodePtr;
  int length;

  nodePtr=&P->NodeList[*n];

  (*n)++;

  //  if(treeString[*ctr]==','){
  //  if(P->verb>5) printf(", skip\n");
  //  (*ctr)++;
  // }

  switch(treeString[*ctr]){
    
    case'(':
      (*ctr)++;
    nodePtr->left=ReadNodes(treeString,ctr,n,P);
    if(P->verb>5)
      printf("set l node %d %d\n",nodePtr->key,nodePtr->left->key);
    break;
    
   

  case ')': 
    fprintf(stderr,"Error in tree at char %d\n",*ctr);
    exit(1);
    break;
    
 
    // character is a number
  default:{
    if(P->dim==2)
      sscanf(treeString+*ctr,"%d:%d:%d:%d",&nodePtr->key,&nodePtr->x0,&nodePtr->y0,&length);
    else if(P->dim==1)
      sscanf(treeString+*ctr,"%d:%d:%d",&nodePtr->key,&nodePtr->x0,&length);
    else{
      fprintf(stderr,"Dimension must equal 1 or 2.\n");
      exit(1);
    }
    if(P->verb>5){
      printf("Reading data: %d %d %d %d\n",nodePtr->key,nodePtr->x0,nodePtr->y0,length);
    }
    if(nodePtr->key<P->nTips){
      nodePtr->tB=0;
      nodePtr->nDescendants=1;
    } else {
      nodePtr->tB=nodePtr->left->tD;    
      nodePtr->nDescendants+=nodePtr->left->nDescendants;
      nodePtr->nDescendants+=nodePtr->right->nDescendants;
    }
    nodePtr->tD=nodePtr->tB+length;
    return(nodePtr);
  }
  }
    
  // Advance to ',' or ')'
  while(treeString[*ctr]!=',' && treeString[*ctr]!=')')
    (*ctr)++;

  if(treeString[*ctr]==','){
    (*ctr)++;
    nodePtr->right=ReadNodes(treeString,ctr,n,P);
    if(P->verb>5)
      printf("set r node %d %d\n",nodePtr->key,nodePtr->right->key);
  } 

  // Advance to ')'
  while(treeString[*ctr]!=')')
    (*ctr)++;

  (*ctr)++;
  if(P->dim==2)
      sscanf(treeString+*ctr,"%d:%d:%d:%d",&nodePtr->key,&nodePtr->x0,&nodePtr->y0,&length);
  else if(P->dim==1)
    sscanf(treeString+*ctr,"%d:%d:%d",&nodePtr->key,&nodePtr->x0,&length);
  else{
    fprintf(stderr,"Dimension must equal 1 or 2.\n");
    exit(1);
  }
  if(P->verb>5){
    printf("Reading data: %d %d %d %d\n",nodePtr->key,nodePtr->x0,nodePtr->y0,length);
  }
  if(nodePtr->key<P->nTips){
    nodePtr->tB=0;
    nodePtr->nDescendants=1;
  }
  else {
    nodePtr->tB=nodePtr->left->tD;   
    nodePtr->nDescendants+=nodePtr->left->nDescendants;
    nodePtr->nDescendants+=nodePtr->right->nDescendants;
/*     if(nodePtr->left->key<P->nTips) */
/*       nodePtr->nDescendants+=nodePtr->left->nDescendants+1; */
/*     if(nodePtr->right->key<P->nTips) */
/*       nodePtr->nDescendants+=nodePtr->right->nDescendants+1; */
  }
  nodePtr->tD=nodePtr->tB+length;
    
  
  return(nodePtr);
} 



struct node* DropMutation(struct node* root, struct Params * P){

  // For implementing mutation
  double *MutBranchProbs;
  double r;
  int i,j;
  struct node * MutantRoot;

  MutBranchProbs=(double *)calloc((size_t)(P->nNodes),sizeof(double));
  if(!MutBranchProbs){
    fprintf(stderr,"Error declaring memory.  Too large of a sample?");
    exit(1);
  }

  // First find sum of branch lengths, i.e. normalization constant
  r=0;
  for(i=0;i<P->nNodes;i++)
    r+=P->NodeList[i].tD-P->NodeList[i].tB;

  // Set up vector of probabilities
  for(i=0;i<P->nNodes;i++){
    MutBranchProbs[i]=(double)(P->NodeList[i].tD-P->NodeList[i].tB)/r;
  }
    
  // Choose branch
  j=IntFromProbsRV(MutBranchProbs,0,P->nNodes);
  
  if(P->verb>1){
    printf("\nPotential locations for mutation to occur:\n");
    for(i=0;i<P->nNodes;i++){
      printf("%d %d %f %d %f\n",P->NodeList[i].key, P->NodeList[i].nDescendants,MutBranchProbs[i],P->NodeList[i].tD-P->NodeList[i].tB,(P->NodeList[i].tD+P->NodeList[i].tB)/2.0);
    }
    printf("Chosen branch: %d\n",P->NodeList[j].key);
  }
    
  MutantRoot=&P->NodeList[j];    
    
  return(MutantRoot);
}

void PrintTipData(struct node *p,struct Params *P){
  if(p->left)
    PrintTipData(p->left,P);
  if(p->right)
    PrintTipData(p->right,P);
  if(!p->right && !p->left){
    if(P->dim==2){
      if(P->verb>1)
	printf("%d %d ",p->x0, p->y0);
      fprintf(P->outfile,"%d %d ",p->x0, p->y0);}
    else{
      if(P->verb>1)
	printf("%d ",p->x0);
      fprintf(P->outfile,"%d ",p->x0);
    }
  }
}




void PrintData(struct node *p,int dim){
  
  // if a tip
  if(p->tB==0){

    if(dim==1){
      printf("%d:%d:%d",p->key, p->x0,p->tD-p->tB);
    }else{
      printf("%d:%d:%d:%d",p->key, p->x0,p->y0,p->tD-p->tB);
    }
    
  } else{
    printf("(");

    PrintData(p->left,dim);

    printf(",");

    PrintData(p->right,dim);

    if(dim==1){
      printf(")");
      printf("%d:%d:%d",p->key, p->x0,p->tD-p->tB);

    }else{
      printf(")");
      printf("%d:%d:%d:%d",p->key, p->x0,p->y0,p->tD-p->tB);
    }
  }
  
  
}
