#include <stdio.h>
#include <stdlib.h>
#include "JN_TreeUtils.h"


int gTreeUtilsDebug;

/* PrintbTreeWithData: Print a tree to outfile but use PrintbNodeData function to control how node data is printed */

void PrintbTreeWithData(FILE * outfile,struct bNode *p,void (*PrintbNodeData)(FILE *,struct bNode)){

  // if a tip
  if(p->lChild==0){
    PrintbNodeData(outfile,*p);

    // if a root
  } else if(p->Anc==0){
    fprintf(outfile,"(");
    PrintbTreeWithData(outfile,p->lChild,PrintbNodeData);
    fprintf(outfile,")");
    //printf("-");
    PrintbNodeData(outfile,*p);

  }else {
    fprintf(outfile,"(");

    PrintbTreeWithData(outfile,p->lChild,PrintbNodeData);

    fprintf(outfile,",");

    PrintbTreeWithData(outfile,p->rChild,PrintbNodeData);

    fprintf(outfile,")");
    PrintbNodeData(outfile,*p);


  }

}

/* PrintbNodeDataSimple: Print a nodes key and branch length */

void PrintbNodeDataSimple(FILE * outfile,struct bNode p){

  if(p.Tip)
    fprintf(outfile,"%d",p.key);
  fprintf(outfile,":%g",p.v);
}




/* ReadbTree:  Read a tree from treeString.  Assumes NodeList is allocated.  Use Read bNodeData to control how data is input.  */


void ReadbTree(char *treeString,struct bNode * NodeList, void (*ReadbNodeData)(char *,struct bNode *)){

  struct bNode* ReadbNodes(char *treeString,struct bNode * NodeList,int *ctr,int *n,void (*ReadbNodeData)(char *,struct bNode *));

  int n=0;
  int ctr=0;
  struct bNode *BaseNode;


  BaseNode=ReadbNodes(treeString,NodeList,&ctr,&n,ReadbNodeData);
  (void)BaseNode;

  LinkAncestors(&NodeList[0],0);

  SortNodeList(NodeList+1,n-1);
  if(gTreeUtilsDebug>0)
    PrintbTreeWithData(stdout,&NodeList[0],PrintbNodeDataSimple);

  return;
}



/* ReadbNodes: Read the nodes in a tree.  Used by ReadbTree.*/

struct bNode* ReadbNodes(char *treeString,struct bNode * NodeList,int *ctr,int *n,void(*ReadbNodeData)(char *,struct bNode *)){

  struct bNode *nodePtr;

  nodePtr=&NodeList[*n];

  (*n)++;


  switch(treeString[*ctr]){

    case'(':
      (*ctr)++;
    nodePtr->lChild=ReadbNodes(treeString,NodeList,ctr,n,ReadbNodeData);
      if(gTreeUtilsDebug>5)
      printf("set l node %d %d\n",nodePtr->key,nodePtr->lChild->key);
    break;


  case ')':
    fprintf(stderr,"Error in tree at char %d\n",*ctr);
    exit(1);
    break;


    // character is a number
  default:{
    nodePtr->Tip=1;
    ReadbNodeData(treeString+*ctr,nodePtr);
    if(gTreeUtilsDebug>5)
      printf("Reading tip data: %d %g\n",nodePtr->key,nodePtr->v);
    nodePtr->t=0;

    return(nodePtr);
  }
  }

  // Advance to ',' or ')'
  while(treeString[*ctr]!=',' && treeString[*ctr]!=')')
    (*ctr)++;

  if(treeString[*ctr]==','){
    (*ctr)++;
    nodePtr->rChild=ReadbNodes(treeString,NodeList,ctr,n,ReadbNodeData);
    if(gTreeUtilsDebug>5)
      printf("set r node %d %d\n",nodePtr->key,nodePtr->rChild->key);
  }

  // Advance to ')'
  while(treeString[*ctr]!=')')
    (*ctr)++;

  (*ctr)++;
  nodePtr->Tip=0;
  ReadbNodeData(treeString+*ctr,nodePtr);
  nodePtr->t=nodePtr->lChild->t+nodePtr->lChild->v;
  if(gTreeUtilsDebug>5)
    printf("Reading interior data: %d %g\n",nodePtr->key,nodePtr->v);


  return(nodePtr);
}


/* ReadbNodeDataSimple: Read's key and then branch length.  Assumes tips have keys < n*/
void ReadbNodeDataSimple(char * S, struct bNode *p){

  if(p->Tip)
    sscanf(S,"%d:%lf",&p->key,&p->v);
  else
    sscanf(S,":%lf",&p->v);

}

void SortNodeList(struct bNode *NodeList,int ListLength){

  int i,j;
  int CompareNodeTimes(void const *a, void const *b);
  struct bNode temp;
  struct bNode *tempAnc1;
  int tempAnc1Child=0;
  int tempAnc2Child=0;
  struct bNode *tempAnc2;
  int MRCA;

  MRCA=0;
  for(i=0;i<ListLength-1;i++){
    for(j=ListLength-1;j>i;j--){
      if(NodeList[j].t<NodeList[j-1].t){
        /*printf("\nSwap %d %d key1: %d key2: %d t1 %g t2 %g\n",j,j-1,NodeList[j].key,NodeList[j-1].key,NodeList[j].t,NodeList[j-1].t);
        for(k=0;k<ListLength;k++){
          printf("k: %d key: %d t: %g s: %d l: %d r: %d Anc: %d\n",k, NodeList[k].key,NodeList[k].t,&NodeList[k],NodeList[k].lChild,NodeList[k].rChild,NodeList[k].Anc);
        }*/

        if(j-1==MRCA)
          MRCA++;
               /* Swap, but also be careful to swap ancestors pointers */
        tempAnc1=NodeList[j].Anc;
        if(tempAnc1->lChild==&NodeList[j])
          tempAnc1Child=1;
        else if(tempAnc1->rChild==&NodeList[j])
          tempAnc1Child=2;

        tempAnc2=NodeList[j-1].Anc;
        if(tempAnc2->lChild==&NodeList[j-1])
          tempAnc2Child=1;
        if(tempAnc2->rChild==&NodeList[j-1])
          tempAnc2Child=2;

        /* Note: if Ancestor of one node, is the other node then... */
        if(tempAnc1==&NodeList[j-1])
          tempAnc1=&NodeList[j];
        if(tempAnc2==&NodeList[j])
          tempAnc2=&NodeList[j-1];

        /* Swap nodes in Nodelist */
        temp=NodeList[j];
        NodeList[j]=NodeList[j-1];
        NodeList[j-1]=temp;

        if(tempAnc1Child==1)
          tempAnc1->lChild=&NodeList[j-1];
        else if(tempAnc1Child==2)
          tempAnc1->rChild=&NodeList[j-1];
        if(tempAnc2Child==1)
          tempAnc2->lChild=&NodeList[j];
        else if(tempAnc2Child==2)
          tempAnc2->rChild=&NodeList[j];

        /* Update ancestors pointers */
        if(NodeList[j].lChild)
          NodeList[j].lChild->Anc=&NodeList[j];
        if(NodeList[j].rChild)
          NodeList[j].rChild->Anc=&NodeList[j];
        if(NodeList[j-1].lChild)
          NodeList[j-1].lChild->Anc=&NodeList[j-1];
        if(NodeList[j-1].rChild)
          NodeList[j-1].rChild->Anc=&NodeList[j-1];


        //PrintbTreeWithData(stdout,&NodeList[MRCA],PrintbNodeDataSimple);

      }
    }
  }
  // Sort NodeList by Node times
  //qsort(NodeList,ListLength,sizeof(struct bNode),CompareNodeTimes);
  // Assign keys in order
  for(i=0;i<ListLength;i++){
    NodeList[i].key=i+1;
  }


};


int CompareNodeTimes(void const *a, void const *b){

  double temp = ((struct bNode *)a)->t-((struct bNode*)b)->t;
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;

  return 0;
}



void LinkAncestors(struct bNode * D, struct bNode * A){

  D->Anc=A;
  if(D->lChild)
    LinkAncestors(D->lChild,D);
  if(D->rChild)
    LinkAncestors(D->rChild,D);

}


void GenerateLevelsList(int n, struct bNode * NL, struct bNode ** LL){
  int i,j;
  int L;
  int CoalNodeL;
  // Init LevelList at L=n
  for(i=0;i<n;i++)
    LL[n-1][i]=NL[i+1];


  // Build recursively
  for(L=n-1;L>=2;L--){
    CoalNodeL=2*n-L;



    // Set level list accordingly
    for(i=0,j=0;i<=L;i++){
      // If not part of coalescent at L add to LevelList
      if((&LL[L][i]!=NL[CoalNodeL].lChild)&&(&LL[L][i]!=NL[CoalNodeL].rChild)){
	LL[L-1][j]=LL[L][i];
	j++;
      }
    };
    LL[L-1][j]=NL[CoalNodeL];
  }

}


/*
int AssembleDescendantList(struct bNode *p){


  int DesKey;

 // if a tip
  if(p->lChild==0){
    CalcVarCovarMatrix(*p);
    p->nlDes=0;
    p->nrDes=0;
    return(p->key)
    // if a root
  } else if(p->Anc==0){
    CalcVarCovarMatrix(p->lChild,D,M);
    return(-1);
  }else {

    DesKey=CalcVarCovarMatrix(p->lChild);

    if(DesKey>1) {// Then lChild is a tip
      p->nlDes=1;
      p->lDesList[0]=DesKey;
    }else{  // lChild is interior node with Descendants to copy
      p->nlDes=p->lChild->nlDes+p->lChild->nrDes;
      for(i=0;i<p->lChild->nlDes)
	p->lDesList[i]=p->lChild->lDesList[i];
      for(i=p->lChild->nlDes;i<p->nlDes;i++)
	p->lDesList[i]=p->lChild->rDesList[i];
    }

   DesKey=CalcVarCovarMatrix(p->rChild);

    if(DesKey>1) {// Then rChild is a tip
      p->nrDes=1;
      p->rDesList[0]=DesKey;
    }else{  // rChild is interior node with Descendants to copy
      p->nrDes=p->rChild->nlDes+p->rChild->nrDes;
      for(i=0;i<p->rChild->nlDes)
	p->rDesList[i]=p->rChild->lDesList[i];
      for(i=p->rChild->nrDes;i<p->nrDes;i++)
	p->rDesList[i]=p->rChild->rDesList[i];
    }
    return(-1);
  }

}

*/
