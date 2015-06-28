
#ifndef __JN_TREE_UTILS_FLAG__
#define __JN_TREE_UTILS_FLAG__

#ifdef UN_EXTERN
	#define GLOB
#else 
	#define GLOB extern
#endif

#define MAX_NODE_LABEL_LENGTH 5


typedef struct bNode{
  struct bNode *lChild;
  struct bNode *rChild;
  struct bNode *Anc;

  int Tip;   /* Boolean 1 if tip 0 if interior */
  char Label[MAX_NODE_LABEL_LENGTH];
  int key;   /* Generally set so 1-n are tips, interior nodes are n+1..., root=0 */
  double v;  /* branch length below node*/
  double t;  /* absolute time of node */
  struct bNodeData * Data;

  

} bNode;




/* Output functions*/
extern void PrintbTreeWithData(FILE *outfile,struct bNode *p,void (*PrintbNodeData)(FILE *,struct bNode));
extern void PrintbNodeDataSimple(FILE * outfile,struct bNode);

/* Input functions */
extern void ReadbTree(char *treeString,struct bNode * NodeList, void (*ReadbNodeData)(char *,struct bNode *));
extern void ReadbNodeDataSimple(char * S, struct bNode *p);

/* Tree management */
extern void SortNodeList(struct bNode *NodeList,int ListLength);
extern void LinkAncestors(struct bNode * D, struct bNode * A);
extern void GenerateLevelsList(int n, struct bNode * NL, struct bNode ** LL);

/* Wish-list */
/* GenerateRandomTopo */
/* */




#endif
