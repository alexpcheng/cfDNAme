#ifndef METSEG_H
#define METSEG_H
/*
 *
 *	metseg.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 09/10/2014 01:46:32 PM CEST  
 *
 */


typedef struct segm{

  double prob;
  double test;
  double methA;
  double methB;
  double **value; //cpgs
  int *pos; //index
  int n;
  int start;
  int stop;
  double meandiff;
  int length;
  char *chr;
  char init;
  struct segm *next;
  struct segm *parent;
} segment_t;

typedef struct{
  char *chr;
  int start;
  int stop;
  double p;
  double mwu;
  int n;
  double meandiff;
  int length;
  double methA;
  double methB;
} segment_out;

typedef struct{
  segment_out *segment_out;
  int n;
  int i;
  int numberTests;
} list_out;


typedef struct{
  segment_t *seg;
  int n;
  char *chr;
  
  int firststop;
  char *nextchr;
  int nextstart;
  segment_t *head;
  segment_t *tail;
} segmentset_t;

typedef struct{
  int start;
  int stop;
  char *chr;
  double *groupA;
  double *groupB;
  int noA;
  int noB;
  double methA;
  double methB;
} cpg_t;

typedef struct{
  int maxdist;
  int mincpgs;
  int threads;
  int mode;
  char *nameA;
  char *nameB;
  double trend;
  double minFactor; // facor for setting min number of present values
  double minMethDist;
  int minNoA; //min number of present values to fill missing numbers in group A, below discard input line
  int minNoB;
  double valley;
  //only used for threaded segmentation
  char **chr;
  int *pos;
  double **value;
  int n;
  int *grpA;
  int noA;
  int *grpB;
  int noB;
  double ***MWU;
  segment_t *seg; 
  cpg_t *cpg;
  segment_t *List;
  int nList;
  list_out *outputList;
  
  int threadno;

} metseg_t;

void initSegment(segment_t *seg);


typedef struct{
  int a;
  int b;
  int ab1;
  int ab2;
  int child;
  segment_t *max;
  
  double ks11;
  double ks12;
  double ks13;
  
  double ks21;
  double ks22;
  double ks23;
  
  double ks31;
  double ks32;
  double ks33;

  double KS1;
  double KS2;
  double KS3;

} segment_p_t;


#endif
