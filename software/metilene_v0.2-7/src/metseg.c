/*
 *  metseg.c
 *  
 *
 *  @author Frank Juehling and Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 09/09/2014 08:54:52 AM CEST
 *  
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>
#include "basic-types.h"
#include "manopt.h"
#include "info.h"
#include "stringutils.h"
#include "fileio.h"
#include "metseg.h"
#include <string.h>
#include "mathematics.h"
#include "vstack.h"
#include "segmentstack.h"
#include <time.h>
#include <ctype.h>
#include <float.h>

#define MAXN 13
#define MAXM 13

char *version = "0.2-7";
unsigned char mute=0;
pthread_mutex_t updatemtx;
double get_ratio(double *a, int m, double *b, int n);
/*------------------------------- initSegment --------------------------------
 *    
 * @brief initialize new segment
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void
initSegment(segment_t *seg) {
  seg->n =0;
  seg->prob=-1;
  seg->test=-1;
  seg->length = 0;
//  seg->chr="";
  seg->init=0;
  seg->value = NULL;
  seg->pos = NULL;
  seg->next = NULL;
  seg->parent = NULL;
}

void
destructSegment(segment_t *seg) {
    for(int i=0;i<seg->n;i++) {
         FREEMEMORY(NULL, seg->chr);
         FREEMEMORY(NULL, seg->value[i]);
    }
  FREEMEMORY(NULL, seg->chr);
  FREEMEMORY(NULL, seg->pos);
  FREEMEMORY(NULL, seg->value);
  FREEMEMORY(NULL, seg->next);
  FREEMEMORY(NULL, seg->parent);
  FREEMEMORY(NULL, seg);
  
  
}

void
destructCpg(cpg_t *cpg) {
    
  FREEMEMORY(NULL, cpg->groupA);
  FREEMEMORY(NULL, cpg->groupB);
  FREEMEMORY(NULL, cpg->chr);
  FREEMEMORY(NULL, cpg);
}

void
destructSegmentSet(segmentset_t *set) {
    segment_t *tmp = set->head;
    segment_t *next;
    while (set->seg) {
        next = tmp->next;
        destructSegment(tmp);
        tmp = next;
    }
  FREEMEMORY(NULL, set->chr);
  FREEMEMORY(NULL, set->nextchr);
  FREEMEMORY(NULL, set);
}


void
initSegmentSet(segmentset_t *set) {
  set->n=0;
  set->firststop=-1;
  set->nextchr=NULL;
  set->nextstart=-1;
  set->seg=NULL;
  set->head=NULL;
  set->tail=NULL;
  set->chr=NULL;
}

segment_t 
*addNewSegmentToSet(segmentset_t *set) {
  segment_t *seg = NULL;
  seg = ALLOCMEMORY(NULL, NULL, segment_t, 1);
  initSegment(seg);
 
  seg->parent= set->tail;
  if(set->head != NULL) {
      set->tail->next = seg;
  }
  else {
      set->head=seg;
  }
  set->tail=seg;   
  set->n++;
  return seg;
}


segment_t 
*getSegment(segmentset_t *set, int n) {
    segment_t *s = NULL;
    if(n<0 || n>=set->n) return s;
    s = set->head; 
    for(int i=1;i<=n;i++) {
        s = s->next;
    }
    return s;
}


void
removeThisSegmentFromSet(segmentset_t *set,segment_t *seg) {
    if(seg->parent == NULL) {
        if(set->n==1) {
            set->head=NULL;
            set->tail=NULL;
            set->n--;

            seg->next=NULL;
            seg->parent=NULL;
            return;
        }
        else {
            set->head=set->head->next;
            set->head->parent=NULL;
            set->n--;

            seg->next=NULL;
            seg->parent=NULL;
            return;
        }        
    }
//remove tail        
    if(seg->next == NULL) {
        if(set->n==1) {
            set->head=NULL;
            set->tail=NULL;
            set->n--;

            seg->next=NULL;
            seg->parent=NULL;
            return;
        }
        else {
            set->tail = set->tail->parent;
            set->tail->next=NULL;
            set->n--;

            seg->next=NULL;
            seg->parent=NULL;
            return;
        }
    }
//remove an object within the list
    seg->parent->next=seg->next;
    seg->next->parent=seg->parent;
    set->n--;

    seg->next=NULL;
    seg->parent=NULL;
    return;
}
    

void
removeSegmentFromSet(segmentset_t *set, int n) {
    //no object to remove
    if(set->n == 0 || n>= set->n) return;
    
    //remove head
    if(n == 0) {
        if(set->n==1) {
            segment_t *seg = set->head;
            set->head=NULL;
            set->tail=NULL;
            set->n--;
            
            seg->next=NULL;
            seg->parent=NULL;
            return;
        }
        else {
            segment_t *seg = set->head;
            set->head=set->head->next;
            set->head->parent=NULL;
            set->n--;
            
            seg->next=NULL;
            seg->parent=NULL;
            return;
        }        
    }
    else 
//remove tail        
        if(n == set->n-1) {
            if(set->n==1) {
                segment_t *seg = set->tail;
                set->head=NULL;
                set->tail=NULL;
                set->n--;
            
                seg->next=NULL;
                seg->parent=NULL;
                return;
            }
            else {
                segment_t *seg = set->tail;
                set->tail = set->tail->parent;
                set->tail->next=NULL;
                set->n--;
             
                seg->next=NULL;
                seg->parent=NULL;
                return;
           }
        }
    //remove an object within the list
        else {
            segment_t *seg = set->head; 
            for(int i=1;i<=n;i++) {
                seg = seg->next;
            }
            seg->parent->next=seg->next;
            seg->next->parent=seg->parent;
            set->n--;
            
            seg->next=NULL;
            seg->parent=NULL;
            return;
        }
}





/*-------------------------------- setSegment --------------------------------
 *    
 * @brief set values in a segment
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void 
setSegment(segment_t *seg, char *chr, int start, int stop, 
    double prob, double meandiff, double test) {

  seg->chr = chr;
  seg->start = start;
  seg->stop = stop;
  seg->prob = prob;
  seg->test = test;
  seg->meandiff = meandiff;
}
/*--------------------------------- get_ratio ---------------------------------
 *    
 * @brief compute ratio of two group means
 * @author Stephan Bernhart
 *   
 */

double 
get_ratio(double *a, int m, double *b, int n) {
  int i;
  double mean1, mean2, ratio;
  mean1=0;
  mean2=0;
  for(i=0; i<m; i++) {
    mean1+=a[i];
  }
  mean1/=(double) m;
  for(i=0; i<n; i++) {
    mean2+=b[i];
  }
  mean2/=(double) n;
  ratio=mean1/mean2;
  return ratio;
}

/*--------------------------------- get_meandiff ---------------------------------
 *    
 * @brief compute ratio of two group means
 * @author Frank Juehling
 *   
 */

double 
get_meandiff(cpg_t *cpg,double *a, int m, double *b, int n) {
  int i;
  double mean1, mean2, meandiff;
  mean1=0;
  mean2=0;
  for(i=0; i<m; i++) {
    mean1+=a[i];
  }
  mean1/=(double) m;
  for(i=0; i<n; i++) {
    mean2+=b[i];
  }
  mean2/=(double) n;
  meandiff=mean1-mean2;
//  fprintf(stdout,"meandiff %f\n",meandiff);
  cpg->methA=mean1;
  cpg->methB=mean2;
  return meandiff;
}

/*--------------------------------- calcMax ----------------------------------
 *    
 * @brief get the maximum quadrant
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

double 
calcMax(double *l1, double *l2, double c1, double c2, int a, int b, 
    double m, double n) {

  double d = 0;

  d=MAX(fabs(((c1)/m)-((c2)/n)),fabs(((c1+l1[0])/m)-((c2+l2[0])/n)));
  d=MAX(fabs(((c1+l1[a])/m)-((c2+l2[a])/n)),d);
  d=MAX(fabs(((c1+l1[b])/m)-((c2+l2[b])/n)),d);
  d=MAX(fabs(((c1+l1[0]+l1[a])/m)-((c2+l2[0]+l2[a])/n)),d);
  d=MAX(fabs(((c1+l1[0]+l1[b])/m)-((c2+l2[0]+l2[b])/n)),d);
  d=MAX(fabs(((c1+l1[b]+l1[a])/m)-((c2+l2[b]+l2[a])/n)),d);
  d=MAX(fabs(((c1+l1[b]+l1[a]+l1[0])/m)-((c2+l2[b]+l2[a]+l2[0])/n)),d);

  return d; 
}

/*--------------------------------- counter ----------------------------------
 *    
 * @brief count the data points in the four quadrants with center (x,y)
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

double 
counter(double x, double y, double *x0, double *x1, int m, 
    double *y0, double *y1, int n) {

  double **c, **l;
  int i;

  c = ALLOCMEMORY(NULL, NULL, double*, 2);
  l = ALLOCMEMORY(NULL, NULL, double*, 2);

  c[0] = ALLOCMEMORY(NULL, NULL, double, 4);
  c[1] = ALLOCMEMORY(NULL, NULL, double, 4);
  l[0] = ALLOCMEMORY(NULL, NULL, double, 5);
  l[1] = ALLOCMEMORY(NULL, NULL, double, 5);
  memset(c[0], 0, sizeof(double)*4);
  memset(c[1], 0, sizeof(double)*4);
  memset(l[0], 0, sizeof(double)*5);
  memset(l[1], 0, sizeof(double)*5);


  for(i=0; i < m; i++) {
    //central point
    if(x0[i] == x && x1[i] == y) {l[0][0]++; continue;}  
    //   if(lx[i] == x || ly[i] == y) {continue;}  
    //points on borders of quadrants
    if(x0[i] == x && x1[i] < y) {l[0][1]++; continue;}  
    if(x0[i] == x && x1[i] > y) {l[0][2]++; continue;}  
    if(x0[i] < x && x1[i] == y) {l[0][3]++; continue;}  
    if(x0[i] > x && x1[i] == y) {l[0][4]++; continue;}  

    if(x0[i] > x) {
      if(x1[i] > y)
        c[0][0]++;
      else
        c[0][1]++;
    }
    else {
      if(x1[i] > y)
        c[0][2]++;
      else
        c[0][3]++;
    }
  }

  for(i=0; i < n; i++) {
    //central point
    if(y0[i] == x && y1[i] == y) {l[1][0]++; continue;}  
    //  if(lx[i] == x || ly[i] == y) {continue;}  
    //points on borders of quadrants
    if(y0[i] == x && y1[i] < y) {l[1][1]++; continue;}  
    if(y0[i] == x && y1[i] > y) {l[1][2]++; continue;}  
    if(y0[i] < x && y1[i] == y) {l[1][3]++; continue;}  
    if(y0[i] > x && y1[i] == y) {l[1][4]++; continue;}  

    if(y0[i] > x) {
      if(y1[i] > y)
        c[1][0]++;
      else
        c[1][1]++;
    }
    else {
      if(y1[i] > y)
        c[1][2]++;
      else
        c[1][3]++;
    }
  }

  double d[4];

  d[0]=calcMax(l[0],l[1],c[0][0],c[1][0],2,4,m,n);
  d[1]=calcMax(l[0],l[1],c[0][1],c[1][1],1,4,m,n);
  d[2]=calcMax(l[0],l[1],c[0][2],c[1][2],2,3,m,n);
  d[3]=calcMax(l[0],l[1],c[0][3],c[1][3],1,3,m,n);

  //fprintf(stdout, "d[0]:%f, d[1]:%f, d[2]:%f, d[3]:%f, d[4]:%f\n", d[0], d[1], d[2], d[3], d[4]);

  double D = MAX(MAX(MAX(d[0],d[1]),d[2]),d[3]);

  FREEMEMORY(NULL, c[0]);
  FREEMEMORY(NULL, c[1]);
  FREEMEMORY(NULL, l[0]);
  FREEMEMORY(NULL, l[1]);
  FREEMEMORY(NULL, c);
  FREEMEMORY(NULL, l);

  return D;  
}

/*--------------------------------- kstest2d ---------------------------------
 *    
 * @brief two-dimensional ks-test
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void 
kstest2d(double *x0, double *x1, int m, double *y0, double *y1, int n, 
    double *kstest) {

  double cor[] = {0.0,0.0};
  double ks[] = {0.0,0.0};
  double d[] = {0.0,0.0};
  double dl1, dl2, s;

  for(int j=0; j < m; j++) {
    d[0] = MAX(d[0],counter(x0[j],x1[j],x0,x1,m,y0,y1,n));
  }

  for(int j=0; j< n; j++) {
    d[1]=  MAX(d[1],counter(y0[j],y1[j],x0,x1,m,y0,y1,n));
  }

  ks[1]=(d[0]+d[1])*0.5;

  cor[0] = rho(NULL, x0, x1, m);
  cor[1] = rho(NULL, y0, y1, n);

  dl1 = m;
  dl2 = n;
  s = sqrt(dl1*dl2/(dl1+dl2));
  kstest[0]  = kscdf(ks[1]*s/(1.0+sqrt(1.0-0.5*(cor[0]*cor[0]+cor[1]*cor[1]))*(0.25-0.75/s)));

  //fprintf(stdout, "s:%f, cor[0]:%f, cor[1]:%f, ks[1]:%f, test:%f\n", s, 
  //cor[0], cor[1], ks[1], kstest[0]);

  return;
}

/*---------------------------------- kstest ----------------------------------
 *    
 * @brief calculated the ks test
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void kstest(segment_t *seg , int a, int b, char mindiff, char mincpgs, char test, 
    double *ks, int *grpA, int noA, int *grpB, int noB, metseg_t* nfo){

  int i,j;
  int l = b-a+1;
  double **la;
  double **lb;
  double *x;
  double *y;
  double mean1=0;
  double mean2=0;
  double meandiff = 0;
  double dl1=noA;
  double dl2=noB;
  double dl = l;
  double p = 2;
  double faskstest[2] = {0,0};

  //mincpgs set to false by default so this condition is never fullfilled

  mindiff = 0;
  mincpgs = 0;
  if(l< nfo->mincpgs && mincpgs) {
    ks[0]=2;ks[1]=0;ks[2]=2;
  }
  //debug
  //  fprintf(stderr, "kstest for: %d\n",l);
  la = ALLOCMEMORY(NULL, NULL, double*, 2);
  lb = ALLOCMEMORY(NULL, NULL, double*, 2);

  la[0] = ALLOCMEMORY(NULL, NULL, double, l*noA);
  la[1] = ALLOCMEMORY(NULL, NULL, double, l*noA);
  lb[0] = ALLOCMEMORY(NULL, NULL, double, l*noB);
  lb[1] = ALLOCMEMORY(NULL, NULL, double, l*noB);
  memset(la[0], 0, sizeof(double)*l*noA);
  memset(la[1], 0, sizeof(double)*l*noA);
  memset(lb[0], 0, sizeof(double)*l*noB);
  memset(lb[1], 0, sizeof(double)*l*noB);

  x = ALLOCMEMORY(NULL, NULL, double, l);
  y = ALLOCMEMORY(NULL, NULL, double, l);

  ks[0]=0;ks[1]=0;ks[2]=0;

  for(i=a; i<=b ; i++) {
    for(j=0; j < noA; j++) {
      la[0][((i-a)*noA)+j]=seg->value[i][grpA[j]];
      la[1][((i-a)*noA)+j]=seg->pos[i];
      mean1+=seg->value[i][grpA[j]];
      x[i-a] += seg->value[i][grpA[j]];
    }
    x[i-a]/=dl1;
  }
  for(i=a; i<=b; i++) {
    for(j=0;j<noB;j++) {
      lb[0][((i-a)*noB)+j]=seg->value[i][grpB[j]];
      lb[1][((i-a)*noB)+j]=seg->pos[i];
      mean2+=seg->value[i][grpB[j]];
      y[i-a] += seg->value[i][grpB[j]];
    }
    y[i-a]/=dl2;
  }
  int u = mannwhitney( la[0], l*noA , lb[0], l*noB);
  p = mannwhitneyPvalue(u, l*noA, l*noB, nfo->MWU, MAXM, MAXN);
 // fprintf(stdout,"pVALUE %f\n",p);
  

  mean1/=dl1*dl;
  mean2/=dl2*dl;
  meandiff = mean1-mean2;
  seg->methA=-1.0;
  seg->methB=-1.0;
  
  kstest2d(la[0], la[1], l*noA, lb[0],lb[1], l*noB, faskstest);

  //fprintf(stdout, "ks: %f\n", faskstest[0]);
  ks[0]=faskstest[0];
  ks[1]=meandiff;
  ks[2]=p;

  FREEMEMORY(NULL, la[0]);
  FREEMEMORY(NULL, la[1]);
  FREEMEMORY(NULL, lb[0]);
  FREEMEMORY(NULL, lb[1]);
  FREEMEMORY(NULL, la);
  FREEMEMORY(NULL, lb);
  FREEMEMORY(NULL, x);
  FREEMEMORY(NULL, y);

}

/*---------------------------------- kstest ----------------------------------
 *    
 * @brief calculated the means
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
void means(segment_t *seg , int a, int b, int *grpA, int noA, int *grpB, int noB, double *means){

  int i,j;
  double mean1=0;
  double mean2=0;
  double dl1=0;
  double dl2=0;
  
   for(i=a; i<=b ; i++) {
    for(j=0; j < noA; j++) {
     mean1+=seg->value[i][grpA[j]];
     dl1+=1;
    }
  }
  for(i=a; i<=b; i++) {
    for(j=0;j<noB;j++) {
      mean2+=seg->value[i][grpB[j]];
      dl2+=1;
    }
  }
  
  mean1/=dl1;
  mean2/=dl2;
  means[0]=mean1;
  means[1]=mean2;
  
  

}


/*---------------------------- calcSingleDiffSum -----------------------------
 *    
 * @brief calculate single difference sum
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

double **
calcSingleDiffSum(segment_t *s, int *grpA, int noA, int *grpB, int noB){

  int i, j;
  double **S, s1, s2;

  S = ALLOCMEMORY(NULL, NULL, double*, s->n);

  for(i=0; i < s->n; i++) {
    S[i] = ALLOCMEMORY(NULL, NULL, double, 3);
    memset(S[i], 0, sizeof(double)*3);
  }

  for(i=0; i < s->n; i++) {
    s1 = 0;
    for(int j=0;j<noA;j++)
      s1+=s->value[i][grpA[j]];
    s1/=noA;

    s2 = 0;
    for(j=0; j < noB; j++)
      s2+=s->value[i][grpB[j]];
    s2/=noB;

    if(i==0) {

      S[i][1]=s1-s2;
      S[i][0]=fabs(S[i][1]);
      if(s1-s2==0) {
        S[i][2] = 0;

      }
      else{
        if(s1-s2>0)
          S[i][2] = 1;
        else
          S[i][2] = -1;
      }

    } else {

      S[i][1]=S[i-1][1]+s1-s2;
      S[i][0]=fabs(S[i][1]);
      if(s1-s2==0)
        S[i][2] = 0;
      else{
        if(s1-s2>0)
          S[i][2] = S[i-1][2]+1;
        else
          S[i][2] = S[i-1][2]-1;
      }                    
    }
  }

  return S;
}

/*----------------------------- calcSingleTrend ------------------------------
 *    
 * @brief calculate single segment trend
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

double 
calcSingleTrendAbs(double **S,int s, int t) {
  double ds = s;
  double dt=t;
  double trend; 
  if(s==0)
    trend = fabs(S[t][2])/dt;
  else
    trend = fabs(S[t][2]-S[s-1][2])/(dt-ds+1);

  return trend;
}

/*--------------------------------- noValley ---------------------------------
 *    
 * @brief check for local valley
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

char 
noValley(double **S, int s, int t, metseg_t *nfo){

  if(t-s+1< nfo->mincpgs) return 1;

  int i;
  double ds = s;
  double dt=t;
  double Sst=S[t][1];
  double Si, mean, imean;
  int minlength=(nfo->mincpgs>10)?nfo->mincpgs:10;
  
  if(s>0) {
    Sst-=S[s-1][1];
  }

  mean=fabs(Sst/(dt-ds+1));
//check if mean of windowsize=mincpgs is < valleyfactor*mean
  for(i=s; i+minlength-1 <= t; i++) {
    Si=S[i+minlength-1][1];
    if(i>0) Si-=S[i-1][1];
    imean=fabs(Si/(minlength));

//    if(fabs(imean)<mean*0.7) {
    if(fabs(imean)<mean* nfo->valley) {
      return 0;
    }    
  }
  return 1;
}



/*--------------------------------- findMaxZ ---------------------------------
 *    
 * @brief find maximum Z score
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

double 
findMaxZ(double **Z, int s, int t, int *ab, metseg_t *nfo) {
  int a, b;
  double max=0;

  //fprintf(stderr, "FindMaxN\t%d\t%d\n",s,t); 
  for(a=s; a<=t; a++)
    for(b=a+nfo->mincpgs-1; b<=t; b++) {
      if(((a!=s)||(b!=t) ) && Z[a][b]>max && (IsFiniteNumber(Z[a][b]))){
        max=Z[a][b];
        ab[0]=a;
        ab[1]=b;
      }  
    }
  //fprintf(stderr, "FindMaxNnew\t%d\t%d\n",ab[0],ab[1]); 
  return max;
}


/*---------------------------- calcSingleDiffZabs ----------------------------
 *    
 * @brief calculate single Z score
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void 
calcSingleDiffZabs(double **S, int Ssize, double **Z, int s, int t, 
    char zero, metseg_t *nfo) {

  int a, b;
  double ds = s;
  double dt=t;
  double da, db, Sst, Sab, u;

  for(a=s; a<=t; a++){ 
    for(b=a+nfo->mincpgs-1; b<=t; b++) {
      da =a;
      db=b;
      Sst=S[t][1];
      if(s>0) { 
        Sst-=S[s-1][1];
      }
      Sab=S[b][1];
      if(a>0) { 
        Sab-=S[a-1][1];
      }
      Sab=fabs(Sab);
      Sst=fabs(Sst);

      if(db-da-nfo->mincpgs == 0 || (a==s && b==t) || a==b) { 
        Z[a-s][b-s]=0;
      } else {
        if(zero) {
          u = Sab;
        } else {
          u = Sab-(((db-da+1)*(Sst))/(dt-ds+1));
        }
        u*=u;
        u/=(db-da+1)*((1-((db-da+1)/(dt-ds+1))));
        Z[a-s][b-s]=u;
      }
    }
  }

  return;
}    


/*------------------------------ pushSegment_p -------------------------------
 *    
 * @brief a helper function to push segments to a stack
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
 
void
pushSegment_p (Segmentstack *stack, int a, int b, int *ab, int child,
    double *ks1, double *ks2, double *ks3, double *KS, segment_t *max)
{

    segment_p_t segment;

    segment.a = a;
    segment.b = b;
    if(child) { 
    segment.ab1 = ab[0];
    segment.ab2 = ab[1];
    segment.child = child;

    segment.ks11 =ks1[0];
    segment.ks12 =ks1[1];
    segment.ks13 =ks1[2];
 
    segment.ks21 =ks2[0];
    segment.ks22 =ks2[1];
    segment.ks23 =ks2[2];

    segment.ks31 =ks3[0];
    segment.ks32 =ks3[1];
    segment.ks33 =ks3[2];

    segment.KS1 =KS[0];
    segment.KS2 =KS[1];
    segment.KS3 =KS[2];
    }
    if(max) { 
      segment_t *copy = ALLOCMEMORY(NULL, NULL, segment_t, 1);
      memmove(copy, max, sizeof(segment_t));
      segment.max = copy;
    }
    bl_segmentstackPush(stack, &segment);

	return ;
}

/*------------------------------- popSegment_p -------------------------------
 *    
 * @brief a helper function to pop segments from a stack
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
 
void
popSegment_p (Segmentstack *stack, int *a, int *b, int *ab, int *child, 
    double *ks1, double *ks2, double *ks3, double *KS, segment_t **max)
{
  segment_p_t *segment;

  segment = bl_segmentstackPop(stack);
     
  *a = segment->a;
  *b = segment->b;

  if(child) { 
    *child = segment->child;

    ab[0] = segment->ab1;
    ab[1] = segment->ab2;

    ks1[0] = segment->ks11;
    ks1[1] = segment->ks12;
    ks1[2] = segment->ks13;

    ks2[0] = segment->ks21;
    ks2[1] = segment->ks22;
    ks2[2] = segment->ks23;

    ks3[0] = segment->ks31;
    ks3[1] = segment->ks32;
    ks3[2] = segment->ks33;

    KS[0] = segment->KS1;
    KS[1] = segment->KS2;
    KS[2] = segment->KS3;
  }

  if(max) { 
    *max = segment->max;
  }

  
  return;
}


/*------------------------------ stackSegment_p ------------------------------
 *    
 * @brief a helper function to initalize a segment stack
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
 
  Segmentstack*
stackSegment_p(Segmentstack *stack)
{
    bl_segmentstackInit(stack, 300);

	return stack;
}

/*-------------------------------- segment_p ---------------------------------
 *    
 * @brief calculate segment probability
 * @author Frank Juehling and Steve Hoffmann  
 *   
 */


segment_t*
segment_pSTKopt(segment_t *seg, segment_t *breaks, int *nbreaks, double **XS, 
    int a, int b, double *KS, char first, 
    int *grpA, int noA, int *grpB, int noB, metseg_t *nfo) {

  int i, n, m;
  int dimZ, child = 0;
  int ab[2] = {-1,0};
  double ks1[] = {2,0,2};
  double ks2[] = {2,0,2};
  double ks3[] = {2,0,2};
  double init[] = {2,0,2};

  double **XZ;
  double newp; // prob, meandiff;

  Segmentstack stack;
  stackSegment_p(&stack);

  while(!bl_segmentstackIsEmpty(&stack) || a != -1) { 

    if(a != -1 && child <= 2) { 
      if(ab[0] == -1) {  

        ab[0] = 0;
        ab[1] = 0;

        memmove(ks1, init, sizeof(double)*3);
        memmove(ks2, init, sizeof(double)*3);
        memmove(ks3, init, sizeof(double)*3);

        dimZ = b - a + 1;
        XZ = ALLOCMEMORY(NULL, NULL, double*, dimZ);

        for(i=0; i < dimZ; i++) {
          XZ[i] = ALLOCMEMORY(NULL, NULL, double, dimZ);
          memset(XZ[i], 0, sizeof(double)*dimZ);
        }

        //calculated the Z scores and find the maximum interval
        calcSingleDiffZabs(XS, seg->n, XZ , a, b, 0, nfo);
        findMaxZ(XZ, 0, dimZ-1, ab, nfo);
        ab[0]+=a; 
        ab[1]+=a;

        //check the left side of the maximum interval with ks
        n=a; m=ab[0]-1;
        if(ab[0] > 0 && m-n+1 >= nfo->mincpgs 
            && calcSingleTrendAbs(XS,a,ab[0]-1) > nfo->trend 
            && noValley(XS, a, ab[0]-1, nfo)) {

          kstest(seg, a, ab[0]-1, 0, 1, 1, ks1, grpA, noA, grpB, noB, nfo);
        }

        //check the maximum interval interval with ks
        n=ab[0];m=ab[1];
        if(m-n+1 >= nfo->mincpgs 
            && calcSingleTrendAbs(XS,ab[0],ab[1]) > nfo->trend 
            && noValley(XS, ab[0], ab[1], nfo) ) {

          kstest(seg, ab[0], ab[1], 0, 1, 1, ks2, grpA, noA, grpB, noB, nfo);
        }

        //check the right side of the maximum interval with ks
        n=ab[1]+1;m=b;
        if(m-n+1 >= nfo->mincpgs 
            && calcSingleTrendAbs(XS,ab[1]+1,b)> nfo->trend 
            && noValley(XS, ab[1]+1, b, nfo)) {

          kstest(seg, ab[1]+1, b, 0, 1, 1, ks3, grpA, noA, grpB, noB, nfo);
        }

        for(i=0; i < dimZ; i++) {
          FREEMEMORY(NULL, XZ[i]);
        }

        FREEMEMORY(NULL, XZ);
      }

      //if one of the children has a good ks we check all
      newp = MIN(ks1[0],MIN(ks2[0],ks3[0]));

      //if the p-value of one of the three intervals is better than the original
      //recurse down to find the best subinterval
      if((newp<KS[0] || (newp>1 && KS[0]>1)) && (b-a >= nfo->mincpgs)) {

        pushSegment_p (&stack, a, b, ab, child+1, ks1, ks2, ks3, KS, NULL);

        //left interval child
        if(child == 0) {
          n = a;
          m = ab[0]-1;
          a = -1;
          if(ab[0] > 0 && n <= m) {
            if(m-n >= nfo->mincpgs) { 
              a = n;
              b = m;
              memmove(KS, ks1, sizeof(double)*3);
              child = 0;
              ab[0] = -1;
            } else { 
              breaks = ALLOCMEMORY(NULL, breaks, segment_t, (*nbreaks)+1);
              setSegment(&breaks[(*nbreaks)], seg->chr, n, m, ks1[0], ks1[1], ks1[2]);
              (*nbreaks)+=1;
            }
          }
        }

        //middle interval child
        if(child == 1) {
          n = ab[0];
          m = ab[1];
          a = -1;
          if(n <= m) {
            if(m-n>= nfo->mincpgs) { 
              a = n;
              b = m;
              memmove(KS, ks2, sizeof(double)*3);
              child = 0;
              ab[0] =-1;
            } else {
              breaks = ALLOCMEMORY(NULL, breaks, segment_t, (*nbreaks)+1);
              setSegment(&breaks[(*nbreaks)], seg->chr, n, m, ks2[0], ks2[1], ks2[2]);
              (*nbreaks)+=1;
            }  
          }
        }

        //right interval child
        if(child == 2) {
          n = ab[1]+1;
          m = b;
          a = -1;
          if(n<=m) {
            if (m-n>= nfo->mincpgs) { 
              a = n;
              b = m;
              memmove(KS, ks3, sizeof(double)*3);
              child = 0;
              ab[0] = -1;
            } else {
              breaks = ALLOCMEMORY(NULL, breaks, segment_t, (*nbreaks)+1);
              setSegment(&breaks[(*nbreaks)], seg->chr, n, m, ks3[0], ks3[1], ks3[2]);
              (*nbreaks)+=1;
            } 
          }
        }
      } else {
        if(child ==0) { 
          breaks = ALLOCMEMORY(NULL, breaks, segment_t, (*nbreaks)+1);
          setSegment(&breaks[(*nbreaks)], seg->chr, a, b, KS[0], KS[1], KS[2]);
          (*nbreaks)+=1;
        }

        a = -1;
        ab[0] = -1;
      }
    }  else {

      popSegment_p (&stack, &a, &b, ab, &child, ks1, ks2, ks3, KS, NULL);

      if(child == 3) a = -1;

    }
  }

  bl_segmentstackDestruct(&stack);

  return breaks;
}
/*-------------------------------- segmenter ---------------------------------
 *    
 * @brief first segmenter function
 * @author Frank Juehling and Steve Hoffmann 
 *    
 */
segment_t*
segmenterSTK(segment_t *seg, segment_t *globalbreaks, int *nglobal, double **XS, 
    int a, int b, double *KS, int *grpA, int noA, int *grpB, int noB, 
    metseg_t *nfo) {

  int nbreaks=0, i, n, m; //, *s, *t;
  double trend;
  segment_t *breaks=NULL, *max;

  Segmentstack stack; 
  stackSegment_p(&stack);
  //inorder traversal of tree of intervals
  //while(!bl_vstackIsEmpty(stack) || a != -1) {
  while(!bl_segmentstackIsEmpty(&stack) || a != -1) {

    if(a != -1) { 
      //push current interval node
      nbreaks = 0;
      breaks = NULL;
      breaks = segment_pSTKopt(seg, breaks, &nbreaks, XS, a, b, KS, 0, 
          grpA, noA, grpB, noB, nfo);

      max = &breaks[0];

      for(i=0; i < nbreaks; i++) { 
        if(max->prob > breaks[i].prob) { 
          max = &breaks[i];
        }
      }

      pushSegment_p (&stack, a, b, NULL, 0, NULL, NULL, NULL, NULL, max);

      //set next
      n = a;
      m = max->start-1;

      if(max->start > 0 && n <= m) {

        trend = calcSingleTrendAbs(XS, n, m);
        double ks[] = {2,0,2};

        if(m-n+1 >= nfo->mincpgs && trend > nfo->trend && noValley(XS, n, m, nfo)) { 
          kstest(seg, n, m, 0, 1, 1, ks, grpA, noA, grpB, noB, nfo);
        }

        a = n;
        b = m;
        KS = ks;
      } else {
        a = -1;
      }

      FREEMEMORY(NULL, breaks);

    } else {

       
      popSegment_p (&stack, &a, &b, NULL, 0, NULL, NULL, NULL, NULL, &max);

      //here the segment is registered!
      globalbreaks = ALLOCMEMORY(NULL, globalbreaks, segment_t, (*nglobal)+1);
      memmove(&globalbreaks[(*nglobal)], max, sizeof(segment_t)); 
      (*nglobal)+=1;

      n = max->stop+1;
      m = b;
      if(n<=m) {
        trend = calcSingleTrendAbs(XS,n,m);
        double ks[] = {2,0,2};
        if(m-n+1 >= nfo->mincpgs && trend > nfo->trend && noValley(XS, n, m, nfo)) { 
          kstest(seg, n, m, 0, 1, 1, ks, grpA, noA, grpB, noB, nfo);
        }

        a = n;
        b = m;
        KS = ks;
      } else {
        a = -1;
      }

      FREEMEMORY(NULL, max);
    }
  }

  bl_segmentstackDestruct(&stack); 
  return globalbreaks;
}
/*------------------------------- calcMeandiff -------------------------------
 *    
 * @brief calculate mean difference between grpA and grpB
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
double 
calcMeandiff(segment_t *s , int a, int b, int *grpA, int noA, 
    int *grpB, int noB, metseg_t* nfo){

  int i, j;
  double l = b-a+1;
  double meandiff=0, mean1=0, mean2=0;
  double dl1=noA, dl2=noB;

  for(i=a; i<=b; i++) { 
    for(j=0; j<noA; j++) {
      mean1+=s->value[i][grpA[j]];
    }
  }

  for(i=a; i<=b; i++){ 
    for(j=0; j<noB; j++) {
      mean2+=s->value[i][grpB[j]];
    }
  }

  mean1/=dl1*l;
  mean2/=dl2*l;
  meandiff = mean1-mean2;
  return meandiff;
}
/*---------------------------------- output ----------------------------------
 *    
 * @brief output routine
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
void
outputO(segment_t *seg, segment_t *breaks, int nglobal, double **XS, 
    int *grpA, int noA, int *grpB, int noB, metseg_t *nfo) {

  int i;
  double meandiff, trend;
  segment_t *b, *tmp=NULL;

  for(i=0; i<nglobal; i++) {

    b = &breaks[i];
    

    if(b->prob > 1) {
      if(tmp == NULL) {

        tmp = ALLOCMEMORY(NULL, NULL, segment_t, 1);
        tmp->chr = b->chr;
        tmp->start=b->start;
        tmp->stop=b->stop;
        tmp->prob=b->prob;
        tmp->meandiff=b->meandiff;
        tmp->test=b->test;

      } else {
        tmp->stop=b->stop;
      }

    } else {

        
  // ks[0]=faskstest[0];
 // ks[1]=meandiff;
 // ks[2]=p;
     
        
        
        
      if(tmp != NULL) {

        trend = calcSingleTrendAbs(XS,tmp->start,tmp->stop);
        double ks[] = {2,0,2};

        if(tmp->stop-tmp->start + 1 >= nfo->mincpgs && trend>nfo->trend 
            && noValley(XS, tmp->start, tmp->stop, nfo)) {

          kstest(seg, tmp->start,tmp->stop, 0, 1, 1, ks, 
              grpA, noA, grpB, noB, nfo);
        }

        if(ks[0]<2) {

          fprintf(stdout, "%s\t%d\t%d\t%.2g\t%f\t%d\t%.2g\n", 
              seg->chr, seg->pos[tmp->start]-1, 
              seg->pos[tmp->stop], ks[0], ks[1],
              (tmp->stop-tmp->start+1), ks[2]);
	  fflush(stdout); 

        } else {
  
          meandiff = calcMeandiff(seg, tmp->start, tmp->stop, 
              grpA, noA, grpB, noB, nfo);

	   fprintf(stdout, "%s\t%d\t%d\tNA\t%f\t%d\tNA\n", 
              seg->chr,seg->pos[tmp->start]-1, seg->pos[tmp->stop], meandiff, 
              (tmp->stop-tmp->start+1));
        }

        FREEMEMORY(NULL, tmp);
        tmp=NULL;
      }

      fprintf(stdout, "%s\t%d\t%d\t%.2g\t%f\t%d\t%.2g\n", seg->chr, 
          seg->pos[b->start]-1, seg->pos[b->stop], b->prob, b->meandiff,
          (b->stop-b->start+1), b->test);
      fflush(stdout); 
    }
  }


  if(tmp != NULL) {

    trend = calcSingleTrendAbs(XS,tmp->start,tmp->stop);
    double ks[] = {2,0,2};

    if(tmp->stop-tmp->start + 1 >= nfo->mincpgs && trend > nfo->trend 
        && noValley(XS, tmp->start, tmp->stop, nfo)) {
      kstest(seg,tmp->start,tmp->stop,0, 1, 1, ks, grpA, noA, grpB, noB, nfo);
    }

    if(ks[0]<2) {
      fprintf(stdout, "%s\t%d\t%d\t%.2g\t%f\t%d\t%.2g\n", seg->chr, 
          seg->pos[tmp->start]-1, seg->pos[tmp->stop], ks[0], ks[1], 
          (tmp->stop-tmp->start+1), ks[2]);
      fflush(stdout); 
    } else {
  
      meandiff = calcMeandiff(seg, tmp->start, tmp->stop, 
          grpA, noA, grpB, noB, nfo);

          fprintf(stdout, "%s\t%d\t%d\tNA\t%f\t%d\tNA\n", seg->chr, 
          seg->pos[tmp->start]-1, seg->pos[tmp->stop], meandiff,(
	  tmp->stop-tmp->start+1));
    }


    FREEMEMORY(NULL, tmp);
    tmp=NULL;
  }

  return;
}
/*---------------------------------- output ----------------------------------
 *    
 * @brief output routine
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
void
output(segment_t *seg, segment_t *breaks, int nglobal, double **XS, 
    int *grpA, int noA, int *grpB, int noB, metseg_t *nfo) {
  
    
  if(nfo->outputList->i >= nfo->outputList->n) {
        nfo->outputList->n+=1000000;
        nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
  }
  
        
  int i;
  double trend;
  segment_t *b, *tmp=NULL;
  nfo->outputList->numberTests+=nglobal;
  for(i=0; i<nglobal; i++) {

    b = &breaks[i];
    

    if(b->prob > 1) {
      if(tmp == NULL) {

        tmp = ALLOCMEMORY(NULL, NULL, segment_t, 1);
        tmp->chr = b->chr;
        tmp->start=b->start;
        tmp->stop=b->stop;
        tmp->prob=b->prob;
        tmp->meandiff=b->meandiff;
        tmp->test=b->test;
        
        double me[] = {-2,-2};
        means(seg, tmp->start,tmp->stop, grpA, noA, grpB, noB, me);
        
        tmp->methA=me[0];
        tmp->methB=me[1];
        

      } else {
        tmp->stop=b->stop;
      }

    } else {

        
  // ks[0]=faskstest[0];
 // ks[1]=meandiff;
 // ks[2]=p;
     
        
        
        
      if(tmp != NULL) {

        trend = calcSingleTrendAbs(XS,tmp->start,tmp->stop);
        double ks[] = {2,0,2};

        if(tmp->stop-tmp->start + 1 >= nfo->mincpgs && trend>nfo->trend 
            && noValley(XS, tmp->start, tmp->stop, nfo)) {

          kstest(seg, tmp->start,tmp->stop, 0, 1, 1, ks, 
              grpA, noA, grpB, noB, nfo);
        }

        if(ks[0]<2) {
            nfo->outputList->segment_out[nfo->outputList->i].chr = ALLOCMEMORY(NULL, NULL, char, strlen(seg->chr)+1);
            nfo->outputList->segment_out[nfo->outputList->i].chr = strcpy(nfo->outputList->segment_out[nfo->outputList->i].chr,seg->chr);
            nfo->outputList->segment_out[nfo->outputList->i].start = seg->pos[tmp->start]-1;
            nfo->outputList->segment_out[nfo->outputList->i].stop = seg->pos[tmp->stop];
            nfo->outputList->segment_out[nfo->outputList->i].p = ks[0];
            nfo->outputList->segment_out[nfo->outputList->i].meandiff = ks[1];
            nfo->outputList->segment_out[nfo->outputList->i].mwu = ks[2];
            nfo->outputList->segment_out[nfo->outputList->i].length = (tmp->stop-tmp->start+1);
            
            double me[] = {-2,-2};
            means(seg,tmp->start,tmp->stop, grpA, noA, grpB, noB, me);
            nfo->outputList->segment_out[nfo->outputList->i].methA = me[0];
            nfo->outputList->segment_out[nfo->outputList->i].methB = me[1];
            
            nfo->outputList->i+=1;
            if(nfo->outputList->i >= nfo->outputList->n) {
                nfo->outputList->n+=1000000;
                nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
            }

        } 

        FREEMEMORY(NULL, tmp);
        tmp=NULL;
      }
       nfo->outputList->segment_out[nfo->outputList->i].chr = ALLOCMEMORY(NULL, NULL, char, strlen(seg->chr)+1);
        nfo->outputList->segment_out[nfo->outputList->i].chr = strcpy(nfo->outputList->segment_out[nfo->outputList->i].chr,seg->chr);
        nfo->outputList->segment_out[nfo->outputList->i].start = seg->pos[b->start]-1;
        nfo->outputList->segment_out[nfo->outputList->i].stop = seg->pos[b->stop];
        nfo->outputList->segment_out[nfo->outputList->i].p = b->prob;
        nfo->outputList->segment_out[nfo->outputList->i].meandiff = b->meandiff;
        nfo->outputList->segment_out[nfo->outputList->i].mwu = b->test;
        nfo->outputList->segment_out[nfo->outputList->i].length = (b->stop-b->start+1);
        
        double me[] = {-2,-2};
        means(seg, b->start,b->stop, grpA, noA, grpB, noB, me);
        nfo->outputList->segment_out[nfo->outputList->i].methA = me[0];
        nfo->outputList->segment_out[nfo->outputList->i].methB = me[1];
           
        
        nfo->outputList->i+=1;
        if(nfo->outputList->i >= nfo->outputList->n) {
            nfo->outputList->n+=1000000;
            nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
        }
    }
  }


  if(tmp != NULL) {

    trend = calcSingleTrendAbs(XS,tmp->start,tmp->stop);
    double ks[] = {2,0,2};

    if(tmp->stop-tmp->start + 1 >= nfo->mincpgs && trend > nfo->trend 
        && noValley(XS, tmp->start, tmp->stop, nfo)) {
      kstest(seg,tmp->start,tmp->stop,0, 1, 1, ks, grpA, noA, grpB, noB, nfo);
    }

    if(ks[0]<2) {
        nfo->outputList->segment_out[nfo->outputList->i].chr = ALLOCMEMORY(NULL, NULL, char, strlen(seg->chr)+1);
        nfo->outputList->segment_out[nfo->outputList->i].chr = strcpy(nfo->outputList->segment_out[nfo->outputList->i].chr,seg->chr);
        nfo->outputList->segment_out[nfo->outputList->i].start = seg->pos[tmp->start]-1;
        nfo->outputList->segment_out[nfo->outputList->i].stop = seg->pos[tmp->stop];
        nfo->outputList->segment_out[nfo->outputList->i].p = ks[0];
        nfo->outputList->segment_out[nfo->outputList->i].meandiff = ks[1];
        nfo->outputList->segment_out[nfo->outputList->i].mwu = ks[2];
        nfo->outputList->segment_out[nfo->outputList->i].length = (tmp->stop-tmp->start+1);
        
        double me[] = {-2,-2};
        means(seg, tmp->start,tmp->stop, grpA, noA, grpB, noB, me);
        nfo->outputList->segment_out[nfo->outputList->i].methA = me[0];
        nfo->outputList->segment_out[nfo->outputList->i].methB = me[1];
        
        nfo->outputList->i+=1;
        if(nfo->outputList->i >= nfo->outputList->n) {
            nfo->outputList->n+=1000000;
            nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
        }
    }
    FREEMEMORY(NULL, tmp);
    tmp=NULL;
  }

  return;
}
//these variables need to be volatile to stop the compiler from optimizing
//as they are manipulated by the threads!
static volatile unsigned int idle;
volatile char *schedule = NULL;
static pthread_mutex_t out;
static pthread_mutex_t cnt;

/*------------------------------- segmentation -------------------------------
 *    
 * @brief do the segmentation
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */


int 
segmentation(char **chr, int *pos, double **value, int n, 
    int *grpA, int noA, int *grpB, int noB, metseg_t *nfo) {

  double **S, trend;
  double ks[] = {2,0,2};
  char novalley;
  segment_t *seg, *global=NULL;
  int i, nglobal = 0;

  seg = ALLOCMEMORY(NULL, NULL, segment_t, 1);

  initSegment(seg);
  seg->n = n;
  seg->chr = chr[0];
  seg->pos = pos;
  seg->value = value;

  S = calcSingleDiffSum(seg, grpA, noA, grpB, noB);
  trend = calcSingleTrendAbs(S, 0, n-1);
  novalley = noValley(S, 0, n-1, nfo);

  if(seg->n-1 >= nfo->mincpgs && trend > nfo->trend && novalley) {
    kstest(seg , 0, n-1, 0, 1, 1, ks, grpA, noA, grpB, noB, nfo);
  }

  global = segmenterSTK(seg, global, &nglobal, S, 0, n-1, ks, 
      grpA, noA, grpB, noB, nfo);

  //set lock if necessary 
  if(nfo->threads >1) { 
    pthread_mutex_lock(&out);
  }
  //output here   
  output(seg, global, nglobal, S, grpA, noA, grpB, noB, nfo);
 
  //unlock if necessary
  if(nfo->threads > 1) {
    pthread_mutex_unlock(&out);
  }

  for(i=0; i < seg->n; i++) {
    FREEMEMORY(NULL, S[i]);
  }


  FREEMEMORY(NULL, global);
  FREEMEMORY(NULL, seg);
  FREEMEMORY(NULL, S);

  return 0;
}



/*-------------------------------- segworker ---------------------------------
 *    
 * @brief for threaded segmentation
 * @author Frank Juehling and Steve Hoffma2nn 
 *   
 */
 
void*
segworker (void *args)
{
  int i;
  metseg_t *t;
  t = (metseg_t*) args;
   
  segmentation(t->chr, t->pos, t->value, t->n, t->grpA, t->noA, t->grpB, t->noB, t);
  
  //cleanup own data
  for(i=0; i < t->n; i++) {
    FREEMEMORY(NULL, t->chr[i]);
    FREEMEMORY(NULL, t->value[i]);
  }
  FREEMEMORY(NULL, t->chr);
  FREEMEMORY(NULL, t->pos);
  FREEMEMORY(NULL, t->value);


  pthread_mutex_lock(&cnt);
  schedule[t->threadno] = 0;
  //decrement of idle at creation time
  idle = idle+1;
  pthread_mutex_unlock(&cnt);

  pthread_exit(args);
}

/*------------------------------- regionTest -------------------------------
 *    
 * @brief test a region
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void 
regionTest(segment_t *seg, int *grpA, int noA, int *grpB, int noB, metseg_t *nfo) {
    double ks[] = {2,0,2};
    if(seg->n>0) {
        kstest(seg , 0, seg->n-1, 0, 0, 1, ks, grpA, noA, grpB, noB, nfo);
       
    }
    double me[] = {-2,-2};
    means(seg, 0, seg->n-1,grpA, noA, grpB, noB,me);
    
//    void kstest(segment_t *seg , int a, int b, char mindiff, char mincpgs, char test, 
//  (segment_t *seg , int a, int b, char mindiff, char mincpgs, char test, 
//    double *ks, int *grpA, int noA, int *grpB, int noB, metseg_t* nfo){
    
    
    if(nfo->threads >1) { 
      pthread_mutex_lock(&out);
    }
 
/*    
  //output here   
    fprintf(stdout, "#%s\t%d\t%d\t%.2g\t%f\t%d\t%.2g\n", 
	    seg->chr, seg->start-1, 
	    seg->stop, ks[0], ks[1],
	    seg->n, ks[2]);
    fflush(stdout); 
  */  
    
    if(nfo->outputList->i >= nfo->outputList->n-1) {
      nfo->outputList->n+=1000000;
      nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
    }
    
    nfo->outputList->segment_out[nfo->outputList->i].chr = ALLOCMEMORY(NULL, NULL, char, strlen(seg->chr)+1);
    nfo->outputList->segment_out[nfo->outputList->i].chr = strcpy(nfo->outputList->segment_out[nfo->outputList->i].chr,seg->chr);
    nfo->outputList->segment_out[nfo->outputList->i].start = seg->start-1;
    nfo->outputList->segment_out[nfo->outputList->i].stop = seg->stop;
    nfo->outputList->segment_out[nfo->outputList->i].p = ks[0];
    nfo->outputList->segment_out[nfo->outputList->i].meandiff = ks[1];
    nfo->outputList->segment_out[nfo->outputList->i].mwu = ks[2];
    nfo->outputList->segment_out[nfo->outputList->i].length = seg->n;
    nfo->outputList->segment_out[nfo->outputList->i].methA = me[0];
    nfo->outputList->segment_out[nfo->outputList->i].methB = me[1];
    
    nfo->outputList->i+=1;
    nfo->outputList->numberTests+=1;
    
    
    
    
    
    //unlock if necessary
    if(nfo->threads > 1) {
    pthread_mutex_unlock(&out);
  }
    
    
    
    

    destructSegment(seg);
    return;
}
void 
cpgTest(char *chr, int start, int stop, double ratio, double p, metseg_t *nfo, double methA, double methB) {
//cpg->chr, cpg->start,cpg->stop,ratio,p);    
    
    if(nfo->threads >1) { 
    pthread_mutex_lock(&out);
  }
    if(nfo->outputList->i >= nfo->outputList->n-1) {
      nfo->outputList->n+=1000000;
      nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
    }
    nfo->outputList->segment_out[nfo->outputList->i].chr = ALLOCMEMORY(NULL, NULL, char, strlen(chr)+1);
    nfo->outputList->segment_out[nfo->outputList->i].chr = strcpy(nfo->outputList->segment_out[nfo->outputList->i].chr,chr);
    nfo->outputList->segment_out[nfo->outputList->i].start = start-1;
    nfo->outputList->segment_out[nfo->outputList->i].stop = stop;
    //nfo->outputList->segment_out[nfo->outputList->i].p = p;
    nfo->outputList->segment_out[nfo->outputList->i].meandiff = ratio;
    nfo->outputList->segment_out[nfo->outputList->i].mwu = p;
    nfo->outputList->segment_out[nfo->outputList->i].length = 1;
    nfo->outputList->segment_out[nfo->outputList->i].methA = methA;
    nfo->outputList->segment_out[nfo->outputList->i].methB = methB;
    
    nfo->outputList->i+=1;
    nfo->outputList->numberTests+=1;
   //unlock if necessary
    if(nfo->threads > 1) {
    pthread_mutex_unlock(&out);
    }
    return;
}
/*-------------------------------- segworker_CpG -----------------------------
 *    
 * @brief for threaded segmentation in single CpG mode
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
void*
segworker_CpG (void *args)
{
    
  metseg_t *t;
  t = (metseg_t*) args;
  cpg_t *cpg = t->cpg;
  
  
  int ua = mannwhitney (cpg->groupA, cpg->noA, cpg->groupB, cpg->noB);
  double p= mannwhitneyPvalue(ua, cpg->noA, cpg->noB, t->MWU, MAXM, MAXN);
  //double p = mannwhitney (cpg->groupA, cpg->noA, cpg->groupB, cpg->noB);
  double ratio = get_meandiff(cpg, cpg->groupA, cpg->noA, cpg->groupB, cpg->noB);
  
  
  
  
  /*
  //output here   
  fprintf(stdout, "#%s\t%d\t%d\t.\t1\t%f\n", 
              cpg->chr, cpg->start-1, 
              cpg->stop,p);
  */
  cpgTest(cpg->chr, cpg->start,cpg->stop,ratio,p,t, cpg->methA,cpg->methB);
  
  
  destructCpg(t->cpg);
  
  pthread_mutex_lock(&cnt);
  schedule[t->threadno] = 0;
  //decrement of idle at creation time
  idle = idle+1;
  pthread_mutex_unlock(&cnt);

  pthread_exit(args);
}
/*-------------------------------- segworker_region --------------------------
 *    
 * @brief for threaded segmentation in region mode
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
void*
segworker_region (void *args)
{
  metseg_t *t;
  t = (metseg_t*) args;
  regionTest(t->seg, t->grpA, t->noA, t->grpB, t->noB, t);
          
  pthread_mutex_lock(&cnt);
  schedule[t->threadno] = 0;
  //decrement of idle at creation time
  idle = idle+1;
  pthread_mutex_unlock(&cnt);

  pthread_exit(args);
}



/*-------------------------------- checkSetNAN ---------------------------------
 *    
 * @brief for checking NaNs in input data
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
 
int
checkSetNAN(stringset_t **csv, double *values){
//check line for NANs and set . to NAN
    int nan=0;
    values[0] = -1;
    values[1] = -1;
    for(int k=2; k < csv[0]->noofstrings; k++) { 
        char *s = csv[0]->strings[k].str;
        values[k] = atof(s);
        if(strcmp(".", s) == 0 || strcmp("-", s) == 0 || strcmp("", s) == 0 || strcmp(" ", s) == 0){
            values[k] = NAN;
        }
        int in=0;
        while (values[k] == values[k] && s[in]) {
          if (isalpha(s[in])) {
              values[k] = NAN;
              break;
          }
//          else printf ("character %c is not alphabetic\n",s[in]);
          in++;
        }  
        if(values[k] != values[k] ) {
            nan+=1;
        }
    }
    return nan;
}
        

   
/*-------------------------------- fillNAN ---------------------------------
 *    
 * @brief for replacing NaNs with betaDist
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
 
int 
fillNAN(double *values, int *grpA, int noA, int *grpB, int noB, metseg_t *nfo) {
  //    fprintf(stderr,"#Ueberhaupt schaetzen\n");
    int na=0;
    int nb=0;
    double *groupA = ALLOCMEMORY(NULL, NULL, double, noA);
    double *groupB = ALLOCMEMORY(NULL, NULL, double, noB);
    double varA;
    double varB;
    double meanA = 0.0;
    double meanB = 0.0;
    int j;

//    fprintf(stdout,"\ngroupA\t");
    j=0;
    for(int i=0; i<noA; i++) 
        if(values[grpA[i]+2] == values[grpA[i]+2]) {
          groupA[j] = values[grpA[i]+2];
 //         fprintf(stdout,"%f\t",groupA[j]);
          na+=1;
          meanA+=groupA[j];
          j++;
    }
 //   fprintf(stdout,"\ngroupB\t");
   
    j=0;
    for(int i=0; i<noB; i++) 
        if(values[grpB[i]+2] == values[grpB[i]+2]) {
          groupB[j] = values[grpB[i]+2];
//          fprintf(stdout,"%f\t",groupB[j]);
          nb+=1;
          meanB+=groupB[j];
          j++;
    }
//    fprintf(stdout,"\n");
   
    if(na<1 || nb<1 || na<nfo->minNoA || nb<nfo->minNoB) {
        FREEMEMORY(NULL, groupA);
        FREEMEMORY(NULL, groupB);
	//        fprintf(stderr,"#REMOVING POSITION CUTOFF\n");
        return 1;
    }
    //    fprintf(stderr,"#NOT REMOVING POSITION CUTOFF\n");
    meanA/=(double)na;
    meanB/=(double)nb;
    if(na == 1) {
        varA=0.000001;
    }
    else {
        varA = var(groupA,na);
    }
        
    if(nb == 1) {
        varB=0.000001;
    }
    else {
        varB = var(groupB,nb);
    }
 //   fprintf(stderr,"#meanA: %f\tmeanB: %f\tvarA: %f\tvarB: %f\n",meanA,meanB,varA,varB);
    
    
    if(meanA < 0.000001)
        meanA = 0.000001;
    if(meanB < 0.000001)
        meanB = 0.000001;
   
 //   fprintf(stdout,"#new As:\n");
    for(int i=0; i<noA; i++) {
      if(isnan(values[grpA[i]+2])) {
          values[grpA[i]+2] = rbeta_mv(meanA, varA);
          while(isnan(values[grpA[i]+2])) {
//            values[grpA[i]+2] = meanA;
            values[grpA[i]+2] = rbeta_mv(meanA, varA);
	    //            fprintf(stderr,"#betaAnot %d\n",i);
          }
          //else {
            //  fprintf(stdout,"#betaA %d\n",i);
       //   }
        }
//        fprintf(stdout,"%f\t",values[grpA[i]+2]);
    }
  //  fprintf(stdout,"#new Bs:\n");
    
    for(int i=0; i<noB; i++) {
      if(isnan(values[grpB[i]+2])) {
          values[grpB[i]+2] = rbeta_mv(meanB, varB);
          while(isnan(values[grpB[i]+2])) {
//            values[grpB[i]+2] = meanB;
            values[grpB[i]+2] = rbeta_mv(meanB, varB);
//	      fprintf(stderr,"#betaBnot %d\n",i);
          }
         // else {
        //      fprintf(stdout,"#betaB %d\n",i);
        //  }
        }
//        fprintf(stdout,"%f\t",values[grpB[i]+2]);
    }
//    fprintf(stdout,"\n");

//   fprintf(stdout,"#A:\t%f / %f\t%f / %f\n",meanA,varA,meanB,varB);
    
    FREEMEMORY(NULL, groupA);
    FREEMEMORY(NULL, groupB);
    return 0;
    
 //   double
//var (double *x, int n)
    
    
}








/*---------------------------- initProgramParams -----------------------------
 *    
 * @brief initialize program parameters
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void
initProgramParams (metseg_t *nfo)
{
  nfo->maxdist=300;
  nfo->mincpgs=10;
  nfo->threads=1;
  nfo->mode=1;
  nfo->nameA = "g1";
  nfo->nameB = "g2";
  nfo->trend = 0.6;
  nfo->minNoA = -1;
  nfo->minNoB = -1;
  nfo->minFactor = 0.8;
  nfo->valley = 0.7;
  nfo->minMethDist = 0.1;
  nfo->MWU = generateMannWhitneyCDFMatrix(MAXM, MAXN);
  //  testMannWhitneyApprox (MAXM, MAXN, nfo->MWU);
  return ;
}

/*---------------------------- initProgramParams2 -----------------------------
 *    
 * @brief initialize program parameters that depend on group sizes
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void
initProgramParams2 (metseg_t *nfo, int noA, int noB)
{
  if(nfo->minNoA<0) {
        nfo->minNoA = ceil((double)noA * nfo->minFactor);
  }
  
  if(nfo->minNoB<0) {
        nfo->minNoB = ceil((double)noB * nfo->minFactor);
  }
  
  return ;
}

/*-------------------------------- variance ---------------------------------
 *    
 * @brief calculate the variance
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
 
double
variance (double *x, int n,double m)
{   
    int i;
    double r, sum=0;

    
    for (i=0; i < n; i++) {
      r = x[i]-m;
      sum += (r*r);
    }

	return sum/n;
}



/*----------------------------------- main -----------------------------------
 *    
 * @brief the main routine
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */


int main(int argc, char** argv) {

  manopt_optionset optset;
  manopt_arg *args;
  metseg_t nfo; 
  metseg_t *th_nfo;
  stringset_t **csv, **bedcsv;
  fileiterator_t *fi, *bedfi;
  unsigned int i, j, k, ln, bedln;
  pthread_t *threads = NULL;
  pthread_attr_t tattr;
  char *bedfile = NULL;


  char **chr = NULL;
  int *grpA = NULL, noA=0;
  int *grpB = NULL, noB=0;
  double **val = NULL;
  int *pos = NULL;
  srand (26061981);
  
  initProgramParams(&nfo);
  
   //we want to detach the threads to automatically have their resources freed
  pthread_attr_init(&tattr); 
  pthread_attr_setdetachstate(&tattr,PTHREAD_CREATE_DETACHED);


  pthread_mutex_init(&out, NULL);
  pthread_mutex_init(&cnt, NULL);

// Options  
  manopt_initoptionset(&optset, argv[0], NULL, 
      "metilene - a tool for fast and sensitive detection of differential DNA methylation\n\nDataInputFile\t\tneeds to be SORTED for chromosomes and genomic positions",
      "Implemented by Frank Juehling and Steve Hoffmann\n  2015-2016 Bioinformatik Leipzig\n",
      version,
      "Please report bugs to [frank,steve]@bioinf.uni-leipzig.de");

  manopt(&optset, REQUINTOPT, 0, 'M', "maxdist", 
      "maximum distance", "<n>", NULL, &nfo.maxdist);
  manopt(&optset, REQUINTOPT, 0, 'm', "mincpgs", 
      "minimum cpgs", "<n>", NULL, &nfo.mincpgs);
  manopt(&optset, REQDBLOPT, 0, 'd', "minMethDiff", 
      "minimum mean methylation difference", "<n>", NULL, &nfo.minMethDist);
  manopt(&optset, REQUINTOPT, 0, 't', "threads", 
      "number of threads", "<n>", NULL, &nfo.threads);
  manopt(&optset, REQUINTOPT, 0, 'f', "mode", 
      "number of method: 1: de-novo, 2: pre-defined regions, 3: DMCs", "<n>", NULL, &nfo.mode);
  manopt(&optset, REQSTRINGOPT, 0, 'a', "groupA", 
      "name of group A", "<string>", NULL, &nfo.nameA);
  manopt(&optset, REQSTRINGOPT, 0, 'b', "groupB", 
      "name of group B", "<string>", NULL, &nfo.nameB);
  manopt(&optset, REQSTRINGOPT, 0, 'B', "bed", 
      "bed-file for mode 2 containing pre-defined regions; needs to be SORTED equally to the DataInputFile", "<string>", NULL, &bedfile);
  manopt(&optset, REQUINTOPT, 0, 'X', "minNoA", 
      "minimal number of values in group A", "<n>", NULL, &nfo.minNoA);
  manopt(&optset, REQUINTOPT, 0, 'Y', "minNoB", 
      "minimal number of values in group B", "<n>", NULL, &nfo.minNoB);
  manopt(&optset, REQDBLOPT, 0, 'v', "valley", 
      "valley filter (0.0 - 1.0)", "<n>", NULL, &nfo.valley);

  args = manopt_getopts(&optset, argc, argv);
  if(args->noofvalues == 1) {
    manopt_help(&optset, "no source file provided.\n");
  }


  fi = initFileIterator(NULL, args->values[1]);

  //handle groups START
  ln = readcsvlines(NULL, fi, '\t', 1, &csv);
  for(k=2; k < csv[0]->noofstrings; k++) {
    if(strstr(csv[0]->strings[k].str, nfo.nameA)) {
      grpA = ALLOCMEMORY(NULL, grpA, int, noA+1);
      grpA[noA] = k-2;
      noA++;
    } else 
      if (strstr(csv[0]->strings[k].str, nfo.nameB)) {
        grpB = ALLOCMEMORY(NULL, grpB, int, noB+1);
        grpB[noB] = k-2;
        noB++;
      }
  }
  initProgramParams2(&nfo, noA, noB);

  fprintf(stderr, "the following ids belong to group A (n=%d):\n", noA);
  for(i=0; i < noA; i++) {
    fprintf(stderr, "%d: column: %d, name:%s\n", i, grpA[i], 
        csv[0]->strings[grpA[i]+2].str);
  }

  fprintf(stderr, "the following ids belong to group B (n=%d):\n", noB);
  for(i=0; i < noB; i++) {
    fprintf(stderr, "%d: column: %d, name:%s\n", i, grpB[i], 
        csv[0]->strings[grpB[i]+2].str);
  }

  destructStringset(NULL, csv[0]);
  FREEMEMORY(NULL, csv);
  
  //init schedules and idle counter
  idle = nfo.threads;
  schedule = ALLOCMEMORY(NULL, NULL, volatile char, nfo.threads);
  for(i=0; i < nfo.threads; i++) {
    schedule[i]=0;
  }
  
  threads = ALLOCMEMORY(space, NULL, pthread_t, nfo.threads);
  th_nfo = ALLOCMEMORY(space, NULL, metseg_t, nfo.threads);
  
  for(i=0; i < nfo.threads; i++) {
    memmove(&th_nfo[i], &nfo, sizeof(metseg_t));
  }
 
  
  
  
  
  
 /* //####################################### REMOVE
  fprintf(stdout,"minFactors\t%d\t%d\n",nfo.minNoA,nfo.minNoB);
  int x = 11;
  int a = ceil(0.8*(double)x);
  fprintf(stdout,"@%d\n",a);
  
        ln = readcsvlines(NULL, fi, '\t', 1, &csv);
        j = 0;
        while(ln) { 
            double *values = ALLOCMEMORY(NULL, NULL, double, csv[0]->noofstrings);
            int nan = checkSetNAN(csv, values);
            if(nan>0) {
                fprintf(stdout,"call fillNAN");
                nan = fillNAN(values, grpA, noA, grpB, noB, &nfo);
            }
            if(nan>0) {
                destructStringset(NULL, csv[0]);
                csv[0] = NULL;
                FREEMEMORY(NULL, csv);
                csv = NULL;
                ln = readcsvlines(NULL, fi, '\t', 1, &csv);
                FREEMEMORY(NULL, values); 
                continue;
            }
            
            
            
            
            fprintf(stdout,"###@@##################@@###\n");    
            
            
            
            destructStringset(NULL, csv[0]);
            csv[0] = NULL;
            FREEMEMORY(NULL, csv);
            csv = NULL;
            ln = readcsvlines(NULL, fi, '\t', 1, &csv);
            FREEMEMORY(NULL, values); 
        }
 nfo.mode=4;
  //####################################### TILL HERE */
  
 








 
  
//###################### SINGLE CpG mode #########################
  if(nfo.mode == 3) {
      
      nfo.outputList = ALLOCMEMORY(NULL, NULL, list_out, 1);
      nfo.outputList->segment_out = ALLOCMEMORY(NULL, NULL, segment_out, 1);
      nfo.outputList->n=1000000;
      nfo.outputList->i=0;
      nfo.outputList->numberTests=0;
      nfo.outputList->segment_out = ALLOCMEMORY(NULL, NULL, segment_out, nfo.outputList->n);
    
      
      
        ln = readcsvlines(NULL, fi, '\t', 1, &csv);
        j = 0;
        while(ln) {
//check missing numbers            
            double *values = ALLOCMEMORY(NULL, NULL, double, csv[0]->noofstrings);
            int nan = checkSetNAN(csv, values);
            if(nan>0) {
          //      fprintf(stderr,"call fillNAN");
                nan = fillNAN(values, grpA, noA, grpB, noB, &nfo);
         //       fprintf(stderr,"...done\n");
                
            }
            if(nan>0) {
                destructStringset(NULL, csv[0]);
                csv[0] = NULL;
                FREEMEMORY(NULL, csv);
                csv = NULL;
                ln = readcsvlines(NULL, fi, '\t', 1, &csv);
                FREEMEMORY(NULL, values); 
                continue;
            }

            cpg_t *cpg = ALLOCMEMORY(NULL, NULL, cpg_t, 1);
            double *groupA = ALLOCMEMORY(NULL, NULL, double, noA);
            double *groupB = ALLOCMEMORY(NULL, NULL, double, noB);
            
            for(i=0; i<noA; i++) { 
//                  groupA[i] = atof(csv[0]->strings[grpA[i]+2].str);
                    groupA[i] = values[grpA[i]+2];
            }
            for(i=0; i<noB; i++) { 
//                  groupB[i] = atof(csv[0]->strings[grpB[i]+2].str);
                    groupB[i] = values[grpB[i]+2];
            }
            cpg->groupA=groupA;
            cpg->groupB=groupB;
            cpg->chr = ALLOCMEMORY(NULL, NULL, char, csv[0]->strings[0].len+1);
            strcpy(cpg->chr, csv[0]->strings[0].str);
            cpg->noA=noA;
            cpg->noB=noB;
            cpg->start=atoi(csv[0]->strings[1].str);
            cpg->stop=atoi(csv[0]->strings[1].str);
            
            if(nfo.threads > 1) { 
                //wait for a free thread
                while(idle == 0);
                //look for the free slot
                for(i=0; i < nfo.threads; i++) {
                  if(schedule[i] == 0) break;
                }
                //this must always hold because of idle variable
                assert(i < nfo.threads);
     //           fprintf(stderr, "starting thread %d\n", i);
                //decrement the idle and set the schedule
                pthread_mutex_lock(&cnt);
                idle--;
                schedule[i]=1;
                pthread_mutex_unlock(&cnt);
                //assign task
                th_nfo[i].threadno = i;
                th_nfo[i].cpg=cpg;
                th_nfo[i].outputList = nfo.outputList;
                fprintf(stderr, "CpG testing %s-[%d]\n", cpg->chr, cpg->start);
                //create the thread (detached!)
                pthread_create(&threads[i], &tattr, segworker_CpG, &th_nfo[i]);
                //now we must make sure that each thread keeps his own chunk
                //of the input data, thus the three arrays are simply set to NULL
                //the thread is going to take care of the deallocation
                cpg = NULL;

              } else { 
                fprintf(stderr, "CpG testing %s-[%d]\n", cpg->chr, cpg->start);
                int ua = mannwhitney (cpg->groupA, noA, cpg->groupB, noB);
		double p= mannwhitneyPvalue(ua, noA, noB, nfo.MWU, MAXM, MAXN);
		double ratio = get_meandiff(cpg, cpg->groupA, noA, cpg->groupB, noB);
            //    fprintf(stdout, "#%s\t%d\t%d\t.\t%.2g\t%f\n",cpg->chr, cpg->start,cpg->stop,ratio,p);
                
                
                cpgTest(cpg->chr, cpg->start,cpg->stop,ratio,p,&nfo, cpg->methA, cpg->methB);
                destructCpg(cpg);
              }
            
            
            
            destructStringset(NULL, csv[0]);
            csv[0] = NULL;
            FREEMEMORY(NULL, csv);
            csv = NULL;
            FREEMEMORY(NULL, values); 
           ln =readcsvlines(NULL, fi, '\t', 1, &csv);
    }
        
        
        

        
    
  }


  


  
  
  
  



  
//###################### DEFINED REGIONS mode#####################
  if(nfo.mode == 2) {
      fprintf(stderr, "Mode 2 -- pre-defined regions\n");
      nfo.outputList = ALLOCMEMORY(NULL, NULL, list_out, 1);
      nfo.outputList->segment_out = ALLOCMEMORY(NULL, NULL, segment_out, 1);
      nfo.outputList->n=1000000;
      nfo.outputList->i=0;
      nfo.outputList->numberTests=0;
      nfo.outputList->segment_out = ALLOCMEMORY(NULL, NULL, segment_out, nfo.outputList->n);
    
      
      
      
    //  exit(EXIT_SUCCESS);
//init segmentset of active regions
      bedfi = initFileIterator(NULL, bedfile);
      segmentset_t *set = ALLOCMEMORY(NULL, NULL, segmentset_t, 1);
      initSegmentSet(set);
//first region      
      bedln = readcsvlines(NULL, bedfi, '\t', 1, &bedcsv);
      set->nextchr = ALLOCMEMORY(NULL, NULL, char, bedcsv[0]->strings[0].len+1);
      strcpy(set->nextchr, bedcsv[0]->strings[0].str);
      set->nextstart = atoi(bedcsv[0]->strings[1].str)+1;
      
      set->chr = ALLOCMEMORY(NULL, NULL, char, bedcsv[0]->strings[0].len+1);
      strcpy(set->chr, bedcsv[0]->strings[0].str);
      
      ln = readcsvlines(NULL, fi, '\t', 1, &csv);
      int l=-1;
      while(ln) {
          l++;
//add missing values          
          double *values = ALLOCMEMORY(NULL, NULL, double, csv[0]->noofstrings);
            int nan = checkSetNAN(csv, values);
            if(nan>0) {
  //              fprintf(stdout,"call fillNAN");
                nan = fillNAN(values, grpA, noA, grpB, noB, &nfo);
            }
            if(nan>0) {
                destructStringset(NULL, csv[0]);
                csv[0] = NULL;
                FREEMEMORY(NULL, csv);
                csv = NULL;
                ln = readcsvlines(NULL, fi, '\t', 1, &csv);
                FREEMEMORY(NULL, values); 
                continue;
            }
         
          
          int pos = atoi(csv[0]->strings[1].str);
          
          
      //   fprintf(stdout,"########## %s %d (currChrom %s ,firststop %d, nextChr %s, nextStart %d)\n",csv[0]->strings[0].str,pos,set->chr,set->firststop, set->nextchr,set->nextstart);
          
//remove filled segments if
//Size of set >0 AND (currentChrNotNULL  OR  curr.Positon > FirstStopInSet OR  ChromosomeChangeForCpGInput)
          if(set->n>0 && (set->chr == NULL || pos > set->firststop || (strcmp(set->chr, csv[0]->strings[0].str) != 0))) {

   //           if((strcmp(set->chr, csv[0]->strings[0].str) != 0))
   //               fprintf(stdout,"CHROMCHANGE\n");
              
              segment_t *seg = set->head;
              set->firststop=-1;
              while(seg) {
                  segment_t *tmp = seg->next;
                  if(set->n>0 && ( set->chr == NULL || (strcmp(set->chr, csv[0]->strings[0].str) != 0) || seg->stop<pos)) {
       //                 fprintf(stdout,"@@@@@@@@@@Removing seg %s:%d-%d next%d parent%d\n",seg->chr,seg->start,seg->stop,seg->next == NULL,seg->parent == NULL);
                        removeThisSegmentFromSet(set,seg);
                        if(nfo.threads > 1) { 
                            //wait for a free thread
                            while(idle == 0);
                            //look for the free slot
                            for(i=0; i < nfo.threads; i++) {
                              if(schedule[i] == 0) break;
                            }
                            //this must always hold because of idle variable
                            assert(i < nfo.threads);
//                            fprintf(stderr, "starting thread %d\n", i);
                            //decrement the idle and set the schedule
                            pthread_mutex_lock(&cnt);
                            idle--;
                            schedule[i]=1;
                            pthread_mutex_unlock(&cnt);
                            //assign task
                            th_nfo[i].seg = seg;
                            th_nfo[i].grpA = grpA;
                            th_nfo[i].grpB = grpB;
                            th_nfo[i].noA = noA;
                            th_nfo[i].noB = noB;
                            th_nfo[i].threadno = i;
                            th_nfo[i].outputList = nfo.outputList;
          
                            fprintf(stderr, "region testing %s-[%d,%d]\n", seg->chr, seg->start, seg->stop);
                            //create the thread (detached!)
                            pthread_create(&threads[i], &tattr, segworker_region, &th_nfo[i]);
                            //now we must make sure that each thread keeps his own chunk
                            //of the input data, thus the three arrays are simply set to NULL
                            //the thread is going to take care of the deallocation
                            seg = NULL;

                          } else { 
                            fprintf(stderr, "region testing %s-[%d,%d]\n", seg->chr, seg->start, seg->stop);
                            regionTest(seg, grpA, noA, grpB, noB, &nfo);
                              
                          }
                  }
                  else {
                      if(set->firststop == -1) { set->firststop = seg->stop; }
                      else      { set->firststop = MIN(seg->stop, set->firststop); }
                  }
                  seg = tmp;
              }
              if(set->n < 1) {
                  FREEMEMORY(NULL, set->chr); 
                  set->chr=NULL;
              }
          }
          
            if((!set->chr) && set->nextchr && (strcmp(set->nextchr, csv[0]->strings[0].str) != 0)) {
                 destructStringset(NULL, csv[0]);
                 csv[0] = NULL;
                 FREEMEMORY(NULL, csv);
                 csv = NULL;
                 ln = readcsvlines(NULL, fi, '\t', 1, &csv);
                 FREEMEMORY(NULL, values); 
                 continue;
             }
          
          
          
//add new segments 
       //   if(set->nextchr)
       //           fprintf(stdout,"while %d\t%d\t%d\n",bedln, pos>=set->nextstart ,(strcmp(set->chr,set->nextchr) == 0));
          while(bedln && pos>=set->nextstart && (!set->nextchr  || !set->chr || (strcmp(set->chr,set->nextchr) == 0)) ) {
              segment_t *tmp = addNewSegmentToSet(set);
              initSegment(tmp);
              tmp->chr = ALLOCMEMORY(NULL, NULL, char, bedcsv[0]->strings[0].len+1);
              strcpy(tmp->chr, bedcsv[0]->strings[0].str);
              tmp->start = atoi(bedcsv[0]->strings[1].str)+1;
              tmp->stop = atoi(bedcsv[0]->strings[2].str);
              
              if(tmp->stop < set->firststop) {
                  set->firststop = tmp->stop;
              }
              
              FREEMEMORY(NULL, set->chr); 
              if(set->n > 0) {
                  set->chr = ALLOCMEMORY(NULL, NULL, char, bedcsv[0]->strings[0].len+1);
                  strcpy(set->chr, bedcsv[0]->strings[0].str);
                  
                  
              }
              
              
         //read in next region
              destructStringset(NULL, bedcsv[0]);
              bedcsv[0] = NULL;
              FREEMEMORY(NULL, bedcsv);
              bedcsv = NULL;
              
              bedln = readcsvlines(NULL, bedfi, '\t', 1, &bedcsv); 

       //       if(bedln)
       //               fprintf(stdout,"seg %s:%s-%s\n",bedcsv[0]->strings[0].str,bedcsv[0]->strings[1].str,bedcsv[0]->strings[2].str);

              FREEMEMORY(NULL, set->nextchr);
              if(bedln) {
                set->nextchr = ALLOCMEMORY(NULL, NULL, char, bedcsv[0]->strings[0].len+1);
                strcpy(set->nextchr, bedcsv[0]->strings[0].str);
                set->nextstart = atoi(bedcsv[0]->strings[1].str)+1;
              }
              else {
                  set->nextstart = -1;
              }
              
          }
      
          
//add CpGs to regions          
        //add CpGs to regions
         
          segment_t *seg = set->head;
        //  if(seg && set->firststop==-1)
        //        { set->firststop = seg->stop; }
        //  int Notbreaking=1;
/*
          if(seg && seg->chr && csv[0] && (strcmp(seg->chr,csv[0]->strings[0].str) || ( (!strcmp(seg->chr,csv[0]->strings[0].str)) && seg->start > atoi(csv[0]->strings[1].str)))) {
  //          if(seg && seg->chr && csv[0] && (strcmp(seg->chr,csv[0]->strings[0].str))) {
              
              
              int z=0;
              while(csv[0] && seg && seg->chr && (strcmp(seg->chr,csv[0]->strings[0].str) || ( (!strcmp(seg->chr,csv[0]->strings[0].str)) && seg->start > atoi(csv[0]->strings[1].str)))){
    //          while(csv[0] && (strcmp(seg->chr,csv[0]->strings[0].str ))){
                  z++;
                destructStringset(NULL, csv[0]);
                csv[0] = NULL;
                FREEMEMORY(NULL, csv);
                csv = NULL;
                FREEMEMORY(NULL, values); 
                ln = readcsvlines(NULL, fi, '\t', 1, &csv);
              }
              fprintf(stdout,"removed %d lines\n",z);
              fprintf(stdout,"continue: segChr=%s csvChr=%s segStart=%d csvStart=%d\n",seg->chr,csv[0]->strings[0].str,seg->start,atoi(csv[0]->strings[1].str)+1);
              Notbreaking = 0;
          }
          */
        //  pos = atoi(csv[0]->strings[1].str)+1;
       //   fprintf(stdout,"########## %s %d (currChrom %s ,firststop %d, nextChr %s, nextStart %d)\n",csv[0]->strings[0].str,pos,set->chr,set->firststop, set->nextchr,set->nextstart);
          
          
          
          while(seg)  {
//              fprintf(stdout,"@@@@@@@@@@@@@@@@@@@@@adding seqs now %s %s\n",seg->chr,csv[0]->strings[0].str);
              if(!strcmp(seg->chr,csv[0]->strings[0].str)) {
//                fprintf(stderr,"#Adding CpG %s:%s to region %s:%d-%d\n",csv[0]->strings[0].str,csv[0]->strings[1].str,seg->chr,seg->start,seg->stop);
                seg->pos = ALLOCMEMORY(NULL, seg->pos, int,    seg->n+1); //index
                seg->value = ALLOCMEMORY(NULL, seg->value, double*, seg->n+1); //cpgs  

                seg->pos[seg->n] = atoi(csv[0]->strings[1].str);
                seg->value[seg->n] = ALLOCMEMORY(NULL, NULL, double, csv[0]->noofstrings);
                for(k=2; k < csv[0]->noofstrings; k++) { 
//                  seg->value[seg->n][k-2] = atof(csv[0]->strings[k].str);
                    seg->value[seg->n][k-2] = values[k];
                }
                seg->n++;
              }
          seg = seg->next;    
          }
       //   if(Notbreaking==1){
            destructStringset(NULL, csv[0]);
            csv[0] = NULL;
            FREEMEMORY(NULL, csv);
            csv = NULL;
            FREEMEMORY(NULL, values); 
            ln = readcsvlines(NULL, fi, '\t', 1, &csv);
    //      }      
      
       }
      
//remaining regions to test 
      segment_t *seg = set->head;
      while(seg) {
                  segment_t *tmp = seg->next;
        //          fprintf(stderr,"@@@@@@@@@@Removing seg %s:%d-%d\n",seg->chr,seg->start,seg->stop);
                  if(nfo.threads > 1) { 
                            //wait for a free thread
                            while(idle == 0);
                            //look for the free slot
                            for(i=0; i < nfo.threads; i++) {
                              if(schedule[i] == 0) break;
                            }
                            //this must always hold because of idle variable
                            assert(i < nfo.threads);
                            fprintf(stderr, "starting thread %d\n", i);
                            //decrement the idle and set the schedule
                            pthread_mutex_lock(&cnt);
                            idle--;
                            schedule[i]=1;
                            pthread_mutex_unlock(&cnt);
                            //assign task
                            th_nfo[i].seg = seg;
                            th_nfo[i].grpA = grpA;
                            th_nfo[i].grpB = grpB;
                            th_nfo[i].noA = noA;
                            th_nfo[i].noB = noB;
                            th_nfo[i].threadno = i;
                            
                            fprintf(stderr, "region testing %s-[%d,%d]\n", seg->chr, seg->start, seg->stop);
                            //create the thread (detached!)
                            pthread_create(&threads[i], &tattr, segworker_region, &th_nfo[i]);
                            //now we must make sure that each thread keeps his own chunk
                            //of the input data, thus the three arrays are simply set to NULL
                            //the thread is going to take care of the deallocation
                            seg = NULL;

                          } else { 
    //                        fprintf(stderr, "region testing %s-[%d,%d]\n", seg->chr, seg->start, seg->stop);
                            regionTest(seg, grpA, noA, grpB, noB, &nfo);
                              
                          }
                  seg = tmp;
              }
    //  if(seg)
              destructSegmentSet(set);
              
          
              
  }
  

  
  
  
  
  
  
  
  
  
  

  
 















//###################### SEGMENTER (main) mode ###########################
  if(nfo.mode == 1) {
   //   fprintf(stderr,"#MODE2\n");
    nfo.outputList = ALLOCMEMORY(NULL, NULL, list_out, 1);
    nfo.outputList->segment_out = ALLOCMEMORY(NULL, NULL, segment_out, 1);
    nfo.outputList->n=1000000;
    nfo.outputList->i=0;
    nfo.outputList->numberTests=0;
    nfo.outputList->segment_out = ALLOCMEMORY(NULL, NULL, segment_out, nfo.outputList->n);
    //fprintf(stderr,"output->n: %d\n",nfo.outputList->n);
      
      
    ln = readcsvlines(NULL, fi, '\t', 1, &csv);
    j = 0;
    while(ln) { 
        //fprintf(stderr,"#new LINE\n");
//check missing numbers            
        double *values = ALLOCMEMORY(NULL, NULL, double, csv[0]->noofstrings);
        int nan = checkSetNAN(csv, values);
        if(nan>0) {
    //     fprintf(stderr,"#call fillNAN\n");
            nan = fillNAN(values, grpA, noA, grpB, noB, &nfo);
    //      fprintf(stderr,"#...done\n");
        }
     //   fprintf(stdout,"#LINES INPUT\n");
        if(nan>0) {
     //       fprintf(stdout,"#REMOVING LINE\n");
            destructStringset(NULL, csv[0]);
            csv[0] = NULL;
            FREEMEMORY(NULL, csv);
            csv = NULL;
            ln = readcsvlines(NULL, fi, '\t', 1, &csv);
            FREEMEMORY(NULL, values); 
            continue;
        }
        else {
    //            fprintf(stdout,"#LINE OKAY \n");
        }
      char *x = my_strdup(csv[0]->strings[0].str);
      int y = atoi(csv[0]->strings[1].str);

      if(j > 0 && (strcmp(x, chr[j-1]) || y > pos[j-1] + nfo.maxdist)) {

        if(nfo.threads > 1) { 
          //wait for a free thread
          while(idle == 0);
          //look for the free slot
          for(i=0; i < nfo.threads; i++) {
            if(schedule[i] == 0) break;
          }
          //this must always hold because of idle variable
          assert(i < nfo.threads);
  //        fprintf(stderr, "starting thread %d\n", i);
          //decrement the idle and set the schedule
          pthread_mutex_lock(&cnt);
          idle--;
          schedule[i]=1;
          pthread_mutex_unlock(&cnt);
          //assign task
          th_nfo[i].chr = chr;
          th_nfo[i].pos = pos;
          th_nfo[i].value = val;
          th_nfo[i].n = j;
          th_nfo[i].grpA = grpA;
          th_nfo[i].grpB = grpB;
          th_nfo[i].noA = noA;
          th_nfo[i].noB = noB;
          th_nfo[i].threadno = i;
          th_nfo[i].outputList = nfo.outputList;
          
          
          fprintf(stderr, "segmenting %s-[%d,%d], %d CpGs\n", chr[0], pos[0], pos[j-1],j);
          //create the thread (detached!)
          pthread_create(&threads[i], &tattr, segworker, &th_nfo[i]);
          //now we must make sure that each thread keeps his own chunk
          //of the input data, thus the three arrays are simply set to NULL
          //the thread is going to take care of the deallocation
          chr = NULL;
          pos = NULL;
          val = NULL;

        } else { 
          fprintf(stderr, "Segmenting %s-[%d,%d], %d CpGs\n", chr[0], pos[0], pos[j-1],j);
          segmentation(chr, pos, val, j, grpA, noA, grpB, noB, &nfo);
          for(i=0; i < j; i++) { 
            FREEMEMORY(NULL, chr[i]);
            FREEMEMORY(NULL, val[i]);
          }   
        }

        j = 0;
      } 

      chr = ALLOCMEMORY(NULL, chr, char*,   j+1); //chr
      pos = ALLOCMEMORY(NULL, pos, int,    j+1); //index
      val = ALLOCMEMORY(NULL, val, double*, j+1); //cpgs  

      chr[j] = x;
      pos[j] = y;
      val[j] = ALLOCMEMORY(NULL, NULL, double, csv[0]->noofstrings);

      for(k=2; k < csv[0]->noofstrings; k++) { 
//        val[j][k-2] = atof(csv[0]->strings[k].str);
        val[j][k-2] = values[k];
      }


      j+=1;
      destructStringset(NULL, csv[0]);
      csv[0] = NULL;
      FREEMEMORY(NULL, csv);
      csv = NULL;
      FREEMEMORY(NULL, values); 
      ln =readcsvlines(NULL, fi, '\t', 1, &csv);
    } 
  
    fprintf(stderr, "segmenting %s-[%d,%d], %d CpGs \n", chr[0], pos[0], pos[j-1],j);
    segmentation(chr, pos, val, j, grpA, noA, grpB, noB, &nfo);
    for(i=0; i < j; i++) { 
        FREEMEMORY(NULL, chr[i]);
        FREEMEMORY(NULL, val[i]);
  }
}
  
  
  
  
  
  
  
  
  
  
  //wait for all threads to terminate
  while(idle != nfo.threads);

  if(nfo.mode == 1 || nfo.mode == 2) {
    double nps = 0+nfo.outputList->numberTests;
//    double nps = 0+nfo.outputList->i;
//    fprintf(stderr, "NUMBER TESTS: %d\n",nfo.outputList->numberTests);
    for(int i=0;i<nfo.outputList->i;i++){
        double q = nfo.outputList->segment_out[i].p*nps;
        if(q>1) { q=1.0; }
        if(nfo.outputList->segment_out[i].meandiff >= nfo.minMethDist || nfo.outputList->segment_out[i].meandiff <= -1* nfo.minMethDist) {
        fprintf(stdout, "%s\t%d\t%d\t%.2g\t%f\t%d\t%.2g\t%.2g\t%.5g\t%.5g\n", 
                nfo.outputList->segment_out[i].chr,
                nfo.outputList->segment_out[i].start,
                nfo.outputList->segment_out[i].stop,
//                nfo.outputList->segment_out[i].p,
                q,
                nfo.outputList->segment_out[i].meandiff,
                nfo.outputList->segment_out[i].length,
                nfo.outputList->segment_out[i].mwu,
//                q);
                nfo.outputList->segment_out[i].p,
                nfo.outputList->segment_out[i].methA,
                nfo.outputList->segment_out[i].methB);
        }
    }
  }
    
  if(nfo.mode == 3) {
              double nps = 0+nfo.outputList->numberTests;         
  //  fprintf(stderr, "NUMBER TESTS: %d\n",nfo.outputList->numberTests);
    for(int i=0;i<nfo.outputList->i;i++){
        double q = nfo.outputList->segment_out[i].mwu*nps;
        if(q>1) { q=1.0; }
        if(nfo.outputList->segment_out[i].meandiff >= nfo.minMethDist || nfo.outputList->segment_out[i].meandiff <= -1* nfo.minMethDist) {
        fprintf(stdout, "%s\t%d\t%d\t%.2g\t%f\t%d\t%.2g\t.\t%.5g\t%.5g\n", 
                nfo.outputList->segment_out[i].chr,
                nfo.outputList->segment_out[i].start,
                nfo.outputList->segment_out[i].stop,
//                nfo.outputList->segment_out[i].p,
                q,
                nfo.outputList->segment_out[i].meandiff,
                nfo.outputList->segment_out[i].length,
                nfo.outputList->segment_out[i].mwu,
                nfo.outputList->segment_out[i].methA,
                nfo.outputList->segment_out[i].methB
//                q);
 //               nfo.outputList->segment_out[i].p
                );
        }
    }
  }
 
  fflush(stdout); 
  
  
  FREEMEMORY(NULL, chr);
  FREEMEMORY(NULL, pos);
  FREEMEMORY(NULL, val);
  FREEMEMORY(NULL, grpA);
  FREEMEMORY(NULL, grpB);
  
  closeFileIterator(NULL, fi);
  FREEMEMORY(NULL, fi);

  if(csv && csv[0]) destructStringset(NULL, csv[0]);
  if(csv) FREEMEMORY(NULL, csv);

  pthread_attr_destroy(&tattr);  
  FREEMEMORY(NULL, threads);
  FREEMEMORY(NULL, th_nfo);
  free((void*)schedule);

  destructMannWhitneyCDFMatrix(nfo.MWU, MAXM, MAXN);
  manopt_destructoptionset(&optset);
  manopt_destructarg(args);
  FREEMEMORY(NULL, args);

  exit(EXIT_SUCCESS);
}

