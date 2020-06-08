/***********************
*Generic EM Model Header
************************/
#ifndef _GEM_GEMCORE_H
#define _GEM_GEMCORE_H

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <limits>         // for numeric_limits<flt>::infinity()
#include <tr1/unordered_map>
#include "mydebug.hpp"
#include "mmio.h"

using namespace std;
//Forward Declaration for std::tr1::unordered_map, for compatibility with SWIG, see tTR in test
//swig notes:
//swig has problem wrapping nested typenames it doesn't recognoize, adhoc solved by obj.disown()
//a forward declaration is needed for swig to know the type
namespace std{
	namespace tr1{
		template <typename T1, typename T2, typename T3, typename T4, typename T5> class unordered_map;
	}
}

/* TYPES */
//in most architectures, short=16bit, long=32bit, longlong=64bit, float=32bit, double=64bit
//Backend Matrix Type is Preselected by tBE.h in testing suite
//Best Performance Currently is std::tr1::unordered_map
//Real Backend Used sometimes is the inherited with interface: B.at(i).at(j) to read and B[i][j]= to assign
//Customizd Numerical Types and Constants
typedef unsigned long iint;		   // for sparse matrix indecies, up to 2^32 = 4294967296L, or namely index int
typedef float flt; 				       // to save memory use float, otherwise use double
typedef numeric_limits<flt> fnl;
extern const flt dft_rt_pflt;                                // to use with P containers
extern const flt dft_rt_lflt;                                // to use with P containers
//extern flt dft_rt_pflt;                                // to use with P containers
//extern flt dft_rt_lflt;                                // to use with P containers
//flt dft_rt_lflt = -numeric_limits<flt>::infinity();

//flt dft_rt_lflt = -fnl::infinity();           // to use with LogP containers
//Derived Complex Types
typedef vector<flt> VF;
typedef vector<VF> MF;
typedef vector<iint> VII;
//Extended tr1::unordered_map template
//TODO sep tcdc and tcdf: template <typename K, typename V, V& > class tr1_unordered_map, see tHC.C; 
template <typename K, typename V, V const& dft_rt_value> 
class tr1_unordered_map : public tr1::unordered_map<K, V> {
public:
  //const static V dft_rt_v = dft_rt_v;  //more abstraction of a sparse matrix behavior to be summarized
  //static const V& dft_rt_v() { return dft_rt_value; };
  //V dft_rt_v; 
  V const& at( K k ) const {
    if ( tr1::unordered_map<K,V>::find(k) != tr1::unordered_map<K,V>::end() )
      return tr1::unordered_map<K,V>::find(k)->second;
    else
      return dft_rt_value;
  };
  //tr1_unordered_map(): tr1::unordered_map<K, V>(){ cout<<"tr1_unordered_map init"<<endl; };
}; //Extends tr1::unordered_map to fullfil interface needs
typedef tr1_unordered_map<iint, flt, dft_rt_lflt> TR1_UM_L;
extern const TR1_UM_L dft_rt_TR1_UM_L;
typedef tr1_unordered_map<iint, TR1_UM_L, dft_rt_TR1_UM_L> TR1_UUM_L; //candidate
typedef const tr1_unordered_map<iint, flt, dft_rt_pflt> TR1_UM_P;
extern TR1_UM_P dft_rt_TR1_UM_P;
typedef tr1_unordered_map<iint, TR1_UM_P, dft_rt_TR1_UM_P> TR1_UUM_P; //candidate

/* FUNCTIONS */
flt vdif( VF&, VF& ); //fdc of vdif( VF&, VF& )
flt dmax( VF& a, VF& b); //fdc of dmax( VF&, VF& )

/* TEMPLATES */
//template <class V, class Vit> void showvec( V& ); //tfdc of showvec<T,Tit>()

//TODO: template <class M, class Mit> void showmat( M& ); //tfdc of showmat<M,Mit>()
//template <class M1> class MEM; //tcdc of MEM<M1>, base class for Mixture Expectatin Maximization(MEM)
template <class M1> // M1: only accepts FSM based on hash map
class MEM {
  public:
	  iint N;          // # of reads
    iint M;          // # of genomes
	  flt R;           // current residual, in P scale
    flt logL;        // loglikelihood, logP scale, base e
    flt plogL;       // previous log likelihood, logP scale, base e
    VF V;            // current parameters, always in P scale
    VF pV;           // previous parameters
    VII Sample;      // bootstrap sequence of samples, original [1,N], change by .bootstrap_sample(), reset by .original_sample()
    M1 P;            // posterior probability matrix, always in logP scale
    MEM( iint n, iint m, VF v ) : N(n), M(m), R(fnl::infinity()), logL(-fnl::infinity()), plogL(-fnl::infinity()), V(v), pV(v){
      Sample = VII(N,0);
      original_sample();
      cerr<<"MEM init in original order!"<<endl; 
    }; 
    ~MEM() { };
    void get_P() { P=M1(); cerr<<"P for hash_based matrix"<<endl; }
    void clr_P() { /*cout<<"check mem now"; cin.get();*/ P.clear(); /*cerr<<"P cleared, for both hash and vector based matrix, check mem again"; cin.get();*/ }
    void get_P( M1* p ) { P=*p; cerr<<"P for matrix need external allocation"<<endl; }
    virtual void mstep(){ cerr<<"virtual mstep"<<endl; }; // delegate for model-dependent m step
    virtual void estep(){ cerr<<"virtual estep"<<endl; }; // delegate for model-dependent e step
    virtual void init_V( char i_mode ) { cerr<<"virtual init_V"<<endl; }; //delegate for model-dependent parameter initialization
    virtual void loglikelihood(){ cerr<<"virtual loglikelihood"<<endl; }; //delegate for model-dependnet likelihood function
    void set_R( flt r ) { R = r; }; // set residual to R
    int solve( flt c, int n, char c_mode ); // solve at most n times to error c
    MF bootstrap( int b, flt c, int n, char i_mode, char c_mode ); // call solve B times with B bootstrap samples
    void bootstrap_sample() { for (iint i=0; i<N; i++) Sample[i] = int( (double(rand())/RAND_MAX)*N ); }; // from TCPL RandInt
    void original_sample() { for (iint i=0; i<N; i++) Sample[i] = i; };
}; //tcdf of MEM<M1>, base class for Mixture Expectatin Maximization

//template <class M1> extern MF MEM<M1>::bootstrap( int,  ); //tfdc of MEM<M1>::bootstrap(...)
template <class M1>
MF MEM<M1>::bootstrap( int b, flt c, int n, char i_mode = 'M', char c_mode = 'U' ){
  MF bootstrap_estimates = MF();
  for( int i=0; i<b; i++ ){
    bootstrap_sample();
    //showvec<VII>( Sample, N ); cerr<<endl;
    init_V( i_mode );
    //showvec<VF>( V, M ); cerr<<endl;
    solve( c, n, c_mode );
    //showvec<VF>( V, M ); cerr<<endl;
    bootstrap_estimates.push_back( V );
  }
  original_sample();
  return bootstrap_estimates;
};
    
//template <class M1> extern int MEM<M1>::solve( flt, int, char ); //tfdc of MEM<M1>::solve(...)
template <class M1> 
int MEM<M1>::solve( flt c, int n, char c_mode = 'L' ){
  int s(1);
  set_R(fnl::infinity());     // reset residue
  clr_P();                    // clear P matrix
  while( (R > c) && (s <= n) ){
    pV=V; estep(); mstep(); // execute a round of EM
    if ( c_mode == 'L' ){ // L for likelihood convergence
      plogL=logL; loglikelihood(); R=abs(exp(plogL)-exp(logL));
    }
    if ( c_mode == 'P' ) // P for parameter convergence
      R=vdif(pV,V);
    if ( c_mode == 'U' ) // U for uniform convergence
      R=dmax(pV,V);
    s++;
    //cerr<<"control: c="<<c<<",s="<<s<<",R="<<R;
    //if( R>c && s<=n ) cerr<<" => continue"<<endl;
    //else cerr<<" => stop"<<endl;
  }
  if ( c_mode != 'L' ) loglikelihood();  //calculate likelihood at last step
  //cerr<<"final test_mem.P="; showmat<TR1_UUM_L>(P,N,M); cerr<<endl;
  return s;
}; //tfdf of MEM<M1>::solve( flt, int, char ), the solver for MEM

//template <class M1, class M2> class MONO_MEM : MEM<M1>; //tcdc of inherited MONO_MEM<M1,M2>
template <class M1, class M2> //M2: ideally accepts NSM, HSM and FSM
class MONO_MEM : public MEM<M1> {
public:
  iint EN;         // # of effective reads
  iint AN;         // # of ambiguous reads
  iint MN;         // # of mapped reads
	M2 Prg;  // Prg[i][j] = prob of r(i) from g(j), always in logP scale
	MONO_MEM( iint n, iint m, VF v ) : MEM<M1>(n, m, v), EN(0), AN(0), MN(0) { cerr<<"MONO_MEM init"<<endl; };
	~MONO_MEM() { };
  void get_Prg() { Prg=M2(); /*cerr<<"Prg for hash_based matrix"<<endl;*/ }
  void clr_Prg() { Prg.clear(); /*cerr<<"Prg cleared, for both hash and vector based matrix"<<endl;*/ }
  void get_Prg( M2* prg ) { Prg=*prg; cerr<<"Prg for matrix need external allocation"<<endl; }
	void loglikelihood(); // default likelihood calculation function
	void init_Prg( iint, string );  // initiate a Prg from existent matrix file
  void init_ENANMN();
  void rand_V();                  //init_V by 'R'andom
  void moment_V();                //init_V by 'M'oment
	void init_V( char);             //init_V initialization function
  void estep();
  void mstep();
}; //tcdf of inherited MONO_MEM<M1,M2>

//template <class M1, class M2> extern void MONO_MEM<M1, M2>::loglikelihood() //tfdc of MONO_MEM<M1,M2>::loglikelihood()
template <class M1, class M2>
void MONO_MEM<M1, M2>::init_V( char i_mode = 'M' ){ 
  switch (i_mode) {
    case 'M': moment_V(); break;
    case 'R': rand_V(); break;
    default: moment_V();
  }
};

//template <class M1, class M2> extern void MONO_MEM<M1, M2>::loglikelihood() //tfdc of MONO_MEM<M1,M2>::loglikelihood()
template <class M1, class M2>
void MONO_MEM<M1, M2>::loglikelihood(){ 
  //cerr<<"calculating loglikelihood..."<<endl;
  double zij, logp;
  MEM<M1>::logL = 0.;
  for (iint i=0; i< MEM<M1>::N; i++)
    for (iint j=0; j< MEM<M1>::M; j++){
      zij = exp( MEM<M1>::P.at(i).at(j) ); 
      logp = Prg.at(MEM<M1>::Sample[i]).at(j);
      if ( logp != -fnl::infinity() && MEM<M1>::V[j] != 0 ) // if Prg[i][j] = 0 or V[j] = 0, implies zij=0, as show in estep
        MEM<M1>::logL += zij * ( log( MEM<M1>::V[j] ) + logp ); 
    }
}; //tcdf of MONO_MEM<M1,M2>::loglikelihood()

//template <class M1, class M2> extern void MONO_MEM<M1, M2>::estep() //tfdc of MONO_MEM<M1, M2>::estep()
template <class M1, class M2>
void MONO_MEM<M1, M2>::estep(){
  //cerr<<"in estep..."<<endl; cin.get();
  for(iint i = 0; i < MEM<M1>::N; i++){ //there is an implicit cast from int to ulint
    double zij, sum_j_pij = 0.;  // sum of each row of P[i][j]
    for(iint j = 0; j < MEM<M1>::M; j++){
    	sum_j_pij += exp( Prg.at(MEM<M1>::Sample[i]).at(j) ) * MEM<M1>::V[j];
    }
    for(iint j = 0; j < MEM<M1>::M; j++) {
	    if(sum_j_pij != 0.){
        zij = log( exp( Prg.at(MEM<M1>::Sample[i]).at(j) ) * MEM<M1>::V[j] / sum_j_pij );
        if ( zij != -fnl::infinity() )
		      MEM<M1>::P[i][j] = zij; 
      }
      //else { //at 1st round, imply every Prg[i][j] is zero, read mapped no where
        //MEM<M1>::P[i][j] = -fnl::infinity(); //keep it mapped now where, so effective_N is still effective, noop needed here
      //}
    }
  }
}; //tfdf of MONO_MEM<M1, M2>::estep()
//NOTE sum_j_pij may underflow, how to deal?

//template <class M1, class M2> extern void MONO_MEM<M1, M2>::mstep() //tfdc of MONO_MEM<M1, M2>::mstep()
template <class M1, class M2>
void MONO_MEM<M1, M2>::mstep(){
  //cerr<<"in mstep..."<<endl; //cin.get();
  flt sum_pij = 0.;
	for(iint j = 0; j < MEM<M1>::M; j++){
    double sum_i_pij = 0.;
    for( iint i = 0; i < MEM<M1>::N; i++){
      sum_i_pij += exp(MEM<M1>::P.at(i).at(j));
    }
    MEM<M1>::V[j] = sum_i_pij/double(MEM<M1>::N);
    sum_pij += sum_i_pij;
  }
  assert( round(sum_pij) == MEM<M1>::N );
  //cerr<<"verifying sum_pij="<<sum_pij<<"  N="<<MEM<M1>::N<<endl;
}; //tfdf of MONO_MEM<M1, M2>::mstep()

//tfdf of MONO_MEM<M1, M2>::rand_V()
template <class M1, class M2>
void MONO_MEM<M1, M2>::rand_V() {
  double sum_j = 0.;
  for(iint j=0; j < MEM<M1>::M; j++){
    MEM<M1>::V[j] = (double)rand()/(RAND_MAX);
    sum_j += MEM<M1>::V[j];
  }
  for(iint j=0; j < MEM<M1>::M; j++)
    MEM<M1>::V[j] = MEM<M1>::V[j]/sum_j;
} //tfdf of MONO_MEM<M1, M2>::rand_V()

//template <class M1, class M2> extern void MONO_MEM<M1, M2>::moment_V() //tfdc of MONO_MEM<M1,M2>::moment_V()
template <class M1, class M2>
void MONO_MEM<M1, M2>::moment_V() { //init using Moment Estimates
  double max_v = -fnl::infinity();
  iint max_j = 0; float co_max = 1;
  VF V_moments = VF(MEM<M1>::M,0);
  for( iint i = 0; i < MEM<M1>::N; i++ ){
    max_v = -fnl::infinity(); max_j = 0; co_max = 1;
    for (iint j = 0; j < MEM<M1>::M; j++ ){
      if ( max_v < Prg.at(MEM<M1>::Sample[i]).at(j) ){ max_v = Prg.at(MEM<M1>::Sample[i]).at(j); max_j = j; co_max=1; }
      if ( max_v == Prg.at(MEM<M1>::Sample[i]).at(j) ){ co_max += 1; max_j = ( double(rand())/RAND_MAX>1/co_max ) ? j : max_j ; }
    }
    assert( max_j >= 0 && max_j < MEM<M1>::M ); //find at least one max element within seq
    V_moments[ max_j ] += 1;
  }
  for(iint j = 0; j < MEM<M1>::M; j++){
    MEM<M1>::V[j] = V_moments[j]/MEM<M1>::N;
  }
}; //tfdf of MONO_MEM<M1,M2>::moment_V()

//template <class M1, class M2> extern void MONO_MEM<M1, M2>::init_Prg( string mmf ) //tfdc of MONO_MEM<M1, M2>::init_Prg(), adapted from SparseLib++ example
template <class M1, class M2> 
void MONO_MEM<M1, M2>::init_Prg( iint roffset, string mmf ){
  int ret_code;
  int N, M, nnz;
  MM_typecode matcode;
  FILE *f;
  int I, J;
  double val;
  
  if ( (f = fopen( mmf.c_str(), "r"))== NULL ) { cerr<<"Can't open mmf..."<<endl; exit(1); }  
  if ( mm_read_banner(f, &matcode) != 0 ) { cerr<<"Could not process Matrix Market banner."<<endl; exit(1); }
  if ( (ret_code = mm_read_mtx_crd_size(f, &N, &M, &nnz)) !=0 ){ cerr<<"Matrix size incoherent..."<<endl; exit(1); }
  
  assert( N == MEM<M1>::N && M == MEM<M1>::M ); // must be full rank matrix 
  for (int i=0; i<nnz; i++){
    fscanf(f, "%d %d %lg\n", &I, &J, &val);
    I--;  // adjust from 1-based to 0-based 
    J--;
    Prg[I+roffset][J] = log(val);  //log Prg
  }
  //cerr<<"final test_mem.Prg="; showmat<TR1_UUM_L>(Prg,roffset+N,M); cerr<<endl;
}; //tfdf of MONO_MEM<M1, M2>::init_Prg, adapted from SparseLib++ example
//typedef MONO_MEM<TR1_UUM_L,TR1_UUM_L> UUML_UUML_MONO;

//tfdf for MONO_MEM<M1,M2>::init_EM()
template <class M1, class M2>
void MONO_MEM<M1, M2>::init_ENANMN(){
  //showmat( Prg, MEM<M1>::N, MEM<M1>::M );
  EN = 0;   // Reads Number, should be N
  AN = 0;   // Ambiguous Reads Number, to deal with ambiguous reds
  MN = 0;   // Mapped Reads Number, mapped reads
  int unmapped = 0; int effective = 0;
  for( iint i = 0; i < MEM<M1>::N; i++ ){
    unmapped = 0; effective = 0;
    for( iint j = 0; j < MEM<M1>::M; j++ ) {
      if ( Prg.at(i).at(j) != -fnl::infinity() ) {
        effective += 1;
        if ( j == MEM<M1>::M - 1 ) unmapped += 1;
      }
    }
    assert( unmapped * effective <= 1 ); //ensure umapped reads is only mapped to unknown
    if(effective>0) EN += 1;
    if(effective>1) AN += 1;
    if(unmapped<1) MN += 1;
  }
  cerr<<"Total Reads Number = "<<MEM<M1>::N<<endl;
  cerr<<"Effective Reads Number = "<<EN<<endl;
  cerr<<"Mapped Reads Number = "<<MN<<endl;
  cerr<<"Ambiguous Reads Number = "<<AN<<endl;
  assert( EN == MEM<M1>::N );  // every read is mapped somewhere
} //tfdf for MONO_MEM<M1,M2>::init_EM()

#endif // _GEM_GEMCORE_H
