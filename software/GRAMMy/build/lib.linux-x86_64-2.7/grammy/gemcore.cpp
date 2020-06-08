/***********************
*Generic EM Model Source
************************/
#ifndef _GEM_GEMCORE_C
#define _GEM_GEMCORE_C

#include "gemcore.hpp"

/* CONSTANTS */
// constants defined here and declared extern in .hpp to avoid redefiniton errors from swig
const flt dft_rt_pflt = 0;
const flt dft_rt_lflt = -fnl::infinity();
//flt dft_rt_pflt = 0;
//flt dft_rt_lflt = -fnl::infinity();
const TR1_UM_P dft_rt_TR1_UM_P = TR1_UM_P();
const TR1_UM_L dft_rt_TR1_UM_L = TR1_UM_L();

/* FUNCTIONS */

//fdf of vdif( VF&, VF& )
flt vdif( VF& a, VF& b){
  flt dif = 0.;
  for ( unsigned i=0; i<a.size(); i++ )
    dif+=pow(a.at(i)-b.at(i), (flt)2.);
  return sqrt(dif); // differ as euclidean distance
}; //fdf of vdif( VF&, VF& )

//fdf of dmax( VF&, VF& )
flt dmax( VF& a, VF& b){
  flt dmax = 0;
  for ( unsigned i=0; i<a.size(); i++ )
    //dmax = ( dmax <= (abs(a.at(i)-b.at(i))/a.at(i)) ) ? ( abs(a.at(i)-b.at(i))/a.at(i) ) : dmax;
    dmax = ( dmax <= abs(a.at(i)-b.at(i)) ) ?  abs(a.at(i)-b.at(i)) : dmax; 
  return dmax; // dmax as maximum of pairwise differ
}; //fdf of dmax( VF&, VF& )

/* TEMPLATES */
//tcdf tr1_unordered_map, extending tr1::unordered_map

//tcdf of MEM<M1>, base class for Mixture Expectatin Maximization

//tfdf of MEM<M1>::solve( flt, int, char ), the solver for MEM

//tcdf of inherited MONO_MEM<M1,M2>:MEM<M1>

//tfdf of MONO_MEM<M1,M2>::loglikelihood()

//tfdf of MONO_MEM<M1,M2>::init_V()

//tfdf of MONO_MEM<M1, M2>::mstep()

//tfdf of MONO_MEM<M1, M2>::estep()

//tfdf of MONO_MEM<M1, M2>::init_Prg, adapted from SparseLib++ example

#endif // _GEM_GEMCORE_C
