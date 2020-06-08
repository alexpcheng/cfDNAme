/**************************
*Kmer Counting Model Source
***************************/

#ifndef _GEM_KMERCORE_C
#define _GEM_KMERCORE_C

#include "kmercore.hpp"

/* FUNCTIONS */
//fdf of base2num( const char& )
kint base2num( const char& b ){
  switch (b){
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
    default:
      return -1; 
  }
}; //fdf of base2num( const char& )

#endif //_GEM_KMERCORE_C
