/**************************
*Kmer Counting Model Header
***************************/
#ifndef _GEM_KMERCORE_H
#define _GEM_KMERCORE_H

#include <string>
#include <vector>

using namespace std;
typedef unsigned long long kint; // for kmer key, up to 2^64 = 18446744073709551616L
typedef vector<kint> VK;


/* CONSTANTS */

const kint p4[24] = { 	      1ull, 4ull, 16ull, 64ull, 256ull, 1024ull, 4096ull, 16384ull, 65536ull, 262144ull, //10
	 					                  1048576ull, 4194304ull, 16777216ull, 67108864ull, 268435456ull, 1073741824ull,     //16
                      	      4294967296ull, 17179869184ull, 68719476736ull, 274877906944ull, 1099511627776ull,  //21
 						                  4398046511104ull, 17592186044416ull, 70368744177664ull };                          //24
const VK p4tab( p4, p4 + sizeof(p4)/sizeof(*p4) ); //support modulus based kmer counting

/* FUNCTIONS */
extern kint base2num( const char& ); //fdc of base2num( const char& )

/* TEMPLATES */
//template <class M1, class M1_IT> M1 kmer2map( string, int, bool);
template <class C1, typename flt>  // assume default return is 0
void kmer2cont( C1& v, string seq, int k=6, bool n=true ){ // kmer counter by numeric method, for k<=16
  kint f = 0; int t = 0; //f is which feature, t is total kmers
  for(kint i=-1; i<(kint)seq.size(); ){
    //cerr<<"i="<<i<<endl;
    int b = i<0 ? -1 : base2num(seq[i]);
    if(b==-1){ // either it is from start or need to restart, we have to read k consecutive
      int c=0; f=0; // counter for the number we have read
      while( c<k ){
        i++;
        if( i>=(kint)seq.size() ) break; // no more to be read
        int b = base2num(seq[i]); // try to read one more
        if(b==-1) { f=0; c=0; continue; } // again fall back to jump ahead k letters and restart
        else { c++; f=f*4+b; } 
      }
      if( i<( kint)seq.size() ) { v[f]+=1; t+=1; i++; }
    } 
    else { f=(f%p4tab.at(k-1))*4+b; v[f]+=1; t+=1; i++; } 
  }
  flt dt = (flt)t;
  if (n) {
    for( kint j=0; j<p4tab[k]; j++; ) { 
      if (v.at(f) != 0)  v[f]=v.at(f)/dt;  
    };
  }
}; //fdf of kmer2map()

#endif //_GEM_KMERCORE_H
