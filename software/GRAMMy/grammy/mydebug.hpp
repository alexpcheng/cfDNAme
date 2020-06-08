/**************************
*Debug Helper Model Header
***************************/
#ifndef _GEM_MYDEBUG_H
#define _GEM_MYDEBUG_H

#include <vector>
#include <iostream>

using namespace std;
typedef unsigned long int iint;
//typedef vector<float> VF;
//typedef VF::iterator VF_it;

//tfdf of showvec<T,Tit>()
template <class V> 
void showvec( V& v, iint m ){
  cerr<<"[";
  for( iint i=0; i<m; i++ ) 
    cerr<<v.at(i)<<",";
  cerr<<" ]";
}; //tfdf of showvec<T,Tit>()
//template<> void showvec<VF, VF_it>( VF& );

template <class M>
void showmat( M& mat, iint n, iint m ){
  for (iint i=0; i<n; i++){
    if (i==0) cerr<<'['; 
    cerr<<'[';
    for (iint j=0; j<m; j++)
      cerr<<mat.at(i).at(j)<<',';
    cerr<<"],";
    if (i==n-1) cerr<<']';
    cerr<<endl;
  }
}

#endif //_GEM_MYDEBUG_H
