%module gemcore
%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"
%include "typemaps.i"
%include "std_pair.i"
%include "mmio.h"
%{ 
  //this part literally included in gemcore_wrap.cxx
  using namespace std;
  #include "gemcore.hpp"
  #include "mydebug.hpp"
%}
%include "gemcore.hpp"
%include "mydebug.hpp"

namespace std{
  %template(VectorFloat) vector<float>;
  %template(VectorInt) vector<int>;
  %template(VectorDouble) vector<double>;
  %template(MatrixFloat) vector< vector<float> >;
};
%template(cMEM) MEM<TR1_UUM_L>;
%template(cMONO) MONO_MEM<TR1_UUM_L, TR1_UUM_L>;
%template(showuum) showmat<TR1_UUM_L>;
