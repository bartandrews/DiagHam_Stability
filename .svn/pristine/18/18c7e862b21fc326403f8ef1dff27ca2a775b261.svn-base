////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          global configuration file                         //
//                                                                            //
//                        last modification : 18/01/2001                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
using std::cout;
using std::endl;

#include "config.h"


int main()
{

// all options

// byte ordering
#ifdef HAVE_LITTLEENDIAN 
  cout <<"__LITTLEENDIAN__ defined"<<endl;
#else
  cout <<"__BIGENDIAN__ defined"<<endl;
#endif

// use sstream instead of strstream
#ifdef __SSTREAM_STYLE__ 
  cout <<" __SSTREAM_STYLE__ defined"<<endl;
#endif

// machine precision
  cout <<"MACHINE_PRECISION 1e-14 defined"<<endl;

// SMP flag
#ifdef __SMP__
  cout <<"__SMP__ defined"<<endl;
#endif
  
// debug flag
#ifdef __DEBUG__
  cout <<"__DEBUG__ defined"<<endl;
#endif

// 64 bits architecture
#ifdef __64_BITS__
  cout <<"__64_BITS__ defined"<<endl;
#endif

// MPI flag
#ifdef HAVE_MPI
  cout <<"__MPI__ defined"<<endl;
#endif


// LAPACK flag
#ifdef HAVE_LAPACK
  cout <<"__LAPACK__ defined"<<endl;
#endif

// use LAPACK routines instead of DiagHam routines
#ifdef HAVE_LAPACK_ONLY
  cout <<"__LAPACKONLY__ defined"<<endl;
#endif


// architecture dependant options
//

// DEC CXX specific options

# if defined __DECC || defined __DECCXX

  cout << "DEC system recognized"<<endl;
// enable cxx options
  cout <<"__USE_STD_IOSTREAM defined"<<endl;

// 64 bits architecture
  cout <<"__64_BITS__ defined"<<endl;

#endif



// xlC and AIX specific options (assume 64bits compilation)

# if defined __TOS_AIX__ && __xlC__
  cout <<"AIX system recognized"<<endl;
// 64 bits architecture
  cout <<"__64_BITS__ defined"<<endl;

#endif



// gcc and x86_64 specific options (assume 64bits compilation)

#ifdef __x86_64__

  cout <<"x86 Architecture recognized"<<endl;
// 64 bits architecture

  cout <<"__64_BITS__ defined"<<endl;

#endif

// pathscale 64 bit options - test for 64 bit compilation with pathscale suite

#ifdef __LP64__

  cout << "Pathscale compiler found"<<endl,

// 64 bits architecture
  cout <<"__64_BITS__ defined"<<endl;

#endif

#ifdef __INTEL_COMPILER

  cout << "Intel compiler found"<<endl;

#endif

// define long long types (128 bits on x86_64 architecture, 64 bits elsewhere)

#ifdef __x86_64__


#ifdef __INTEL_COMPILER

  cout <<"__128_BIT_LONGLONG__ defined"<<endl;
  cout <<"LONGLONG __m128i defined"<<endl;
  cout <<"ULONGLONG __m128i defined"<<endl;

#else

 cout <<"__128_BIT_LONGLONG__ defined"<<endl;
 cout <<"LONGLONG int128_t defined"<<endl;
 cout <<"ULONGLONG uint128_t defined"<<endl;

#endif

#else

 cout <<"LONGLONG long long defined"<<endl;
 cout <<"ULONGLONG unsigned long long defined"<<endl;

#endif

// package option
//

// provide output package
#ifdef USE_OUTPUT
 cout <<"USE_OUTPUT defined"<<endl;
#endif



// provide polynomial package
#ifdef USE_POLYNOMIAL
 cout <<"USE_POLYNOMIAL defined"<<endl;
#endif

// provide use of the cluster architecture package
#ifdef USE_CLUSTER_ARCHITECTURE
 cout <<"USE_CLUSTER_ARCHITECTURE defined"<<endl;
#endif

// provide use of the generic Hilbert space package
#ifdef USE_HILBERT_SPACE
 cout <<"USE_HILBERT_SPACE defined"<<endl;
#endif

}
