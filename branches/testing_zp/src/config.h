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


#ifndef CONFIG_H
#define CONFIG_H


#include "config_ac.h"

// all options

// byte ordering
#ifdef HAVE_LITTLEENDIAN 
#define __LITTLEENDIAN__
#else
#define __BIGENDIAN__
#endif

// use sstream instead of strstream
#define  __SSTREAM_STYLE__

// machine precision
#define MACHINE_PRECISION 1e-14

// SMP flag 
#define __SMP__

// debug flag
#define __DEBUG__

// 64 bits architecture
//#define __64_BITS__

// MPI flag
#ifdef HAVE_MPI
#define __MPI__
#endif


// LAPACK flag
#ifdef HAVE_LAPACK
#define __LAPACK__
// some LAPACK type definitions
typedef double doublereal;
typedef struct { doublereal r, i; } doublecomplex;
#ifdef __STDC__
#define FORTRAN_NAME(n_)	n_ ## _
#else
#define FORTRAN_NAME(n_)	n_/**/_
#endif
#endif

// use LAPACK routines instead of DiagHam routines
#ifdef HAVE_LAPACK_ONLY
#define __LAPACKONLY__
#endif


// bz2 flag
#ifdef HAVE_BZ2
#define __BZ2LIB__
#endif

// architecture dependant options
//

// DEC CXX specific options

# if defined __DECC || defined __DECCXX

// enable cxx options
#define __USE_STD_IOSTREAM

// 64 bits architecture
#ifndef __64_BITS__
#define __64_BITS__
#endif

#endif



// xlC and AIX specific options (assume 64bits compilation)

# if defined __TOS_AIX__ && __xlC__

// 64 bits architecture
#ifndef __64_BITS__
#define __64_BITS__
#endif

#endif



// gcc and x86_64 specific options (assume 64bits compilation)

#ifdef __x86_64__

// 64 bits architecture
#ifndef __64_BITS__
#define __64_BITS__
#endif

#endif

// pathscale 64 bit options - test for 64 bit compilation with pathscale suite

#ifdef __LP64__

// 64 bits architecture
#ifndef __64_BITS__
#define __64_BITS__
#endif

#endif

// define long long types (128 bits on x86_64 architecture, 64 bits elsewhere)

#define LONGLONG long long
#define ULONGLONG unsigned long long

// package option
//

// provide output package
#define USE_OUTPUT

// provide polynomial package
//#define USE_POLYNOMIAL

// provide use of the cluster architecture package
#define USE_CLUSTER_ARCHITECTURE

// provide use of the generic Hilbert space package
#define USE_HILBERT_SPACE

#endif
