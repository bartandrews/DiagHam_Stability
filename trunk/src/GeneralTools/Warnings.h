////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2022 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         some macros to generate warnings and compile-time switches         //
//                                                                            //
//                        last modification : 30/04/2022                      //
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

#ifndef WARNINGS_H
#define WARNINGS_H
#include <stdexcept>
#include "config.h"



/// Custom overload macro.
#undef NoOverload
#undef __NoOverload
#ifdef NDEBUG
  #define	NoOverload()	((void)0)
#else
  #ifndef NOFUNC
    #define NoOverload()  \
        __NoOverload (__func__, __FILE__, __LINE__)
    #define __NoOverload(function, file, line) \
        ((void) fprintf(stderr, "Function not overloaded: function %s, file %s, line %u.\n", function, file, line), throw std::runtime_error("Function not overloaded"))
  #else
    // If compiler does not support __func__
    #define NoOverload()  \
        __NoOverload (__FILE__, __LINE__)
    #define __NoOverload(file, line) \
        ((void) fprintf(stderr, "Function not overloaded: file %s, line %u.\n", file, line), throw std::runtime_error("Function not overloaded"))
  #endif
#endif


#endif // end header
