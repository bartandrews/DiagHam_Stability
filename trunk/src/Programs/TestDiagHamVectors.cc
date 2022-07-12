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

#include "config.h"

#include <iostream>
using std::cout;
using std::endl;

#include "Vector/RealVector.h"
#include "Vector/RealPtrVector.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"

#include <sys/time.h>



int main()
{

 // test some basic vector routines

  cout << "Testing vector routines"<<endl;


#ifdef __INTEL_COMPILER

  cout << "Using Intel compiler"<<endl;

#endif

 NumRecRandomGenerator MyGenerator;

  int Length=1000000;
  RealVector V1(Length);
  RealVector V2(Length);

  for (int i=0; i<Length; ++i)
    {
      V1[i] = MyGenerator.GetRealRandomNumber();
      V2[i] = MyGenerator.GetRealRandomNumber();
    }

  RealPtrVector VPt2(Length);
  for (int i=0; i<Length; ++i)
    VPt2[i] = &V2[i];

  
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  int Attempts=10000;
  double prod, Dt;

  V1.HaveBlas(true);
  
  gettimeofday (&(TotalStartingTime), 0);
  for (int i=0; i<Attempts; ++i)
    prod = V1*V2;
  gettimeofday (&(TotalEndingTime), 0);
  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1.0e6);		      
  cout << "Took " << Dt << " s for regular scalar product : "<<prod <<endl;


  gettimeofday (&(TotalStartingTime), 0);
  for (int i=0; i<Attempts; ++i)
    prod = V1*VPt2;
  gettimeofday (&(TotalEndingTime), 0);
  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1.0e6);		      
  cout << "Took " << Dt << " s for pointered scalar product with pointers in contiguos order: "<<prod <<endl;  
  
  for (int i=0; i<Length; ++i)
    VPt2[i] = &V2[(int)(MyGenerator.GetRealRandomNumber()*Length)];

  gettimeofday (&(TotalStartingTime), 0);
  for (int i=0; i<Attempts; ++i)
    prod = V1*VPt2;
  gettimeofday (&(TotalEndingTime), 0);
  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1.0e6);		      
  cout << "Took " << Dt << " s for pointered scalar product with pointers in random order: "<<prod <<endl;  
  
  exit(1);


}
