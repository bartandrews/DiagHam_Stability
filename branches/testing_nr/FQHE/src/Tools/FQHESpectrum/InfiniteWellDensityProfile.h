////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class for calculation of Pseudopotential Coefficients  //
//                                                                            //
//                        last modification : 19/11/2007                      //
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


#ifndef INFINITEWELLDENSITYPROFILE_H
#define INFINITEWELLDENSITYPROFILE_H


#include "config.h"

#include "AbstractZDensityProfile.h"

class InfiniteWellDensityProfile : public AbstractZDensityProfile
{
 private:

  // width of the well
  double Width;

  // the norm of the (square) wavefunction
  double SqNorm;

  // the band index (starting to number at zero)
  int Band;

 public:

  // the standard constructor -> well with 1 magnetic length width
  InfiniteWellDensityProfile();

  // constructor
  // width = width of the well
  // band = band index (allows for excited states to be obtained)
  // 
  InfiniteWellDensityProfile(double width, int band = 0);
  
  // virtual destructor
  virtual ~InfiniteWellDensityProfile();

  // get minimum and maximum value of density profile where the probability density is larger than precision
  // min = minimum value of z offset
  // max = maximum value of z offset
  // precision = requested precision
  virtual void GetSupport(double &min, double &max, double precision=1e-10);

  // evaluate the density for a given offset
  // z = offset of distribution
  virtual double GetValue(double z);

  

  // get type of the density profile
  virtual int GetType();

};

// get type of the density profile
inline int InfiniteWellDensityProfile::GetType()
{
  return AbstractZDensityProfile::InfiniteWell;
}


#endif
