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


#ifndef INFINITEWELLBILAYERDENSITYPROFILE_H
#define INFINITEWELLBILAYERDENSITYPROFILE_H


#include "config.h"

#include "AbstractZDensityProfile.h"

class InfiniteWellBilayerDensityProfile : public AbstractZDensityProfile
{
 private:

  // width of the well
  double Width;

  // the norm of the (square) wavefunction
  double SqNorm;

  // variable indicating whether we are in the left (-1) or right (+1) layer
  double Symmetry;

 public:

  // the standard constructor -> well with 1 magnetic length width
  InfiniteWellBilayerDensityProfile();

  // constructor
  // width = width of the well
  // band = band index (allows for excited states to be obtained)
  // 
  InfiniteWellBilayerDensityProfile(double width, int band = 0);
  
  // virtual destructor
  virtual ~InfiniteWellBilayerDensityProfile();

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
inline int InfiniteWellBilayerDensityProfile::GetType()
{
  if (this->Symmetry<0.0)
    return AbstractZDensityProfile::InfiniteWellBilayerLeft;
  else
    return AbstractZDensityProfile::InfiniteWellBilayerRight;
}


#endif
