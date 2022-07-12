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


#ifndef TABULATEDDENSITYPROFILE_H
#define TABULATEDDENSITYPROFILE_H


#include "config.h"

#include "AbstractZDensityProfile.h"

#ifdef HAVE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#endif


class TabulatedDensityProfile : public AbstractZDensityProfile
{
 private:

  // minimum z where function is defined
  double Zmin;

  // maximum z where function is defined
  double Zmax;

  // the norm of the (square) wavefunction
  double SqNorm;

  // number of data-points
  int NbrDataPoints;

  // abscissa of the data-points
  double *ZValues;

  // Density according to abscissa
  double *Density;

  // flag indicating whether interpolation is present
  bool SplineInterpolate;

#ifdef HAVE_GSL
  // gsl spline object
  gsl_spline *Spline;
  // gsl spline accelerator object
  gsl_interp_accel *SplineAccelerator;
#endif

 public:
  
  // the standard constructor 
  TabulatedDensityProfile();

  // constructor
  // dataFile = file name of text-file to be opened
  // scale = scale-factor in the z-direction
  // squareValues = option to square the values of individual points (wavefunction provided)
  // splineInterpolate = flag indicating whether a spline interpolation should be used  
  // checkNorm = flag indicating whether the norm should be recalculated
  TabulatedDensityProfile(char *dataFile, double scale, bool squareValues=false, bool checkNorm=false);
  
  // virtual destructor
  virtual ~TabulatedDensityProfile();

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
inline int TabulatedDensityProfile::GetType()
{
  return AbstractZDensityProfile::TabulatedProfile;
}


#endif
