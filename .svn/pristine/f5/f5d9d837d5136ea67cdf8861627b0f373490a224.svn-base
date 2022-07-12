////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
// class for the eigenfunctions of the triangular well confinement potential  //
//                                                                            //
//                        last modification : 04/10/2012                      //
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

#ifndef TRIANGULARWELLEIGENFUNCTION_H
#define TRIANGULARWELLEIGENFUNCTION_H

#include "config.h"

//#include <cstring>
//#include <stdio.h>
//#include <cstdlib>
#include <iostream>

#include <cmath>
#ifdef __GSL__
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#endif

class TriangularWellEigenfunction
{
 private:
  // width of the well
  double w;
  
  // bias of the well, i.e. height of the 'triangularity'
  double V;
  
  // The two firsts energies
  double E[2];
  
  // Norms of the two firsts eigenfunctions
  double norm[2];
  
  // Number of steps for the determination of eigenrenergies via root finding of the eigenfunctions
  int NewtonMaxStepsNbr;
  double NewtonEpsRelative;
  
  // absolute and relative errors allowed for the computation of the norm
  double integration_absolute_error;
  double integration_relative_error;

  //constants used by methods
  double c1,c2,c3[2];
 
 public:
  // Constructor
  //
  // width = thickness of the quantum well
  // bias = bias of the triangular quantum well in energy per length
  //
  TriangularWellEigenfunction(double width, double bias);

  // Destructor
  //
  ~TriangularWellEigenfunction();
  
  // Returns the energy of the nth eigenstate of the triangular quantum well
  //
  // n = number of the eigenstate, n=0 (ground state) or n=1 (first excited state)
  //
  // return value = nth eigenenergy 
  //
  double GetEnergy(int n);
  
  // General solution of the triangular-well potential Schroedinger equation with energy E
  //
  // x = position of the confinement potential where the function is evaluated
  // E = energy 
  //
  // return value = value of the function  at x
  //
  double Psi(double x,double E);

  // Derivative with respect to E of Psi (defined above) at x=0
  //
  // x = position of the confinement potential where the derivative is evaluated
  // E = energy 
  //
  // return value = value of the derivative at E
  //
  double PsiPrime(double E);

  // nth eigenfunction of the triangular well confinement potential
  //
  // x = position of the confinement potential where the function is evaluated
  // n = number of the eigenstate (can be either 0 for the ground state or 1)
  //
  // return value = value of the eigenfunction  at x
  //
  double Eigenfunction(double x, int n);

  // Define a static GSL_function to integrate the square of the eigenfunction Psi, and then get its norm
  //
  // x = position
  // p = pointer, which adress points toward a structure
  //
  // return value = Psi^2 in a GSL_function form
  //
  static double GSL_f(double f, void *params);

  // Norm of the nth eigenfunction, norm=1/sqrt(\int_0^w Psi^2(x,E) dx)
  //
  // n = number of the eigenstate (can be either 0 for the ground state or 1)
  //
  // return value = norm of the nth eigenstate
  //
  double Norm(int n);
  
  // Energy of the nth level, given by a WKB-type approximation, used for root finding of psi(x) as the initial value, which gives accurates energies
  //
  // n = number of the eigenstate (can be either 0 for the ground state or 1)
  //
  // return value = wkb approximation of the nth eigenenergy
  //
  double wkb(int n);

  // Energy of the nth level
  //
  // n = number of the eigenstate (can be either 0 for the ground state or 1)
  //
  // return value = nth eigenenergy
  //
  double Eigenenergy(int n);

  // GSL functions used for root-finding
  static double GSL_Newton_f(double E, void *params);
  static double GSL_Newton_df(double E, void *params);
  static void GSL_Newton_fdf(double E, void *params, double *y, double *dy);
};

#endif
