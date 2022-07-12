////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//         Copyright (C) 2001-2002 Nicolas Regnault and Gunnar Moeller        //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 16/09/2004                      //
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


#ifndef ABSTRACTPARTICLECOLLECTIONONSPHERE_H
#define ABSTRACTPARTICLECOLLECTIONONSPHERE_H
#include "config.h"
#include "AbstractParticleCollection.h"
#include "MathTools/Complex.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "GeneralTools/GarbageFlag.h"
#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"

class AbstractParticleCollectionOnSphere : public AbstractParticleCollection {
 protected:
  
 public:

  //destructor
  virtual ~AbstractParticleCollectionOnSphere();
  
  // rotate the coordinates of all particles with a common random vector (cannot be reverted for now)
  virtual void RotateAll() = 0;

  // rotate the coordinates of all particles by the rotation that brings the north pole to angles theta and phi
  virtual void RotateAll(double theta, double phi) = 0;

  // get previous coordinates of last particle that was moved
  virtual void GetPreviousPos(Complex &lastU, Complex &lastV) = 0;

  // get pointers to spinor coordinates
  virtual void GetSpinorCoordinates(Complex* &U, Complex* &V) = 0;

  // stretch the default steplength by a factor
  virtual void MultiplyStepLength(double multiplier) = 0;

  // access particle positions of a single particle 
  virtual double Theta(int nbrParticle) = 0;
  virtual double Phi(int nbrParticle) = 0;

  // set new particle position
  virtual void SetPosition(int nbrParticle, double theta, double phi) = 0;

  // get absolute values of all relative distances
  // distances = matrix in which to return the distances
  virtual void GetDistances(RealSymmetricMatrix &distances) = 0;

  // get reference to internal matrix with all relative distances
  // return = matrix with distances
  virtual const RealSymmetricMatrix& GetDistances() = 0;


  // toggle positions of first N/2 particles with the remaining N/2 positions
  //
  virtual void ToggleHalfHalf() = 0;
  
};


#endif
