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


#ifndef PARTICLEONSPHERECOLLECTION_H
#define PARTICLEONSPHERECOLLECTION_H
#include "config.h"
#include "MathTools/Complex.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "GeneralTools/GarbageFlag.h"
#include "Vector/RealVector.h"
#include "AbstractParticleCollection.h"
#include "Matrix/RealSymmetricMatrix.h"

class ParticleOnSphereCollection : public AbstractParticleCollection {
 private:
  int NbrParticles;
  double Theta0;
  Complex LastU;
  Complex LastV;
  double LastTheta;
  double LastPhi;
  Complex *SpinorUCoordinates;
  Complex *SpinorVCoordinates;      
  AbstractRandomNumberGenerator *Generator;
  bool ExternalGenerator;
  GarbageFlag Flag;
  RealVector ThetaPhi;
  
 public:

  // constructors:
  ParticleOnSphereCollection();
  ParticleOnSphereCollection(int N, long seed =-1);
  ParticleOnSphereCollection(int N, AbstractRandomNumberGenerator *generator);
  ParticleOnSphereCollection(const ParticleOnSphereCollection &tocopy);

  //destructor
  ~ParticleOnSphereCollection();
  
  // randomly moves particle number nbrParticle
  void Move(int nbrParticle);

  // randomly select a particle and move it
  // return value = number of particle that was moved
  int Move();

  // rotate the coordinates of all particles with a common random vector (cannot be reverted for now)
  void RotateAll();

  // rotate the coordinates of all particles by the rotation that brings the north pole to angles theta and phi
  void RotateAll(double theta, double phi);

  // get previous coordinates of last particle that was moved
  void GetPreviousPos(Complex &lastU, Complex &lastV)
  { lastU=this->LastU; lastV=this->LastV; }

  // get pointers to spinor coordinates
  void GetSpinorCoordinates(Complex* &U, Complex* &V)
  { U=this->SpinorUCoordinates; V=this->SpinorVCoordinates; }

  // get number of particles
  int GetNbrParticles(){ return NbrParticles; }

  // restore last move
  void RestoreMove();

  // stretch the default steplength by a factor
  void MultiplyStepLength(double multiplier);

  // access particle positions of a single particle 
  double Theta(int nbrParticle);
  double Phi(int nbrParticle);

  // set new particle position
  void SetPosition(int nbrParticle, double theta, double phi);

  // get all particle positions
  RealVector& GetPositions();

  // allow access to internal Random number generator:
  double GetRandomNumber();

  // randomize particle positions
  void Randomize();

  // get type of collection
  virtual int GetCollectionType();

  // get absolute values of all relative distances
  // distances = matrix in which to return the distances
  void GetDistances(RealSymmetricMatrix &distances);

  // toggle positions of first N/2 particles with the remaining N/2 positions
  //
  void ToggleHalfHalf();
  
};


inline RealVector& ParticleOnSphereCollection::GetPositions()
{
  return this->ThetaPhi;
}

// get type of collection
inline int ParticleOnSphereCollection::GetCollectionType()
{
  return AbstractParticleCollection::OnSphereCollection;
}
  


#endif
