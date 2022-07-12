////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//         Copyright (C) 2001-2002 Nicolas Regnault and Gunnar Moeller        //
//                                                                            //
//                                                                            //
//           class of particles with complex coordinates on a disk            //
//                                                                            //
//                        last modification : 16/07/2016                      //
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


#ifndef PARTICLEONDISKCOLLECTION_H
#define PARTICLEONDISKCOLLECTION_H

#include "config.h"
#include "MathTools/Complex.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "GeneralTools/GarbageFlag.h"
#include "Vector/RealVector.h"
#include "AbstractParticleCollectionOnDisk.h"
#include "Matrix/RealSymmetricMatrix.h"

class ParticleOnDiskCollection : public AbstractParticleCollectionOnDisk {
 protected:
  // number of particles
  int NbrParticles;
  // typical radius for moving particles
  double R0;

  // target filling factor
  double Nu;

  Complex LastZ;
  /* double LastR; */
  /* double LastPhi; */
  Complex *CoordinatesZ;
  AbstractRandomNumberGenerator *Generator;
  bool ExternalGenerator;
  GarbageFlag Flag;
  RealVector Positions;
  
 public:

  // constructors:
  ParticleOnDiskCollection();
  ParticleOnDiskCollection(int N, double nu, long seed =-1);
  ParticleOnDiskCollection(int N, double nu, AbstractRandomNumberGenerator *generator);
  ParticleOnDiskCollection(const ParticleOnDiskCollection &tocopy);

  //destructor
  virtual ~ParticleOnDiskCollection();
  
  // randomly moves particle number nbrParticle
  virtual void Move(int nbrParticle);

  // randomly select a particle and move it
  // return value = number of particle that was moved
  virtual int Move();

  // get previous coordinates of last particle that was moved
  virtual void GetPreviousPos(Complex &lastZ)
  { lastZ=this->LastZ; }

  // get pointers to spinor coordinates
  virtual void GetCoordinates(Complex* &Z)
  { Z=this->CoordinatesZ; }

  // get number of particles
  virtual int GetNbrParticles(){ return NbrParticles; }

  // restore last move
  virtual void RestoreMove();

  // stretch the default steplength by a factor
  virtual void MultiplyStepLength(double multiplier);

  // access particle positions of a single particle 
  virtual double GetR(int nbrParticle);
  virtual double Phi(int nbrParticle);

  // set new particle position
  virtual void SetPosition(int nbrParticle, double x, double y);

  // get all particle positions
  virtual RealVector& GetPositions();

  // allow access to internal Random number generator:
  virtual double GetRandomNumber();

  // randomize particle positions
  virtual void Randomize();

  // get type of collection
  virtual int GetCollectionType();

  // get absolute values of all relative distances
  // distances = matrix in which to return the distances
  virtual void GetDistances(RealSymmetricMatrix &distances);

  // toggle positions of first N/2 particles with the remaining N/2 positions
  //
  virtual void ToggleHalfHalf();

  // print positions (mainly for testing)
  ostream &PrintPositions(ostream &Str);
};


inline RealVector& ParticleOnDiskCollection::GetPositions()
{
  return this->Positions;
}

// get type of collection
inline int ParticleOnDiskCollection::GetCollectionType()
{
  return AbstractParticleCollection::OnDiskCollection;
}
  


#endif
