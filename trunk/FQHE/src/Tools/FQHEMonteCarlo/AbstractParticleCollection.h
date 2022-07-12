////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
// abstract class for collection of particles used in a Monte Carlo algorithm // 
//                                                                            //
//                     last modification : 18/02/2008                         //
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


#ifndef ABSTRACTPARTICLECOLLECTION_H
#define ABSTRACTPARTICLECOLLECTION_H

#include "config.h"

#include "Vector/RealVector.h"

class AbstractParticleCollection
{
 public:
  enum Types
    {
      Other = 0x0,
      OnSphereCollection = 0x01,
      OnDiskCollection = 0x02
    };

 protected:
  

  // index of last moved particle
  int LastMoved;

 public:
  // destructor
  virtual ~AbstractParticleCollection();
  
  // randomly moves particle number nbrParticle
  virtual void Move(int nbrParticle) = 0;

  // randomly select a particle and move it
  // return value = number of particle that was moved
  virtual int Move() = 0;

  // get number of last particle that was moved
  int GetMovedNbr() { return LastMoved; }

  // get number of particles
  virtual int GetNbrParticles() = 0;
 
  // restore last move
  virtual void RestoreMove() = 0;

  // get all particle positions
  virtual RealVector& GetPositions() = 0;

  // allow access to internal Random number generator:
  virtual double GetRandomNumber() = 0;

  // randomize particle positions
  virtual void Randomize() = 0;

  // get type of collection
  virtual int GetCollectionType() = 0;
  
};

#endif
