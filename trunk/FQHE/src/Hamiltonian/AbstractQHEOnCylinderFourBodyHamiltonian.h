////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                         to particles on a torus with                       //
//                                                                            //
//                        last modification : 27/06/2003                      //
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


#ifndef ABSTRACTQHEONCYLINDERFOURBODYHAMILTONIAN_H
#define ABSTRACTQHEONCYLINDERFOURBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"

#include <iostream>


using std::ostream;


class AbstractArchitecture;


class AbstractQHEOnCylinderFourBodyHamiltonian : public AbstractQHEHamiltonian
{

  friend class QHEParticlePrecalculationOperation;
  
 protected:
  
  // Hilbert space associated to the system
  ParticleOnSphere* Particles;

  // number of particles
  int NbrParticles;

  // global energy shift (can be used to store the energy of the Wigner crystal)
  double EnergyShift;

  // ratio between the width in the x direction and the width in the y direction
  double Ratio;
  // ratio between the width in the y direction and the width in the x direction
  double InvRatio;

  // maximum momentum value reached by a particle in the state
  int MaxMomentum;
  // number of Lz values in a state
  int NbrLzValue;

  //amplitude of the quadratic confinement potential
  double Confinement;

  //parameter for the electric field
  double ElectricField;

  //parameter for the magnetic field (needed when electric field is nonzero)
  double MagneticField;

  // array containing all interaction factors 
  Complex* InteractionFactors;
  // number of interaction factors
  int NbrInteractionFactors;
  // arrays for indices attached to each interaction factor
  unsigned* CreationIndices;
  unsigned* AnnihilationIndices;


  // array that contains all one-body interaction factors
  Complex* OneBodyInteractionFactors;
  
  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  // step between each precalculated index
  int FastMultiplicationStep;
  // number of non-null term in the hamiltonian for each state
  int* NbrInteractionPerComponent;
  // index of the state obtained for each term of the hamiltonian when applying on a given state
  int** InteractionPerComponentIndex;
  // multiplicative coefficient obtained for each term of the hamiltonian when applying on a given state and with a given destination state
  Complex** InteractionPerComponentCoefficient;

  // flag for implementation of hermitian symmetry
  bool HermitianSymmetryFlag;

 public:

  // constructor
  //
  AbstractQHEOnCylinderFourBodyHamiltonian();

  // destructor
  //
  virtual ~AbstractQHEOnCylinderFourBodyHamiltonian() = 0;

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  virtual AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  virtual int GetHilbertSpaceDimension ();
  
  // ask if Hamiltonian implements hermitian symmetry operations
  // 
  virtual bool IsHermitian();

  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  virtual void ShiftHamiltonian (double shift);

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
					  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				     int firstComponent, int nbrComponent);
 
  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  virtual List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  virtual List<Matrix*> RightInteractionOperators();  

  // save precalculations in a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be stored
  // return value = true if no error occurs
  virtual bool SavePrecalculation (char* fileName);

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

  // test the amount of memory needed for fast multiplication algorithm
  //
  // allowedMemory = amount of memory that cam be allocated for fast multiplication
  // return value = amount of memory needed
  virtual long FastMultiplicationMemory(long allowedMemory);

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm
  //
  virtual void EnableFastMultiplication();

  // enable fast multiplication algorithm (partial evaluation)
  //
  // jobIndex = index of the job that proceeds part of the fast multiplication evaluation
  // nbrJob = number of jobs that proceed the fast multiplication evaluation
  virtual void PartialEnableFastMultiplication(int jobIndex, int nbrJob);

  // load precalculations from a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be read
  // return value = true if no error occurs
  virtual bool LoadPrecalculation (char* fileName);

};

#endif
