////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2006 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                               delta interaction                            //
//                                                                            //
//                        last modification : 16/09/2002                      //
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


#ifndef PARTICLEONSPHEREDELTAHAMILTONIAN_H
#define PARTICLEONSPHEREDELTAHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian : public AbstractQHEHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:

  // Hilbert space associated to the system -> redefined with respect to inherited field! Check if properly used!
  ParticleOnSphereWithSpin* Particles;

    // number of particles
  int NbrParticles;

  // maximum Lz value reached by a particle in the state
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;

  // Interaction coefficients in Delta and LaplacianDelta channels:
  double V0;
  double V1;

  // array containing all interaction factors linking the same spin species:
  double* UUInteractionFactors;
  // number of interaction factors
  int UUNbrInteractionFactors;
  // arrays for indices attached to each interaction factor
  int* UUM1Value;
  int* UUM2Value;
  int* UUM3Value;
  

  // array containing all interaction factors linking different spin species:
  double* UDInteractionFactors;
  // number of interaction factors
  int UDNbrInteractionFactors;
  // arrays for indices attached to each interaction factor
  int* UDM1Value;
  int* UDM2Value;
  int* UDM3Value;


    // shift to apply to go from precalculation index to the corresponding index in the HilbertSpace
  int PrecalculationShift;

  // amount of memory (in bytes) that can be used to store precalculated matrix elements
  long Memory; 
  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  // step between each precalculated index (main part: start at 0, FastMultiplicationStep, 2*FastMultiplicationStep, ...)
  int FastMultiplicationStep;
  // step between each precalculated index (optional part: start at 1, FastMultiplicationSubStep, 2 * FastMultiplicationSubStep, ...)
  int FastMultiplicationSubStep;
  // indicate the poistion of the data relative to the sub step precalculations in NbrInteractionPerComponent, InteractionPerComponentIndex, and InteractionPerComponentCoefficient
  int FastMultiplicationSubStepPosition;
  // number of non-null term in the hamiltonian for each state
  int* NbrInteractionPerComponent;
  // index of the state obtained for each term of the hamiltonian when applying on a given state
  int** InteractionPerComponentIndex;
  // multiplicative coefficient obtained for each term of the hamiltonian when applying on a given state and with a given destination state
  double** InteractionPerComponentCoefficient;

  // flag to indicate if a hamiltonian is temporary stored on disk
  bool DiskStorageFlag;
  // name of the file that contains hamiltonian matrix elements
  char* DiskStorageFileName;
  // index of the first row that appears in the on-disk hamiltonian
  int DiskStorageStart;
  // maximum number of non-null terms in the hamiltonian for each state
  int MaxNbrInteractionPerComponent;
  // size of the in-memory temporary buffer
  long BufferSize;

  // shift to apply to the Hamiltonian diagonal elements
  double HamiltonianShift;

  


 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax,
						   double v0, double v1, AbstractArchitecture* architecture, 
				   long memory = -1, bool onDiskCacheFlag = false, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // reset interaction coefficients
  // v0 interaction in s-wave Delta-Interaction channel
  // v1 interaction in p-wave Laplacian-Delta-Interaction channel
  virtual void SetInteraction(double v0, double v1);

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

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();

  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);


  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  List<Matrix*> RightInteractionOperators();


  /* Following: Functions that require new definition due to the spin of particles */
  /* List taken from ParticleOnTorusCoulombWithSpin.cc     */

    // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of idinces 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
			       int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
				  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				     int firstComponent, int nbrComponent);

  

  

  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian& H);

  // save precalculations in a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be stored
  // return value = true if no error occurs
  virtual bool SavePrecalculation (char* fileName);


 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  // test the amount of memory needed for fast multiplication algorithm
  //
  // return value = amount of memory needed
  long FastMultiplicationMemory(long allowedMemory);
  
  // enable fast multiplication algorithm
  //
  void EnableFastMultiplication();

  // enable fast multiplication algorithm using on disk cache 
  //
  // fileName = prefix of the name of the file where temporary matrix elements will be stored
  
  void EnableFastMultiplicationWithDiskStorage(char* fileName);

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element

  long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored

  RealVector& LowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
					     int firstComponent, int nbrComponent);


  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored

  virtual RealVector& LowLevelAddMultiplyPartialFastMultiply(RealVector& vSources, RealVector& vDestinations, 
							     int firstComponent, int nbrComponent);


  // load precalculations from a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be read
  // return value = true if no error occurs
  virtual bool LoadPrecalculation (char* fileName);
  
};

#endif
