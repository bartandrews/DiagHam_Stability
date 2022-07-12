////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//                 class of quatum Hall hamiltonian associated                //
//                 to particles on a lattice in magnetic field                //
//                                                                            //
//                        last modification : 12/02/2008                      //
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


#ifndef ABSTRACTQHEONLATTICEHAMILTONIAN_H
#define ABSTRACTQHEONLATTICEHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"
#include "GeneralTools/RealUniqueArray.h"
#include "GeneralTools/ComplexUniqueArray.h"

#include <iostream>


using std::ostream;


class AbstractArchitecture;


class AbstractQHEOnLatticeHamiltonian : public AbstractQHEHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:
  
  // Hilbert space associated to the system
  ParticleOnLattice* Particles;

  // number of particles
  int NbrParticles;

  // lattice dimension in the x-direction
  int Lx;
  // lattice dimension in the y-direction
  int Ly;
  // number of sublattice sites per unit cell
  int SubLattices;
  // ky-symmetry-flag
  bool HaveKySymmetry;
  // maximum ky-momentum
  int KyMax;  
  // number of sites
  int NbrSites;  
  // number of unit cells
  int NbrCells;
  // number of flux through lattice
  int NbrFluxQuanta;
  // flux through a unit cell
  double FluxDensity;


  

  // array containing kinetic energy (hopping) terms
  Complex* HoppingTerms;
  // number of hopping terms
  int NbrHoppingTerms;
  // arrays for indices attached to hopping terms (i nitial and f inal)
  int* KineticQi;
  int* KineticQf;

  // fields for general four-fermion interactions:
  // array containing interaction factors
  Complex* InteractionFactors;
  // number of interaction factors
  int NbrInteractionFactors;
  // arrays for quantum numbers attached to each interaction factor
  int* Q1Value;
  int* Q2Value;
  int* Q3Value;
  int* Q4Value;
  
  // alternative method used to store indices attached to each interaction factor
  // Q1Value and Q2Value still store the NbrQ12Indices different combinations of Q1/Q2
  // InteractionFactors still hold the matrix elements in a completely linearized array
  int NbrQ12Indices;
  // number of different combinations of Q3/Q4 for each pair of Q1/Q2
  int* NbrQ34Values;
  // and the corresponding values of Q3/Q4 themselves
  int** Q3PerQ12;
  int** Q4PerQ12;

  // separate out completely diagonal contact-interactions:
  double* DiagonalInteractionFactors;
  // number of interaction factors (usually = number of sites)
  int NbrDiagonalInteractionFactors;
  // quantum number attached to each "self-interaction" factor
  int* DiagonalQValues;
  

  // shift to apply to go from precalculation index to the corresponding index in the HilbertSpace
  int PrecalculationShift;

  
  // amount of memory (in bytes) that can be used to store precalculated matrix elements
  long AllowedMemory;
  // flag indicating whether fast multiplication data was loaded from a file
  bool LoadedPrecalculation;
  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  // step between each precalculated index (main part: start at 0, FastMultiplicationStep, 2*FastMultiplicationStep, ...)
  int FastMultiplicationStep;
  // number of non-null real terms in the hamiltonian for each state (typically a small number)
  unsigned short* NbrRealInteractionPerComponent;
  // number of non-null complex terms in the hamiltonian for each state
  unsigned short* NbrComplexInteractionPerComponent;
  // index of the state obtained for each term of the hamiltonian when applying on a given state
  // holding indices for both real (1st) and complex matrix elements
  int** InteractionPerComponentIndex;
  // index of real (first in enumeration) and complex (following) multiplicative coefficients obtained for each term of the hamiltonian when applying on a given state and with a given destination state
  unsigned short** InteractionPerComponentCoefficientIndex;
  

  // array storing a single copy of each real matrix element value
  RealUniqueArray RealInteractionCoefficients;

  // array storing a single copy of each real matrix element value
  ComplexUniqueArray ComplexInteractionCoefficients;
  

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

  // default constructor
  //
  AbstractQHEOnLatticeHamiltonian();

  // destructor
  //
  virtual ~AbstractQHEOnLatticeHamiltonian() = 0;

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // set flux density in units of flux quanta through the lattice
  //
  // nbrFluxQuanta = flux quantua piercing the lattice
  virtual void SetNbrFluxQuanta(int nbrFluxQuanta);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  virtual AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  virtual int GetHilbertSpaceDimension ();
  
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

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of idinces 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and store result in another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors where result has to be stored
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
					       int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
					  int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent);

  
  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied<
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
 
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							     int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
						     int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
								     int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  ComplexVector* LowLevelMultipleAddMultiplyDiskStorage(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						     int firstComponent, int nbrComponent);

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
  // nbrComponent  = number of components that have to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int nbrComponent);

  // enable fast multiplication algorithm
  //
  virtual void EnableFastMultiplication();

  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int lastComponent);
  
  // enable fast multiplication algorithm using on disk cache 
  //
  // fileName = prefix of the name of the file where temporary matrix elements will be stored
  virtual void EnableFastMultiplicationWithDiskStorage(char* fileName);

  // load precalculations from a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be read
  // return value = true if no error occurs
  virtual bool LoadPrecalculation (char* fileName);

};

#endif
