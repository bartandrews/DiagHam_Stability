////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of NDMAP spin chain hamiltonian                   //
//                                                                            //
//                        last modification : 23/10/2003                      //
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


#ifndef NDMAPSPINCHAINHAMILTONIAN_H
#define NDMAPSPINCHAINHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class NDMAPSpinChainHamiltonian : public AbstractHamiltonian
{

  friend class NDMAPPrecalculationOperation;

 protected:

  // pointer to the Hilbert space describing the spin chain 
  AbstractSpinChainWithTranslations* Chain;

  // architecture used for precalculation
  AbstractArchitecture* Architecture;

  // coupling constant between spin in the +- and -+ direction 
  double J;
  // coupling constant between spin in the zz direction 
  double Jz;
  // half coupling constant between spin in the +- and -+ direction 
  double HalfJ;
  // magnetic field value times the coupling constant with the magnetic field (g mu_b) in the direction parallel to the chain
  double ParallelMagneticField;
  // magnetic field value times the coupling constant with the magnetic field (g mu_b) in the direction perpendicular to the chain
  double PerpendicularMagneticField;
   // half the value of the magnetic field times the coupling constant with the magnetic field (g mu_b) in the direction perpendicular to the chain
  double HalfPerpendicularMagneticField;
  // single ion anisotropy constant
  double D;
  // single ion in-plane anisotropy constant
  double E;
  // half the single ion in-plane anisotropy constant
  double HalfE;

  // number of spin 
  int NbrSpin;

  // array conating all matrix diagonal elements
  double* SzSzContributions;

  //array containing all the cosinus that are needed when computing matrix elements
  double* CosinusTable;
  //array containing all the sinus that are needed when computing matrix elements
  double* SinusTable;

  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  // step between each precalculated index
  int FastMultiplicationStep;
  // number of non-null term in the hamiltonian for each state
  int* NbrInteractionPerComponent;
  // index of the state obtained for each term of the hamiltonian when applying on a given state
  int** InteractionPerComponentIndex;
  // multiplicative coefficient obtained for each term of the hamiltonian when applying on a given state and with a given destination state
  double** InteractionPerComponentCoefficient;
  // number of translations obtained for each term of the hamiltonian when applying on a given state and with a given destination state
  int** InteractionPerComponentNbrTranslations;

 public:

  // constructor from default datas
  //
  // chain = pointer to Hilbert space of the associated system
  // nbrSpin = number of spin
  // j = coupling constants between spins
  // jz = coupling constants between spins in the z direction
  // parallelMagneticField = magnetic field value times the coupling constant with the magnetic field (g mu_b) in the direction parallel to the chain
  // perpendicularMagneticField = magnetic field value times the coupling constant with the magnetic field (g mu_b) in the direction perpendicular to the chain
  // d = single ion anisotropy constant
  // e = single ion in-plane anisotropy constant
  // architecture = architecture to use for precalculation
  // memory = amount of memory (in bytes) that can allocated for precalcutation
  NDMAPSpinChainHamiltonian(AbstractSpinChainWithTranslations* chain, int nbrSpin, double j, double jz, double parallelMagneticField, 
			    double perpendicularMagneticField, double d = 0.0, double e = 0.0,  AbstractArchitecture* architecture = 0, unsigned long memory = 100000000);

  // destructor
  //
  ~NDMAPSpinChainHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set chain
  // 
  // chain = pointer on Hilbert space of the associated system
  // return value = reference on current Hamiltonian
  NDMAPSpinChainHamiltonian& SetChain(AbstractSpinChainWithTranslations* chain);

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

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

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
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				     int firstComponent, int nbrComponent);
 
  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  List<Matrix*> RightInteractionOperators();  

  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, NDMAPSpinChainHamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, NDMAPSpinChainHamiltonian& H);

 private:
 
  // evaluate all cosinus/sinus that are needed when computing matrix elements
  //
  void EvaluateCosinusTable();

  // evaluate all matrix elements
  //   
  void EvaluateDiagonalMatrixElements();

  // test the amount of memory needed for fast multiplication algorithm
  //
  // allowedMemory = amount of memory that cam be allocated for fast multiplication
  // return value = amount of memory needed
  long FastMultiplicationMemory(long allowedMemory);

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element
  long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm
  //
  void EnableFastMultiplication();

  // enable fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  void PartialEnableFastMultiplication(int firstComponent, int lastComponent);

};

#endif
