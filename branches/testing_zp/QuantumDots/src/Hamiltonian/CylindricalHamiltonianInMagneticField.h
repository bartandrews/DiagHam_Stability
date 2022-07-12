////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2003 Duc-Phuong Nguyen                    //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 25/03/2003                        //
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


#ifndef CYLINDRICALHAMILTONIANINMAGNETICFIELD_H
#define CYLINDRICALHAMILTONIANINMAGNETICFIELD_H


#include "config.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "HilbertSpace/PlanarRotationSymmetryZPeriodicOneParticle.h"

#include <iostream>

using std::ostream;

class MathematicaOutput;
class ThreeDConstantCylinderPotential;

class CylindricalHamiltonianInMagneticField : public AbstractHamiltonian
{

 protected:

  // Hilbert space associated to the system
  PlanarRotationSymmetryZPeriodicOneParticle* Space;

  // quantum number of kinetic momentum in Z direction
  int NumberM;

  // number of Landau states in plane
  int NbrStateR;

  // wave function basis dimension in the z direction
  int NbrStateZ;
  int LowerImpulsionZ;

  // system dimension in the z direction (in Angstrom unit)
  double ZSize;

  // effective mass in the plane
  double Mur;
  // effective mass in the z direction (in electron mass unit)
  double Muz;

  // partial diagonal terms (z kinetic energy and magnetic energy)
  double* PartialDiagonalElement;

  // array containing all the Hamiltonian elements with symmetry exploitation
  double*** RealHamiltonian;
  double*** ImaginaryHamiltonian;

 public:
  
  // constructor from data
  //
  // space = Hilbert space
  // mur = effective mass in plane
  // muz = effective mass in Z direction
  // bz = Z magnetic field component
  // waveVectorZ = wave vector of Bloch function in Z direction
  // PotentialInput = pointer to a 3D potential with constant value in a cell
  CylindricalHamiltonianInMagneticField(PlanarRotationSymmetryZPeriodicOneParticle* space, double mur, double muz, double bz, double waveVectorZ, ThreeDConstantCylinderPotential* PotentialInput);

  // copy constructor (without duplicating datas)
  //
  // hamiltonian = reference on hamiltonian to copy  
  CylindricalHamiltonianInMagneticField(const CylindricalHamiltonianInMagneticField& hamiltonian);

  // destructor
  //  
  ~CylindricalHamiltonianInMagneticField();
  
  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  
  AbstractHamiltonian* Clone ();

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

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of idinces 
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
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent);

  // determine the maximal value of partial diagonal array
  //
  // return = the wanted value
  double MaxPartialDiagonalElement();
 
 private:

  // evaluate all interaction factors
  //
  // Bz = magnetic field component in Z direction
  // waveVectorZ = wave vector of Bloch function in Z direction
  // potential = pointer to the potential
  void EvaluateInteractionFactors(double Bz, double waveVectorZ, ThreeDConstantCylinderPotential* &potential);

  // evaluate the plane wave function overlap
  //
  // potential = pointer to the potential
  // nbrState = number of states chosen for this direction
  // realArray = 2D array containing the real elements of the overlap
  // imaginaryArray = 2D array containing the imaginary elements of the overlap
  bool EvaluatePlaneWaveFunctionOverlap(ThreeDConstantCylinderPotential* &potential, int nbrState, double** &realArray, double** &imaginaryArray);

};

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

inline AbstractHilbertSpace* CylindricalHamiltonianInMagneticField::GetHilbertSpace ()
{
  return this->Space;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

inline int CylindricalHamiltonianInMagneticField::GetHilbertSpaceDimension ()
{
  return this->Space->GetHilbertSpaceDimension ();
}

#endif
