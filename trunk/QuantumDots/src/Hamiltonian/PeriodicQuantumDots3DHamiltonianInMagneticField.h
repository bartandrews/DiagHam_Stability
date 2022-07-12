////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2004 Duc-Phuong Nguyen                    //
//                                                                            //
//                                                                            //
//     class of hamiltonian associated quantum dots in a magnetic field       //
//                                                                            //
//                      last modification : 19/04/2004                        //
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


#ifndef PERIODICQUANTUMDOTS3DHAMILTONIANINMAGNETICFIELD_H
#define PERIODICQUANTUMDOTS3DHAMILTONIANINMAGNETICFIELD_H


#include "config.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "HilbertSpace/PeriodicThreeDOneParticle.h"
#include <iostream>


using std::ostream;


class MathematicaOutput;
class ThreeDConstantCellPotential;

class PeriodicQuantumDots3DHamiltonianInMagneticField : public AbstractHamiltonian
{

 protected:

  // Hilbert space associated to the system
  PeriodicThreeDOneParticle* Space;

  // wave function basis dimension in the x direction
  int NbrStateX;  

  int LowerImpulsionX; 

  // wave function basis dimension in the y direction
  int NbrStateY;

  int LowerImpulsionY; 

  // wave function basis dimension in the z direction
  int NbrStateZ;

  int LowerImpulsionZ;  

  // number of cells in the x direction
  int NbrCellX;
  // number of cells in the y direction
  int NbrCellY;
  // number of cells in the z direction
  int NbrCellZ;

  // system dimension in the x direction (in Angstrom unit)
  double XSize;
  // system dimension in the y direction (in Angstrom unit)
  double YSize;
  // system dimension in the z direction (in Angstrom unit)
  double ZSize;

  // effective mass in the x direction (in electron mass unit)
  double Mux;
  // effective mass in the y direction (in electron mass unit)
  double Muy;
  // effective mass in the z direction (in electron mass unit)
  double Muz;

  // magnetic field components (in Tesla unit)
  //double Bx;
  //double By;
  double Bz;

  // cache for hamiltonian diagonal elements
  double* KineticElements;

  // flag to indicate how many dimensions are precalculated for hamiltonian interaction elements
  int NbrPrecalculatedDimension;

  // tridimensionnal array containing all cell interaction factors
  double*** InteractionFactors;
  
  // paramagnetic terms
  double** RealParamagneticTermPxY;
  double** ImaginaryParamagneticTermPxY;
  double** RealParamagneticTermPyX;
  double** ImaginaryParamagneticTermPyX;

  // tridimensionnal array to store hamiltonian elements
  double*** RealPrecalculatedHamiltonian;
  double*** ImaginaryPrecalculatedHamiltonian;

 public:

  // constructor from data
  //
  // space = Hilbert space
  // xSize = the sample length in X direction
  // ySize = the sample length in Y direction
  // zSize = the sample length in Z direction
  // mux = effective mass in X direction
  // muy = effective mass in Y direction
  // muz = effective mass in Z direction
  // bx = X magnetic field component
  // by = Y magnetic field component
  // bz = Z magnetic field component
  // PotentialInput = pointer to a 3D potential with constant value in a cell
  // waveVectorZ = wave vector of Bloch function in Z direction
  PeriodicQuantumDots3DHamiltonianInMagneticField(PeriodicThreeDOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, double bx, double by, double bz, ThreeDConstantCellPotential* PotentialInput, double waveVectorZ = 0.0);

  // copy constructor (without duplicating datas)
  //
  // hamiltonian = reference on hamiltonian to copy  
  PeriodicQuantumDots3DHamiltonianInMagneticField(const PeriodicQuantumDots3DHamiltonianInMagneticField& hamiltonian);

  // destructor
  //
  ~PeriodicQuantumDots3DHamiltonianInMagneticField();
  
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
  
  // evaluate confinement potential factors
  //
  // waveVectorZ = wave vector of Bloch function in Z direction 
  void EvaluateConfinementPotentialFactors(double waveVectorZ);

  // evaluate magnetic field factors
  //
  void EvaluateMagneticFieldFactors();
  
  // evaluate the wave function overlap
  //
  // nbrStep = number of steps in the given direction
  // nbrState = number of states chosen for this direction
  // realArray = 2D array containing the real elements of the overlap
  // imaginaryArray = 2D array containing the imaginary elements of the overlap
  bool EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray);

  // evaluate the mean position operator in a given direction
  //
  // size = size of sample in the given direction
  // nbrState = number of states wanted
  // real = reference to the array of mean position's real component
  // imaginary = reference to the array of mean position's imaginary component
  // realSquared = reference to the array of X²'s real component
  // imaginarySquared = reference to the array of X²'s imaginary component
  // return = true if successful, otherwise false  
  bool EvaluateMeanPositionOperator(double size, int nbrState, double* &real, double* &imaginary, double* &realSquared, double* &imaginarySquared);

};

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

inline AbstractHilbertSpace* PeriodicQuantumDots3DHamiltonianInMagneticField::GetHilbertSpace ()
{
  return this->Space;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

inline int PeriodicQuantumDots3DHamiltonianInMagneticField::GetHilbertSpaceDimension ()
{
  return this->Space->GetHilbertSpaceDimension ();
}

#endif
