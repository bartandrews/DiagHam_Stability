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


#ifndef REFLEXIONSYMMETRICPERIODIC3DHAMILTONIAN_H
#define REFLEXIONSYMMETRICPERIODIC3DHAMILTONIAN_H

#include "config.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "HilbertSpace/PeriodicXReflexionYZPeriodicThreeDOneParticle.h"

#include <iostream>

using std::ostream;

class MathematicaOutput;
class ThreeDConstantCellPotential;

class ReflexionSymmetricPeriodic3DHamiltonian : public AbstractHamiltonian
{

 protected:

  // Hilbert space associated to the system
  PeriodicXReflexionYZPeriodicThreeDOneParticle* Space;

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

  // cache for hamiltonian diagonal elements
  double* KineticElements;

  // tridimensionnal array containing all cell interaction factors
  double*** InteractionFactors;

  // wave function overlaps on a cell in a the x direction (with symmetric access type for the first two indices)
  double*** WaveFunctionOverlapX;

  // wave function overlaps on a cell in a the y direction (with symmetric access type for the first two indices)
  double** RealWaveFunctionOverlapY;
  double** ImaginaryWaveFunctionOverlapY;
 
  // wave function overlaps on a cell in a the z direction (with symmetric access type for the first two indices)
  double** RealWaveFunctionOverlapZ;
  double** ImaginaryWaveFunctionOverlapZ;

  // tridimensionnal array (with symmetric access type i > j for the two first indices) to store partial calculation to construct hamiltonian
  // elements (integration over two dimensions, first and second indices are of the form n1 * dim + m1)
  double**** RealHamiltonian;
  double**** ImaginaryHamiltonian;

 public:
  
  // constructor from data
  //
  // space = Hilbert space
  // pairX = whether basis is pair in X direction, if not impair
  // xSize = the sample length in X direction
  // ySize = the sample length in Y direction
  // zSize = the sample length in Z direction
  // mux = effective mass in X direction
  // muy = effective mass in Y direction
  // muz = effective mass in Z direction
  // nbrCellX = number of steps in X direction
  // nbrCellY = number of steps in Y direction
  // nbrCellZ = number of steps in Z direction
  // PotentielInput = pointer to a 3D potential with constant value in a cell
  // waveVectorY = wave vector of Bloch function in Y direction
  // waveVectorZ = wave vector of Bloch function in Z direction
  ReflexionSymmetricPeriodic3DHamiltonian(PeriodicXReflexionYZPeriodicThreeDOneParticle* space, bool pairX, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDConstantCellPotential* PotentialInput, double waveVectorY = 0.0, double waveVectorZ = 0.0);

  // copy constructor (without duplicating datas)
  //
  // hamiltonian = reference on hamiltonian to copy  
  ReflexionSymmetricPeriodic3DHamiltonian(const ReflexionSymmetricPeriodic3DHamiltonian& hamiltonian);

  // destructor
  //  
  ~ReflexionSymmetricPeriodic3DHamiltonian();
  
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
  // pairX = whether basis is pair in X direction, if not impair
  // waveVectorY = wave vector of Bloch function in Y direction
  // waveVectorZ = wave vector of Bloch function in Z direction
  void EvaluateInteractionFactors(bool pairX, double waveVectorY, double waveVectorZ);

  // evaluate sinus wave function overlaps on a cell in a given direction
  //
  // size = system length in the choosen direction
  // nbrStep = number of subdivision in the choosen direction
  // nbrState = number of state in the choosen direction
  // return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)
  double*** EvaluateSinusWaveFunctionOverlap(double size, int nbrStep, int nbrState);

  // evaluate cosinus wave function overlaps on a cell in a given direction
  //
  // size = system length in the choosen direction
  // nbrStep = number of subdivision in the choosen direction
  // nbrState = number of state in the choosen direction
  // return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)
  double*** EvaluateCosinusWaveFunctionOverlap(double size, int nbrStep, int nbrState);

  // evaluate the plane wave function overlap
  //
  // nbrStep = number of steps in the given direction
  // nbrState = number of states chosen for this direction
  // realArray = 2D array containing the real elements of the overlap
  // imaginaryArray = 2D array containing the imaginary elements of the overlap
  bool EvaluatePlaneWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray);

};

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

inline AbstractHilbertSpace* ReflexionSymmetricPeriodic3DHamiltonian::GetHilbertSpace ()
{
  return this->Space;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

inline int ReflexionSymmetricPeriodic3DHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Space->GetHilbertSpaceDimension ();
}

#endif
