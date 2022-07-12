////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 26/02/2003                        //
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


#ifndef QUANTUMDOTS3DHAMILTONIAN_H
#define QUANTUMDOTS3DHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "HilbertSpace/Confined3DOneParticle.h"
#include "Potential/ThreeDPotential.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class QuantumDots3DHamiltonian : public AbstractHamiltonian
{

 protected:

  // Hilbert space associated to the system
  Confined3DOneParticle* Space;

  // wave function basis dimension in the x direction
  int NbrStateX;
  // wave function basis dimension in the y direction
  int NbrStateY;
  // wave function basis dimension in the z direction
  int NbrStateZ;

  // number of cells in the x direction
  int NbrCellX;
  // number of cells in the y direction
  int NbrCellY;
  // number of cells in the z direction
  int NbrCellZ;
  // total number of cells
  int TotalNbrCells;

  // system dimension in the x direction (in Angstrom unit)
  double XSize;
  // system dimension in the y direction (in Angstrom unit)
  double YSize;
  // system dimension in the z direction (in Angstrom unit)
  double ZSize;

  // number of layers before the active layers
  int LeftNumber;
  // region size in the z direction where potential is constant in the plane (regions before active layers)
  double* PreConstantRegionSize;
  // value of the potential in the region before the active layers
  double* PreConstantRegionPotential;
  // number of layers after the active layers
  int RightNumber;
  // region size in the z direction where potential is constant in every direction (regions after active layers)
  double* PostConstantRegionSize;
  // value of the potential in the region after the active layers
  double* PostConstantRegionPotential;

  // effective mass in the x direction (in electron mass unit)
  double Mux;
  // effective mass in the y direction (in electron mass unit)
  double Muy;
  // effective mass in the z direction (in electron mass unit)
  double Muz;

  // cache for hamiltonian diagonal elements
  double* DiagonalElements;

  // flag to indicate how many dimensions are precalculated for hamiltonian interaction elements
  int NbrPrecalculatedDimension;

  // tridimensionnal array containing all cell interaction factors
  double*** InteractionFactors;
  //ThreeDPotential* InteractionFactors;

  // wave function overlaps on a cell in a the x direction (with symmetric access type for the first two indices)
  double*** WaveFunctionOverlapX;
  // wave function overlaps on a cell in a the y direction (with symmetric access type for the first two indices)
  double*** WaveFunctionOverlapY;
  // wave function overlaps on a cell in a the z direction (with symmetric access type for the first two indices)
  double*** WaveFunctionOverlapZ;

  // bidimensionnal array (with symmetric access type i > j) to store hamiltonian when all matrix elements can be stored in memory
  double** FullPrecalculatedHamiltonian;

  // tridimensionnal array (with symmetric access type i > j for the two first indices) to store partial calculation to construct hamiltonian
  // elements (integration over two dimensions, first and second indices are of the form n1 * dim + m1)
  double*** Partial2DPrecalculatedHamiltonian;
  // bidimensionnal array to store partial diagonal calculation to construct hamiltonian
  // elements (integration over two dimensions, first and second indices are of the form n1 * dim + m1)
  double** Partial2DDiagonalPrecalculatedHamiltonian;
  // elements (intergration over z, the first end second indices are of p1, p2)
  double ** PartialZPrecalculatedHamiltonian;

 public:

  // constructor from default data
  //
  // space = Hilbert space associated to the system
  // xSize = system dimension in the x direction (in Angstrom unit)
  // ySize = system dimension in the y direction (in Angstrom unit)
  // zSize = system dimension in the z direction (in Angstrom unit)
  // preConstantRegionSize = region size in the z direction where potential is constant in every direction (region before gradiant zone)
  // postConstantRegionSize = region size in the z direction where potential is constant in every direction (region after gradiant zone)
  // postConstantRegionPotential = value of the potential in the region after the gradiant zone
  // mux = effective mass in the x direction (in electron mass unit)
  // muy = effective mass in the y direction (in electron mass unit)
  // muz = effective mass in the z direction (in electron mass unit)
  // nbrCellX = number of cells in the x direction
  // nbrCellY = number of cells in the y direction
  // nbrCellZ = number of cells in the z direction
  // overlapingFactors = tridimensionnal array where overlaping factors are stored
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  QuantumDots3DHamiltonian(Confined3DOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDPotential* PotentialInput, int memory = -1);

  // copy constructor (without duplicating datas)
  //
  // hamiltonian = reference on hamiltonian to copy
  QuantumDots3DHamiltonian(const QuantumDots3DHamiltonian& hamiltonian);

  // destructor
  //
  ~QuantumDots3DHamiltonian();

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
  friend ostream& operator << (ostream& Str, QuantumDots3DHamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, QuantumDots3DHamiltonian& H);

  // determine the maximal value of partial diagonal array
  //
  // return = the wanted value
  double MaxPartialDiagonalElement();

 private:
 
  // evaluate all interaction factors
  //   
  // memory = amount of memory available to store precalculated values
  void EvaluateInteractionFactors(int memory);

  // evaluate wave function overlaps on a cell in a given direction
  //
  // size = system length in the choosen direction
  // nbrStep = number of subdivision in the choosen direction
  // nbrState = number of state in the choosen direction
  // memory = reference on current memory usage (will be increment with memory used to store evaluated overlap)
  // return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)
  double*** EvaluateWaveFunctionOverlap(double size, int nbrStep, int nbrState, int& memory);

  // evaluate wave function overlaps on a cell in the z direction
  //
  // memory = reference on current memory usage (will be increment with memory used to store evaluated overlap)
  // return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)
  double*** EvaluateWaveFunctionZOverlap(int& memory);

};

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

inline AbstractHilbertSpace* QuantumDots3DHamiltonian::GetHilbertSpace ()
{
  return this->Space;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

inline int QuantumDots3DHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Space->GetHilbertSpaceDimension ();
}
  
#endif
