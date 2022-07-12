////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 10/13/2005                        //
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


#ifndef QUANTUMWELLDHAMILTONIANINMAGNETICFIELD_H
#define QUANTUMWELLDHAMILTONIANINMAGNETICFIELD_H


#include "config.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/HermitianMatrix.h"
#include "MathTools/Complex.h"
#include "Tools/Potential/BinaryThreeDConstantCellPotential.h"


#include <iostream>

using std::ostream;

class MathematicaOutput;



class QuantumWellHamiltonianInMagneticField : public AbstractHamiltonian
{

 protected:

  // hamiltonian representation 
  HermitianMatrix Hamiltonian;

  // wave function basis dimension in the x direction
  int NbrStateX;
  // wave function basis dimension in the y direction
  int NbrStateY;
  // wave function basis dimension in the z direction
  int NbrStateZ;

  // number of cells in the x direction
  int NbrXCells;
  // number of cells in the y direction
  int NbrYCells;
  // number of cells in the z direction
  int NbrZCells;
  // total number of cells
  int TotalNbrCells;

  // system dimension in the x direction (in Angstrom unit)
  double XSize;
  // system dimension in the y direction (in Angstrom unit)
  double YSize;
  // system dimension in the z direction (in Angstrom unit)
  double ZSize;

  // effective mass  (in electron mass unit)
  double Mass;

  // B field value (in Tesla)
  double BField;
  // value of the magnetic length (in Angstrom unit)
  double MagneticLength;
  // value of the cyclotron energy (in meV unit)
  double CyclotronEnergy;
  // degeneracy of each Landau level
  int LandauDegeneracy;

  // z confinement in the first subband
  double ZEnergy1;
  // z confinement in the second subband
  double ZEnergy2;
  // landauIndex1 = Landau index of the first subband
  int LandauIndex1;
  // landauIndex2 = Landau index of the second subband
  int LandauIndex2;

  
  double MailleParameter;
  // conduction band offset between GaAs and InAs
  double BandOffset;
  // In/Ga dopage ratio (=x with Ga_(1-x) In_x As)
  double InDopage;
  // number of crystal elementary cell in the qunatum well
  int NbrCells;

  // description of the potential
  ThreeDConstantCellPotential* Potential;
  
  double* GaXPosition;
  double* GaYPosition;
  double* GaZPosition;
  double* InXPosition;
  double* InYPosition;
  double* InZPosition;
  

 public:

  // constructor from default data
  //
  // xSize = system dimension in the x direction (in Angstrom unit)
  // ySize = system dimension in the y direction (in Angstrom unit)
  // zSize = system dimension in the z direction (in Angstrom unit)
  // mass = effective mass in the x direction (in electron mass unit)
  // bField = B field value (in Tesla)
  // zEnergy1 = z confinement in the first subband
  // zEnergy2 = z confinement in the second subband
  // landauIndex1 = Landau index of the first subband
  // landauIndex2 = Landau index of the second subband
  // mailleParameter =
  // bandOffset = conduction band offset between GaAs and InAs
  // inDopage = In/Ga dopage ratio (=x with Ga_(1-x) In_x As)
  // potentialDescription = name of the file that contains the potential description (null if the potential has to be evaluated)
  QuantumWellHamiltonianInMagneticField(double xSize, double ySize, double zSize, double mass, double bField, double zEnergy1, double zEnergy2,
					int landauIndex1, int landauIndex2, double mailleParameter, double bandOffset, double inDopage, char* potentialDescription = 0);

  // copy constructor (without duplicating datas)
  //
  // hamiltonian = reference on hamiltonian to copy
  QuantumWellHamiltonianInMagneticField(const QuantumWellHamiltonianInMagneticField& hamiltonian);

  // destructor
  //
  ~QuantumWellHamiltonianInMagneticField();

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

  // store Hamiltonian into an hermitian matrix
  //
  // M = reference on matrix where Hamiltonian has to be stored
  // return value = reference on  corresponding hermitian matrix
  HermitianMatrix& GetHamiltonian (HermitianMatrix& M);

  // save potential on disk
  // 
  // filename = name of the file (with path) where potential has to be saved
  // return value = true if no error occured
  bool SavePotential(char* filename);
  
 private:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

};

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

inline AbstractHilbertSpace* QuantumWellHamiltonianInMagneticField::GetHilbertSpace ()
{
  return 0;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

inline int QuantumWellHamiltonianInMagneticField::GetHilbertSpaceDimension ()
{
  return (this->LandauDegeneracy * 2);
}
  
#endif
