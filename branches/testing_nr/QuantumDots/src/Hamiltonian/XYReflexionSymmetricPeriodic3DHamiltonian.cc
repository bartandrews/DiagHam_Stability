////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2004 Duc-Phuong Nguyen                    //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 25/03/2004                        //
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


#include "config.h"
#include "HilbertSpace/PeriodicReflexionSymmetryOneDOneParticle.h"
#include "Hamiltonian/XYReflexionSymmetricPeriodic3DHamiltonian.h"
#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"
#include "Tools/Potential/ThreeDConstantCellPotential.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>

using std::ostream;
using std::cout;
using std::endl;


#define PERIODIC_HAMILTONIAN_FACTOR 150.4
#define BLOCH_FACTOR 7.644


// constructor from data
//
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
// waveVectorZ = wave vector of Bloch function in Z direction

XYReflexionSymmetricPeriodic3DHamiltonian::XYReflexionSymmetricPeriodic3DHamiltonian(PeriodicXYReflexionZPeriodicThreeDOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDConstantCellPotential* PotentialInput, double waveVectorZ)
{
  this->Space = space;
  this->XSize = xSize;
  this->YSize = ySize;
  this->ZSize = zSize;
  this->Mux = mux;
  this->Muy = muy;
  this->Muz = muz;
  this->NbrCellX = nbrCellX;
  this->NbrCellY = nbrCellY;
  this->NbrCellZ = nbrCellZ;

  this->NbrStateX = this->Space->GetNbrStateX();
  this->LowerImpulsionX = this->Space->GetLowerImpulsionX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->LowerImpulsionY = this->Space->GetLowerImpulsionY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->LowerImpulsionZ = this->Space->GetLowerImpulsionZ();

  this->InteractionFactors = new double** [this->NbrCellZ];
  int CenterX = this->NbrCellX / 2; int CenterY = this->NbrCellY / 2;
  for (int k = 0; k < this->NbrCellZ; ++k)
    {
      this->InteractionFactors[k] = new double* [this->NbrCellY];
      for (int j = 0; j < this->NbrCellY; ++j)
	{
	  this->InteractionFactors[k][j] = new double [this->NbrCellX];
	  for (int i = 0; i < this->NbrCellX; ++i)
	    {
	      this->InteractionFactors[k][j][i] = PotentialInput->GetPotential(i, j, k);

	      if ((i > CenterX) && (j > CenterY))
	        if ((this->InteractionFactors[k][j][i] != this->InteractionFactors[k][this->NbrCellY - 1 - j][i]) || (this->InteractionFactors[k][j][i] != this->InteractionFactors[k][j][this->NbrCellX - 1 - i]) ||
		    (this->InteractionFactors[k][j][i] != this->InteractionFactors[k][this->NbrCellY - 1 - j][this->NbrCellX - 1 - i]))
		  {
		    cout << "The potential is not reflexion symmetric in X or Y direction. Exit now!" << endl;
		    exit(0);
		  }
	      
	    }
	}	        
    }
  bool pairX = ((PeriodicReflexionSymmetryOneDOneParticle*) this->Space->GetStateX())->IsEven ();
  bool pairY = ((PeriodicReflexionSymmetryOneDOneParticle*) this->Space->GetStateY())->IsEven ();

  cout << "Hamiltonian dimension: " << this->Space->GetHilbertSpaceDimension () << endl;
  cout << "Evaluation of Hamiltionian elements ..." << endl;
  this->EvaluateInteractionFactors(pairX, pairY, waveVectorZ);
  cout << "Evaluation finished ..." << endl;
}

// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

XYReflexionSymmetricPeriodic3DHamiltonian::XYReflexionSymmetricPeriodic3DHamiltonian(const XYReflexionSymmetricPeriodic3DHamiltonian& hamiltonian)
{
  this->Space = hamiltonian.Space;
  this->XSize = hamiltonian.XSize;
  this->YSize = hamiltonian.YSize;
  this->ZSize = hamiltonian.ZSize;
  this->Mux = hamiltonian.Mux;
  this->Muy = hamiltonian.Muy;
  this->Muz = hamiltonian.Muz;
  this->NbrCellX = hamiltonian.NbrCellX;
  this->NbrCellY = hamiltonian.NbrCellY;
  this->NbrCellZ = hamiltonian.NbrCellZ;
  this->NbrStateX = this->Space->GetNbrStateX();
  this->NbrStateY = this->Space->GetNbrStateY();
  this->NbrStateZ = this->Space->GetNbrStateZ();
  this->LowerImpulsionX = this->Space->GetLowerImpulsionX();
  this->LowerImpulsionY = this->Space->GetLowerImpulsionY();
  this->LowerImpulsionZ = this->Space->GetLowerImpulsionZ();
  this->KineticElements = hamiltonian.KineticElements;
  this->InteractionFactors = hamiltonian.InteractionFactors;
  this->WaveFunctionOverlapX = hamiltonian.WaveFunctionOverlapX;
  this->WaveFunctionOverlapY = hamiltonian.WaveFunctionOverlapY;
  this->RealWaveFunctionOverlapZ = hamiltonian.RealWaveFunctionOverlapZ;
  this->ImaginaryWaveFunctionOverlapZ = hamiltonian.ImaginaryWaveFunctionOverlapZ;
  this->RealHamiltonian =  hamiltonian.RealHamiltonian;
  this->ImaginaryHamiltonian = hamiltonian.ImaginaryHamiltonian;
}

// destructor
//

XYReflexionSymmetricPeriodic3DHamiltonian::~ XYReflexionSymmetricPeriodic3DHamiltonian()
{  
  delete[] this->KineticElements;
  delete[] this->InteractionFactors;
  delete[] this->WaveFunctionOverlapX;
  delete[] this->WaveFunctionOverlapY;
  delete[] this->RealWaveFunctionOverlapZ;
  delete[] this->ImaginaryWaveFunctionOverlapZ;
  delete[] this->RealHamiltonian;
  delete[] this->ImaginaryHamiltonian;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* XYReflexionSymmetricPeriodic3DHamiltonian::Clone ()
{
  return new XYReflexionSymmetricPeriodic3DHamiltonian(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void XYReflexionSymmetricPeriodic3DHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void XYReflexionSymmetricPeriodic3DHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Space->GetHilbertSpaceDimension (); ++i)
    this->KineticElements[i] += shift;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex XYReflexionSymmetricPeriodic3DHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
{
  double x = 0.0;
  int dim = this->Space->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex XYReflexionSymmetricPeriodic3DHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& XYReflexionSymmetricPeriodic3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Space->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of idinces 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& XYReflexionSymmetricPeriodic3DHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,
						       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination.Re(i) = 0.0;
      vDestination.Im(i) = 0.0;
    }
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& XYReflexionSymmetricPeriodic3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  int tmpDim = this->NbrStateY * this->NbrStateZ;
  int** RI = new int* [this->NbrStateX];
  int m1 = 0, n1 = 0, p1 = 0, m2 = 0, n2 = 0;
  for (m1 = 0; m1 < this->NbrStateX; ++m1)
    {
      RI[m1] = new int [this->NbrStateY];
      for (n1 = 0; n1 < this->NbrStateY; ++n1)
	RI[m1][n1] = m1 * tmpDim + n1 * this->NbrStateZ;
    }

  int Index1 = 0, Index2 = 0, Index3 = 0, Index4 = 0, Index5 = 0, Index6 = 0, Index7 = 0, Index8 = 0;
  double TmpRe = 0.0, TmpIm = 0.0;
  int IndexZ = 0, tmpIndexZ = 0;
  double* TmpRealHamiltonian; double* TmpImaginaryHamiltonian;
  double VSourceRe1 = 0.0, VSourceIm1 = 0.0, VSourceRe3 = 0.0, VSourceIm3 = 0.0, VSourceRe5 = 0.0, VSourceIm5 = 0.0, VSourceRe7 = 0.0, VSourceIm7 = 0.0;
  double TmpTotalRe1 = 0.0, TmpTotalIm1 = 0.0, TmpTotalRe3 = 0.0, TmpTotalIm3 = 0.0, TmpTotalRe5 = 0.0, TmpTotalIm5 = 0.0, TmpTotalRe7 = 0.0, TmpTotalIm7 = 0.0;
  // Index1 = (m1 * N + n1) * H + p1; Index3 = (m1 * N + n2) * H + p1; Index5 = (m2 * N + n1) * H + p1; Index7 = (m2 * N + n2) * H + p1
  // Index2 = (m2 * N + n2) * H + p2; Index4 = (m2 * N + n1) * H + p2; Index6 = (m1 * N + n2) * H + p2; Index8 = (m1 * N + n1) * H + p2
  for (m1 = 0; m1 < this->NbrStateX; ++m1)
    {
      for (n1 = 0; n1 < this->NbrStateY; ++n1)
	{
	  for (m2 = 0; m2 < m1; ++m2)
	    {
	      for (n2 = 0; n2 < n1; ++n2)
		{
		  Index1 = RI[m1][n1]; // Index1 = (m1 * N + n1) * H
		  Index3 = RI[m1][n2]; // Index3 = (m1 * N + n2) * H
		  Index5 = RI[m2][n1]; // Index5 = (m2 * N + n1) * H
		  Index7 = RI[m2][n2]; // Index7 = (m2 * N + n2) * H
		  TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m2][n2];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m2][n2];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      Index2 = RI[m2][n2]; // Index2 = (m2 * N + n2) * H
		      Index4 = RI[m2][n1]; // Index4 = (m2 * N + n1) * H		      
		      Index6 = RI[m1][n2]; // Index6 = (m1 * N + n2) * H
		      Index8 = RI[m1][n1]; // Index8 = (m1 * N + n1) * H

		      VSourceRe1 = vSource.Re(Index1); VSourceIm1 = vSource.Im(Index1);
		      VSourceRe3 = vSource.Re(Index3); VSourceIm3 = vSource.Im(Index3);
		      VSourceRe5 = vSource.Re(Index5); VSourceIm5 = vSource.Im(Index5);
		      VSourceRe7 = vSource.Re(Index7); VSourceIm7 = vSource.Im(Index7);

		      TmpTotalRe1 = 0.0; TmpTotalIm1 = 0.0; TmpTotalRe3 = 0.0; TmpTotalIm3 = 0.0;
		      TmpTotalRe5 = 0.0; TmpTotalIm5 = 0.0; TmpTotalRe7 = 0.0; TmpTotalIm7 = 0.0;

		      for (IndexZ = p1; IndexZ > 0; --IndexZ)
			{
			  TmpRe = TmpRealHamiltonian[IndexZ];
			  TmpIm = TmpImaginaryHamiltonian[IndexZ];

			  vDestination.Re(Index2) += (TmpRe * VSourceRe1 - TmpIm * VSourceIm1);
			  vDestination.Im(Index2) += (TmpRe * VSourceIm1 + TmpIm * VSourceRe1);
			  vDestination.Re(Index4) += (TmpRe * VSourceRe3 - TmpIm * VSourceIm3);
			  vDestination.Im(Index4) += (TmpRe * VSourceIm3 + TmpIm * VSourceRe3);
			  vDestination.Re(Index6) += (TmpRe * VSourceRe5 - TmpIm * VSourceIm5);
			  vDestination.Im(Index6) += (TmpRe * VSourceIm5 + TmpIm * VSourceRe5);
			  vDestination.Re(Index8) += (TmpRe * VSourceRe7 - TmpIm * VSourceIm7);
			  vDestination.Im(Index8) += (TmpRe * VSourceIm7 + TmpIm * VSourceRe7);

			  TmpTotalRe1 += (TmpRe * vSource.Re(Index2) + TmpIm * vSource.Im(Index2));
			  TmpTotalIm1 += (TmpRe * vSource.Im(Index2) - TmpIm * vSource.Re(Index2));
			  TmpTotalRe3 += (TmpRe * vSource.Re(Index4) + TmpIm * vSource.Im(Index4));
			  TmpTotalIm3 += (TmpRe * vSource.Im(Index4) - TmpIm * vSource.Re(Index4));
			  TmpTotalRe5 += (TmpRe * vSource.Re(Index6) + TmpIm * vSource.Im(Index6));
			  TmpTotalIm5 += (TmpRe * vSource.Im(Index6) - TmpIm * vSource.Re(Index6));
			  TmpTotalRe7 += (TmpRe * vSource.Re(Index8) + TmpIm * vSource.Im(Index8));
			  TmpTotalIm7 += (TmpRe * vSource.Im(Index8) - TmpIm * vSource.Re(Index8));
			  
			  ++Index2; ++Index4; ++Index6; ++Index8;
			}
		      // case p2 = p1
		      TmpRe = TmpRealHamiltonian[IndexZ];
		      TmpIm = TmpImaginaryHamiltonian[IndexZ];
		      
		      vDestination.Re(Index2) += (TmpRe * VSourceRe1 + TmpIm * VSourceIm1);
		      vDestination.Im(Index2) += (TmpRe * VSourceIm1 - TmpIm * VSourceRe1);
		      vDestination.Re(Index4) += (TmpRe * VSourceRe3 + TmpIm * VSourceIm3);
		      vDestination.Im(Index4) += (TmpRe * VSourceIm3 - TmpIm * VSourceRe3);

		      vDestination.Re(Index1) += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
		      vDestination.Im(Index1) += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));
		      vDestination.Re(Index3) += (TmpRe * vSource.Re(Index4) - TmpIm * vSource.Im(Index4));
		      vDestination.Im(Index3) += (TmpRe * vSource.Im(Index4) + TmpIm * vSource.Re(Index4));    

		      vDestination.Re(Index1) += TmpTotalRe1; vDestination.Im(Index1) += TmpTotalIm1;
		      vDestination.Re(Index3) += TmpTotalRe3; vDestination.Im(Index3) += TmpTotalIm3;
		      vDestination.Re(Index5) += TmpTotalRe5; vDestination.Im(Index5) += TmpTotalIm5;
		      vDestination.Re(Index7) += TmpTotalRe7; vDestination.Im(Index7) += TmpTotalIm7;

		      ++Index1; ++Index3; ++Index5; ++Index7;
		    }		  
		}
	      // case n2 = n1
	      Index1 = RI[m1][n1]; // Index1 = (m1 * N + n1) * H	      
	      Index5 = RI[m2][n1]; // Index5 = (m2 * N + n1) * H
	      TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m2][n1];
	      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m2][n1];
	      for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		{
		  Index2 = RI[m2][n1]; // Index2 = (m2 * N + n1) * H
		  Index6 = RI[m1][n1]; // Index6 = (m1 * N + n1) * H

		  VSourceRe1 = vSource.Re(Index1); VSourceIm1 = vSource.Im(Index1);
		  VSourceRe5 = vSource.Re(Index5); VSourceIm5 = vSource.Im(Index5);

		  TmpTotalRe1 = 0.0; TmpTotalIm1 = 0.0; TmpTotalRe5 = 0.0; TmpTotalIm5 = 0.0;

		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {
		      TmpRe = TmpRealHamiltonian[IndexZ];
		      TmpIm = TmpImaginaryHamiltonian[IndexZ];

		      vDestination.Re(Index2) += (TmpRe * VSourceRe1 - TmpIm * VSourceIm1);
		      vDestination.Im(Index2) += (TmpRe * VSourceIm1 + TmpIm * VSourceRe1);
		      vDestination.Re(Index6) += (TmpRe * VSourceRe5 - TmpIm * VSourceIm5);
		      vDestination.Im(Index6) += (TmpRe * VSourceIm5 + TmpIm * VSourceRe5);

		      TmpTotalRe1 += (TmpRe * vSource.Re(Index2) + TmpIm * vSource.Im(Index2));
		      TmpTotalIm1 += (TmpRe * vSource.Im(Index2) - TmpIm * vSource.Re(Index2));
		      TmpTotalRe5 += (TmpRe * vSource.Re(Index6) + TmpIm * vSource.Im(Index6));
		      TmpTotalIm5 += (TmpRe * vSource.Im(Index6) - TmpIm * vSource.Re(Index6));
		      
		      ++Index2; ++Index6;
		    }
		  // case n2 = n1 && p2 = p1
		  TmpRe = TmpRealHamiltonian[IndexZ];
		  TmpIm = TmpImaginaryHamiltonian[IndexZ];		  
 
		  vDestination.Re(Index2) += (TmpRe * VSourceRe1 + TmpIm * VSourceIm1);
		  vDestination.Im(Index2) += (TmpRe * VSourceIm1 - TmpIm * VSourceRe1);

		  vDestination.Re(Index1) += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
		  vDestination.Im(Index1) += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2));   
		  
		  vDestination.Re(Index1) += TmpTotalRe1; vDestination.Im(Index1) += TmpTotalIm1;
		  vDestination.Re(Index5) += TmpTotalRe5; vDestination.Im(Index5) += TmpTotalIm5;

		  ++Index1; ++Index5;
		}
	    }
	  //case m2 = m1
	  for (n2 = 0; n2 < n1; ++n2)
	    {			
	      Index1 = RI[m1][n1]; // Index1 = (m1 * N + n1) * H
	      Index3 = RI[m1][n2]; // Index3 = (m1 * N + n2) * H		  		  
	      TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m1][n2];
	      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m1][n2];
	      for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		{
		  Index2 = RI[m1][n2]; // Index2 = (m1 * N + n2) * H
		  Index4 = RI[m1][n1]; // Index4 = (m1 * N + n1) * H

		  VSourceRe1 = vSource.Re(Index1); VSourceIm1 = vSource.Im(Index1);
		  VSourceRe3 = vSource.Re(Index3); VSourceIm3 = vSource.Im(Index3);

		  TmpTotalRe1 = 0.0; TmpTotalIm1 = 0.0; TmpTotalRe3 = 0.0; TmpTotalIm3 = 0.0;

		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {
		      TmpRe = TmpRealHamiltonian[IndexZ];
		      TmpIm = TmpImaginaryHamiltonian[IndexZ];

		      vDestination.Re(Index2) += (TmpRe * VSourceRe1 - TmpIm * VSourceIm1);
		      vDestination.Im(Index2) += (TmpRe * VSourceIm1 + TmpIm * VSourceRe1);
		      vDestination.Re(Index4) += (TmpRe * VSourceRe3 - TmpIm * VSourceIm3);
		      vDestination.Im(Index4) += (TmpRe * VSourceIm3 + TmpIm * VSourceRe3);

		      TmpTotalRe1 += (TmpRe * vSource.Re(Index2) + TmpIm * vSource.Im(Index2));
		      TmpTotalIm1 += (TmpRe * vSource.Im(Index2) - TmpIm * vSource.Re(Index2));
		      TmpTotalRe3 += (TmpRe * vSource.Re(Index4) + TmpIm * vSource.Im(Index4));
		      TmpTotalIm3 += (TmpRe * vSource.Im(Index4) - TmpIm * vSource.Re(Index4));

		      ++Index2; ++Index4;
		    }
		  // case m2 = m1 && p2 = p1
		  TmpRe = TmpRealHamiltonian[IndexZ];
		  TmpIm = TmpImaginaryHamiltonian[IndexZ];

		  vDestination.Re(Index2) += (TmpRe * VSourceRe1 + TmpIm * VSourceIm1);
		  vDestination.Im(Index2) += (TmpRe * VSourceIm1 - TmpIm * VSourceRe1);

		  vDestination.Re(Index1) += (TmpRe * vSource.Re(Index2) - TmpIm * vSource.Im(Index2));
		  vDestination.Im(Index1) += (TmpRe * vSource.Im(Index2) + TmpIm * vSource.Re(Index2)); 

		  vDestination.Re(Index1) += TmpTotalRe1; vDestination.Im(Index1) += TmpTotalIm1;
		  vDestination.Re(Index3) += TmpTotalRe3; vDestination.Im(Index3) += TmpTotalIm3;

		  ++Index1; ++Index3;
		}
	    }
	  // case m2 = m1 && n2 = n1
	  Index1 = RI[m1][n1]; // Index1 = (m1 * N + n1) * H
	  TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m1][n1];
	  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m1][n1];
	  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
	    {
	      Index2 = RI[m1][n1]; // Index2 = (m1 * N + n1) * H

	      VSourceRe1 = vSource.Re(Index1); VSourceIm1 = vSource.Im(Index1);

	      TmpTotalRe1 = 0.0; TmpTotalIm1 = 0.0; TmpTotalRe5 = 0.0; TmpTotalIm5 = 0.0;
	      for (IndexZ = p1; IndexZ > 0; --IndexZ)
		{
		  TmpRe = TmpRealHamiltonian[IndexZ];
		  TmpIm = TmpImaginaryHamiltonian[IndexZ];

		  vDestination.Re(Index2) += (TmpRe * VSourceRe1 - TmpIm * VSourceIm1);
		  vDestination.Im(Index2) += (TmpRe * VSourceIm1 + TmpIm * VSourceRe1);

		  TmpTotalRe1 += (TmpRe * vSource.Re(Index2) + TmpIm * vSource.Im(Index2));
		  TmpTotalIm1 += (TmpRe * vSource.Im(Index2) - TmpIm * vSource.Re(Index2));

		  ++Index2;
		}
	      // case m2 = m1 && n2 = n1 && p2 = p1
	      TmpRe = TmpRealHamiltonian[IndexZ];
	      TmpIm = TmpImaginaryHamiltonian[IndexZ];
	      
	      vDestination.Re(Index2) += (TmpRe * VSourceRe1 + TmpIm * VSourceIm1);
	      vDestination.Im(Index2) += (TmpRe * VSourceIm1 - TmpIm * VSourceRe1);
	      // kinetic elements added
	      vDestination.Re(Index2) += (KineticElements[Index1] * VSourceRe1);
	      vDestination.Im(Index2) += (KineticElements[Index1] * VSourceIm1);

	      ++Index2;

	      vDestination.Re(Index1) += TmpTotalRe1; vDestination.Im(Index1) += TmpTotalIm1;

	      ++Index1;
	    }
	}
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& XYReflexionSymmetricPeriodic3DHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  if ((firstComponent == 0) && (nbrComponent == this->Space->GetHilbertSpaceDimension()))
    return this->LowLevelAddMultiply(vSource, vDestination);
  else
    {
      int lastComponent = firstComponent + nbrComponent;
      int m1, m2, n1, n2, p1;
      int IndexZ;
      double* TmpRealHamiltonian;
      double* TmpImaginaryHamiltonian;
      double TmpRe = 0.0; double TmpIm = 0.0;
      
      int Index1 = firstComponent; int Index2 = 0;
      int ReducedIndex1 = firstComponent / this->NbrStateZ;
      int p1Begin = firstComponent - ReducedIndex1 * this->NbrStateZ;
      int m1Begin = ReducedIndex1 / this->NbrStateY;
      int n1Begin = ReducedIndex1 - m1Begin * this->NbrStateY;

      ReducedIndex1 = (lastComponent - 1) / this->NbrStateZ;
      int p1Limit = (lastComponent - 1) - ReducedIndex1 * this->NbrStateZ;
      int m1Limit = ReducedIndex1 / this->NbrStateY;
      int n1Limit = ReducedIndex1 - m1Limit * this->NbrStateY;

      int** TotalIndex = new int* [this->NbrStateX]; 
      for (int m = 0; m < this->NbrStateX; ++m)
        {
          TotalIndex[m] = new int [this->NbrStateY];
          for (int n = 0; n < this->NbrStateY; ++n)
            TotalIndex[m][n] = (m * this->NbrStateY + n) * this->NbrStateZ;
        }
      int* TmpTotalIndex1; int* TmpTotalIndex2;

      int LimitZ = 0; 
      for (int Index = firstComponent; Index < lastComponent; ++Index)
	{
	  vDestination.Re(Index) += vSource.Re(Index) * this->KineticElements[Index];
	  vDestination.Im(Index) += vSource.Im(Index) * this->KineticElements[Index];
	}
      // p1 : p1Begin -> NbrStateZ
      m1 = m1Begin; n1 = n1Begin;
      TmpTotalIndex1 = TotalIndex[m1];
      for (m2 = 0; m2 < m1; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];
	  for (n2 = 0; n2 < n1; ++n2)
	    {
	      Index1 = firstComponent;
	      TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m2][n2];
	      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m2][n2];
	      for (p1 = p1Begin; p1 < this->NbrStateZ; ++p1)
		{
		  Index2 = TmpTotalIndex2[n2];
		  TmpRe = 0.0; TmpIm = 0.0;
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  LimitZ = this->NbrStateZ - p1;
		  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;		  
		    }
		  vDestination.Re(Index1) += TmpRe;
		  vDestination.Im(Index1) += TmpIm;
		  ++Index1;
		}
	    }
	  for (n2 = n1; n2 < this->NbrStateY; ++n2)
	    {
	      Index1 = firstComponent;
	      TmpRealHamiltonian = this->RealHamiltonian[m1][n2][m2][n1];
	      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n2][m2][n1];
	      for (p1 = p1Begin; p1 < this->NbrStateZ; ++p1)
		{
		  Index2 = TmpTotalIndex2[n2];
		  TmpRe = 0.0; TmpIm = 0.0;
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  LimitZ = this->NbrStateZ - p1;
		  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;		  
		    }
		  vDestination.Re(Index1) += TmpRe;
		  vDestination.Im(Index1) += TmpIm;
		  ++Index1;
		}	      
	    }
	}
      for (m2 = m1; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];
	  for (n2 = 0; n2 < n1; ++n2)
	    {
	      Index1 = firstComponent;
	      TmpRealHamiltonian = this->RealHamiltonian[m2][n1][m1][n2];
	      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n1][m1][n2];
	      for (p1 = p1Begin; p1 < this->NbrStateZ; ++p1)
		{
		  Index2 = TmpTotalIndex2[n2];
		  TmpRe = 0.0; TmpIm = 0.0;
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  LimitZ = this->NbrStateZ - p1;
		  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;		  
		    }
		  vDestination.Re(Index1) += TmpRe;
		  vDestination.Im(Index1) += TmpIm;
		  ++Index1;
		}	      
	    }
	  for (n2 = n1; n2 < this->NbrStateY; ++n2)
	    {
	      Index1 = firstComponent;
	      TmpRealHamiltonian = this->RealHamiltonian[m2][n2][m1][n1];
	      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n2][m1][n1];
	      for (p1 = p1Begin; p1 < this->NbrStateZ; ++p1)
		{
		  Index2 = TmpTotalIndex2[n2];
		  TmpRe = 0.0; TmpIm = 0.0;
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  LimitZ = this->NbrStateZ - p1;
		  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;		  
		    }
		  vDestination.Re(Index1) += TmpRe;
		  vDestination.Im(Index1) += TmpIm;
		  ++Index1;
		}	      
	    }
	}
      
      // n1 : n1Begin -> NbrStateY
      m1 = m1Begin;
      TmpTotalIndex1 = TotalIndex[m1];
      for (n1 = n1Begin + 1; n1 < this->NbrStateY; ++n1)
	{
	  for (m2 = 0; m2 < m1; ++m2)
	    {
	      TmpTotalIndex2 = TotalIndex[m2];
	      for (n2 = 0; n2 < n1; ++n2)
		{
		  Index1 = TmpTotalIndex1[n1];
		  TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m2][n2];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m2][n2];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      Index2 = TmpTotalIndex2[n2];
		      TmpRe = 0.0; TmpIm = 0.0;
		      for (IndexZ = p1; IndexZ > 0; --IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;
			}
		      LimitZ = this->NbrStateZ - p1;
		      for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;		  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		    }	      
		}
	      for (n2 = n1; n2 < this->NbrStateY; ++n2)
		{
		  Index1 = TmpTotalIndex1[n1];
		  TmpRealHamiltonian = this->RealHamiltonian[m1][n2][m2][n1];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n2][m2][n1];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      Index2 = TmpTotalIndex2[n2];
		      TmpRe = 0.0; TmpIm = 0.0;
		      for (IndexZ = p1; IndexZ > 0; --IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;
			}
		      LimitZ = this->NbrStateZ - p1;
		      for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;		  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		  }	      
		}
	    }
	  for (m2 = m1; m2 < this->NbrStateX; ++m2)
	    {
	      TmpTotalIndex2 = TotalIndex[m2];
	      for (n2 = 0; n2 < n1; ++n2)
		{
		  Index1 = TmpTotalIndex1[n1];
		  TmpRealHamiltonian = this->RealHamiltonian[m2][n1][m1][n2];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n1][m1][n2];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      Index2 = TmpTotalIndex2[n2];
		      TmpRe = 0.0; TmpIm = 0.0;
		      for (IndexZ = p1; IndexZ > 0; --IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;
			}
		      LimitZ = this->NbrStateZ - p1;
		      for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;		  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		  }	      
		}
	      for (n2 = n1; n2 < this->NbrStateY; ++n2)
		{
		  Index1 = TmpTotalIndex1[n1];
		  TmpRealHamiltonian = this->RealHamiltonian[m2][n2][m1][n1];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n2][m1][n1];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      Index2 = TmpTotalIndex2[n2];
		      TmpRe = 0.0; TmpIm = 0.0;
		      for (IndexZ = p1; IndexZ > 0; --IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;
			}
		      LimitZ = this->NbrStateZ - p1;
		      for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;		  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		  }	      
		}
	    }
	}         
      // m1 : m1Begin + 1 -> m1Limit - 1
      for (m1 = m1Begin + 1; m1 < m1Limit; ++m1)
	{
	  TmpTotalIndex1 = TotalIndex[m1];
	  for (n1 = 0; n1 < this->NbrStateY; ++n1)
	    {
	      for (m2 = 0; m2 < m1; ++m2)
		{
		  TmpTotalIndex2 = TotalIndex[m2];
		  for (n2 = 0; n2 < n1; ++n2)
		    {
		      Index1 = TmpTotalIndex1[n1];
		      TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m2][n2];
		      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m2][n2];
		      for (p1 = 0; p1 < this->NbrStateZ; ++p1)
			{
			  Index2 = TmpTotalIndex2[n2];
			  TmpRe = 0.0; TmpIm = 0.0;
			  for (IndexZ = p1; IndexZ > 0; --IndexZ)
			    {
			      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			      ++Index2;
			    }
			  LimitZ = this->NbrStateZ - p1;
			  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			    {
			      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			      ++Index2;		  
			    }
			  vDestination.Re(Index1) += TmpRe;
			  vDestination.Im(Index1) += TmpIm;
			  ++Index1;
			}	      
		    }
		  for (n2 = n1; n2 < this->NbrStateY; ++n2)
		    {
		      Index1 = TmpTotalIndex1[n1];
		      TmpRealHamiltonian = this->RealHamiltonian[m1][n2][m2][n1];
		      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n2][m2][n1];
		      for (p1 = 0; p1 < this->NbrStateZ; ++p1)
			{
			  Index2 = TmpTotalIndex2[n2];
			  TmpRe = 0.0; TmpIm = 0.0;
			  for (IndexZ = p1; IndexZ > 0; --IndexZ)
			    {
			      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			      ++Index2;
			    }
			  LimitZ = this->NbrStateZ - p1;
			  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			    {
			      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			      ++Index2;		  
			    }
			  vDestination.Re(Index1) += TmpRe;
			  vDestination.Im(Index1) += TmpIm;
			  ++Index1;
			}	      
		    }
		}
	      for (m2 = m1; m2 < this->NbrStateX; ++m2)
		{
		  TmpTotalIndex2 = TotalIndex[m2];
		  for (n2 = 0; n2 < n1; ++n2)
		    {
		      Index1 = TmpTotalIndex1[n1];
		      TmpRealHamiltonian = this->RealHamiltonian[m2][n1][m1][n2];
		      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n1][m1][n2];
		      for (p1 = 0; p1 < this->NbrStateZ; ++p1)
			{
			  Index2 = TmpTotalIndex2[n2];
			  TmpRe = 0.0; TmpIm = 0.0;
			  for (IndexZ = p1; IndexZ > 0; --IndexZ)
			    {
			      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			      ++Index2;
			    }
			  LimitZ = this->NbrStateZ - p1;
			  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			    {
			      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			      ++Index2;		  
			    }
			  vDestination.Re(Index1) += TmpRe;
			  vDestination.Im(Index1) += TmpIm;
			  ++Index1;
			}	      
		    }
		  for (n2 = n1; n2 < this->NbrStateY; ++n2)
		    {
		      Index1 = TmpTotalIndex1[n1];
		      TmpRealHamiltonian = this->RealHamiltonian[m2][n2][m1][n1];
		      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n2][m1][n1];
		      for (p1 = 0; p1 < this->NbrStateZ; ++p1)
			{
			  Index2 = TmpTotalIndex2[n2];
			  TmpRe = 0.0; TmpIm = 0.0;
			  for (IndexZ = p1; IndexZ > 0; --IndexZ)
			    {
			      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			      ++Index2;
			    }
			  LimitZ = this->NbrStateZ - p1;
			  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			    {
			      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			      ++Index2;		  
			    }
			  vDestination.Re(Index1) += TmpRe;
			  vDestination.Im(Index1) += TmpIm;
			  ++Index1;
			}	      
		    }
		}
	    }     
	}
      // n1 : 0 -> n1Limit - 1
      m1 = m1Limit;
      TmpTotalIndex1 = TotalIndex[m1];
      for (n1 = 0; n1 < n1Limit; ++n1)
	{
	  for (m2 = 0; m2 < m1; ++m2)
	    {
	      TmpTotalIndex2 = TotalIndex[m2];
	      for (n2 = 0; n2 < n1; ++n2)
		{
		  Index1 = TmpTotalIndex1[n1];
		  TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m2][n2];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m2][n2];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      Index2 = TmpTotalIndex2[n2];
		      TmpRe = 0.0; TmpIm = 0.0;
		      for (IndexZ = p1; IndexZ > 0; --IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;
			}
		      LimitZ = this->NbrStateZ - p1;
		      for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;		  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		  }	      
		}
	      for (n2 = n1; n2 < this->NbrStateY; ++n2)
		{
		  Index1 = TmpTotalIndex1[n1];
		  TmpRealHamiltonian = this->RealHamiltonian[m1][n2][m2][n1];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n2][m2][n1];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      Index2 = TmpTotalIndex2[n2];
		      TmpRe = 0.0; TmpIm = 0.0;
		      for (IndexZ = p1; IndexZ > 0; --IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;
			}
		      LimitZ = this->NbrStateZ - p1;
		      for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;		  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		    }	      
		}
	    }
	  for (m2 = m1; m2 < this->NbrStateX; ++m2)
	    {
	      TmpTotalIndex2 = TotalIndex[m2];
	      for (n2 = 0; n2 < n1; ++n2)
		{
		  Index1 = TmpTotalIndex1[n1];
		  TmpRealHamiltonian = this->RealHamiltonian[m2][n1][m1][n2];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n1][m1][n2];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      Index2 = TmpTotalIndex2[n2];
		      TmpRe = 0.0; TmpIm = 0.0;
		      for (IndexZ = p1; IndexZ > 0; --IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;
			}
		      LimitZ = this->NbrStateZ - p1;
		      for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;		  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		  }	      
		}
	      for (n2 = n1; n2 < this->NbrStateY; ++n2)
		{
		  Index1 = TmpTotalIndex1[n1];
		  TmpRealHamiltonian = this->RealHamiltonian[m2][n2][m1][n1];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n2][m1][n1];
		  for (p1 = 0; p1 < this->NbrStateZ; ++p1)
		    {
		      Index2 = TmpTotalIndex2[n2];
		      TmpRe = 0.0; TmpIm = 0.0;
		      for (IndexZ = p1; IndexZ > 0; --IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;
			}
		      LimitZ = this->NbrStateZ - p1;
		      for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
			{
			  TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
			  TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
			  ++Index2;		  
			}
		      vDestination.Re(Index1) += TmpRe;
		      vDestination.Im(Index1) += TmpIm;
		      ++Index1;
		    }	      
		}
	    }
	}
      // p1 : 0 ->p1Limit
      m1 = m1Limit; n1 = n1Limit;
      TmpTotalIndex1 = TotalIndex[m1];
      for (m2 = 0; m2 < m1; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];
	  for (n2 = 0; n2 < n1; ++n2)
	    {
	      Index1 = TmpTotalIndex1[n1];;
	      TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m2][n2];
	      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m2][n2];
	      for (p1 = 0; p1 <= p1Limit; ++p1)
		{
		  Index2 = TmpTotalIndex2[n2];
		  TmpRe = 0.0; TmpIm = 0.0;
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  LimitZ = this->NbrStateZ - p1;
		  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;		  
		    }
		  vDestination.Re(Index1) += TmpRe;
		  vDestination.Im(Index1) += TmpIm;
		  ++Index1;
		}	      
	    }
	  for (n2 = n1; n2 < this->NbrStateY; ++n2)
	    {
	      Index1 = TmpTotalIndex1[n1];;
	      TmpRealHamiltonian = this->RealHamiltonian[m1][n2][m2][n1];
	      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n2][m2][n1];
	      for (p1 = 0; p1 <= p1Limit; ++p1)
		{
		  Index2 = TmpTotalIndex2[n2];
		  TmpRe = 0.0; TmpIm = 0.0;
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  LimitZ = this->NbrStateZ - p1;
		  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;		  
		    }
		  vDestination.Re(Index1) += TmpRe;
		  vDestination.Im(Index1) += TmpIm;
		  ++Index1;
		}	      
	    }
	}
      for (m2 = m1; m2 < this->NbrStateX; ++m2)
	{
	  TmpTotalIndex2 = TotalIndex[m2];
	  for (n2 = 0; n2 < n1; ++n2)
	    {
	      Index1 = TmpTotalIndex1[n1];;
	      TmpRealHamiltonian = this->RealHamiltonian[m2][n1][m1][n2];
	      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n1][m1][n2];
	      for (p1 = 0; p1 <= p1Limit; ++p1)
		{
		  Index2 = TmpTotalIndex2[n2];
		  TmpRe = 0.0; TmpIm = 0.0;
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  LimitZ = this->NbrStateZ - p1;
		  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;		  
		    }
		  vDestination.Re(Index1) += TmpRe;
		  vDestination.Im(Index1) += TmpIm;
		  ++Index1;
		}	      
	    }
	  for (n2 = n1; n2 < this->NbrStateY; ++n2)
	    {
	      Index1 = TmpTotalIndex1[n1];;
	      TmpRealHamiltonian = this->RealHamiltonian[m2][n2][m1][n1];
	      TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m2][n2][m1][n1];
	      for (p1 = 0; p1 <= p1Limit; ++p1)
		{
		  Index2 = TmpTotalIndex2[n2];
		  TmpRe = 0.0; TmpIm = 0.0;
		  for (IndexZ = p1; IndexZ > 0; --IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;
		    }
		  LimitZ = this->NbrStateZ - p1;
		  for (IndexZ = 0; IndexZ < LimitZ; ++IndexZ)
		    {
		      TmpRe += (TmpRealHamiltonian[IndexZ] * vSource.Re(Index2) - TmpImaginaryHamiltonian[IndexZ] * vSource.Im(Index2));
		      TmpIm += (TmpRealHamiltonian[IndexZ] * vSource.Im(Index2) + TmpImaginaryHamiltonian[IndexZ] * vSource.Re(Index2));
		      ++Index2;		  
		    }
		  vDestination.Re(Index1) += TmpRe;
		  vDestination.Im(Index1) += TmpIm;
		  ++Index1;
		}	      
	    }
	}
      delete[] TotalIndex;  
      return vDestination;
	}
}

// evaluate all interaction factors
// 
// pairX = whether basis is pair in X direction, if not impair
// pairY = whether basis is pair in Y direction, if not impair
// waveVectorZ = wave vector of Bloch function in Z direction

void XYReflexionSymmetricPeriodic3DHamiltonian::EvaluateInteractionFactors(bool pairX, bool pairY, double waveVectorZ)
{
  if (pairX)    
    this->WaveFunctionOverlapX = this->EvaluateCosinusWaveFunctionOverlap(this->XSize, this->NbrCellX, this->NbrStateX);
  else
    this->WaveFunctionOverlapX = this->EvaluateSinusWaveFunctionOverlap(this->XSize, this->NbrCellX, this->NbrStateX);
  if (pairY)
    this->WaveFunctionOverlapY = this->EvaluateCosinusWaveFunctionOverlap(this->YSize, this->NbrCellY, this->NbrStateY);
  else
    this->WaveFunctionOverlapY = this->EvaluateSinusWaveFunctionOverlap(this->YSize, this->NbrCellY, this->NbrStateY);

  if (!this->EvaluatePlaneWaveFunctionOverlap(this->NbrCellZ, this->NbrStateZ, this->RealWaveFunctionOverlapZ, this->ImaginaryWaveFunctionOverlapZ))
    {
      cout << "Error in evaluation of function overlap in Z direction. Stop!" << endl;
      exit(0);
    }

  double InvXFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Mux * this->XSize * this->XSize);
  double InvYFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muy * this->YSize * this->YSize);
  double InvZFactor = PERIODIC_HAMILTONIAN_FACTOR / (this->Muz * this->ZSize * this->ZSize);
  double ShiftSquareKz = BLOCH_FACTOR * waveVectorZ * waveVectorZ / (2.0 * this->Muz);
  double ShiftKz = BLOCH_FACTOR * waveVectorZ * 2.0 * M_PI/ (this->Muz * this->ZSize); 
  this->KineticElements = new double[this->Space->GetHilbertSpaceDimension()];
  // this->NbrStateX * this->NbrStateY * this->NbrStateZ this->Space->GetHilbertSpaceDimension()
  double FactorX = 0.0, FactorY = 0.0;
  int TotalIndex = 0;
  for (int i = 0; i < this->NbrStateX; ++i)
    {
      FactorX = double((i + this->LowerImpulsionX) * (i + this->LowerImpulsionX)) * InvXFactor;
      for (int j = 0; j < this->NbrStateY; ++j)
	{	  
	  FactorY = double((j + this->LowerImpulsionY) * (j + this->LowerImpulsionY)) * InvYFactor + FactorX;
	  for (int k = 0; k < this->NbrStateZ; ++k)
	    {
	      this->KineticElements[TotalIndex] = FactorY + double((k + this->LowerImpulsionZ) * (k + this->LowerImpulsionZ)) * InvZFactor;
	      this->KineticElements[TotalIndex] += (ShiftSquareKz + ShiftKz * double(k + this->LowerImpulsionZ));
	      ++TotalIndex;
	    }
	}
    }

  double* TmpRealHamiltonian;
  double* TmpImaginaryHamiltonian;
  double* TmpWaveFunctionOverlapX;
  double* TmpWaveFunctionOverlapY;
  double* TmpRealWaveFunctionOverlapZ;
  double* TmpImaginaryWaveFunctionOverlapZ;
  int m1 = 0, m2 = 0, n1 = 0, n2 = 0, p = 0;
  double TmpRe = 0.0, TmpIm = 0.0, TmpXY = 0.0, TmpX = 0.0, TmpY = 0.0;
  int CellX = 0, CellY = 0, CellZ = 0;

  this->RealHamiltonian = new double**** [this->NbrStateX];
  this->ImaginaryHamiltonian = new double**** [this->NbrStateX];
  double* TmpInter = new double [this->NbrCellZ];
  double* TmpInter2; double* TmpInteraction;
  for (m1 = 0; m1 < this->NbrStateX; ++m1)
    {	  
      this->RealHamiltonian[m1] = new double*** [this->NbrStateY];	
      this->ImaginaryHamiltonian[m1] = new double*** [this->NbrStateY];	
      double*** Inter2 = new double** [m1 + 1];
      for (m2 = 0; m2 <= m1; ++m2)
	{     
	  Inter2[m2] = new double* [this->NbrCellZ];
	  TmpWaveFunctionOverlapX = this->WaveFunctionOverlapX[m1][m2]; 
	  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
	    {
	      Inter2[m2][CellZ] = new double [this->NbrCellY];
	      for (CellY = 0; CellY < this->NbrCellY; ++CellY)
		{
		  TmpX = 0.0;
		  TmpInteraction = this->InteractionFactors[CellZ][CellY];
		  for (CellX = 0; CellX < this->NbrCellX; ++CellX)			     
		    TmpX += TmpInteraction[CellX] * TmpWaveFunctionOverlapX[CellX];
		  Inter2[m2][CellZ][CellY] = TmpX;
		}
	    }
	}
      
      for (n1 = 0; n1 < this->NbrStateY; ++n1)
	{
	  this->RealHamiltonian[m1][n1] = new double** [m1 + 1];
	  this->ImaginaryHamiltonian[m1][n1] = new double** [m1 + 1];	      	      
	  for (m2 = 0; m2 <= m1; ++m2)
	    {	      
	      this->RealHamiltonian[m1][n1][m2] = new double* [n1 + 1];	
 	      this->ImaginaryHamiltonian[m1][n1][m2] = new double* [n1 + 1];	      
	      for (n2 = 0; n2 <= n1; ++n2)
		{
		  TmpWaveFunctionOverlapY = this->WaveFunctionOverlapY[n1][n2];	  
		  for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)
		    {
		      TmpXY = 0.0;
		      TmpInter2 = Inter2[m2][CellZ];
		      for (CellY = 0; CellY < this->NbrCellY; ++CellY)			
			TmpXY += (TmpInter2[CellY] * TmpWaveFunctionOverlapY[CellY] );     		      
		      TmpInter[CellZ] = TmpXY;
		    }

		  this->RealHamiltonian[m1][n1][m2][n2] = new double [this->NbrStateZ];
		  this->ImaginaryHamiltonian[m1][n1][m2][n2] = new double [this->NbrStateZ];
		  TmpRealHamiltonian = this->RealHamiltonian[m1][n1][m2][n2];
		  TmpImaginaryHamiltonian = this->ImaginaryHamiltonian[m1][n1][m2][n2];	
		  
		  for (p = 0; p < this->NbrStateZ; ++p)
		    {
		      TmpRealWaveFunctionOverlapZ = RealWaveFunctionOverlapZ[p];
		      TmpImaginaryWaveFunctionOverlapZ = this->ImaginaryWaveFunctionOverlapZ[p];
		      
		      TmpRe = 0.0; TmpIm = 0.0;
		      for (CellZ = 0; CellZ < this->NbrCellZ; ++CellZ)		    
			{			  			  
			  TmpRe += (TmpInter[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ]);
			  TmpIm += (TmpInter[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ]);				    
			}
		      TmpRealHamiltonian[p] = TmpRe;
		      TmpImaginaryHamiltonian[p] = TmpIm;			      
		    }
		}		  	      
	    }
	}
      delete[] Inter2;
    }  
  delete[] TmpInter;
}

// evaluate sinus wave function overlaps on a cell in a given direction
//
// size = system length in the choosen direction
// nbrStep = number of subdivision in the choosen direction
// nbrState = number of state in the choosen direction
// memory = reference on current memory usage (will be increment with memory used to store evaluated overlap)
// return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)

double*** XYReflexionSymmetricPeriodic3DHamiltonian::EvaluateSinusWaveFunctionOverlap(double size, int nbrStep, int nbrState)
{
  double*** TmpArray = new double** [nbrState];
  double StepInc = 1.0 / ((double) nbrStep);
  double Tmp;
  double Diff;
  for (int i = 0; i < nbrState; ++i)
    {
      TmpArray[i] = new double* [i + 1];
      for (int j = 0; j < i; ++j)
	{
	  TmpArray[i][j] = new double [nbrStep];
	  for (int k = 0; k < nbrStep; ++k)
	    {
	      Diff = (double) 2 * (i - j);
	      Tmp = M_PI * Diff * StepInc;	      
	      TmpArray[i][j][k] = M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / Diff;
	      Diff = (double) 2 * (i + j + 2);
	      Tmp = M_PI * Diff * StepInc;	      
	      TmpArray[i][j][k] -= M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / Diff;
	    }
	}
      TmpArray[i][i] = new double [nbrStep];
      for (int k = 0; k < nbrStep; ++k)
	{
	  Tmp = M_PI * (double) (4 * i + 4) * StepInc;	      
	  TmpArray[i][i][k] = StepInc - M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / ((double) (4 * i + 4));
	}     
    }
  return TmpArray;
}

// evaluate cosinus wave function overlaps on a cell in a given direction
//
// size = system length in the choosen direction
// nbrStep = number of subdivision in the choosen direction
// nbrState = number of state in the choosen direction
// memory = reference on current memory usage (will be increment with memory used to store evaluated overlap)
// return value = tridimensionnal array containg all matrix elements for all cells (first two indices using symmetric storage)

double*** XYReflexionSymmetricPeriodic3DHamiltonian::EvaluateCosinusWaveFunctionOverlap(double size, int nbrStep, int nbrState)
{
  double*** TmpArray = new double** [nbrState];
  double StepInc = 1.0 / ((double) nbrStep);
  double Tmp;
  double Diff;
  
  TmpArray[0] = new double* [1];
  TmpArray[0][0] = new double [nbrStep];
  for (int k = 0; k < nbrStep; ++k)
    TmpArray[0][0][k] = StepInc;

  for (int i = 1; i < nbrState; ++i)
    {
      TmpArray[i] = new double* [i + 1];

      TmpArray[i][0] = new double [nbrStep];
      for (int k = 0; k < nbrStep; ++k)
	{
	  Diff = (double) 2 * i;
	  Tmp = M_PI * Diff * StepInc;	  
	  TmpArray[i][0][k] = M_SQRT2 * M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / Diff;	  
	}     
      for (int j = 1; j < i; ++j)
	{
	  TmpArray[i][j] = new double [nbrStep];
	  for (int k = 0; k < nbrStep; ++k)
	    {
	      Diff = (double) 2 * (i - j);
	      Tmp = M_PI * Diff * StepInc;	      
	      TmpArray[i][j][k] = M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / Diff;
	      Diff = (double) 2 * (i + j);
	      Tmp = M_PI * Diff * StepInc;
	      TmpArray[i][j][k] += M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / Diff;
	    }
	}
      TmpArray[i][i] = new double [nbrStep];
      for (int k = 0; k < nbrStep; ++k)
	{
	  Tmp = M_PI * (double) (4 * i) * StepInc;	      
	  TmpArray[i][i][k] = StepInc + M_1_PI * (sin (Tmp * ((double) (k + 1))) - sin (Tmp * ((double) (k)))) / ((double) (4 * i));
	}     
    }
  return TmpArray;
}

// evaluate the plane wave function overlap
//
// nbrStep = number of steps in the given direction
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool XYReflexionSymmetricPeriodic3DHamiltonian::EvaluatePlaneWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray)
{
  double Diff = 0.0;
  double Tmp = 0.0;
  double Tmp1 = 1.0 / double (nbrStep); 
  realArray = new double* [nbrState];
  imaginaryArray = new double* [nbrState];    
  for (int delta = 0; delta < nbrState; ++delta)
    {
      realArray[delta] = new double [nbrStep];
      imaginaryArray[delta] = new double [nbrStep];
      if (delta != 0)
	{
	  Diff = 2.0 * M_PI * double(delta);
	  Tmp = Diff / nbrStep;	
	  Diff = 1.0 / Diff;	
	  for (int i = 0; i < nbrStep; ++i)
	    {
	      realArray[delta][i] = Diff * (sin(Tmp * (i + 1)) - sin(Tmp * i));
	      imaginaryArray[delta][i] = Diff * (-cos(Tmp * (i + 1)) + cos(Tmp * i));
	    }
	}
      else
	for (int i = 0; i < nbrStep; ++i)
	  {
	    realArray[delta][i] = Tmp1;
	    imaginaryArray[delta][i] = 0.0;
	  }	
    }
  return true;
}

// determine the maximal value of partial diagonal array
//
// return = the wanted value

double XYReflexionSymmetricPeriodic3DHamiltonian::MaxPartialDiagonalElement()
{
  double tmp = this->KineticElements[0];
  for (int i = 1; i < this->Space->GetHilbertSpaceDimension(); ++i)
    if (tmp < this->KineticElements[i])
      tmp = this->KineticElements[i];
  return tmp;
}
