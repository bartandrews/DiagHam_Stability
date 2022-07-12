////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                    class author: Gunnar Möller                             //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                to particles on a torus with magnetic translations          //
//                                                                            //
//                        last modification : 28/11/2007                      //
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
#include "Hamiltonian/AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>
#include <fstream>

// testing only:
// #include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslations.h"

using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::hex;
using std::dec;
using std::endl;
using std::ostream;


// destructor
//

AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::~AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian()
{
  
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  // delete existing interaction factors and recalculate
  delete[] InteractionFactorsUpUp;
  delete[] InteractionFactorsDownDown;
  delete[] InteractionFactorsUpDown;  
  delete[] OneBodyInteractionFactorsUpUp;
  delete[] OneBodyInteractionFactorsDownDown;    
  this->Particles = (ParticleOnTorusWithSpinAndMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::ShiftHamiltonian (double shift)
{
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  return Complex();
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  /*
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Cosinus;
  double Sinus;
  int NbrTranslation;
  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  Complex TmpZ;
  Complex Z (0.0, 0.0);
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
  for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
    {
      m1 = this->M1Value[j];
      m2 = this->M2Value[j];
      m3 = this->M3Value[j];
      m4 = this->M4Value[j];
      TmpInteraction = this->InteractionFactors[j];
      for (int i = 0; i < Dim; ++i)
	{
	  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslation);
	  if (Index < Dim)
	    {
	      Coefficient *= TmpInteraction;
	      Cosinus = Coefficient * this->CosinusTable[NbrTranslation];
	      Sinus = Coefficient * this->SinusTable[NbrTranslation];
	      TmpZ.Re = ((V2.Re(i) * Cosinus) - (V2.Im(i) * Sinus));
	      TmpZ.Im += ((V2.Re(i) * Sinus) + (V2.Im(i) * Cosinus));
	      Z += Conj(V1[Index]) * TmpZ;
	    }
	}
    }
  m1 = this->M1Value[ReducedNbrInteractionFactors];
  m2 = this->M2Value[ReducedNbrInteractionFactors];
  m3 = this->M3Value[ReducedNbrInteractionFactors];
  m4 = this->M4Value[ReducedNbrInteractionFactors];
  TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
  for (int i = 0; i < Dim; ++i)
    {
      Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslation);
      if (Index < Dim)
	{
	  Coefficient *= TmpInteraction;
	  Cosinus = Coefficient * this->CosinusTable[NbrTranslation];
	  Sinus = Coefficient * this->SinusTable[NbrTranslation];
	  TmpZ.Re = ((V2.Re(i) * Cosinus) - (V2.Im(i) * Sinus));
	  TmpZ.Im += ((V2.Re(i) * Sinus) + (V2.Im(i) * Cosinus));
	  Z += Conj(V1[Index]) * TmpZ;
	}
      Z += this->EnergyShift * (Conj(V1[i]) * V2[i]);
    }
  return Complex(Z);
  */
  return Complex();
}
  
// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
											       int firstComponent, int nbrComponent) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
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

RealVector& AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
												  int firstComponent, int nbrComponent)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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

ComplexVector& AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
ComplexVector& AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
												 int firstComponent, int nbrComponent)
{
  unsigned L16Mask = (1u<<16)-1;
  unsigned H16Mask = (~0u)^L16Mask;
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Cosinus;
  double Sinus;
  int NbrTranslation;
  if (this->FastMultiplicationFlag == false)
    {
      double Coefficient2;
      double Coefficient3;
      int Index;
      int m1, m2, m3, m4;      
      int SumIndices;
      int TmpNbrM34Values;
      unsigned* TmpM34Values;      
      int ReducedNbrInteractionFactors;
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  ReducedNbrInteractionFactors = 0;
	  for (int m12 = 0; m12 < this->NbrM12IntraIndices; ++m12)
	    {
	      m1 = this->M12IntraValue[m12] & L16Mask;
	      m2 = (this->M12IntraValue[m12] & H16Mask)>>16;
	      Coefficient = Particles->AuAu(i, m1, m2);
	      if (Coefficient != 0.0)
		{
		  SumIndices = m1 + m2;
		  TmpNbrM34Values = this->NbrM34IntraValues[m12];
		  TmpM34Values = this->M34IntraValues[m12];
		  for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
		    {
		      m3 = (TmpM34Values[m34]) & L16Mask;
		      m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
		      Index = Particles->AduAdu(m3, m4, Coefficient2, NbrTranslation);
		      if (Index < Dim)
			{
			  Coefficient3 = Coefficient * Coefficient2 * this->InteractionFactorsUpUp[ReducedNbrInteractionFactors];
			  // int Index2, NbrTranslation2;
// 			  double Coefficientprime;
// 			  Index2 = Particles->AduAduAuAu (i, m3, m4, m1, m2, Coefficientprime, NbrTranslation2);
// 			  if (Index != Index2)
// 			    cout << "Problem with Index in HilbertSpace from AbstractHamiltonian" << endl;
// 			  if (Coefficient*Coefficient2 != Coefficientprime)
// 			    {
// 			      cout << "Problem with Coefficient in HilbertSpace from AbstractHamiltonian: "<<
// 				Coefficient*Coefficient2 << " vs " << Coefficientprime << " for i=" << i
// 				   << " m1= " << m1 << " m2= " << m2 << " m3= " << m3 << " m4= " << m4<< endl;
// // 			      cout << "Coefficient: " << Coefficient << " Coefficient2: " << Coefficient2 << endl;
// // 			      ((FermionOnTorusWithSpinAndMagneticTranslations*)Particles)->AduAduAuAuV (i, m3, m4, m1, m2, Coefficientprime, NbrTranslation2);
// // 			      ((FermionOnTorusWithSpinAndMagneticTranslations*)Particles)->AuAuV(i, m1, m2);
// // 			      ((FermionOnTorusWithSpinAndMagneticTranslations*)Particles)->AduAduV(m3, m4, Coefficient2, NbrTranslation);
// // 			      cout << "Again, Coefficient: " << Coefficient << " Coefficient2: " << Coefficient2 << endl;
// 			    }
// 			  if (NbrTranslation != NbrTranslation2)
// 			    cout << "Problem with NbrTranslation in HilbertSpace from AbstractHamiltonian" << endl;
	       		  Cosinus = Coefficient3 * this->CosinusTable[NbrTranslation];			  
			  Sinus = Coefficient3 * this->SinusTable[NbrTranslation];
			  vDestination.Re(Index) += ((vSource.Re(i) * Cosinus) - (vSource.Im(i) * Sinus));
			  vDestination.Im(Index) += ((vSource.Re(i) * Sinus) + (vSource.Im(i) * Cosinus));
			}
		      ++ReducedNbrInteractionFactors;
		    }
		}
	      else
		ReducedNbrInteractionFactors += this->NbrM34IntraValues[m12];
	    }
	  ReducedNbrInteractionFactors = 0;
	  for (int m12 = 0; m12 < this->NbrM12IntraIndices; ++m12)
	    {
	      m1 = this->M12IntraValue[m12] & L16Mask;
	      m2 = (this->M12IntraValue[m12] & H16Mask)>>16;
	      Coefficient = Particles->AdAd(i, m1, m2);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = m1 + m2;
		  TmpNbrM34Values = this->NbrM34IntraValues[m12];
		  TmpM34Values = this->M34IntraValues[m12];
		  for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
		    {
		      m3 = (TmpM34Values[m34]) & L16Mask;
		      m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
		      Index = Particles->AddAdd(m3, m4, Coefficient2, NbrTranslation);
		      if (Index < Dim)
			{
			  Coefficient3 = Coefficient * Coefficient2 * this->InteractionFactorsDownDown[ReducedNbrInteractionFactors];
			  Cosinus = Coefficient3 * this->CosinusTable[NbrTranslation];
			  Sinus = Coefficient3 * this->SinusTable[NbrTranslation];
			  vDestination.Re(Index) += ((vSource.Re(i) * Cosinus) - (vSource.Im(i) * Sinus));
			  vDestination.Im(Index) += ((vSource.Re(i) * Sinus) + (vSource.Im(i) * Cosinus));
			}
		      ++ReducedNbrInteractionFactors;
		    }
		}
	      else
		ReducedNbrInteractionFactors += this->NbrM34IntraValues[m12];
	    }
	  
	  ReducedNbrInteractionFactors = 0;
	  for (int m12 = 0; m12 < this->NbrM12InterIndices; ++m12)
	    {
	      m1 = this->M12InterValue[m12] & L16Mask;
	      m2 = (this->M12InterValue[m12] & H16Mask)>>16;
	      Coefficient = Particles->AuAd(i, m1, m2);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = m1 + m2;
		  TmpNbrM34Values = this->NbrM34InterValues[m12];
		  TmpM34Values = this->M34InterValues[m12];
		  for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
		    {
		      m3 = (TmpM34Values[m34]) & L16Mask;
		      m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
		      Index = Particles->AduAdd(m3, m4, Coefficient2, NbrTranslation);
		      if (Index < Dim)
			{
			  Coefficient3 = Coefficient * Coefficient2 * this->InteractionFactorsUpDown[ReducedNbrInteractionFactors];
			  Cosinus = Coefficient3 * this->CosinusTable[NbrTranslation];
			  Sinus = Coefficient3 * this->SinusTable[NbrTranslation];
			  vDestination.Re(Index) += ((vSource.Re(i) * Cosinus) - (vSource.Im(i) * Sinus));
			  vDestination.Im(Index) += ((vSource.Re(i) * Sinus) + (vSource.Im(i) * Cosinus));
			}		      ++ReducedNbrInteractionFactors;
		    }
		}
	      else
		ReducedNbrInteractionFactors += this->NbrM34InterValues[m12];
	    }
	}  

      if (this->OneBodyInteractionFactorsUpUp != 0) 
	if (this->OneBodyInteractionFactorsDownDown != 0)
	  {
	    double TmpDiagonal = 0.0;
	    for (int i = firstComponent; i < LastComponent; ++i)
	      { 
		TmpDiagonal = 0.0;
		for (int j = 0; j <= this->MaxMomentum; ++j) 
		  {
		    TmpDiagonal += this->OneBodyInteractionFactorsUpUp[j] * Particles->AduAu(i, j);
		    TmpDiagonal += this->OneBodyInteractionFactorsDownDown[j] * Particles->AddAd(i, j);
		  }
		vDestination.Re(i) += (this->EnergyShift + TmpDiagonal) * vSource.Re(i);
		vDestination.Im(i) += (this->EnergyShift + TmpDiagonal) * vSource.Im(i);
	      }
	  }
	else
	  {
	    double TmpDiagonal = 0.0;
	    for (int i = firstComponent; i < LastComponent; ++i)
	      { 
		TmpDiagonal = 0.0;
		for (int j = 0; j <= this->MaxMomentum; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactorsUpUp[j] * Particles->AduAu(i, j);
		vDestination.Re(i) += (this->EnergyShift + TmpDiagonal) * vSource.Re(i);
		vDestination.Im(i) += (this->EnergyShift + TmpDiagonal) * vSource.Im(i);
	      }
	  }
      else
	if (this->OneBodyInteractionFactorsDownDown != 0)
	  {
	    double TmpDiagonal = 0.0;
	    for (int i = firstComponent; i < LastComponent; ++i)
	      { 
		TmpDiagonal = 0.0;
		for (int j = 0; j <= this->MaxMomentum; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactorsDownDown[j] * Particles->AddAd(i, j);
		vDestination.Re(i) += (this->EnergyShift + TmpDiagonal) * vSource.Re(i);
		vDestination.Im(i) += (this->EnergyShift + TmpDiagonal) * vSource.Im(i);
	      }
	  }	
	else
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      vDestination.Re(i) += (this->EnergyShift) * vSource.Re(i);
	      vDestination.Im(i) += (this->EnergyShift) * vSource.Im(i);
	    }
    }
  else // fast multiplication
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int* TmpNbrTranslationArray;
	  int j;
	  int TmpNbrInteraction;
	  double TmpRe;
	  double TmpIm;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslation[i];
	      TmpRe = vSource.Re(i);
	      TmpIm = vSource.Im(i);
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Cosinus = TmpCoefficientArray[j];
		  NbrTranslation = TmpNbrTranslationArray[j];
		  Sinus = Cosinus * this->SinusTable[NbrTranslation];
		  Cosinus *= this->CosinusTable[NbrTranslation];
		  vDestination.Re(TmpIndexArray[j]) += ((Cosinus * TmpRe) - (Sinus * TmpIm));
		  vDestination.Im(TmpIndexArray[j]) += ((Sinus * TmpRe) + (Cosinus * TmpIm));
		}
	      vDestination.Re(i) += EnergyShift * TmpRe;
	      vDestination.Im(i) += EnergyShift * TmpIm;
	    }
	}
      else
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int* TmpNbrTranslationArray;
	  int j;
	  int TmpNbrInteraction;
	  int Pos = firstComponent / this->FastMultiplicationStep; 
	  int PosMod = firstComponent % this->FastMultiplicationStep;
	  double TmpRe;
	  double TmpIm;
	  if (PosMod != 0)
	    {
	      ++Pos;
	      PosMod = this->FastMultiplicationStep - PosMod;
	    }
	  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
	      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
	      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslation[Pos];
	      TmpRe = vSource.Re(i);
	      TmpIm = vSource.Im(i);
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Cosinus = TmpCoefficientArray[j];
		  NbrTranslation = TmpNbrTranslationArray[j];
		  Sinus = Cosinus * this->SinusTable[NbrTranslation];
		  Cosinus *= this->CosinusTable[NbrTranslation];
		  vDestination.Re(TmpIndexArray[j]) += ((Cosinus * TmpRe) - (Sinus * TmpIm));
		  vDestination.Im(TmpIndexArray[j]) += ((Sinus * TmpRe) + (Cosinus * TmpIm));
		}
	      vDestination.Re(i) += EnergyShift * TmpRe;
	      vDestination.Im(i) += EnergyShift * TmpIm;
	      ++Pos;
	    }

	  double Coefficient2;
	  double Coefficient3;
	  int Index;
	  int m1, m2, m3, m4;      
	  int SumIndices;
	  int TmpNbrM34Values;
	  unsigned* TmpM34Values;      
	  int ReducedNbrInteractionFactors;
	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    ReducedNbrInteractionFactors = 0;
		    for (int m12 = 0; m12 < this->NbrM12IntraIndices; ++m12)
		      {
			m1 = this->M12IntraValue[m12] & L16Mask;
			m2 = (this->M12IntraValue[m12] & H16Mask)>>16;
			Coefficient = Particles->AuAu(i, m1, m2);
			if (Coefficient != 0.0)
			  {
			    SumIndices = m1 + m2;
			    TmpNbrM34Values = this->NbrM34IntraValues[m12];
			    TmpM34Values = this->M34IntraValues[m12];
			    for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
			      {
				m3 = (TmpM34Values[m34]) & L16Mask;
				m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
				Index = Particles->AduAdu(m3, m4, Coefficient2, NbrTranslation);
				if (Index < Dim)
				  {
				    Coefficient3 = Coefficient * Coefficient2 * this->InteractionFactorsUpUp[ReducedNbrInteractionFactors];
				    Cosinus = Coefficient3 * this->CosinusTable[NbrTranslation];
				    Sinus = Coefficient3 * this->SinusTable[NbrTranslation];
				    vDestination.Re(Index) += ((vSource.Re(i) * Cosinus) - (vSource.Im(i) * Sinus));
				    vDestination.Im(Index) += ((vSource.Re(i) * Sinus) + (vSource.Im(i) * Cosinus));
				  }
				++ReducedNbrInteractionFactors;
			      }
			  }
			else
			  ReducedNbrInteractionFactors += this->NbrM34IntraValues[m12];
		      }
		    ReducedNbrInteractionFactors = 0;
		    for (int m12 = 0; m12 < this->NbrM12IntraIndices; ++m12)
		      {
			m1 = this->M12IntraValue[m12] & L16Mask;
			m2 = (this->M12IntraValue[m12] & H16Mask)>>16;
			Coefficient = Particles->AdAd(i, m1, m2);	  
			if (Coefficient != 0.0)
			  {
			    SumIndices = m1 + m2;
			    TmpNbrM34Values = this->NbrM34IntraValues[m12];
			    TmpM34Values = this->M34IntraValues[m12];
			    for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
			      {
				m3 = (TmpM34Values[m34]) & L16Mask;
				m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
				Index = Particles->AddAdd(m3, m4, Coefficient2, NbrTranslation);
				if (Index < Dim)
				  {
				    Coefficient3 = Coefficient * Coefficient2 * this->InteractionFactorsDownDown[ReducedNbrInteractionFactors];
				    Cosinus = Coefficient3 * this->CosinusTable[NbrTranslation];
				    Sinus = Coefficient3 * this->SinusTable[NbrTranslation];
				    vDestination.Re(Index) += ((vSource.Re(i) * Cosinus) - (vSource.Im(i) * Sinus));
				    vDestination.Im(Index) += ((vSource.Re(i) * Sinus) + (vSource.Im(i) * Cosinus));
				  }
				++ReducedNbrInteractionFactors;
			      }
			  }
			else
			  ReducedNbrInteractionFactors += this->NbrM34IntraValues[m12];
		      }
		    
		    ReducedNbrInteractionFactors = 0;
		    for (int m12 = 0; m12 < this->NbrM12InterIndices; ++m12)
		      {
			m1 = this->M12InterValue[m12] & L16Mask;
			m2 = (this->M12InterValue[m12] & H16Mask)>>16;
			Coefficient = Particles->AuAd(i, m1, m2);	  
			if (Coefficient != 0.0)
			  {
			    SumIndices = m1 + m2;
			    TmpNbrM34Values = this->NbrM34InterValues[m12];
			    TmpM34Values = this->M34InterValues[m12];
			    for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
			      {
				m3 = (TmpM34Values[m34]) & L16Mask;
				m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
				Index = Particles->AduAdd(m3, m4, Coefficient2, NbrTranslation);
				if (Index < Dim)
				  {
				    Coefficient3 = Coefficient * Coefficient2 * this->InteractionFactorsUpDown[ReducedNbrInteractionFactors];
				    Cosinus = Coefficient3 * this->CosinusTable[NbrTranslation];
				    Sinus = Coefficient3 * this->SinusTable[NbrTranslation];
				    vDestination.Re(Index) += ((vSource.Re(i) * Cosinus) - (vSource.Im(i) * Sinus));
				    vDestination.Im(Index) += ((vSource.Re(i) * Sinus) + (vSource.Im(i) * Cosinus));
				  }		      ++ReducedNbrInteractionFactors;
			      }
			  }
			else
			  ReducedNbrInteractionFactors += this->NbrM34InterValues[m12];
		      }
		  }  
		
		if (this->OneBodyInteractionFactorsUpUp != 0) 
		  if (this->OneBodyInteractionFactorsDownDown != 0)
		    {
		      double TmpDiagonal = 0.0;
		      for (int i = firstComponent; i < LastComponent; ++i)
			{ 
			  TmpDiagonal = 0.0;
			  for (int j = 0; j <= this->MaxMomentum; ++j) 
			    {
			      TmpDiagonal += this->OneBodyInteractionFactorsUpUp[j] * Particles->AduAu(i, j);
			      TmpDiagonal += this->OneBodyInteractionFactorsDownDown[j] * Particles->AddAd(i, j);
			    }
			  vDestination.Re(i) += (this->EnergyShift + TmpDiagonal) * vSource.Re(i);
			  vDestination.Im(i) += (this->EnergyShift + TmpDiagonal) * vSource.Im(i);
			}
		    }
		  else
		    {
		      double TmpDiagonal = 0.0;
		      for (int i = firstComponent; i < LastComponent; ++i)
			{ 
			  TmpDiagonal = 0.0;
			  for (int j = 0; j <= this->MaxMomentum; ++j) 
			    TmpDiagonal += this->OneBodyInteractionFactorsUpUp[j] * Particles->AduAu(i, j);
			  vDestination.Re(i) += (this->EnergyShift + TmpDiagonal) * vSource.Re(i);
			  vDestination.Im(i) += (this->EnergyShift + TmpDiagonal) * vSource.Im(i);
			}
		    }
		else
		  if (this->OneBodyInteractionFactorsDownDown != 0)
		    {
		      double TmpDiagonal = 0.0;
		      for (int i = firstComponent; i < LastComponent; ++i)
			{ 
			  TmpDiagonal = 0.0;
			  for (int j = 0; j <= this->MaxMomentum; ++j) 
			    TmpDiagonal += this->OneBodyInteractionFactorsDownDown[j] * Particles->AddAd(i, j);
			  vDestination.Re(i) += (this->EnergyShift + TmpDiagonal) * vSource.Re(i);
			  vDestination.Im(i) += (this->EnergyShift + TmpDiagonal) * vSource.Im(i);
			}
		    }	
		  else
		    for (int i = firstComponent; i < LastComponent; ++i)
		      {
			vDestination.Re(i) += this->EnergyShift * vSource.Re(i);
			vDestination.Im(i) += this->EnergyShift * vSource.Im(i);
		      }
	      }
	}
    }
      
  return vDestination;
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::FastMultiplicationMemory(long allowedMemory)
{  
  this->NbrInteractionPerComponent = new int [this->Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
    this->NbrInteractionPerComponent[i] = 0;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start memory" << endl;

  QHEParticlePrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);

  long Memory = 0;
  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
    Memory += this->NbrInteractionPerComponent[i];

  cout << "nbr interaction = " << Memory << endl;
  long TmpMemory = allowedMemory - (2*sizeof (int*) + sizeof (int) + sizeof(double*)) * this->Particles->GetHilbertSpaceDimension();
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (2*sizeof (int) + sizeof(double)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (2*sizeof (int) + sizeof(double)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  TmpMemory = allowedMemory - (2*sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); i += this->FastMultiplicationStep)
	    Memory += this->NbrInteractionPerComponent[i];
	}
      int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
      for (int i = 0; i < ReducedSpaceDimension; ++i)
	TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
      delete[] this->NbrInteractionPerComponent;
      this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
      Memory = ((2*sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension) + (Memory * ((2 * sizeof (int)) + sizeof(double)));
    }
  else
    {
      Memory = ((2*sizeof (int*) + sizeof (int) + sizeof(double*)) * this->Particles->GetHilbertSpaceDimension()) + (Memory * ((2 * sizeof (int)) + sizeof(double)));
      this->FastMultiplicationStep = 1;
    }

  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  return Memory;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that have to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
  long Memory = 0;  
  unsigned L16Mask = (1u<<16)-1;
  unsigned H16Mask = (~0u)^L16Mask;
  int Index;
  double Coefficient;
  int NbrTranslation;
  int Dim = Particles->GetHilbertSpaceDimension();
  int LastComponent = nbrComponent + firstComponent;
  
  double Coefficient2;
  int m1, m2, m3, m4;      
  int SumIndices;
  int TmpNbrM34Values;
  unsigned* TmpM34Values;      
  int ReducedNbrInteractionFactors;
  for (int i = firstComponent; i < LastComponent; ++i)
    {    
      ReducedNbrInteractionFactors = 0;
      for (int m12 = 0; m12 < this->NbrM12IntraIndices; ++m12)
	{
	  m1 = this->M12IntraValue[m12] & L16Mask;
	  m2 = (this->M12IntraValue[m12] & H16Mask)>>16;
	  Coefficient = Particles->AuAu(i, m1, m2);
	  if (Coefficient != 0.0)
	    {
	      SumIndices = m1 + m2;
	      TmpNbrM34Values = this->NbrM34IntraValues[m12];
	      TmpM34Values = this->M34IntraValues[m12];
	      for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
		{
		  m3 = (TmpM34Values[m34]) & L16Mask;
		  m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
		  Index = Particles->AduAdu(m3, m4, Coefficient2, NbrTranslation);
		  if (Index < Dim)
		    {
		      ++Memory;
		      ++this->NbrInteractionPerComponent[i];
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM34IntraValues[m12];
	}
      ReducedNbrInteractionFactors = 0;
      for (int m12 = 0; m12 < this->NbrM12IntraIndices; ++m12)
	{
	  m1 = this->M12IntraValue[m12] & L16Mask;
	  m2 = (this->M12IntraValue[m12] & H16Mask)>>16;
	  Coefficient = Particles->AdAd(i, m1, m2);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = m1 + m2;
	      TmpNbrM34Values = this->NbrM34IntraValues[m12];
	      TmpM34Values = this->M34IntraValues[m12];
	      for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
		{
		  m3 = (TmpM34Values[m34]) & L16Mask;
		  m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
		  Index = Particles->AddAdd(m3, m4, Coefficient2, NbrTranslation);
		  if (Index < Dim)
		    {
		      ++Memory;
		      ++this->NbrInteractionPerComponent[i];
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM34IntraValues[m12];
	}
      
      ReducedNbrInteractionFactors = 0;
      for (int m12 = 0; m12 < this->NbrM12InterIndices; ++m12)
	{
	  m1 = this->M12InterValue[m12] & L16Mask;
	  m2 = (this->M12InterValue[m12] & H16Mask)>>16;
	  Coefficient = Particles->AuAd(i, m1, m2);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = m1 + m2;
	      TmpNbrM34Values = this->NbrM34InterValues[m12];
	      TmpM34Values = this->M34InterValues[m12];
	      for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
		{
		  m3 = (TmpM34Values[m34]) & L16Mask;
		  m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
		  Index = Particles->AduAdd(m3, m4, Coefficient2, NbrTranslation);
		  if (Index < Dim)
		    {
		      ++Memory;
		      ++this->NbrInteractionPerComponent[i];
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM34InterValues[m12];
	}

      // one-particle terms:
      if (this->OneBodyInteractionFactorsUpUp != 0) 
	if (this->OneBodyInteractionFactorsDownDown != 0)
	  {
	    double TmpDiagonal = 0.0;
	    TmpDiagonal = 0.0;
	    for (int j = 0; j <= this->MaxMomentum; ++j) 
	      {
		TmpDiagonal += this->OneBodyInteractionFactorsUpUp[j] * Particles->AduAu(i, j);
		TmpDiagonal += this->OneBodyInteractionFactorsDownDown[j] * Particles->AddAd(i, j);
	      }
	    if (TmpDiagonal!=0.0)
	      {
		++Memory;
		++this->NbrInteractionPerComponent[i];
	      }
	  }
	else
	  {
	    double TmpDiagonal = 0.0;
	    TmpDiagonal = 0.0;
	    for (int j = 0; j <= this->MaxMomentum; ++j) 
	      TmpDiagonal += this->OneBodyInteractionFactorsUpUp[j] * Particles->AduAu(i, j);
	    if (TmpDiagonal!=0.0)
	      {
		++Memory;
		++this->NbrInteractionPerComponent[i];
	      }
	  }
      else
	if (this->OneBodyInteractionFactorsDownDown != 0)
	  {
	    double TmpDiagonal = 0.0;
	    TmpDiagonal = 0.0;
	    for (int j = 0; j <= this->MaxMomentum; ++j) 
	      TmpDiagonal += this->OneBodyInteractionFactorsDownDown[j] * Particles->AddAd(i, j);
	    if (TmpDiagonal!=0.0)
	      {
		++Memory;
		++this->NbrInteractionPerComponent[i];
	      }
	  }	
    }
  Memory = ((2*sizeof (int*) + sizeof (int) + sizeof(double*)) * this->Particles->GetHilbertSpaceDimension() + 
	    Memory *  (sizeof (int) + sizeof(double) + sizeof(int)));
  return Memory;  
}

// enable fast multiplication algorithm
//

void AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::EnableFastMultiplication()
{
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start enable" << endl;
  int ReducedSpaceDimension = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new double* [ReducedSpaceDimension];
  this->InteractionPerComponentNbrTranslation = new int* [ReducedSpaceDimension];
  
  QHEParticlePrecalculationOperation Operation(this, false);
  Operation.ApplyOperation(this->Architecture);
  
  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  
  
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that have to be precalcualted

void AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{
  unsigned L16Mask = (1u<<16)-1;
  unsigned H16Mask = (~0u)^L16Mask;
  int Index;
  double Coefficient;
  int NbrTranslation;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int* TmpNbrTranslationArray;
  int Dim = Particles->GetHilbertSpaceDimension();
  int LastComponent = firstComponent + nbrComponent;
  double Coefficient2;
  int m1, m2, m3, m4;      
  int SumIndices;
  int TmpNbrM34Values;
  unsigned* TmpM34Values;
  int Pos;
  int ReducedNbrInteractionFactors;
  int PosIndex = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++PosIndex;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)  
    {
      this->InteractionPerComponentIndex[PosIndex] = new int [this->NbrInteractionPerComponent[i]];
      this->InteractionPerComponentCoefficient[PosIndex] = new double [this->NbrInteractionPerComponent[i]];      
      this->InteractionPerComponentNbrTranslation[PosIndex] = new int [this->NbrInteractionPerComponent[i]];
      TmpIndexArray = this->InteractionPerComponentIndex[PosIndex];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[PosIndex];
      TmpNbrTranslationArray = this->InteractionPerComponentNbrTranslation[PosIndex];
      ++PosIndex;

      Pos = 0;
      
      ReducedNbrInteractionFactors = 0;
      for (int m12 = 0; m12 < this->NbrM12IntraIndices; ++m12)
	{
	  m1 = this->M12IntraValue[m12] & L16Mask;
	  m2 = (this->M12IntraValue[m12] & H16Mask)>>16;
	  Coefficient = Particles->AuAu(i, m1, m2);
	  if (Coefficient != 0.0)
	    {
	      SumIndices = m1 + m2;
	      TmpNbrM34Values = this->NbrM34IntraValues[m12];
	      TmpM34Values = this->M34IntraValues[m12];
	      for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
		{
		  m3 = (TmpM34Values[m34]) & L16Mask;
		  m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
		  Index = Particles->AduAdu(m3, m4, Coefficient2, NbrTranslation);
		  if (Index < Dim)
		    {
		      TmpIndexArray[Pos] = Index;
		      TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * this->InteractionFactorsUpUp[ReducedNbrInteractionFactors];
		      TmpNbrTranslationArray[Pos] = NbrTranslation;
		      ++Pos;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM34IntraValues[m12];
	}
      ReducedNbrInteractionFactors = 0;
      for (int m12 = 0; m12 < this->NbrM12IntraIndices; ++m12)
	{
	  m1 = this->M12IntraValue[m12] & L16Mask;
	  m2 = (this->M12IntraValue[m12] & H16Mask)>>16;
	  Coefficient = Particles->AdAd(i, m1, m2);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = m1 + m2;
	      TmpNbrM34Values = this->NbrM34IntraValues[m12];
	      TmpM34Values = this->M34IntraValues[m12];
	      for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
		{
		  m3 = (TmpM34Values[m34]) & L16Mask;
		  m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
		  Index = Particles->AddAdd(m3, m4, Coefficient2, NbrTranslation);
		  if (Index < Dim)
		    {
		      TmpIndexArray[Pos] = Index;
		      TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * this->InteractionFactorsDownDown[ReducedNbrInteractionFactors];
		      TmpNbrTranslationArray[Pos] = NbrTranslation;
		      ++Pos;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM34IntraValues[m12];
	}
      
      ReducedNbrInteractionFactors = 0;
      for (int m12 = 0; m12 < this->NbrM12InterIndices; ++m12)
	{
	  m1 = this->M12InterValue[m12] & L16Mask;
	  m2 = (this->M12InterValue[m12] & H16Mask)>>16;
	  Coefficient = Particles->AuAd(i, m1, m2);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = m1 + m2;
	      TmpNbrM34Values = this->NbrM34InterValues[m12];
	      TmpM34Values = this->M34InterValues[m12];
	      for (int m34 = 0; m34 < TmpNbrM34Values; ++m34)
		{
		  m3 = (TmpM34Values[m34]) & L16Mask;
		  m4 = ((TmpM34Values[m34]) & H16Mask)>>16;
		  Index = Particles->AduAdd(m3, m4, Coefficient2, NbrTranslation);
		  if (Index < Dim)
		    {
		      TmpIndexArray[Pos] = Index;
		      TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * this->InteractionFactorsUpDown[ReducedNbrInteractionFactors];
		      TmpNbrTranslationArray[Pos] = NbrTranslation;
		      ++Pos;
		    }
		  ++ReducedNbrInteractionFactors;
		}
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM34InterValues[m12];
	}
      
      // one-particle terms:
      if (this->OneBodyInteractionFactorsUpUp != 0) 
	if (this->OneBodyInteractionFactorsDownDown != 0)
	  {
	    double TmpDiagonal = 0.0;
	    TmpDiagonal = 0.0;
	    for (int j = 0; j <= this->MaxMomentum; ++j) 
	      {
		TmpDiagonal += this->OneBodyInteractionFactorsUpUp[j] * Particles->AduAu(i, j);
		TmpDiagonal += this->OneBodyInteractionFactorsDownDown[j] * Particles->AddAd(i, j);
	      }
	    if (TmpDiagonal!=0.0)
	      {
		TmpIndexArray[Pos] = i;
		TmpCoefficientArray[Pos] = (TmpDiagonal);
		TmpNbrTranslationArray[Pos] = 0;
		++Pos;
	      }
	  }
	else
	  {
	    double TmpDiagonal = 0.0;
	    TmpDiagonal = 0.0;
	    for (int j = 0; j <= this->MaxMomentum; ++j) 
	      TmpDiagonal += this->OneBodyInteractionFactorsUpUp[j] * Particles->AduAu(i, j);
	    if (TmpDiagonal!=0.0)
	      {
		TmpIndexArray[Pos] = i;
		TmpCoefficientArray[Pos] = (TmpDiagonal);
		TmpNbrTranslationArray[Pos] = 0;
		++Pos;
	      }
	  }
      else
	if (this->OneBodyInteractionFactorsDownDown != 0)
	  {
	    double TmpDiagonal = 0.0;
	    TmpDiagonal = 0.0;
	    for (int j = 0; j <= this->MaxMomentum; ++j) 
	      TmpDiagonal += this->OneBodyInteractionFactorsDownDown[j] * Particles->AddAd(i, j);
	    if (TmpDiagonal!=0.0)
	      {
		TmpIndexArray[Pos] = i;
		TmpCoefficientArray[Pos] = (TmpDiagonal);
		TmpNbrTranslationArray[Pos] = 0;
		++Pos;
	      }
	  }
    }
}


// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs

bool AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::SavePrecalculation (char* fileName)
{
  if (this->FastMultiplicationFlag)
    {
      ofstream File;
      File.open(fileName, ios::binary | ios::out);
      int Tmp = this->Particles->GetHilbertSpaceDimension();
      File.write((char*) &(Tmp), sizeof(int));
      File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
      Tmp /= this->FastMultiplicationStep;
      if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
	++Tmp;
      File.write((char*) this->NbrInteractionPerComponent, sizeof(int) * Tmp);
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * this->NbrInteractionPerComponent[i]);	  
	}
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentCoefficient[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
	}
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentNbrTranslation[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
	}
      File.close();
      return true;
    }
  else
    {
      return false;
    }  
}

// load precalculations from a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be read
// return value = true if no error occurs

bool AbstractQHEOnTorusWithSpinAndMagneticTranslationsHamiltonian::LoadPrecalculation (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  int Tmp;
  File.read((char*) &(Tmp), sizeof(int));
  if (Tmp != this->Particles->GetHilbertSpaceDimension())
    {
      File.close();
      return false;
    }
  File.read((char*) &(this->FastMultiplicationStep), sizeof(int));
  Tmp /= this->FastMultiplicationStep;
  if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
    ++Tmp;
  this->NbrInteractionPerComponent = new int [Tmp];
  File.read((char*) this->NbrInteractionPerComponent, sizeof(int) * Tmp);
  this->InteractionPerComponentIndex = new int* [Tmp];
  this->InteractionPerComponentCoefficient = new double* [Tmp];
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * this->NbrInteractionPerComponent[i]);	  
    }
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentCoefficient[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
    }
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentNbrTranslation[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
    }
  File.close();
  this->FastMultiplicationFlag = true;
  return true;
}

