////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                        class author: Gunnar Möller                         //
//                                                                            //
//  class of tight binding model for the square lattice with homogeneous flux //
//                                                                            //
//                        last modification : 08/05/2012                      //
//                                                                            //
//                                                                            //
//   Note: throughout code x and y refer to directions of basis vectors for   //
//   triangular lattice (1,0) and (1/2, sqrt(3)/2)                            //
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
#include "Tools/FTITightBinding/TightBindingModelHofstadterTriangular.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <iostream>
using std::cout;
using std::endl;

#define DEBUG_OUTPUT


// default constructor
//
// nbrCellsX = number of unit cells in the x direction
// nbrCellsY = number of unit cella in the y direction
// unitCellX = number of sites in unit cell in x direction
// unitCellY = number of sites in unit cell in y direction
// nbrFlux = number of flux quanta per unit cell
// axis = direction of Landau gauge within cell ('x' or 'y')
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
TightBindingModelHofstadterTriangular::TightBindingModelHofstadterTriangular(int nbrCellX, int nbrCellY, int unitCellX, int unitCellY, int nbrFlux, char axis,
									     double gammaX, double gammaY, 
									     AbstractArchitecture* architecture, bool storeOneBodyMatrices, bool useEmbedding)
{
  this->NbrSiteX = nbrCellX;
  this->NbrSiteY = nbrCellY;
  this->UnitCellX = unitCellX;
  this->UnitCellY = unitCellY;
  this->TwistAngle=M_PI/3.0;
  this->LatticeVector1.Resize(2);
  this->LatticeVector1[0] = this->UnitCellX;
  this->LatticeVector1[1] = 0.0;
  this->LatticeVector1.Resize(2);
  this->LatticeVector2[0] = 0.5*this->UnitCellY;
  this->LatticeVector2[1] = 0.5*sqrt(3.0)*this->UnitCellY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = this->UnitCellX * this->UnitCellY;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->LandauGaugeAxis=axis;
  this->Architecture = architecture;

  if (storeOneBodyMatrices == true)
    {
      this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
    }
  else
    {
      this->OneBodyBasis = 0;
    }
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    }
  this->SetNbrFluxQuanta(nbrFlux);
  if (useEmbedding == true)
    {
      this->SetNaturalEmbedding();
    }
  else
    {
      this->SetNoEmbedding();
    }
  this->ComputeBandStructure();  
  
  //or else, use explicitly hardcoded embedding
  //this->CoreComputeBandStructureWithEmbedding(0, this->GetNbrStatePerBand());
}

// destructor
//

TightBindingModelHofstadterTriangular::~TightBindingModelHofstadterTriangular()
{
}

// get the position of a sublattice site
//
// position = reference on a vector where the answer is supplied
// sublatticeIndex = index of the sub-lattice position
void TightBindingModelHofstadterTriangular::GetSublatticeVector(RealVector &position, int sublatticeIndex)
{
  int x1, x2;
  this->DecodeSublatticeIndex(sublatticeIndex, x1, x2);
  if (position.GetVectorDimension()!=2)
    position.Resize(2);
  position.Copy(this->LatticeVector1, (double)x1/this->UnitCellX);
  position.AddLinearCombination((double)x2/this->UnitCellY, this->LatticeVector2);
}



// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelHofstadterTriangular::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  double K1;
  double K2;
  //int sublXf, sublYf; // final site sublattice positions
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int Index = this->GetLinearizedMomentumIndex(kx, ky);
	  if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	    {
	      K1 = this->KxFactor*(((double) kx) + this->GammaX);
	      K2 = this->KyFactor*(((double) ky) + this->GammaY);
#ifdef DEBUG_OUTPUT
	      cout << "Sector kx="<<kx<<", ky="<<ky<<" (KX="<<K1<<", KY="<<K2<<")"<<endl;
#endif
	      // construct magnetic unit cell:
	      HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	      
	      Complex TranslationPhase;
	      switch (this->LandauGaugeAxis)
		{
		case 'y': {
		  for (int i=0; i<UnitCellX; ++i)
		    {
		      Complex Phase=Polar(1.0,2.0*M_PI*this->FluxDensity*(double)i);
		      Complex Phase2=Polar(1.0,2.0*M_PI*this->FluxDensity*((double)i-0.5));
		      Complex Phase3=Polar(1.0,-2.0*M_PI*this->FluxDensity*((double)i+0.5));

		      for (int j=0; j<UnitCellY; ++j)
			{
                  
                  
			  
			  int InitialIndex = this->EncodeSublatticeIndex(i, j, K1, K2, TranslationPhase); // initial TranlationPhase always one, so can be discarded
			  double InitialEmbeddingPhase=this->GetEmbeddingPhase(InitialIndex, K1, K2);
			  int FinalIndex = this->EncodeSublatticeIndex(i-1, j, K1, K2, TranslationPhase);
			  double FinalEmbeddingPhase=this->GetEmbeddingPhase(FinalIndex, K1, K2);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Polar(InitialEmbeddingPhase-FinalEmbeddingPhase));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<", embeddingPhase="<<Polar(InitialEmbeddingPhase-FinalEmbeddingPhase)<<endl;
#endif

			  FinalIndex = this->EncodeSublatticeIndex(i+1, j, K1, K2, TranslationPhase);
			  FinalEmbeddingPhase=this->GetEmbeddingPhase(FinalIndex, K1, K2);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Polar(InitialEmbeddingPhase-FinalEmbeddingPhase));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<", embeddingPhase="<<Polar(InitialEmbeddingPhase-FinalEmbeddingPhase)<<endl;
#endif
			  
			  FinalIndex = this->EncodeSublatticeIndex(i, j-1, K1, K2, TranslationPhase);
			  FinalEmbeddingPhase=this->GetEmbeddingPhase(FinalIndex, K1, K2);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Conj(Phase)*Polar(InitialEmbeddingPhase-FinalEmbeddingPhase));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase*Conj(Phase)<<", embeddingPhase="<<Polar(InitialEmbeddingPhase-FinalEmbeddingPhase)<<endl;
#endif

			  FinalIndex = this->EncodeSublatticeIndex(i, j+1, K1, K2, TranslationPhase);
			  FinalEmbeddingPhase=this->GetEmbeddingPhase(FinalIndex, K1, K2);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Phase*Polar(InitialEmbeddingPhase-FinalEmbeddingPhase));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase*Phase<<", embeddingPhase="<<Polar(InitialEmbeddingPhase-FinalEmbeddingPhase)<<endl;
#endif
                
			  FinalIndex = this->EncodeSublatticeIndex(i-1, j+1, K1, K2, TranslationPhase);
			  FinalEmbeddingPhase=this->GetEmbeddingPhase(FinalIndex, K1, K2);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Phase2*Polar(InitialEmbeddingPhase-FinalEmbeddingPhase));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase*Conj(Phase)<<", embeddingPhase="<<Polar(InitialEmbeddingPhase-FinalEmbeddingPhase)<<endl;
#endif
                
			  FinalIndex = this->EncodeSublatticeIndex(i+1, j-1, K1, K2, TranslationPhase);
			  FinalEmbeddingPhase=this->GetEmbeddingPhase(FinalIndex, K1, K2);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Phase3*Polar(InitialEmbeddingPhase-FinalEmbeddingPhase));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase*Phase<<", embeddingPhase="<<Polar(InitialEmbeddingPhase-FinalEmbeddingPhase)<<endl;
#endif
			}
		    }
		  break;
		}

		case 'x': { // needs checking
		  for (int j=0; j<UnitCellY; ++j)
		    {
		      
		      for (int i=0; i<UnitCellX; ++i)
			{
			  Complex Phase=Polar(1.0,-M_PI*this->FluxDensity*((double)j-(double)i)*((double)i+(double)j));

			  int InitialIndex = this->EncodeSublatticeIndex(i, j, K1, K2, TranslationPhase); // TranlationPhase always one, so can be discarded
			  int FinalIndex = this->EncodeSublatticeIndex(i+1, j, K1, K2, TranslationPhase);
			  double InitialEmbeddingPhase=this->GetEmbeddingPhase(InitialIndex, K1, K2);
			  double FinalEmbeddingPhase=this->GetEmbeddingPhase(FinalIndex, K1, K2);

			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -Conj(Phase)*TranslationPhase*Polar(InitialEmbeddingPhase-FinalEmbeddingPhase));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-Conj(Phase)*TranslationPhase<<endl;
#endif

			  FinalIndex = this->EncodeSublatticeIndex(i-1, j, K1, K2, TranslationPhase);
			  FinalEmbeddingPhase=this->GetEmbeddingPhase(FinalIndex, K1, K2);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -Phase*TranslationPhase*Polar(InitialEmbeddingPhase-FinalEmbeddingPhase));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-Phase*TranslationPhase<<endl;
#endif
			  
			  FinalIndex = this->EncodeSublatticeIndex(i, j+1, K1, K2, TranslationPhase);
			  FinalEmbeddingPhase=this->GetEmbeddingPhase(FinalIndex, K1, K2);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Polar(InitialEmbeddingPhase-FinalEmbeddingPhase));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif

			  FinalIndex = this ->EncodeSublatticeIndex(i, j-1, K1, K2, TranslationPhase);
			  FinalEmbeddingPhase=this->GetEmbeddingPhase(FinalIndex, K1, K2);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Polar(InitialEmbeddingPhase-FinalEmbeddingPhase));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif
			}
		    }
		  break;
		}

		default:
 		  cout << "Invalid Landau quantization axis encountered in ParticleOnLatticeDeltaHamiltonian."<<endl;
		  exit(1);
		  break;
		}

#ifdef DEBUG_OUTPUT
	      cout << TmpOneBodyHamiltonian<< endl;
#endif
	    

	      if (this->OneBodyBasis != 0)
		{
		  ComplexMatrix TmpMatrix(this->NbrBands, this->NbrBands, true);
		  TmpMatrix.SetToIdentity();
		  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
		  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
		  this->OneBodyBasis[Index] = TmpMatrix;
		  for (int i = 0; i < this->NbrBands; ++i)
		    this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		}
	      else
		{
		  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
		  TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif
		  for (int i = 0; i < this->NbrBands; ++i)
		    this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		}
	    }
	}
    }
}


// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelHofstadterTriangular::CoreComputeBandStructureWithEmbedding(long minStateIndex, long nbrStates)
{
  /// @note this function is deprecated: CoreComputeBandStructure now uses the embedding set in 2DAbstractTightBindingModel.
  /// @todo calculate single-particle states using Fourier transform with respect to site position
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  double K1;
  double K2;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int Index = this->GetLinearizedMomentumIndex(kx, ky);
	  if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	    {
	      K1 = this->KxFactor*(((double) kx) + this->GammaX);
	      K2 = this->KyFactor*(((double) ky) + this->GammaY);
#ifdef DEBUG_OUTPUT
	      cout << "Sector kx="<<kx<<", ky="<<ky<<" (KX="<<K1<<", KY="<<K2<<")"<<endl;
#endif
              Complex PhaseEmbeddingX =  Phase(K1/((double)this->UnitCellX));
              Complex PhaseEmbeddingY =  Phase(K2/((double)this->UnitCellY));
	      // construct magnetic unit cell:
	      
	      HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	      
	      Complex TranslationPhase;
	      switch (this->LandauGaugeAxis)
		{
		case 'y': {
		  for (int i=0; i<UnitCellX; ++i)
		    {
		      Complex Phase=Polar(1.0,2.0*M_PI*this->FluxDensity*(double)i);
		      for (int j=0; j<UnitCellY; ++j)
			{
			  int InitialIndex = this->EncodeSublatticeIndex(i, j, 0.0, 0.0, TranslationPhase); // TranlationPhase always one, so can be discarded
			  int FinalIndex = this->EncodeSublatticeIndex(i+1, j, 0.0, 0.0, TranslationPhase);

			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*PhaseEmbeddingX);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<", embeddingPhase="<<PhaseEmbeddingX<<endl;
#endif

			  FinalIndex = this->EncodeSublatticeIndex(i-1, j, 0.0,0.0, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Conj(PhaseEmbeddingX));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<", embeddingPhase="<<Conj(PhaseEmbeddingX)<<endl;
#endif
			  
			  FinalIndex = this->EncodeSublatticeIndex(i, j+1, 0.0, 0.0, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Conj(Phase)*PhaseEmbeddingY);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase*Conj(Phase)<<", embeddingPhase="<<PhaseEmbeddingY<<endl;
#endif

			  FinalIndex = this->EncodeSublatticeIndex(i, j-1,0.0, 0.0, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Phase*Conj(PhaseEmbeddingY));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase*Phase<<", embeddingPhase="<<Conj(PhaseEmbeddingY)<<endl;
#endif

			}
		    }
		  break;
		}


		case 'x': {
		  for (int j=0; j<UnitCellY; ++j)
		    {
		      Complex Phase=Polar(1.0,-2.0*M_PI*this->FluxDensity*(double)j);
		      for (int i=0; i<UnitCellX; ++i)
			{

			  int InitialIndex = this->EncodeSublatticeIndex(i, j, K1, K2, TranslationPhase); // TranlationPhase always one, so can be discarded
			  int FinalIndex = this->EncodeSublatticeIndex(i+1, j, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -Conj(Phase)*TranslationPhase*PhaseEmbeddingX);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-Conj(Phase)*TranslationPhase<<endl;
#endif

			  FinalIndex = this->EncodeSublatticeIndex(i-1, j, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -Phase*TranslationPhase*Conj(PhaseEmbeddingX));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-Phase*TranslationPhase<<endl;
#endif
			  
			  FinalIndex = this->EncodeSublatticeIndex(i, j+1, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*PhaseEmbeddingY);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif

			  FinalIndex = this ->EncodeSublatticeIndex(i, j-1, K1, K2, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -TranslationPhase*Conj(PhaseEmbeddingY));
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif

			}
		    }
		  break;
		}
		default:
 		  cout << "Invalid Landau quantization axis encountered in ParticleOnLatticeDeltaHamiltonian."<<endl;
		  exit(1);
		  break;
		}

#ifdef DEBUG_OUTPUT
	      cout << TmpOneBodyHamiltonian<< endl;
#endif
	    

	      if (this->OneBodyBasis != 0)
		{
		  ComplexMatrix TmpMatrix(this->NbrBands, this->NbrBands, true);
		  TmpMatrix.SetToIdentity();
		  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
		  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
		  this->OneBodyBasis[Index] = TmpMatrix;
		  for (int i = 0; i < this->NbrBands; ++i)
		    this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		}
	      else
		{
		  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
		  TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif
		  for (int i = 0; i < this->NbrBands; ++i)
		    this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		}
	    }
	}
    }
}


// set natural embedding, i.e., at positions of a uniform lattice
//
void TightBindingModelHofstadterTriangular::SetNaturalEmbedding()
{
  Complex phase;
  this->EmbeddingX.Resize(UnitCellX*UnitCellY);
  this->EmbeddingY.Resize(UnitCellX*UnitCellY);
  double invX = 1.0/((double)this->UnitCellX);
  double invY = 1.0/((double)this->UnitCellY);
  for (int i=0; i<this->UnitCellX; ++i)
    for (int j=0; j<this->UnitCellY; ++j)
      {
	int sublattice = this->EncodeSublatticeIndex(i,j,0.0,0.0,phase);
	this->EmbeddingX[sublattice] = i*invX;
	this->EmbeddingY[sublattice] = j*invY;
      }
}


// nbrFluxQuanta = number of quanta of flux piercing the unit cell
void TightBindingModelHofstadterTriangular::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->FluxDensity = ((double)NbrFluxQuanta)/this->NbrBands;
#ifdef DEBUG_OUTPUT
  cout <<"this->FluxDensity = "<< this->FluxDensity<<endl;
#endif
  switch (this->LandauGaugeAxis)
    {
    case 'x':
      this->LxTranslationPhase = 1.0;  // no phase for translating in the y-direction in Landau gauge...
      this->LyTranslationPhase = Polar(1.0, 2.0*M_PI*FluxDensity*this->UnitCellY);
          
      cout << "'x-axis' gauge: LyTranslationPhase= exp(I*"<<2.0*M_PI*FluxDensity*this->UnitCellY<<")="<<LyTranslationPhase<<endl;
      break;
    case 'y':
    
      this->LxTranslationPhase = Polar(1.0, 2.0*M_PI*FluxDensity*this->UnitCellX);
      this->LyTranslationPhase = 1.0;  // no phase for translating in the y-direction in Landau gauge...
      cout << "'y-axis' gauge: LxTranslationPhase= exp(I*"<<-2.0*M_PI*FluxDensity*this->UnitCellX<<")="<<LxTranslationPhase<<endl;
      break;
    default:
      cout << "Unknown Quantization axis! Exiting TightBindingModelHofstadterSquare..."<<endl;
      exit(1);
      break;
    }
}
      
 
// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// KX = current momentum in x-direction
// KY = current momentum in y-direction
// translationPhase = phase factor associated with any crossings of unit cell boundary
//
int TightBindingModelHofstadterTriangular::EncodeSublatticeIndex(int posx, int posy, double KX, double KY, Complex &translationPhase)
{
  /// @note this function overloads a virtual function with the same name, but different signature.
  //cout << "Encoding " << posx<<", "<<posy<<": ";
  int numXTranslations=0, numYTranslations=0;  
  while (posx<0)
    {
      posx+=this->UnitCellX;
      ++numXTranslations;      
    }
  while (posx>=this->UnitCellX)
    {
      posx-=this->UnitCellX;
      --numXTranslations;
    }
  while (posy<0)
    {
      posy+=this->UnitCellY;
      ++numYTranslations;      
    }
  while (posy>=this->UnitCellY)
    {
      posy-=this->UnitCellY;
      --numYTranslations;
    }
  int rst = posx + this->UnitCellX*posy;
  // determine phase for shifting site to the simulation cell:
  Complex tmpPhase(1.0,0.0);
  Complex tmpPhase2;
  translationPhase=tmpPhase;
  if (numXTranslations>0)
    tmpPhase2=LxTranslationPhase;
  else
    tmpPhase2=Conj(LxTranslationPhase);
  for (int i=0; i<abs(numXTranslations); ++i)
    tmpPhase*=tmpPhase2;
  //cout<<" tmpPhaseX="<<tmpPhase;
  for (int y=1;y<=posy; ++y)
    translationPhase*=tmpPhase;
  translationPhase*=Polar(1.0, KX*numXTranslations);
  tmpPhase=1.0;
  if (numYTranslations>0)
    tmpPhase2=LyTranslationPhase;
  else
    tmpPhase2=Conj(LyTranslationPhase);
  for (int i=0; i<abs(numYTranslations); ++i)
    tmpPhase*=tmpPhase2;
  //cout<<" tmpPhaseY="<<tmpPhase;
  for (int x=1;x<=posx; ++x)
    translationPhase*=tmpPhase;
  translationPhase*=Polar(1.0, KY*numYTranslations);
  //  cout << "tX="<<numXTranslations<< ", tY="<<numYTranslations<<", translationPhase= " <<translationPhase<<endl;
  return rst;
}




// get the eigenstates in real space, using CoreComputeBandStructureWithEmbedding
// 
// return value = tight binding eigenvectors

ComplexMatrix TightBindingModelHofstadterTriangular::GetRealSpaceTightBindingEigenstates()
{
  double LogTranslationPhaseX= -2.0*M_PI*this->FluxDensity*this->UnitCellX;
  ComplexMatrix EigenStates(this->NbrBands *  this->NbrStatePerBand,this->NbrBands *  this->NbrStatePerBand ,true);
  int Kx;  int Ky;
  int K1;  int K2;
  int OrbitalIndex = 0;
  int PosXUnitCell = 0;
  int PosYUnitCell = 0;
  for(int i = 0; i <this->NbrBands *  this->NbrStatePerBand;i++)
    {
      int BandNumber = i/this->NbrStatePerBand;
      int MomentumIndex = i%this->NbrStatePerBand;
 
      this->GetLinearizedMomentumIndex(MomentumIndex,Kx,Ky);

      K1 = this->KxFactor*(((double) Kx) + this->GammaX);
      K2 = this->KyFactor*(((double) Ky) + this->GammaY);
      for(int j = 0; j <this->NbrBands *  this->NbrStatePerBand;j++) 
	{ 
	  this->GetRealSpaceTightBindingLinearizedIndex(j, PosXUnitCell, PosYUnitCell, OrbitalIndex);
	  int TotalPosY = PosYUnitCell*this->UnitCellY + OrbitalIndex/this->UnitCellX;
	  EigenStates[i][j] = this->OneBodyBasis[MomentumIndex][BandNumber][OrbitalIndex] * Phase(K1*PosXUnitCell + K2*PosYUnitCell)* Phase(PosXUnitCell* LogTranslationPhaseX*TotalPosY) ;
	}
    }
  return EigenStates;
}


/*
// get the eigenstates in real space, using CoreComputeBandStructureWithEmbedding
// 
// return value = tight binding eigenvectors

ComplexMatrix TightBindingModelHofstadterSquare::GetRealSpaceTightBindingEigenstates()
{
ComplexMatrix EigenStates(this->NbrBands *  this->NbrStatePerBand,this->NbrBands *  this->NbrStatePerBand ,true);
int Kx;  int Ky;
int K1;  int K2;
int OrbitalIndex = 0;
int UnitCellsX = 0;
int UnitCellsY = 0;
for(int i = 0; i <this->NbrBands *  this->NbrStatePerBand;i++)
{
int BandNumber = i/this->NbrStatePerBand;
int MomentumIndex = i%this->NbrStatePerBand;
 
this->GetLinearizedMomentumIndex(MomentumIndex,Kx,Ky);

K1 = this->KxFactor*(((double) Kx) + this->GammaX);
K2 = this->KyFactor*(((double) Ky) + this->GammaY);
for(int j = 0; j <this->NbrBands *  this->NbrStatePerBand;j++) 
{ 
this->GetRealSpaceTightBindingLinearizedIndex(j, UnitCellsX, UnitCellsY, OrbitalIndex);
EigenStates[i][j] = this->OneBodyBasis[MomentumIndex][BandNumber][OrbitalIndex] * Phase(K1*UnitCellsX + K2*UnitCellsY);
}
}
return EigenStates;
}

*/

// get the tight binding hamiltonian in real space 
// 
// return value = tight binding hamiltonian

HermitianMatrix TightBindingModelHofstadterTriangular::GetRealSpaceTightBindingHamiltonian()
{
  int* NbrConnectedOrbitals = new int [this->NbrBands];
  int** OrbitalIndices = new int* [this->NbrBands];
  int** SpatialIndices = new int* [this->NbrBands];
  Complex** HoppingAmplitudes = new Complex* [this->NbrBands];
  for(int i = 0; i < this->NbrBands; i++)
    NbrConnectedOrbitals[i] = 4;
  for (int i = 0; i < this->NbrBands; ++i)
    {
      OrbitalIndices[i] = new int[NbrConnectedOrbitals[i]];
      SpatialIndices[i] = new int[2 * NbrConnectedOrbitals[i]];
      HoppingAmplitudes[i] = new Complex[NbrConnectedOrbitals[i]];
    }


  int   NumXTranslations;
  int   NumYTranslations;

  Complex translationPhase=1.0;
  Complex tmpPhase, tmpPhase2;
  switch (this->LandauGaugeAxis)
    {
    case 'y': {
#ifdef DEBUG_OUTPUT
      cout <<" this->LxTranslationPhase  = " << this->LxTranslationPhase<<endl;
#endif

      for (int PosX = 0; PosX <this->UnitCellX; PosX++)
	{
	  Complex PhaseY = Phase(2.0*M_PI*this->FluxDensity*(double)PosX);

	  for (int PosY = 0; PosY <this->UnitCellY; PosY++)
	    {
	      int TmpPosition =  this->EncodeSublatticeIndex(PosX, PosY,NumXTranslations,NumYTranslations,translationPhase);
	      int TmpIndex = 0;
	      OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(PosX+1, PosY,NumXTranslations,NumYTranslations,translationPhase);
	      SpatialIndices[TmpPosition][TmpIndex * 2] = -1*NumXTranslations;
	      SpatialIndices[TmpPosition][(TmpIndex * 2) +1] =  -1*NumYTranslations;
	      HoppingAmplitudes[TmpPosition][TmpIndex] = -1.0* translationPhase;
#ifdef DEBUG_OUTPUT
	      cout <<TmpPosition<< " " << OrbitalIndices[TmpPosition][TmpIndex]<< " " << SpatialIndices[TmpPosition][TmpIndex * 2]<< " " <<SpatialIndices[TmpPosition][(TmpIndex * 2)+1]<<" "<<HoppingAmplitudes[TmpPosition][TmpIndex]<<endl;
#endif
	      ++TmpIndex;
	      OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(PosX, PosY+1,NumXTranslations,NumYTranslations,translationPhase);
	      SpatialIndices[TmpPosition][TmpIndex * 2] = -1*NumXTranslations;
	      SpatialIndices[TmpPosition][(TmpIndex * 2) +1] = -1*NumYTranslations;
	      HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0*Conj(PhaseY)*translationPhase;
#ifdef DEBUG_OUTPUT
	      cout <<TmpPosition<< " " << OrbitalIndices[TmpPosition][TmpIndex]<< " " << SpatialIndices[TmpPosition][TmpIndex * 2]<< " " <<SpatialIndices[TmpPosition][(TmpIndex * 2)+1]<<" "<<HoppingAmplitudes[TmpPosition][TmpIndex]<<endl;
#endif

	      ++TmpIndex;
	      OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(PosX-1, PosY,NumXTranslations,NumYTranslations, translationPhase);
	      SpatialIndices[TmpPosition][TmpIndex * 2] = -1*NumXTranslations;
	      SpatialIndices[TmpPosition][(TmpIndex * 2) +1] = -1*NumYTranslations;
	      HoppingAmplitudes[TmpPosition][TmpIndex] = -1.0* translationPhase;
#ifdef DEBUG_OUTPUT
	      cout <<TmpPosition<< " " << OrbitalIndices[TmpPosition][TmpIndex]<< " " << SpatialIndices[TmpPosition][TmpIndex * 2]<< " " <<SpatialIndices[TmpPosition][(TmpIndex * 2)+1]<<" "<<HoppingAmplitudes[TmpPosition][TmpIndex]<<endl;
#endif
	      ++TmpIndex;
	      OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(PosX, PosY-1,NumXTranslations,NumYTranslations, translationPhase);
	      SpatialIndices[TmpPosition][TmpIndex * 2] = -1*NumXTranslations;
	      SpatialIndices[TmpPosition][(TmpIndex * 2) +1] = -1*NumYTranslations;
	      HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0 * PhaseY* translationPhase;
#ifdef DEBUG_OUTPUT
	      cout <<TmpPosition<< " " << OrbitalIndices[TmpPosition][TmpIndex]<< " " << SpatialIndices[TmpPosition][TmpIndex * 2]<< " " <<SpatialIndices[TmpPosition][(TmpIndex * 2)+1]<<" "<<HoppingAmplitudes[TmpPosition][TmpIndex]<<endl;
#endif
	    }
	}
      break;
    }


    case 'x': {
      for (int PosY=0; PosY<this->UnitCellY; ++PosY)
	{
	  Complex PhaseX = Phase(-2.0*M_PI*this->FluxDensity*(double)PosY);
	  for (int PosX=0; PosX<this->UnitCellX; ++PosX)
	    {

	      int TmpPosition =  this->EncodeSublatticeIndex(PosX, PosY,NumXTranslations,NumYTranslations,translationPhase);
	      int TmpIndex = 0;
	      OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(PosX+1, PosY,NumXTranslations,NumYTranslations, translationPhase);
	      SpatialIndices[TmpPosition][TmpIndex * 2] = NumXTranslations;
	      SpatialIndices[TmpPosition][(TmpIndex * 2) +1] = NumYTranslations;
	      HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0*Conj(PhaseX);
	      ++TmpIndex;
	      OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(PosX, PosY+1,NumXTranslations,NumYTranslations, translationPhase);
	      SpatialIndices[TmpPosition][TmpIndex * 2] = NumXTranslations;
	      SpatialIndices[TmpPosition][(TmpIndex * 2) +1] = NumYTranslations;
	      HoppingAmplitudes[TmpPosition][TmpIndex] = -1.0;
	      ++TmpIndex;
	      OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(PosX-1, PosY,NumXTranslations,NumYTranslations, translationPhase);
	      SpatialIndices[TmpPosition][TmpIndex * 2] = NumXTranslations;
	      SpatialIndices[TmpPosition][(TmpIndex * 2) +1] = NumYTranslations;
	      HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0*PhaseX;
	      ++TmpIndex;
	      OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(PosX, PosY-1,NumXTranslations,NumYTranslations, translationPhase);
	      SpatialIndices[TmpPosition][TmpIndex * 2] = NumXTranslations;
	      SpatialIndices[TmpPosition][(TmpIndex * 2) +1] = NumYTranslations;
	      HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;

	    }
	}
      break;
    }

    default:
      cout << "Invalid Landau quantization axis encountered in TightBindingModelHofstadterSquare."<<endl;
      exit(1);
      break;
    }

  HermitianMatrix TmpMatrix = this->BuildTightBindingHamiltonianRealSpace(NbrConnectedOrbitals, OrbitalIndices, SpatialIndices, HoppingAmplitudes);
  for (int i = 0; i < this->NbrBands; ++i)
    {
      delete[] HoppingAmplitudes[i];
      delete[] SpatialIndices[i];
      delete[] OrbitalIndices[i];
    }
  delete[] HoppingAmplitudes;
  delete[] SpatialIndices;
  delete[] OrbitalIndices;
  delete[] NbrConnectedOrbitals;
  return TmpMatrix;
}




// build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
//
// nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
// orbitalIndices = array that gives the orbital indices of the connected orbitals
// spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
// hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
// return value = tight binding hamiltonian in real space 

HermitianMatrix  TightBindingModelHofstadterTriangular::BuildTightBindingHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes)
{
  HermitianMatrix TmpHamiltonian(this->NbrBands * this->NbrSiteX * this->NbrSiteY, true);
  ComplexMatrix TmpHamiltonian2(this->NbrBands * this->NbrSiteX * this->NbrSiteY, true);
  int   NumXTranslations;
  int   NumYTranslations;
  Complex  tmpPhase, tmpPhase2;
  Complex TmpXPhase = Complex(1.0,0);
  for (int i = 0; i < this->NbrSiteX; ++i)
    {
      for (int j = 0; j < this->NbrSiteY; ++j)
	{
	  for (int k = 0; k < this->NbrBands; ++k)
	    {
	      int InitialIndex = this->GetRealSpaceTightBindingLinearizedIndex(i, j, k);
	      for (int l = 0; l < nbrConnectedOrbitals[k]; ++l)
		{
		  int FinalIndex = this->GetRealSpaceTightBindingLinearizedIndexSafe(spatialIndices[k][l << 1] + i, spatialIndices[k][(l << 1) + 1] + j, orbitalIndices[k][l],NumXTranslations,NumYTranslations);                 
		  if (InitialIndex>=FinalIndex)
		    {
		      tmpPhase = 1.0;
		      int Tmp = orbitalIndices[k][l] - k;
		      if( ( (orbitalIndices[k][l]%this->UnitCellX - k%this->UnitCellX) ==0  ) && (spatialIndices[k][l << 1]==0 ) )
			{
			  if( spatialIndices[k][(l << 1) + 1] >= 0)
			    for (int p=0; p < i; p++)
			      tmpPhase*=Conj(this->LxTranslationPhase);
			  else
			    {
			      for (int p=0; p < i; p++)
				tmpPhase*=this->LxTranslationPhase;
			    }
			}

		      /*                  if( ( (orbitalIndices[k][l]%this->UnitCellX - k%this->UnitCellX) ==0  ) && (spatialIndices[k][l << 1]==0 ) )
					  {
					  if( spatialIndices[k][(l << 1) + 1] >= 0)
					  for (int p=0; p < i; p++)
					  tmpPhase*=this->LxTranslationPhase;
					  else
					  {
					  for (int p=0; p < i; p++)
					  tmpPhase*=Conj(this->LxTranslationPhase);
					  }
					  }

		      */

		      if( NumXTranslations>0)
			tmpPhase*= Phase(-2.0*M_PI*this->FluxDensity*this->NbrSiteX*this->UnitCellX*(orbitalIndices[k][l]/this->UnitCellX));
		      if( NumXTranslations<0)
			tmpPhase*= Phase(2.0*M_PI*this->FluxDensity*this->NbrSiteX*this->UnitCellX*(orbitalIndices[k][l]/this->UnitCellX));

      
		      int TmpX = spatialIndices[k][l << 1] + i;
		      int TmpY = spatialIndices[k][(l << 1)+1] + j;
		      TmpHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, hoppingAmplitudes[k][l]*tmpPhase);


#ifdef DEBUG_OUTPUT
		      cout <<"x = " <<i<< " y = " <<j <<" k = " <<k<< " going to X = " <<  TmpX  << " Y = "<<TmpY<<" index "<< orbitalIndices[k][l]<<" Coefficient" << hoppingAmplitudes[k][l]*tmpPhase<<"NumTranslation X = "<<NumXTranslations <<endl;
#endif
		    }
		}
	    }
	}      
    }
  return TmpHamiltonian;
}

