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
#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
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
// nbrLayers = number of layers to be stacked (independent)
// nbrFlux = number of flux quanta per unit cell
// axis = direction of Landau gauge within cell ('x' or 'y')
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
TightBindingModelHofstadterSquare::TightBindingModelHofstadterSquare(int nbrCellX, int nbrCellY, int unitCellX, int unitCellY, int nbrLayers, int nbrFlux, char axis,
								     double gammaX, double gammaY, 
								     AbstractArchitecture* architecture, bool storeOneBodyMatrices, bool useEmbedding)
{
  this->NbrSiteX = nbrCellX;
  this->NbrSiteY = nbrCellY;
  this->UnitCellX = unitCellX;
  this->UnitCellY = unitCellY;
  this->NbrLayers = nbrLayers;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = this->UnitCellX * this->UnitCellY * this->NbrLayers;
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

TightBindingModelHofstadterSquare::~TightBindingModelHofstadterSquare()
{
}


// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelHofstadterSquare::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  double K1;
  double K2;
  //int sublXf, sublYf; // final site sublattice positions

  // parameters for number of required images of simulation cell
  int images = 50; // number of images of the simulation cell
  int maxRange = std::sqrt(23./(1.0-this->FluxDensity)); // larger exponents will yield numerical zero
  if (this->Range > maxRange)
    {
      cout << "Testing: range not reduced to new maxRange="<<maxRange<<endl;
      // this->Range = maxRange;
    }
  images = this->Range / (UnitCellX < UnitCellY? UnitCellX:UnitCellY) + 2;
  Complex amplitude;
  //loop over initial sites (i,j)
  
  int *qi = new int[this->NbrBands];
  int *qf = new int[this->NbrBands];

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
		      for (int j=0; j<UnitCellY; ++j)
			{
			  for (int s=0; s<this->NbrLayers; ++s)
			    qi[s] = this->EncodeQuantumNumber(i, j, s, 0.0, 0.0, TranslationPhase);
	      	      
			  // have long-range hopping, so sum over all possible final sites (k,l)
			  for (int k=0; k<UnitCellX; ++k)
			    {
			      for (int l=0; l<UnitCellY; ++l)
				{		  
				  for (int s=0; s<this->NbrLayers; ++s)
				    qf[s] = Particles->EncodeQuantumNumber(k, l, s, K1, K2, TranslationPhase);
		  		  
				  ComplexMatrix sumHopping(this->NbrLayers, this->NbrLayers, true);
				  for (int dX = -images; dX <= images; ++dX)
				    for (int dY = -images; dY <= images; ++dY)
				      {
					if ( this->GetDistance(k+dX*UnitCellX - i, l+dY*UnitCellY - j) <= this->Range // allow hard cut-off
					     // && ( ( dX!=0 || dY!=0 ) || ( i!=k || j!=l ) ) // optionally, exclude onsite terms (needed to get zero energy for lowest band)
					     )
					  {
					    Particles->EncodeQuantumNumber(k+dX*UnitCellX, l+dY*UnitCellY, 0, 0.0, 0.0, TranslationPhase); // just need to get new TranslationPhase, here
					    amplitude = -HoppingSign*this->KapitMuellerHopping(k+dX*UnitCellX, l+dY*UnitCellY, i, j) * Conj(TranslationPhase);

					    if (Norm(amplitude) > 1e-13)
					      {
						// check for any branch cuts being crossed
						int shiftModulo = this->EvaluateBranchCrossings(i, j, k+dX*UnitCellX, l+dY*UnitCellY);
						while (shiftModulo<0) shiftModulo+=this->NbrLayers;
						shiftModulo %= this->NbrLayers;		    

						for (int si=0; si<this->NbrLayers; ++si)
						  {
						    int sf = (si+shiftModulo);
						    while (sf<0) sf+=this->NbrLayers;
						    sf %= this->NbrLayers;		    
						    sumHopping.AddToMatrixElement(si, sf, amplitude);
						  }
#ifdef DEBUG_OUTPUT
						//if (TranslationPhase!=1.0)
						cout << "xx ("<<i<<", "<<j<<")->("<<k+dX*UnitCellX<<", "<<l+dY*UnitCellY<<") amplitude="<<amplitude<<", shift="<< shiftModulo <<", dL=("<<dX<<"," <<dY<<") : "
						     <<"Translation ["<<qi[0]<<"->"<<qf[0]<<"]="
						     << TranslationPhase << endl;
#endif
					      }
					  }
				      }
#ifdef DEBUG_OUTPUT		     
				  //if (Norm(sumHopping[0][0])>1e-15)
				  {
				    cout << "H[("<<i<<", "<<j<<")->("<<k<<", "<<l<<")]= ";
				    for (int si=0; si<this->NbrLayers; ++si)
				      for (int sf=0; sf<this->NbrLayers; ++sf)
					cout << sumHopping[si][sf]<<", ";
				    cout<<" tP="<<TranslationPhase<<endl;
				  }
#endif
				  for (int si=0; si<this->NbrLayers; ++si)
				    for (int sf=0; sf<this->NbrLayers; ++sf)
				      {
					// only take into account terms with non-zero magnitude
					sumHopping.GetMatrixElement(si, sf, amplitude);
					if (Norm(amplitude)>1e-15)
					  {			    
					    KineticQi[TmpNumberTerms] = qi[si];
					    KineticQf[TmpNumberTerms] = qf[sf];
					    HoppingTerms[TmpNumberTerms] = amplitude;
					    ++TmpNumberTerms;
					  }
				      }
				}
			    }
			}
		    }
		  delete [] qi;
		  delete [] qf;

		  break;
		}

		case 'x': { // needs checking
		  cout << "Invalid Landau quantization axis x not yet defined in TightBindingModelKapitMueller."<<endl;
		  exit(1);
		  break;
		}

		default:
 		  cout << "Invalid Landau quantization axis encountered in TightBindingModelKapitMueller."<<endl;
		  exit(1);
		  break;
		}

#ifdef DEBUG_OUTPUT
	      if(this->NbrBands<25)
		cout << TmpOneBodyHamiltonian<< endl;
#endif

	      // try to shift the Hamiltonian
	      double shift = 20.0;
	      for (int i=0; i<this->NbrBands; ++i)
		TmpOneBodyHamiltonian.AddToMatrixElement(i,i,shift);
	    

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
		    this->EnergyBandStructure[i][Index] = TmpDiag(i, i)-shift;

#ifdef DEBUG_OUTPUT
		  char buffer[255];
		  sprintf(buffer, "sp_eigenvector_kx_%d_ky_%d.txt", kx, ky);
		  ofstream EVFile(buffer, std::ios::out);
		  EVFile << "#x y layer psi abs(psi) arg(psi)*q/pi\n";
		  for (int j=0; j<UnitCellY; ++j)
		    for (int i=0; i<UnitCellX; ++i)
		      for (int l=0; i<NbrLayers; ++l)
			{
			  int Origin = this->EncodeSublatticeIndex(0, 0, l, K1, K2, TranslationPhase);
			  int SublIndex = this->EncodeSublatticeIndex(i, j, l, K1, K2, TranslationPhase);
			  Complex Gauge=Polar(1.0,-Arg(TmpMatrix[0][0]));
			  TmpMatrix[0]*=Gauge;
			  EVFile << i <<" "<<j<<" "<<l<<" "<<TmpMatrix[0][SublIndex]<<" "<<Norm(TmpMatrix[0][SublIndex])<<" "<<Arg(TmpMatrix[0][SublIndex])*this->NbrBands/(M_PI)<<"\n";
			}
		  EVFile.close();
#endif
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
		    this->EnergyBandStructure[i][Index] = TmpDiag(i, i) - shift;
		}
	    }
	}
    }
}




// set natural embedding, i.e., at positions of a uniform lattice
//
void TightBindingModelHofstadterSquare::SetNaturalEmbedding()
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
void TightBindingModelHofstadterSquare::SetNbrFluxQuanta(int nbrFluxQuanta)
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
      this->LxTranslationPhase = Polar(1.0, -2.0*M_PI*FluxDensity*this->UnitCellX);
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
int TightBindingModelHofstadterSquare::EncodeSublatticeIndex(int posx, int posy, int layer, double KX, double KY, Complex &translationPhase)
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
  int rst = (posx + this->UnitCellX*posy)*this->NbrLayers + layer;
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

ComplexMatrix TightBindingModelHofstadterSquare::GetRealSpaceTightBindingEigenstates()
{
  double LogTranslationPhaseX= -2.0*M_PI*this->FluxDensity*this->UnitCellX;
  ComplexMatrix EigenStates(this->NbrBands *  this->NbrStatePerBand,this->NbrBands *  this->NbrStatePerBand ,true);
  int Kx;  int Ky;
  int K1;  int K2;
  int OrbitalIndex = 0;
  int PosXUnitCell = 0;
  int PosYUnitCell = 0;
  int PosX=0, PosY=0;
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

	  this->DecodeSublatticeIndex(OrbitalIndex, PosX, PosY);

	  // note: looks suspicious: OrbitalIndex/this->UnitCellX;
	  int TotalPosY = PosYUnitCell*this->UnitCellY + OrbitalIndex/this->UnitCellX;

	  if (PosY != OrbitalIndex/this->UnitCellX)
	    std::cerr<<"Discrepancy of y position in GetRealSpaceTightBindingEigenstates"<<endl;

	  EigenStates[i][j] = this->OneBodyBasis[MomentumIndex][BandNumber][OrbitalIndex] * Phase(K1*PosXUnitCell + K2*PosYUnitCell);
	    // * Phase(PosXUnitCell* LogTranslationPhaseX*TotalPosY); // in Landau-gauge, there should be no extra phase!
	}
    }
  return EigenStates;
}


// find the orbitals connected to those located at the origin unit cell
// 
void TightBindingModelHofstadterSquare::FindConnectedOrbitals()
{
  int p, q;
  int numXTranslations, numYTranslations;
  Complex translationPhase; 
	      
#ifdef DEBUG_OUTPUT
  cout << "FindConnectedOrbitals()"<<endl;
#endif


  if (this->NbrConnectedOrbitals == 0)
    {
      this->NbrConnectedOrbitals = new int [this->NbrBands];
      this->ConnectedOrbitalIndices = new int* [this->NbrBands];
      this->ConnectedOrbitalSpatialIndices = new int* [this->NbrBands];
      this->ConnectedOrbitalHoppingAmplitudes = new Complex* [this->NbrBands];
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  this->NbrConnectedOrbitals[i] = 4;  // four nearest neighbors on square lattice 
	  this->ConnectedOrbitalIndices[i] = new int[this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalSpatialIndices[i] = new int[2 * this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalHoppingAmplitudes[i] = new Complex[this->NbrConnectedOrbitals[i]];
	}
      
      int TmpIndex = 0;
      
      Complex TranslationPhase;
      switch (this->LandauGaugeAxis)
	{
	case 'y': {
	  for (int i=0; i<UnitCellX; ++i)
	    {
	      Complex ABPhase=Polar(1.0,2.0*M_PI*this->FluxDensity*(double)i);
	      for (int j=0; j<UnitCellY; ++j)
		{
		  TmpIndex=0;
		  int InitialIndex = this->EncodeSublatticeIndex(i, j, 0.0, 0.0, TranslationPhase); // TranlationPhase always one, so can be discarded
		  
		  // get final Sublattice Index
		  int FinalIndex = this->EncodeSublatticeIndex(i+1, j, numXTranslations, numYTranslations, TranslationPhase);
		  this->GetRealSpaceIndex(-numXTranslations, -numYTranslations, p , q); // re-map number of translations of magnetic unit cells
		  this->ConnectedOrbitalIndices[InitialIndex][TmpIndex] = FinalIndex;
		  this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2] = p;
		  this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2) +1] = q;
		  this->ConnectedOrbitalHoppingAmplitudes[InitialIndex][TmpIndex] = -TranslationPhase * Phase(-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2]) + EmbeddingX[FinalIndex] - EmbeddingX[InitialIndex]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2) +1]) + EmbeddingY[FinalIndex] - EmbeddingY[InitialIndex]));
#ifdef DEBUG_OUTPUT
		  cout <<"subi="<<InitialIndex<< "->subf=" << this->ConnectedOrbitalIndices[InitialIndex][TmpIndex]<< " dRx=" << this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2]<< " dRY=" <<this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2)+1]<<" tij="<<this->ConnectedOrbitalHoppingAmplitudes[InitialIndex][TmpIndex]<<endl;
#endif
		  ++TmpIndex;

		  FinalIndex = this->EncodeSublatticeIndex(i-1, j, numXTranslations, numYTranslations, translationPhase);
		  this->GetRealSpaceIndex(-numXTranslations, -numYTranslations, p , q); // number of translations of magnetic unit cells
		  this->ConnectedOrbitalIndices[InitialIndex][TmpIndex] = FinalIndex;
		  this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2] = p;
		  this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2) +1] = q;
		  this->ConnectedOrbitalHoppingAmplitudes[InitialIndex][TmpIndex] = -TranslationPhase* Phase(-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2]) + EmbeddingX[FinalIndex] - EmbeddingX[InitialIndex]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2) +1]) + EmbeddingY[FinalIndex] - EmbeddingY[InitialIndex]));
#ifdef DEBUG_OUTPUT
		  cout <<"subi="<<InitialIndex<< "->subf=" << this->ConnectedOrbitalIndices[InitialIndex][TmpIndex]<< " dRx=" << this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2]<< " dRY=" <<this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2)+1]<<" tij="<<this->ConnectedOrbitalHoppingAmplitudes[InitialIndex][TmpIndex]<<endl;
#endif
		  ++TmpIndex;

		  FinalIndex = this->EncodeSublatticeIndex(i, j+1, numXTranslations, numYTranslations, translationPhase);
		  this->GetRealSpaceIndex(-numXTranslations, -numYTranslations, p , q); // number of translations of magnetic unit cells
		  this->ConnectedOrbitalIndices[InitialIndex][TmpIndex] = FinalIndex;
		  this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2] = p;
		  this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2) +1] = q;
		  this->ConnectedOrbitalHoppingAmplitudes[InitialIndex][TmpIndex] = -TranslationPhase*Conj(ABPhase)* Phase(-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2]) + EmbeddingX[FinalIndex] - EmbeddingX[InitialIndex]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2) +1]) + EmbeddingY[FinalIndex] - EmbeddingY[InitialIndex]));
#ifdef DEBUG_OUTPUT
		  cout <<"subi="<<InitialIndex<< "->subf=" << this->ConnectedOrbitalIndices[InitialIndex][TmpIndex]<< " dRx=" << this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2]<< " dRY=" <<this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2)+1]<<" tij="<<this->ConnectedOrbitalHoppingAmplitudes[InitialIndex][TmpIndex]<<endl;
#endif
		  ++TmpIndex;

		  FinalIndex = this->EncodeSublatticeIndex(i, j-1, numXTranslations, numYTranslations, translationPhase);
		  this->GetRealSpaceIndex(-numXTranslations, -numYTranslations, p , q); // number of translations of magnetic unit cells
		  this->ConnectedOrbitalIndices[InitialIndex][TmpIndex] = FinalIndex;
		  this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2] = p;
		  this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2) +1] = q;
		  this->ConnectedOrbitalHoppingAmplitudes[InitialIndex][TmpIndex] = -TranslationPhase*ABPhase* Phase(-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2]) + EmbeddingX[FinalIndex] - EmbeddingX[InitialIndex]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2) +1]) + EmbeddingY[FinalIndex] - EmbeddingY[InitialIndex]));
#ifdef DEBUG_OUTPUT
		  cout <<"subi="<<InitialIndex<< "->subf=" << this->ConnectedOrbitalIndices[InitialIndex][TmpIndex]<< " dRx=" << this->ConnectedOrbitalSpatialIndices[InitialIndex][TmpIndex * 2]<< " dRY=" <<this->ConnectedOrbitalSpatialIndices[InitialIndex][(TmpIndex * 2)+1]<<" tij="<<this->ConnectedOrbitalHoppingAmplitudes[InitialIndex][TmpIndex]<<endl;
#endif
		  ++TmpIndex;
		}
	    }
	  break;
	}


	default:
	  cout << "Invalid Landau quantization axis encountered in ParticleOnLatticeDeltaHamiltonian."<<endl;
	  exit(1);
	  break;
	}
    }      
}


// get the tight binding hamiltonian in real space 
// 
// return value = tight binding hamiltonian

HermitianMatrix TightBindingModelHofstadterSquare::GetRealSpaceTightBindingHamiltonian2()
{

  cout << "GetRealSpaceTightBindingHamiltonian2"<<endl;
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

HermitianMatrix  TightBindingModelHofstadterSquare::BuildTightBindingHamiltonianRealSpace2(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes)
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



// compute the description of the density-density interaction for the unit cell at the origin
//
// nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsPotentials = intensity of each density-density term 
// bosonFlag = true if we are dealing with bosons
// uPotential = nearest neighbor (for fermions) or on-site (for bosons) interaction amplitude
// vPotential = next nearest neighbor (for fermions) or nearest neighbor (for bosons) interaction amplitude

void TightBindingModelHofstadterSquare::ComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
							   int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
							   bool bosonFlag, double uPotential, double vPotential)
{
  nbrInteractingOrbitals = new int[this->GetNbrBands()];
  interactingOrbitalsOrbitalIndices = new int*[this->GetNbrBands()];
  interactingOrbitalsSpatialIndices = new int*[this->GetNbrBands()];
  interactingOrbitalsPotentials = new double*[this->GetNbrBands()];
  int p, q;
  int numXTranslations, numYTranslations;
  Complex TranslationPhase;
  if (bosonFlag == false)
    {

      // bosons
      // allocate maximum number of interacting orbitals per site
      for (int s=0; s<this->GetNbrBands(); ++s)       
	nbrInteractingOrbitals[s] = 4; // four NN interactions
      if (vPotential != 0.0)
	for (int s=0; s<this->GetNbrBands(); ++s)       
	  nbrInteractingOrbitals[s] += 4; // four additional terms for NNN interactions
      for (int s=0; s<this->GetNbrBands(); ++s)
	{
	  interactingOrbitalsOrbitalIndices[s] = new int[nbrInteractingOrbitals[s]];
	  interactingOrbitalsSpatialIndices[s] = new int[nbrInteractingOrbitals[s] * 2];
	  interactingOrbitalsPotentials[s] = new double[nbrInteractingOrbitals[s]];

	  nbrInteractingOrbitals[s] = 0;
      
	  // define NN interactions
	  if (uPotential != 0.0)
	    {
	      int i,j;
	      this->DecodeSublatticeIndex(s, i, j);

	      int nbrV = 4;
	      int dX[4] = {1,-1,0,0};
	      int dY[4] = {0,0,1,-1};
	      
	      for (int nn=0; nn<nbrV; ++nn)
		{
		  int s2=this->EncodeSublatticeIndex(i+dX[nn], j+dY[nn], numXTranslations, numYTranslations, TranslationPhase);
		  if (s2>=s)
		    {
		      this->GetRealSpaceIndex(-numXTranslations, -numYTranslations, p , q);
	      
		      interactingOrbitalsOrbitalIndices[s][nbrInteractingOrbitals[s]] = s2;
		      interactingOrbitalsSpatialIndices[s][2 * nbrInteractingOrbitals[s]] = p;
		      interactingOrbitalsSpatialIndices[s][(2 * nbrInteractingOrbitals[s]) + 1] = q;
		      if (s==s2)
			interactingOrbitalsPotentials[s][nbrInteractingOrbitals[s]] = 0.5*uPotential;
		      else
			interactingOrbitalsPotentials[s][nbrInteractingOrbitals[s]] = uPotential;
		      ++nbrInteractingOrbitals[s];
		    }
		}
	    }


	  if (vPotential != 0.0)
	    {
	      int i,j;
	      this->DecodeSublatticeIndex(s, i, j);

	      int nbrV = 4;
	      int dX[4] = {1,-1,0,0};
	      int dY[4] = {0,0,1,-1};
	      
	      for (int nn=0; nn<nbrV; ++nn)
		{
		  int s2=this->EncodeSublatticeIndex(i+dX[nn], j+dY[nn], numXTranslations, numYTranslations, TranslationPhase);
		  if (s2>=s)
		    {
		      this->GetRealSpaceIndex(-numXTranslations, -numYTranslations, p , q);
	      
		      interactingOrbitalsOrbitalIndices[s][nbrInteractingOrbitals[s]] = s2;
		      interactingOrbitalsSpatialIndices[s][2 * nbrInteractingOrbitals[s]] = p;
		      interactingOrbitalsSpatialIndices[s][(2 * nbrInteractingOrbitals[s]) + 1] = q;
		      if (s==s2)
			interactingOrbitalsPotentials[s][nbrInteractingOrbitals[s]] = 0.5*vPotential;
		      else
			interactingOrbitalsPotentials[s][nbrInteractingOrbitals[s]] = vPotential;
		      ++nbrInteractingOrbitals[s];
		    }
		}
	    }
	}
    }
  else
    { // bosons
      // allocate maximum number of interacting orbitals per site
      for (int s=0; s<this->GetNbrBands(); ++s)       
	nbrInteractingOrbitals[s] = 1; // one interaction for onsite term
      if (vPotential != 0.0)
	for (int s=0; s<this->GetNbrBands(); ++s)       
	  nbrInteractingOrbitals[s] += 4; // four terms for NN interactions
      for (int s=0; s<this->GetNbrBands(); ++s)       
	{
	  interactingOrbitalsOrbitalIndices[s] = new int[nbrInteractingOrbitals[s]];
	  interactingOrbitalsSpatialIndices[s] = new int[nbrInteractingOrbitals[s] * 2];
	  interactingOrbitalsPotentials[s] = new double[nbrInteractingOrbitals[s]];

	  nbrInteractingOrbitals[s] = 0;
      
	  // define onsite interactions
	  interactingOrbitalsOrbitalIndices[s][nbrInteractingOrbitals[s]] = s;
	  this->GetRealSpaceIndex(0, 0, p, q);
	  interactingOrbitalsSpatialIndices[s][2 * nbrInteractingOrbitals[s]] = p;
	  interactingOrbitalsSpatialIndices[s][(2 * nbrInteractingOrbitals[s]) + 1] = q;
	  interactingOrbitalsPotentials[s][nbrInteractingOrbitals[s]] = 0.5*uPotential;
	  ++nbrInteractingOrbitals[s];

	  if (vPotential != 0.0)
	    {
	      int i,j;
	      this->DecodeSublatticeIndex(s, i, j);
	      
	      int nbrV = 4;
	      int dX[4] = {1,-1,0,0};
	      int dY[4] = {0,0,1,-1};
	      
	      for (int nn=0; nn<nbrV; ++nn)
		{
		  int s2=this->EncodeSublatticeIndex(i+dX[nn], j+dY[nn], numXTranslations, numYTranslations, TranslationPhase);
		  if (s2>=s)
		    {
		      this->GetRealSpaceIndex(-numXTranslations, -numYTranslations, p , q);
	      
		      interactingOrbitalsOrbitalIndices[s][nbrInteractingOrbitals[s]] = s2;
		      interactingOrbitalsSpatialIndices[s][2 * nbrInteractingOrbitals[s]] = p;
		      interactingOrbitalsSpatialIndices[s][(2 * nbrInteractingOrbitals[s]) + 1] = q;
		      if (s==s2)
			interactingOrbitalsPotentials[s][nbrInteractingOrbitals[s]] = 0.5*vPotential;
		      else
			interactingOrbitalsPotentials[s][nbrInteractingOrbitals[s]] = vPotential;
		      ++nbrInteractingOrbitals[s];
		    }
		}
	    }
	}
    }
}

/// evaluate branch crossings for a link from site i to j (from ParticleOnLatticeKapitMuellerMultiLayerHamiltonian)
/// @param xi, yi, xj, yj coordinates of the initial (i) and final site (j)
/// @return overall shift for crossing any of the known branch cuts
int TightBindingModelHofstadterSquare::EvaluateBranchCrossings(int xi, int yi, int xj, int yj)
{
  RealVector Ri(2), Rj(2);
  Ri[0]=xi;
  Ri[1]=yi;
  Rj[0]=xj;
  Rj[1]=yj;
  RealVector Rij(Rj, true);
  Rij-=Ri;
  double xiA, xiB, xiJ, xiS, etaA, etaB, etaJ;
  int shift=0;
  for (int i=0; i<NbrBranchCuts; ++i)
    {
      RealVector RA(this->BranchCuts[i].A, true);
      RA -= Ri;
      RealVector RB(this->BranchCuts[i].B, true);
      RB -= Ri;
      xiA=this->BranchCuts[i].E1 * RA;
      etaA=this->BranchCuts[i].E2 * RA;
      xiB=this->BranchCuts[i].E1 * RB;
      etaB=this->BranchCuts[i].E2 * RB;
      xiJ=this->BranchCuts[i].E1 * Rij;
      etaJ=this->BranchCuts[i].E2 * Rij;
      if (fabs(etaJ)>1e-13) // Rij not parallel to cut?
	{
	  xiS = xiJ * etaA / etaJ;
	  double etaS = etaA;
	  // crossing point in interval of branch cut?
	  if ( (xiS >= fmin(xiA,xiB)-1e-13) && (xiS <= fmax(xiA,xiB)+1e-13)
	       && (xiS <= fmax(0.0,xiJ)+1e-13) && (xiS >= fmin(0.0,xiJ)-1e-13)
	       && ( etaS <=  fmax(0.0,etaJ)+1e-13) && (etaS >= fmin(0.0,etaJ)-1e-13) )
	    {
	      shift += this->BranchShift[i] * (1 - 2 * (etaJ < 0));
	    }
	}
    }
  return shift;
}
