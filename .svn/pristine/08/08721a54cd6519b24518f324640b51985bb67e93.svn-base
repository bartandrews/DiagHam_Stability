////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of abstract tight binding model                    //
//                                                                            //
//                        last modification : 01/05/2012                      //
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
#include "Tools/FTITightBinding/AbstractTightBindingModel.h"
#include "Architecture/ArchitectureOperation/FTIComputeBandStructureOperation.cc"
#include "GeneralTools/Endian.h"
#include "GeneralTools/OrderedList.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include <sys/time.h>
#include <fstream>


using std::ofstream;
using std::ios;
using std::endl;

// constructor
//

AbstractTightBindingModel::AbstractTightBindingModel()
{
  this->Architecture = 0;
  this->NbrConnectedOrbitals = 0;
  this->ConnectedOrbitalIndices = 0;
  this->ConnectedOrbitalSpatialIndices = 0;
  this->ConnectedOrbitalHoppingAmplitudes = 0;
  this->EnergyBandStructure = 0;
  this->OneBodyBasis = 0;
}


// destructor
//

AbstractTightBindingModel::~AbstractTightBindingModel()
{
  if (this->NbrConnectedOrbitals != 0)
    {
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  delete[] this->ConnectedOrbitalIndices[i];
	  delete[] this->ConnectedOrbitalSpatialIndices[i];
	  delete[] this->ConnectedOrbitalHoppingAmplitudes[i];
	}
      delete[] this->ConnectedOrbitalIndices;
      delete[] this->ConnectedOrbitalSpatialIndices;
      delete[] this->ConnectedOrbitalHoppingAmplitudes;
      delete[] this->NbrConnectedOrbitals;
    }
  if (this->OneBodyBasis != 0)
    {
      delete[] this->OneBodyBasis;
    }
  if (this->EnergyBandStructure != 0)
    {
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  delete[] this->EnergyBandStructure[i];
	}
      delete[] this->EnergyBandStructure;
    }
}

// write an ASCII header that describes the tight binding model
// 
// output = reference on the output stream
// commentChar = optional ASCII character to add in front of each header line
// return value  = reference on the output stream

ostream& AbstractTightBindingModel::WriteASCIIHeader(ostream& output, char commentChar)
{
  return output;
}

// write an header that describes the tight binding model
// 
// output = reference on the output stream
// return value  = reference on the output stream

ofstream& AbstractTightBindingModel::WriteHeader(ofstream& output)
{
  int HeaderSize = 0;
  WriteLittleEndian(output, HeaderSize);
  return output; 
}

// read the header that describes the tight binding model
// 
// return value  = size of header that was read.
int AbstractTightBindingModel::ReadHeader(ifstream& input)
{
  int HeaderSize = -1;
  ReadLittleEndian(input, HeaderSize);
  return HeaderSize;
}


// write the eigenvalues and eigenstates (if available) to a band structure file
// output = reference on the output stream
// return value  = reference on the output stream
ostream& AbstractTightBindingModel::WriteEigensystem(ofstream& output)
{
  for (int idx = 0; idx < this->NbrStatePerBand; ++idx)
    for (int i = 0; i < this->NbrBands; ++i)
      {
	double Tmp = this->GetEnergy(i, idx);
	WriteLittleEndian(output, Tmp);
      }
  if (this->HaveOneBodyBasis() == true)
    {
      for (int idx = 0; idx < this->NbrStatePerBand; ++idx)
      {
//           cout << idx <<  " "  << endl;
	this->GetOneBodyMatrix(idx).WriteMatrix(output);
      }
    }
  return output;
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool AbstractTightBindingModel::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  File.precision(14);
  this->WriteASCIIHeader(File, '#');
  File << "# index";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  File << endl;
  for (int kx = 0; kx < this->NbrStatePerBand; ++kx)
    {
      File << kx; 
      for (int i = 0; i < this->NbrBands; ++i)
	File << " " << this->GetEnergy(i, kx);
      File << endl;
    }
  File.close();
  return true;
}

// compute the spread of a given band
// 
// band = band index
// return value = band spread 

double AbstractTightBindingModel::ComputeBandSpread(int band)
{
  if ((band < 0) || (band > this->NbrBands))
    {
      return -1.0;
    }
  double MinBand = this->GetEnergy(band, 0);
  double MaxBand = this->GetEnergy(band, 0);
  for (int i = 1; i < this->NbrStatePerBand; ++i)
    {
      if (this->GetEnergy(band, i) > MaxBand)
	{
	  MaxBand = this->GetEnergy(band, i);
	}
      else
	{
	  if (this->GetEnergy(band, i) < MinBand)
	    {
	      MinBand = this->GetEnergy(band, i);
	    }
	}
    }
  return (MaxBand - MinBand);
}

// compute the direct gap between two bands
// 
// band1 = index of the lower band 
// band2 = index of the upper band (-1 if it has to be set to band1 + 1)
// return value =  direct band gap

double AbstractTightBindingModel::ComputeDirectBandGap(int band1, int band2)
{
  if ((band1 < 0) || (band1 > this->NbrBands))
    {
      return -1.0;
    }
  if (band2 < 0)
    band2 = band1 + 1;
  if ((band2 < 0) || (band2 > this->NbrBands) || (band2 <= band1))
    {
      return -1.0;
    }
  double Gap = this->GetEnergy(band2, 0) - this->GetEnergy(band1, 0);
  for (int i = 1; i < this->NbrStatePerBand; ++i)
    {
      double TmpGap = this->GetEnergy(band2, i) - this->GetEnergy(band1, i);
      if (TmpGap < Gap)
	Gap = TmpGap;
    }
  return Gap;
}

// compute the direct gap between two bands
// 
// band1 = index of the lower band 
// band2 = index of the upper band (-1 if it has to be set to band1 + 1)
// return value =  direct band gap

double AbstractTightBindingModel::ComputeIndirectBandGap(int band1, int band2)
{
  if ((band1 < 0) || (band1 > this->NbrBands))
    {
      return -1.0;
    }
  if (band2 < 0)
    band2 = band1 + 1;
  if ((band2 < 0) || (band2 > this->NbrBands) || (band2 <= band1))
    {
      return -1.0;
    }
  double Max1=this->GetEnergy(band1, 0);
  double Min2=this->GetEnergy(band2, 0);
  for (int i = 1; i < this->NbrStatePerBand; ++i)
    {
      double Tmp1 = this->GetEnergy(band1, i);
      double Tmp2 = this->GetEnergy(band2, i);
      if (Tmp1 > Max1)
	Max1 = Tmp1;
      if (Tmp2 < Min2)
	Min2 = Tmp2;
    }
  return Min2-Max1;
}


class SingleParticleStateData
{
public:
  SingleParticleStateData();
  SingleParticleStateData(double e, double b, double i);
  ~SingleParticleStateData(){};
  double Energy;
  double BandIndex;
  double StateIndex;

  friend bool operator<(const SingleParticleStateData &d1, const SingleParticleStateData &d2);
  friend bool operator>(const SingleParticleStateData &d1, const SingleParticleStateData &d2);
  friend bool operator==(const SingleParticleStateData &d1, const SingleParticleStateData &d2);

  friend ostream& operator<<(ostream &str, const SingleParticleStateData &d);
  

};

SingleParticleStateData::SingleParticleStateData()
{
  this->Energy = 0.0;
  this->BandIndex = -1;
  this->StateIndex = -1;
}

SingleParticleStateData::SingleParticleStateData(double e, double b, double i)
{
  this->Energy = e;
  this->BandIndex = b;
  this->StateIndex = i;
}

bool operator<(const SingleParticleStateData &d1, const  SingleParticleStateData &d2)
{
  return (d1.Energy < d2.Energy);
}

bool operator>(const SingleParticleStateData &d1, const SingleParticleStateData &d2)
{
  return (d1.Energy > d2.Energy);
}

bool operator==(const SingleParticleStateData &d1, const SingleParticleStateData &d2)
{
  return (d1.Energy == d2.Energy);
}

ostream &operator<<(ostream &str, const SingleParticleStateData &d)
{
  str <<d.BandIndex<<"_"<<d.StateIndex<<" ("<<d.Energy<<")";
  return str;
}



// compute the ground state energy for a number of fermions filling the band 
// 
// nbrFermions = number of particles in the band structure
// bands = number of bands used in groundstate configuration
// return value =  total groundstate energy
double AbstractTightBindingModel::ComputeGroundstateEnergy(int nbrFermions, int &bands, bool verbose)
{
  bands = -1;
  if (nbrFermions<1) return 0.0;

  if (nbrFermions>this->NbrBands*this->NbrStatePerBand) return 0.0;
  OrderedList<SingleParticleStateData> AllStates(false); // ordered list, which does not eliminate duplicates

  for (int b=0; b<this->NbrBands; ++b)
    for (int i=0; i<this->NbrStatePerBand; ++i)
      AllStates += SingleParticleStateData(this->GetEnergy(b, i),b,i);

  bands = 0;
  double Energy=0.0;

  if ((verbose)&&(false))
    {
      cout << "Total number of states: "<<AllStates.GetNbrElement()<<endl;
      for (int n=0; n<AllStates.GetNbrElement(); ++n)
	{
	  SingleParticleStateData D=AllStates[n];
	  cout << "State "<<n<<": "<<D.BandIndex<<"_"<<D.StateIndex<<" (E="<<D.Energy<<")\n";
	}
    }

  for (int n=0; n<nbrFermions; ++n)
    {
      SingleParticleStateData D=AllStates[n];
      if (verbose) cout << n+1 <<"th occupied state: "<<D.BandIndex<<"_"<<D.StateIndex<<" (E="<<D.Energy<<")\n";
      Energy+=D.Energy;
      if (D.BandIndex>bands)
	bands=D.BandIndex;
    }
  
  ++bands;
  return Energy;
}


// return the energy of the lowest energy single-particle state
// 
double AbstractTightBindingModel::SingleParticleGroundstateEnergy()
{
  double Minimum = this->GetEnergy(0,0);
  for (int b=0; b<this->NbrBands; ++b)
    for (int i=0; i<this->NbrStatePerBand; ++i)
      if (this->GetEnergy(b, i) < Minimum)
	Minimum = this->GetEnergy(b, i);
  return Minimum;
}


  
  


// write the full band structure information in a binary file
//
// fileName = name of the output file 
// return value = true if no error occured  

bool AbstractTightBindingModel::WriteBandStructure(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  WriteLittleEndian(File, this->NbrBands);
  WriteLittleEndian(File, this->NbrStatePerBand);
  this->WriteHeader(File);
  this->WriteEigensystem(File);
  File.close();
  return true;
}

// write the full band structure information in an ASCII file
//
// fileName = name of the output file 
// return value = true if no error occured  

bool AbstractTightBindingModel::WriteBandStructureASCII(char* fileName)
{
  ofstream File;
  this->WriteASCIIHeader(File, '#');
  File << "# index";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  for (int i = 0; i < this->NbrBands; ++i)
    for (int j = 0; j < this->NbrBands; ++j)
      File <<  "    U_{" << i << ", " << j << "}";
  File << endl;
  Complex Tmp;
  for (int LinearizedMomentumIndex = 0; LinearizedMomentumIndex < this->NbrStatePerBand; ++LinearizedMomentumIndex)
    {
      File << LinearizedMomentumIndex; 
      for (int i = 0; i < this->NbrBands; ++i)
	File << " " <<  this->GetEnergy(i, LinearizedMomentumIndex);
      for (int i = 0; i < this->NbrBands; ++i)
	for (int j = 0; j < this->NbrBands; ++j)
	  {
	    this->GetOneBodyMatrix(LinearizedMomentumIndex).GetMatrixElement(i, j, Tmp);
	    File <<  "    " << Tmp;
	  }
      File << endl;
    }

  File.close();
  return true;
}

// compute the band structure
//

void AbstractTightBindingModel::ComputeBandStructure()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);
  if (this->Architecture == 0)
    {
      this->CoreComputeBandStructure(0, this->GetNbrStatePerBand());
      return;
    }
  else
    {
      FTIComputeBandStructureOperation Operation (this);
      Operation.ApplyOperation(this->Architecture);
    }
  timeval TotalEndingTime;
  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout << "One-body diagonalization done in " << Dt << " s" << endl;
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void AbstractTightBindingModel::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  this->FindConnectedOrbitals();
  if (this->NbrConnectedOrbitals != 0)
    {
      if (nbrStates == 0l)
	nbrStates = this->NbrStatePerBand;

      HermitianMatrix TmpOneBodyHamiltonian = this->BuildTightBindingHamiltonianRealSpace(this->NbrConnectedOrbitals, this->ConnectedOrbitalIndices,
											  this->ConnectedOrbitalSpatialIndices, this->ConnectedOrbitalHoppingAmplitudes);
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
	  this->OneBodyBasis[0] = TmpMatrix;
	  for (int i = 0; i < this->NbrBands; ++i)
	    {
	      this->EnergyBandStructure[i][0] = TmpDiag(i, i);
	    }
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
	    this->EnergyBandStructure[i][0] = TmpDiag(i, i);
	}
    }
  else
    {
      cout << "error, using dummy AbstractTightBindingModel::CoreComputeBandStructure" << endl;
    }
}

// build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
//
// nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
// orbitalIndices = array that gives the orbital indices of the connected orbitals
// spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
// hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
// return value = tight binding hamiltonian in real space 

HermitianMatrix AbstractTightBindingModel::BuildTightBindingHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes)
{
  HermitianMatrix TmpHamiltonian (this->NbrBands, true);
  for (int i = 0; i < this->NbrBands; ++i)
    {
      for (int j = 0; j < nbrConnectedOrbitals[i]; ++j)
	{
	  TmpHamiltonian.AddToMatrixElement(i, orbitalIndices[i][j], hoppingAmplitudes[i][j]);
	}
    }
  return TmpHamiltonian;
}

// build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions but without assuming its hermiticiy
//
// nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
// orbitalIndices = array that gives the orbital indices of the connected orbitals
// spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
// hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
// return value = tight binding hamiltonian in real space 

ComplexMatrix AbstractTightBindingModel::BuildTightBindingNonHermitianHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes)
{
  ComplexMatrix TmpHamiltonian;
  return TmpHamiltonian;
}

// get the position of a sublattice site
//
// position = reference on a vector where the answer is supplied
// sublatticeIndex = index of the sub-lattice position
void AbstractTightBindingModel::GetSublatticeVector(RealVector &position, int sublatticeIndex)
{
  cout << "Attention, using dummy method AbstractTightBindingModel::GetSublatticeVector";
  position.ClearVector();
  return;
}

// get the lattice vector for translation along the fundamental lattice directions
//
// latticeVector[out] = reference on a vector where the answer is supplied
// numTranslations = vector of the number of translations along each lattice direction, in units of unit cell size
void AbstractTightBindingModel::GetLatticeVector(RealVector &position, RealVector &numTranslations)
{
  cout << "Attention, using dummy method AbstractTightBindingModel::GetSublatticeVector";
  position.Copy(numTranslations);
  return;
}

// get the elementary lattice vector for translation along the n-th fundamental lattice directions
//
// latticeVector[out] = reference on a vector where the answer is supplied
// dimensionIdx = index of lattice dimension, labelled from 0, ..., d-1
void AbstractTightBindingModel::GetLatticeVector(RealVector &position, int dimensionIdx)
{
  int Dim = position.GetVectorDimension();
  if (dimensionIdx>=0 && dimensionIdx<Dim)
    {
      RealVector tmpVector(Dim, true);
      tmpVector[dimensionIdx] = 1.0;
      this->GetLatticeVector(position, tmpVector);
    }
  cout << "dimensionIdx= "<<dimensionIdx<<" is invalid or position has wrong dimension"<<endl;
  position.ClearVector();
}

// get the lattice vector for translation along the fundamental lattice directions
//
// latticeVector[out] = reference on a vector where the answer is supplied
// numTranslations = coordinates in terms of lattice vectors
// sublatticeIndex = index of the sub-lattice position
void AbstractTightBindingModel::GetSitePosition(RealVector &position, RealVector &numTranslations, int sublatticeIndex)
{
  this->GetLatticeVector(position, numTranslations);
  RealVector TmpVector(position.GetVectorDimension());
  this->GetSublatticeVector(TmpVector, sublatticeIndex);
  position+=TmpVector;
}

// get the tight binding hamiltonian in real space 
// 
// return value = tight binding hamiltonian

HermitianMatrix AbstractTightBindingModel::GetRealSpaceTightBindingHamiltonian()
{
  HermitianMatrix TmpMatrix;
  this->FindConnectedOrbitals();
  if (this->NbrConnectedOrbitals != 0)
    {
      TmpMatrix = this->BuildTightBindingHamiltonianRealSpace(this->NbrConnectedOrbitals, this->ConnectedOrbitalIndices,
							      this->ConnectedOrbitalSpatialIndices, this->ConnectedOrbitalHoppingAmplitudes);
    }
  return TmpMatrix;
}

// get the tight binding hamiltonian in real space , without assuming its hermiticity
// 
// return value = tight binding hamiltonian

ComplexMatrix AbstractTightBindingModel::GetRealSpaceTightBindingNonHermitianHamiltonian()
{
  ComplexMatrix TmpMatrix;
  this->FindConnectedOrbitals();
  if (this->NbrConnectedOrbitals != 0)
    {
      TmpMatrix = this->BuildTightBindingNonHermitianHamiltonianRealSpace(this->NbrConnectedOrbitals, this->ConnectedOrbitalIndices,
									  this->ConnectedOrbitalSpatialIndices, this->ConnectedOrbitalHoppingAmplitudes);
      return TmpMatrix;
    }
  return TmpMatrix;
}

// test the tight binding hamiltonian in real space, checking its hermiticity
// 
// return value = true if the tight binding hamiltonian is hermitian

bool AbstractTightBindingModel::TestRealSpaceTightBindingHamiltonian()
{
  ComplexMatrix TmpMatrix;
  this->FindConnectedOrbitals();
  if (this->NbrConnectedOrbitals != 0)
    {
      TmpMatrix = this->BuildTightBindingNonHermitianHamiltonianRealSpace(this->NbrConnectedOrbitals, this->ConnectedOrbitalIndices,
									  this->ConnectedOrbitalSpatialIndices, this->ConnectedOrbitalHoppingAmplitudes);
      return TmpMatrix.TestHermitian();
    }
  return false;
}

// find the orbitals connected to those located at the origin unit cell
// 
  
void AbstractTightBindingModel::FindConnectedOrbitals()
{
   cout << "warning, tight binding model does not have access to connected orbitals" << endl; 
}

// read the eigenvalues and eigenstates from a band structure file
//
// input = input stream from which header was read
// HeaderSize = size of header that was read
// return value  = size of band structure that was

void AbstractTightBindingModel::ReadEigensystem(ifstream& input, int HeaderSize, unsigned long FileSize)
{
  cout << "warning using dummy method AbstractTightBindingModel::ReadEigensystem" << endl;
}

// get the reciprocal lattice vector for the n-th fundamental lattice direction
//
// latticeVector[out] = reference on a vector where the answer is supplied
// dimensionIdx = index of lattice dimension, labeled from 0, ..., d-1

void AbstractTightBindingModel::GetReciprocalLatticeVector(RealVector &position, int dimensionIdx)
{
  cout << "warning using dummy method AbstractTightBindingModel::GetReciprocalLatticeVector" << endl;
}

