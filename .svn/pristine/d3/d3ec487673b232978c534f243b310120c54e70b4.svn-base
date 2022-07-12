////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 2D tight binding model                  //
//                                                                            //
//                        last modification : 13/10/2012                      //
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
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sys/time.h>

using std::ofstream;
using std::endl;
using std::cout;
using std::max;
using std::min;

#ifdef __FFTW__
#include <fftw3.h>
#endif

#ifdef HAVE_LAPACK
// solve general linear system
extern "C" void FORTRAN_NAME(dgesv)(const int* n, const int* nrhs, const double* a, const int* lda, const int* ipiv, const double* b, const int* ldb, const int* info );
#endif


// default constructor
//

Abstract2DTightBindingModel::Abstract2DTightBindingModel()
{
  this->EmbeddingY = RealVector();
  this->LatticeVector1.Resize(2);
  this->LatticeVector1[0] = 1.0;
  this->LatticeVector1[1] = 0.0;
  this->LatticeVector2.Resize(2);
  this->LatticeVector2[0] = 0.0;
  this->LatticeVector2[1] = 1.0;
  this->TwistAngle = M_PI / 2;
  this->Curvature = NULL;
  this->Chern = NULL;
  this->LLLGammaX = NULL;
  this->LLLGammaY = NULL;
  this->Offset = 0;
  this->OffsetReal = 0;
  this->ProjectedMomenta = 0;
  this->Inversion = ComplexMatrix();
}

// destructor
//

Abstract2DTightBindingModel::~Abstract2DTightBindingModel()
{
  delete[] this->Curvature;
  delete[] this->Chern;
  delete[] this->LLLGammaX;
  delete[] this->LLLGammaY;
  if (this->ProjectedMomenta != 0)
    delete[] this->ProjectedMomenta;
}

// get the position of a sublattice site
//
// position = reference on a vector where the answer is supplied
// sublatticeIndex = index of the sub-lattice position
void Abstract2DTightBindingModel::GetSublatticeVector(RealVector &position, int sublatticeIndex)
{
  cout << "Please overload GetSublatticeVector to make sure conventions for positions match the specific implementation"<<endl;
  if (position.GetVectorDimension()!=2)
    position.Resize(2);
  position[0]=EmbeddingX[sublatticeIndex];
  position[1]=EmbeddingY[sublatticeIndex];
}

// get the lattice vector for translation along the fundamental lattice directions
//
// latticeVector[out] = reference on a vector where the answer is supplied
// numTranslations = vector of the number of translations along each lattice direction, in units of unit cell size
void Abstract2DTightBindingModel::GetLatticeVector(RealVector &position, RealVector &numTranslations)
{
  if (position.GetVectorDimension()!=2)
    position.Resize(2);
  position.ClearVector();
  position.AddLinearCombination(numTranslations[0], this->LatticeVector1);
  position.AddLinearCombination(numTranslations[1], this->LatticeVector2);
}

// get the elementary lattice vector for translation along the n-th fundamental lattice directions
//
// latticeVector[out] = reference on a vector where the answer is supplied
// dimensionIdx = index of lattice dimension, labelled from 0, ..., d-1
void Abstract2DTightBindingModel::GetLatticeVector(RealVector &position, int dimensionIdx)
{
  if (dimensionIdx==0)
    position.Copy(this->LatticeVector1);
  else
    {
      if (dimensionIdx==1)
	position.Copy(this->LatticeVector2);
      else
	position.ClearVector();
    }
}

// get the reciprocal lattice vector for the n-th fundamental lattice direction
//
// latticeVector[out] = reference on a vector where the answer is supplied
// dimensionIdx = index of lattice dimension, labeled from 0, ..., d-1
void Abstract2DTightBindingModel::GetReciprocalLatticeVector(RealVector &position, int dimensionIdx)
{
  position.Resize(2);
  double Prefactor = 2.*M_PI/this->GetUnitCellSize();
  if (dimensionIdx==0)
    {
      position[0] = Prefactor * this->LatticeVector2[1];
      position[1] = -Prefactor * this->LatticeVector2[0];
    }
  else
    {
      if (dimensionIdx==1)
	{
	  position[0] = -Prefactor * this->LatticeVector1[1];
	  position[1] = Prefactor * this->LatticeVector1[0];
	}
      else
	position.ClearVector();
    }
}


// get the size (length / area / volume ... ) of the unit cell
//
// return value =  size
double Abstract2DTightBindingModel::GetUnitCellSize()
{
  return std::fabs(this->LatticeVector1[0] * this->LatticeVector2[1] - this->LatticeVector1[1] * this->LatticeVector2[0]);
}



// write an header that describes the tight binding model
// 
// output = reference on the output stream
// return value  = reference on the output stream

ofstream& Abstract2DTightBindingModel::WriteHeader(ofstream& output)
{
  int Dimension = 2;
  int HeaderSize = (((this->NbrBands + 2) * Dimension + 1) * sizeof(double)) + ((Dimension + 1) * sizeof(int));
  if (this->Inversion.GetNbrRow() == 0)
    {
      for (int i = 0; i < this->NbrBands; ++i)
	for (int j = 0; j < this->NbrBands; ++j)
          {
	    Complex TmpInversion = (i == j);
	    WriteLittleEndian(output, TmpInversion);
          }
    }
  else
    {
      for (int i = 0; i < this->NbrBands; ++i)
	for (int j = 0; j < this->NbrBands; ++j)
          {
	    Complex TmpInversion = this->Inversion[i][j];
	    WriteLittleEndian(output, TmpInversion);
          }
    }
  WriteLittleEndian(output, HeaderSize);
  WriteLittleEndian(output, Dimension);
  WriteLittleEndian(output, this->NbrSiteX);
  WriteLittleEndian(output, this->KxFactor);
  WriteLittleEndian(output, this->GammaX);
  if (this->EmbeddingX.GetVectorDimension() != this->NbrBands)
  {
      double Tmp = 0.0;
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, Tmp);
  }
  else
  {
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, this->EmbeddingX[i]);
  }
  WriteLittleEndian(output, this->NbrSiteY);
  WriteLittleEndian(output, this->KyFactor);
  WriteLittleEndian(output, this->GammaY);
  if (this->EmbeddingY.GetVectorDimension() != this->NbrBands)
  {
      double Tmp = 0.0;
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, Tmp);
  }
  else
  {
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, this->EmbeddingY[i]);
  }
  WriteLittleEndian(output, this->TwistAngle);
  return output; 
}

// open a file and read a header that describes the tight binding model
// 
// return value = size of header that was read (negative if unsuccessful)
int Abstract2DTightBindingModel::ReadHeader(ifstream &File)
{
  int HeaderSize = -1;
  ReadLittleEndian(File, HeaderSize);
  int CorrectDimension = 2;
  int CorrectHeaderSize = (((this->NbrBands + 2) * CorrectDimension + 1) * sizeof(double)) + ((CorrectDimension + 1) * sizeof(int));
  if (HeaderSize >= CorrectHeaderSize)
    {
      int TmpDimension = -1;
      ReadLittleEndian(File, TmpDimension);
      HeaderSize -= sizeof(int);
      if (TmpDimension == CorrectDimension)
	{
	  ReadLittleEndian(File, this->NbrSiteX);
	  ReadLittleEndian(File, this->KxFactor);
	  ReadLittleEndian(File, this->GammaX);	 
          this->EmbeddingX.Resize(this->NbrBands);
          for (int i = 0; i < this->NbrBands; ++i)
            {
              double Tmp = 0.0;
              ReadLittleEndian(File, Tmp);
              this->EmbeddingX[i] = Tmp;
            }
	  ReadLittleEndian(File, this->NbrSiteY);
	  ReadLittleEndian(File, this->KyFactor);
	  ReadLittleEndian(File, this->GammaY);	  
          this->EmbeddingY.Resize(this->NbrBands);
          for (int i = 0; i < this->NbrBands; ++i)
            {
              double Tmp = 0.0;
              ReadLittleEndian(File, Tmp);
              this->EmbeddingY[i] = Tmp;
            }
	  ReadLittleEndian(File, this->TwistAngle);
          HeaderSize -= (CorrectHeaderSize - sizeof(int));
	}
      else
	return -1; // no valid header
      if (HeaderSize > 0) 
	File.seekg (HeaderSize, std::ios::cur);
    }
  else
    {
      return -1; // no valid header
    }
  return HeaderSize;
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool Abstract2DTightBindingModel::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  File << endl;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int LinearizedMomentumIndex = this->GetLinearizedMomentumIndex(kx, ky);
	  File << kx << " " << ky; 
	  for (int i = 0; i < this->NbrBands; ++i)
	    File << " " << this->EnergyBandStructure[i][LinearizedMomentumIndex];
	  File << endl;
	}
    }
  File.close();
  return true;
}



// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool Abstract2DTightBindingModel::WriteAsciiSpectrumColumn(char* fileName)
{
  ofstream File;
  File.precision(14);
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky   E" << endl;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int LinearizedMomentumIndex = this->GetLinearizedMomentumIndex(kx, ky);	  
	  for (int i = 0; i < this->NbrBands; ++i)
          File << kx << " " << ky << " " << this->EnergyBandStructure[i][LinearizedMomentumIndex] << endl;
	}
    }
  File.close();
  return true;
}

// write the energy spectrum in an ASCII file, focusing on lines connecting the high symmetry points
//
// fileName = name of the ASCII file 
// nbrSteps = number of steps between two consecutive high symmetry points
// return value = true if no error occured

bool Abstract2DTightBindingModel::WriteAsciiSpectrumAlongHighSymmetryPoints(char* fileName, int nbrSteps)
{
  double** HighSymmetryPointCoordinates = 0;
  char** HighSymmetryPointNames = 0;
  int NbrHighSymmetryPoints = this->GetHighSymmetryPoints(HighSymmetryPointNames, HighSymmetryPointCoordinates);
  if (NbrHighSymmetryPoints > 0)
    {
      if (nbrSteps < NbrHighSymmetryPoints)
	{
	  nbrSteps = NbrHighSymmetryPoints;
	}

      int* HighSymmetryPointIndices = new int[NbrHighSymmetryPoints];
      double* TmpKxStarting = new double [NbrHighSymmetryPoints];
      double* TmpKyStarting = new double [NbrHighSymmetryPoints];
      double* TmpKxEnding = new double [NbrHighSymmetryPoints];
      double* TmpKyEnding = new double [NbrHighSymmetryPoints];
      double* TmpLengths = new double [NbrHighSymmetryPoints];
      for (int i = 0; i < NbrHighSymmetryPoints; ++i)
	{
	  TmpKxStarting[i] = HighSymmetryPointCoordinates[i][0];
	  TmpKyStarting[i] = HighSymmetryPointCoordinates[i][1];
	  TmpKxEnding[i] = HighSymmetryPointCoordinates[(i + 1) % NbrHighSymmetryPoints][0];
	  TmpKyEnding[i] = HighSymmetryPointCoordinates[(i + 1) % NbrHighSymmetryPoints][1];
	  double TmpLength = this->GetDistanceReciprocalSpace(TmpKxStarting[i], TmpKyStarting[i], TmpKxEnding[i], TmpKyEnding[i]);
	  TmpLengths[i] = TmpLength;
	}
      double TotalLength = 0.0;
      for (int i = 0; i < NbrHighSymmetryPoints; ++i)
	{
	  TotalLength += TmpLengths[i];
	}
      double LengthIncrement = TotalLength / ((double) nbrSteps);
      
      HighSymmetryPointIndices[0] = 0;      
      for (int i = 1; i < NbrHighSymmetryPoints; ++i)
	{
	  HighSymmetryPointIndices[i] = HighSymmetryPointIndices[i - 1] + (int) ((TmpLengths[i - 1] / LengthIncrement));
	  if (HighSymmetryPointIndices[i] == HighSymmetryPointIndices[i - 1])
	    HighSymmetryPointIndices[i]++; 
	}
      nbrSteps = HighSymmetryPointIndices[NbrHighSymmetryPoints - 1] + ((int) (TmpLengths[NbrHighSymmetryPoints - 1] / LengthIncrement));
      if (nbrSteps == HighSymmetryPointIndices[NbrHighSymmetryPoints - 1])
	{
	  ++nbrSteps;
	}

      double* TmpKx = new double [nbrSteps];
      double* TmpKy = new double [nbrSteps];
      int Index = 0;
      for (int i = 0; i < (NbrHighSymmetryPoints - 1); ++i)
	{
	  int TmpNbrSteps = HighSymmetryPointIndices[i + 1] - HighSymmetryPointIndices[i];
	  for (int j = 0; j < TmpNbrSteps; ++j)
	    {
	      TmpKx[Index] = (((((double) j) / ((double) TmpNbrSteps)) * TmpKxEnding[i]) 
			      + ((((double) (TmpNbrSteps - j)) / ((double) TmpNbrSteps)) * TmpKxStarting[i]));
	      TmpKy[Index] = (((((double) j) / ((double) TmpNbrSteps)) * TmpKyEnding[i]) 
			      + ((((double) (TmpNbrSteps - j)) / ((double) TmpNbrSteps)) * TmpKyStarting[i]));
	      ++Index;
	    }
	}
      int TmpNbrSteps = nbrSteps - HighSymmetryPointIndices[NbrHighSymmetryPoints - 1];
      for (int j = 0; j < TmpNbrSteps; ++j)
	{
	  TmpKx[Index] = (((((double) j) / ((double) TmpNbrSteps)) * TmpKxEnding[NbrHighSymmetryPoints - 1]) 
			  + ((((double) (TmpNbrSteps - j)) / ((double) TmpNbrSteps)) * TmpKxStarting[NbrHighSymmetryPoints - 1]));
	  TmpKy[Index] = (((((double) j) / ((double) TmpNbrSteps)) * TmpKyEnding[NbrHighSymmetryPoints - 1]) 
			  + ((((double) (TmpNbrSteps - j)) / ((double) TmpNbrSteps)) * TmpKyStarting[NbrHighSymmetryPoints - 1]));
	  ++Index;
	}
      double* TmpEnergies = new double [this->NbrBands];
      ofstream File;
      File.open(fileName);
      for (int i = 0; i < NbrHighSymmetryPoints; ++i)
	File << " # index = " << HighSymmetryPointIndices[i] << " -> " << HighSymmetryPointNames[i] << endl;
      File << "# index    k_x    k_y";
      for (int i = 0; i < this->NbrBands; ++i)
	File <<  "    E_" << i;
      File << endl;
      for (int i = 0; i < nbrSteps; ++i)
	{
	  this->ComputeBandStructureSinglePoint(TmpKx[i], TmpKy[i], TmpEnergies);
	  File << i << " " << TmpKx[i] << " " << TmpKy[i]; 
	  for (int j = 0; j < this->NbrBands; ++j)
	    File << " " << TmpEnergies[j];
	  File << endl;
	}
      File.close();
      return true;
    }
  else
    {
      cout << "high symmetry points are not available for this tight binding model" << endl;
      return false;
    }
}

// write the full band structure information in an ASCII file
//
// fileName = name of the output file 
// return value = true if no error occured  

bool Abstract2DTightBindingModel::WriteBandStructureASCII(char* fileName)
{
  ofstream File;
  File.open(fileName);
  File.precision(14);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  for (int i = 0; i < this->NbrBands; ++i)
    for (int j = 0; j < this->NbrBands; ++j)
      File <<  "    U_{" << i << ", " << j << "}";
  File << endl;
  Complex Tmp;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int LinearizedMomentumIndex = this->GetLinearizedMomentumIndex(kx, ky);
	  File << kx << " " << ky; 
	  for (int i = 0; i < this->NbrBands; ++i)
	    File << " " << this->EnergyBandStructure[i][LinearizedMomentumIndex];
	  for (int i = 0; i < this->NbrBands; ++i)
	    for (int j = 0; j < this->NbrBands; ++j)
	      {
		this->GetOneBodyMatrix(LinearizedMomentumIndex).GetMatrixElement(i, j, Tmp);
		File <<  "    " << Tmp;
	      }
	  File << endl;
	}
    }

  File.close();
  return true;
}

// compute the exponentiated, unitary Abelian connection
//
// kx = momentum along x
// ky = momentum along y
// qx = momentum transfer along x
// qy = momentum transfer along y
// band = band index
// return value = < u(k) | u(k+q) >

Complex Abstract2DTightBindingModel::GetAbelianConnection(int kx, int ky, int qx, int qy, int band)
{
    // [band][orbital] = Conj(<orbital|band>)
    int k1 = this->GetLinearizedMomentumIndexSafe(kx, ky);
    int k2 = this->GetLinearizedMomentumIndexSafe(kx + qx, ky + qy);
    ComplexVector bra = this->GetOneBodyMatrix(k1)[band];
    ComplexVector ket = this->GetOneBodyMatrix(k2)[band];

    Complex inner = 0.0;
    if (this->EmbeddingX.GetVectorDimension() != this->NbrBands)
        inner = ket * bra;
    else
    {
        for (int i = 0; i < this->NbrBands; ++i)
            inner += Phase(- 2.0 * M_PI * (qx * this->EmbeddingX[i] / this->NbrSiteX + qy * this->EmbeddingY[i] / this->NbrSiteY)) * bra[i] * Conj(ket[i]);
    }
    double n = Norm(inner);
    if (n < 1e-13)
    {
        cout << "Cannot make connection unitary due to orthogonality!" << endl;
        exit(1);
    }
    return inner / n;
}

// compute the exponentiated, unitary Abelian connection times the quantum distance
//
// kx = momentum along x
// ky = momentum along y
// qx = momentum transfer along x
// qy = momentum transfer along y
// band = band index
// return value = < u(k) | u(k+q) >

Complex Abstract2DTightBindingModel::GetAbelianConnectionQuantumDistance(int kx, int ky, int qx, int qy, int band)
{
    // [band][orbital] = Conj(<orbital|band>)
    int k1 = this->GetLinearizedMomentumIndexSafe(kx, ky);
    int k2 = this->GetLinearizedMomentumIndexSafe(kx + qx, ky + qy);
    ComplexVector bra = this->GetOneBodyMatrix(k1)[band];
    ComplexVector ket = this->GetOneBodyMatrix(k2)[band];

    Complex inner = 0.0;
    if (this->EmbeddingX.GetVectorDimension() != this->NbrBands)
        inner = ket * bra;
    else
    {
        for (int i = 0; i < this->NbrBands; ++i)
            inner += Phase(- 2.0 * M_PI * (qx * this->EmbeddingX[i] / this->NbrSiteX + qy * this->EmbeddingY[i] / this->NbrSiteY)) * bra[i] * Conj(ket[i]);
    }
    return inner;
}

// compute the unitary Abelian Wilson loop
//
// ky = momentum along y
// band = band index
// return value = value of the Wilson loop

Complex Abstract2DTightBindingModel::GetAbelianWilsonLoopX(int ky, int band)
{
    Complex prod = 1.0;
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        prod *= this->GetAbelianConnection(kx, ky, 1, 0, band);
    return prod;
}

// compute the unitary Abelian Wilson loop
//
// kx = momentum along x
// band = band index
// return value = value of the Wilson loop

Complex Abstract2DTightBindingModel::GetAbelianWilsonLoopY(int kx, int band)
{
    Complex prod = 1.0;
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
        prod *= this->GetAbelianConnection(kx, ky, 0, 1, band);
    return prod;
}

// compute the stream function for the part of Berry connections that accounts for curvature fluctuations
//
// band = band index
// phi = reference to the vector storing the stream function over linearized BZ
// return = 0 if succeed, otherwise fail

int Abstract2DTightBindingModel::ComputeStreamFunction(int band, RealVector& phi, double& vx, double& vy)
{
    if (this->Curvature == NULL)
        this->ComputeCurvature();


#ifdef __FFTW__
    fftw_complex* source = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->NbrStatePerBand);
    fftw_complex* fsource = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->NbrStatePerBand);

    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int k = kx * this->NbrSiteY + ky; // need to enforce row-major layout
            source[k][0] = this->Curvature[band][ky][kx] - ((double) this->Chern[band]) / this->NbrStatePerBand;
            source[k][1] = 0.0;
        }
    fftw_plan plan_source = fftw_plan_dft_2d(this->NbrSiteX, this->NbrSiteY, source, fsource, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_source);
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            fsource[kx * this->NbrSiteY + ky][0] /= this->NbrStatePerBand;
            fsource[kx * this->NbrSiteY + ky][1] /= this->NbrStatePerBand;
        }


    fftw_complex* solution = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->NbrStatePerBand);
    fftw_complex* fsolution = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * this->NbrStatePerBand);

    fsolution[0][0] = 0.0;
    fsolution[0][1] = 0.0;
    for (int nx = 0; nx < this->NbrSiteX; ++nx)
    {
        for (int ny = 0; ny < this->NbrSiteY; ++ny)
        {
            int n = nx * this->NbrSiteY + ny;
            if (n == 0)
                continue;

            double nsq = - 4.0;
            nsq += 2 * cos(2 * M_PI * ((double) nx) / this->NbrSiteX);
            nsq += 2 * cos(2 * M_PI * ((double) ny) / this->NbrSiteY);

            fsolution[n][0] = fsource[n][0] / nsq;
            fsolution[n][1] = fsource[n][1] / nsq;
        }
    }

    fftw_plan plan_solution = fftw_plan_dft_2d(this->NbrSiteX, this->NbrSiteY, fsolution, solution, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_solution);

    double maxerror = 0.0;
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        int kxm = (kx - 1 + this->NbrSiteX) % this->NbrSiteX;
        int kxp = (kx + 1) % this->NbrSiteX;
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int kym = (ky - 1 + this->NbrSiteY) % this->NbrSiteY;
            int kyp = (ky + 1) % this->NbrSiteY;

            int k = kx * this->NbrSiteY + ky;
            int ku = kx * this->NbrSiteY + kyp;
            int kd = kx * this->NbrSiteY + kym;
            int kl = kxm * this->NbrSiteY + ky;
            int kr = kxp * this->NbrSiteY + ky;

            if (fabs(solution[k][1]) > maxerror)
                maxerror = fabs(solution[k][1]);
            double error = fabs(solution[kl][0] + solution[kr][0] + solution[ku][0] + solution[kd][0] - 4 * solution[k][0] - source[k][0]);
            if (error > maxerror)
                maxerror = error;
        }
    }
    cout << "FFTW-Poisson error: " << maxerror << endl;

    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
            phi[this->GetLinearizedMomentumIndexSafe(kx, ky)] = solution[kx * this->NbrSiteY + ky][0];

    fftw_destroy_plan(plan_source);
    fftw_destroy_plan(plan_solution);
    fftw_free(source);
    fftw_free(fsource);
    fftw_free(solution);
    fftw_free(fsolution);

#else
#ifdef HAVE_LAPACK
    // curvature fluctuations
    double* raw = new double[this->NbrStatePerBand];
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
            raw[this->GetLinearizedMomentumIndexSafe(kx, ky)] = this->Curvature[band][ky][kx] - ((double) this->Chern[band]) / this->NbrStatePerBand;

    // laplacian with PBC
    double* laplacian = new double[this->NbrStatePerBand * this->NbrStatePerBand];
    for (int k = 0; k < this->NbrStatePerBand; ++k)
    {
        for (int l = 0; l < this->NbrStatePerBand; ++l)
        {
            laplacian[k + l * this->NbrStatePerBand] = 0.0;
        }
        laplacian[k + k * this->NbrStatePerBand] = -4.0;
    }
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int k = this->GetLinearizedMomentumIndexSafe(kx, ky);
            laplacian[k + this->GetLinearizedMomentumIndexSafe(kx, ky + 1) * this->NbrStatePerBand] = 1.0;
            laplacian[k + this->GetLinearizedMomentumIndexSafe(kx, ky - 1) * this->NbrStatePerBand] = 1.0;
            laplacian[k + this->GetLinearizedMomentumIndexSafe(kx + 1, ky) * this->NbrStatePerBand] = 1.0;
            laplacian[k + this->GetLinearizedMomentumIndexSafe(kx - 1, ky) * this->NbrStatePerBand] = 1.0;
        }
    }

    // solve poisson: \nabla^2 phi = curvature - 1/Nf
    int n = this->NbrStatePerBand;
    int nrhs = 1;
    int* ipiv = new int[this->NbrStatePerBand];
    int info = 42;
    FORTRAN_NAME(dgesv)(&n, &nrhs, laplacian, &n, ipiv, raw, &n, &info);
    delete[] ipiv;
    if (info)
    {
        cout << "Lapack DGESV failed with exit code " << info << endl;
        return -1;
    }
    for (int k = 0; k < this->NbrStatePerBand; ++k)
        phi[k] = raw[k];
    delete[] laplacian;
    delete[] raw;
#else
    cout << "Need FFTW or LAPACK to solve the Poisson equation for stream function." << endl;
    return -1;
#endif
#endif

    // overall shift in Ax and Ay, to reset Wx(0) and Wy(0) to 1
    vx = 0.0;
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        vx += phi[this->GetLinearizedMomentumIndexSafe(kx, 0)] - phi[this->GetLinearizedMomentumIndexSafe(kx, -1)];
    vx /= this->NbrSiteX;

    vy = 0.0;
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
        vy -= phi[this->GetLinearizedMomentumIndexSafe(0, ky)] - phi[this->GetLinearizedMomentumIndexSafe(-1, ky)];
    vy /= this->NbrSiteY;

    return 0;
}

// build the gauge transform such that gauge(k) * |k⟩_lat is in the "Γ"-shaped parallel-transport gauge
// 
// band = band index
// gauge = reference to the gauge transform

void Abstract2DTightBindingModel::BuildParallelTransportGauge(int band, ComplexMatrix& gauge)
{
    if (this->Curvature == NULL)
        this->ComputeCurvature();

    gauge[0][0] = 1.0;

    // p-t along ky
    // we don't evenly spread Wy(0) over the bonds, in order to match with the quasi-periodic LLL gauge
    // Wy(0) is thus hidden in the bond crossing periodic boundary
    for (int ky = 1; ky < this->NbrSiteY; ++ky)
        gauge[ky][0] = gauge[ky - 1][0] / this->GetAbelianConnection(0, ky - 1, 0, 1, band);

    // p-t with rotation along kx
    // need to fit the phase of the LLL connection Ay = e^{- i 2 Pi (ky + gammay) / Nf}
    double* theta = new double[this->NbrSiteY];
    theta[0] = Arg(this->GetAbelianWilsonLoopX(0, band)); // need to be consistent with GetLLLGammaY
    if (this->Chern[band] > 0)
    {
        for (int ky = 1; ky < this->NbrSiteY; ++ky)
        {
            theta[ky] = Arg(this->GetAbelianWilsonLoopX(ky, band));
            while (theta[ky] <= theta[ky - 1])
                theta[ky] += 2 * M_PI;
            while (theta[ky] > theta[ky - 1])
                theta[ky] -= 2 * M_PI;
        }
    }
    else
    {
        for (int ky = 1; ky < this->NbrSiteY; ++ky)
        {
            theta[ky] = Arg(this->GetAbelianWilsonLoopX(ky, band));
            while (theta[ky] >= theta[ky - 1])
                theta[ky] -= 2 * M_PI;
            while (theta[ky] < theta[ky - 1])
                theta[ky] += 2 * M_PI;
        }
    }

    cout << "Theta(ky) / (2 Pi)" << endl;
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
        cout << "  " << theta[ky] / (2 * M_PI)<< endl;
    double total = fabs(theta[this->NbrSiteY - 1] - theta[0]) / (2 * M_PI);
    if ((total > abs(this->Chern[band])) || (total < (abs(this->Chern[band]) - 1)))
        cout << "Bad branch choice!" << endl;
    cout << endl;

    for (int ky = 0; ky < this->NbrSiteY; ++ky)
    {
        Complex lambda = Phase(theta[ky] / this->NbrSiteX);
        for (int kx = 1; kx < this->NbrSiteX; ++kx)
            gauge[ky][kx] = gauge[ky][kx - 1] * lambda / this->GetAbelianConnection(kx - 1, ky, 1, 0, band);
    }
    delete[] theta;
}

// build the gauge transform such that gauge(k) * |k⟩_lat is in the generalized π/2-rotated Landau gauge
//
// band = band index
// gauge = reference to the gauge transform
// return = 0 if succeed, otherwise fail

int Abstract2DTightBindingModel::BuildGeneralizedLandauGauge(int band, ComplexMatrix& gauge)
{
    double vx, vy;
    // matrices are initialized by (nrow, ncol) and accessed by [col][row]
    RealMatrix ax(this->NbrSiteX, this->NbrSiteY, true);
    RealMatrix ay(this->NbrSiteX, this->NbrSiteY, true);
    RealVector phi(this->NbrSiteX * this->NbrSiteY, true);

    // handle curvature fluctuations
    if (this->Curvature == NULL)
        this->ComputeCurvature();
    if (this->ComputeStreamFunction(band, phi, vx, vy) != 0)
        return -1;
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            ax[ky][kx] += vx + phi[this->GetLinearizedMomentumIndexSafe(kx, ky - 1)] - phi[this->GetLinearizedMomentumIndexSafe(kx, ky)];
            ay[ky][kx] += vy + phi[this->GetLinearizedMomentumIndexSafe(kx, ky)] - phi[this->GetLinearizedMomentumIndexSafe(kx - 1, ky)];
        }
    }

    // twisting angles needed to handle Wx(0) and Wy(0)
    double gammaX = this->GetLLLGammaX(band);
    double gammaY = this->GetLLLGammaY(band);

    // handle curvature average
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            ax[ky][kx] += - this->Chern[band] * (ky + gammaY) / this->NbrStatePerBand;
            if (ky == this->NbrSiteY - 1)
                ay[ky][kx] += this->Chern[band] * (kx + gammaX) / this->NbrSiteX;
        }
    }

    // construct gauge transform to the target gauge given by A?Phase, gauge(k) * |k⟩_lat
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
            gauge[ky][kx] = 1.0;
    for (int kx = 1; kx < this->NbrSiteX; ++kx)
        gauge[0][kx] = gauge[0][kx - 1] * Phase(2 * M_PI * ax[0][kx - 1]) / this->GetAbelianConnection(kx - 1, 0, 1, 0, band);
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 1; ky < this->NbrSiteY; ++ky)
            gauge[ky][kx] = gauge[ky - 1][kx] * Phase(2 * M_PI * ay[ky - 1][kx]) / this->GetAbelianConnection(kx, ky - 1, 0, 1, band);
    }

    // check gauge construction
    bool Fail = false;

    // check connections against plaquette Wilson loops
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            double curvature = ax[ky][kx] + ay[ky][(kx + 1) % this->NbrSiteX] - ax[(ky + 1) % this->NbrSiteY][kx] - ay[ky][kx];
            if ((kx == this->NbrSiteX - 1) && (ky == this->NbrSiteY - 1)) // trivial fix. does not affect phase
                curvature += Chern[band];

            double diff = fabs(curvature - this->Curvature[band][ky][kx]);
            if (diff > 1e-13)
            {
                Fail = true;
                cout << "Curvature Mismatch @ k = (" << kx << "," << ky << "), Error = " << diff << endl;
            }
        }
    }

    // check connections against large Wilson loops
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
    {
        double sum = 0;
        for (int kx = 0; kx < this->NbrSiteX; ++kx)
            sum += ax[ky][kx];
        double diff = Arg(Phase(2 * M_PI * sum) / this->GetAbelianWilsonLoopX(ky, band));
        if (diff > 1e-13)
        {
            Fail = true;
            cout << "Wx Mismatch @ ky = " << ky << ", Error = " << diff << endl;
        }
    }
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        double sum = 0;
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
            sum += ay[ky][kx];
        double diff = Arg(Phase(2 * M_PI * sum) / this->GetAbelianWilsonLoopY(kx, band));
        if (diff > 1e-13)
        {
            Fail = true;
            cout << "Wy Mismatch @ kx = " << kx << ", Error = " << diff << endl;
        }
    }

    // check gauge transform
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        int kxm = (kx - 1 + this->NbrSiteX) % this->NbrSiteX;
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int kym = (ky - 1 + this->NbrSiteY) % this->NbrSiteY;
            double diff = Arg((gauge[ky][kx] / gauge[kym][kx]) / (Phase(2 * M_PI * ay[kym][kx]) / this->GetAbelianConnection(kx, kym, 0, 1, band)));
            if (diff > 1e-13)
            {
                Fail = true;
                cout << "Gauge Mismatch @ k = (" << kx << "," << ky << "), along y, Error = " << diff << endl;
            }
            diff = Arg((gauge[ky][kx] / gauge[ky][kxm]) / (Phase(2 * M_PI * ax[ky][kxm]) / this->GetAbelianConnection(kxm, ky, 1, 0, band)));
            if (diff > 1e-13)
            {
                Fail = true;
                cout << "Gauge Mismatch @ k = (" << kx << "," << ky << "), along x, Error = " << diff << endl;
            }
        }
    }

    return int(Fail);
}

// compute the curvature over each plaquette in the BZ, and also Chern number
//
// band = band index
// return = 0 if succeed, otherwise fail

void Abstract2DTightBindingModel::ComputeCurvature()
{
    if (this->Curvature == NULL)
        this->Curvature = new RealMatrix[this->NbrBands];
    if (this->Chern == NULL)
        this->Chern = new int[this->NbrBands];

    for (int b = 0; b < this->NbrBands; ++b)
        this->Curvature[b].ResizeAndClean(this->NbrSiteX, this->NbrSiteY);

    for (int b = 0; b < this->NbrBands; ++b)
    {
        for (int kx = 0; kx < this->NbrSiteX; ++kx)
            for (int ky = 0; ky < this->NbrSiteY; ++ky)
                this->Curvature[b][ky][kx] = this->ComputeCurvatureSinglePlaquette(kx, ky, b);

        double sum = 0.0;
        for (int kx = 0; kx < this->NbrSiteX; ++kx)
            for (int ky = 0; ky < this->NbrSiteY; ++ky)
                sum += this->Curvature[b][ky][kx];

        if (sum == 0.0)
            this->Chern[b] = 0;
        else
        {
            sum += (sum / fabs(sum)) * 0.5;
            this->Chern[b] = int(sum);
        }

        if (fabs(this->Chern[b] - (sum - (sum / fabs(sum)) * 0.5)) > 1e-13) // shouldn't happen unless connections are seriously screwed up
            cout << "non-integer Chern number for band " << b << "?! (curvature sum = " << sum << ")" << endl;
        if (this->Chern[b] == 0)
            cout << "Zero Chern number for band " << b << "?!" << endl;
    }

    if (this->LLLGammaX == NULL)
        this->LLLGammaX = new double[this->NbrBands];
    if (this->LLLGammaY == NULL)
        this->LLLGammaY = new double[this->NbrBands];
    for (int b = 0; b < this->NbrBands; ++b)
    {
        double gammaX = Arg(this->GetAbelianWilsonLoopY(0, b)); // Arg takes value in (-π, π]
        gammaX *= this->NbrSiteX / (2 * M_PI * this->Chern[b]);
        this->LLLGammaX[b] = gammaX;

        double gammaY = Arg(this->GetAbelianWilsonLoopX(0, b)); // Arg takes value in (-π, π]
        gammaY *= - this->NbrSiteY / (2 * M_PI * this->Chern[b]);
        this->LLLGammaY[b] = gammaY;
    }
}


// write the eigenvalues of the D matrix in an ASCII file
//
// fileName = name of the ASCII file 
//nbrOccupiedBands = nbr of occupied bands
// return value = true if no error occured

bool Abstract2DTightBindingModel::WriteAsciiDMatrixEigenValues(char* fileName, int nbrOccupiedBands)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# ky" ;
  for (int i = 0; i < nbrOccupiedBands; ++i)
    File <<  "   DEigenValue_" << i << "    Theta_" << i;
  File << endl;
  
  double distancePlus;
  double distanceMoins;
  double distanceMod2PiPlus;
  double distanceMod2PiMoins;
 
  double theta1;
  double theta2;
  
  Complex** Lambda = this->ComputeDMatrixEigenvalues(nbrOccupiedBands, 0, this->NbrSiteY - 1, this->NbrSiteY); 
  double** Theta = new double*[this->NbrSiteY];
  for (int ky = 0; ky < this->NbrSiteY; ++ky)
  {
   Theta[ky] = new double[2];
   for (int i = 0; i < 2; ++i)
   {
     theta1 = atan2(Lambda[ky][nbrOccupiedBands - 2].Im,Lambda[ky][nbrOccupiedBands - 2].Re);
     theta2 = atan2(Lambda[ky][nbrOccupiedBands - 1].Im,Lambda[ky][nbrOccupiedBands - 1].Re);
     Theta[ky][0] = max(theta1, theta2);
     Theta[ky][1] = min(theta1, theta2); 
   }
  }
  cout << "test" << endl;
  for (int ky = 0; ky < this->NbrSiteY - 1; ++ ky)
  {
    distancePlus = fabs(Theta[ky][0] - Theta[ky + 1][0]);
    distanceMod2PiPlus = fabs(Theta[ky][0] - Theta[ky + 1][1] - 2*M_PI);
    distanceMoins = fabs(Theta[ky][1] - Theta[ky + 1][1]);
    distanceMod2PiMoins = fabs(Theta[ky][1] - Theta[ky + 1][0] + 2*M_PI);
    
    if (distanceMod2PiPlus < distancePlus)
    {
     double Tmp = Theta[ky + 1][0];
     Theta[ky + 1][0] = Theta[ky + 1][1] + 2*M_PI;
     Theta[ky + 1][1] = Tmp;
    }
    
    if (distanceMod2PiMoins < distanceMoins)
    {
     double Tmp = Theta[ky + 1][1];
     Theta[ky + 1][1] = Theta[ky + 1][0] - 2*M_PI;
     Theta[ky + 1][0] = Tmp;
    }
  }
  for (int ky = 0; ky < this->NbrSiteY; ++ky)
    {
      File << ky; 
      File << " " << Lambda[ky][0] << " " << atan2(Lambda[ky][0].Im,Lambda[ky][0].Re) << " " << Lambda[ky][1] << " "<< atan2(Lambda[ky][1].Im,Lambda[ky][1].Re) << " " << Theta[ky][0] << " " << Theta[ky][1] ;
      File << endl;
    }
   
  File.close();
  return true;
}

// compute the Chern number of a given band
//
// band = band index
// return value = Chern number

double Abstract2DTightBindingModel::ComputeChernNumber(int band)
{
  if (this->HaveOneBodyBasis() == false)
    {
      cout << "error, the tight binding model does not provide the one body basis" << endl;
      return 0.0;
    }
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  Complex TmpChernNumber = 0.0;
  Complex Tmp1[4];
  Complex Tmp2[8];
  for (long LinearizedMomentumIndex = 0l; LinearizedMomentumIndex < this->NbrStatePerBand; ++LinearizedMomentumIndex)
    {
      int Kx;
      int Ky;
      this->GetLinearizedMomentumIndex(LinearizedMomentumIndex, Kx, Ky);
      int LinearizedMomentumIndexIncX = this->GetLinearizedMomentumIndex((Kx + 1) % this->NbrSiteX, Ky);
      int LinearizedMomentumIndexDecX;
      if (Kx > 0)
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex((Kx - 1) % this->NbrSiteX, Ky);
      else
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex(this->NbrSiteX - 1, Ky);
      int LinearizedMomentumIndexIncY = this->GetLinearizedMomentumIndex(Kx, (Ky + 1) % this->NbrSiteY);
      int LinearizedMomentumIndexDecY;
      if (Ky > 0)
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, (Ky - 1) % this->NbrSiteY);
      else
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, this->NbrSiteY - 1);

      ComplexMatrix& LocalBasis = this->OneBodyBasis[LinearizedMomentumIndex];
      ComplexMatrix& LocalBasisIncX = this->OneBodyBasis[LinearizedMomentumIndexIncX];
      ComplexMatrix& LocalBasisDecX = this->OneBodyBasis[LinearizedMomentumIndexDecX];
      ComplexMatrix& LocalBasisIncY = this->OneBodyBasis[LinearizedMomentumIndexIncY];
      ComplexMatrix& LocalBasisDecY = this->OneBodyBasis[LinearizedMomentumIndexDecY];  
      Tmp1[0] = 0.0;
      Tmp1[1] = 0.0;
      Tmp1[2] = 0.0;
      Tmp1[3] = 0.0;

      Tmp2[0] = 0.0;
      Tmp2[1] = 0.0;
      Tmp2[2] = 0.0;
      Tmp2[3] = 0.0;
      Tmp2[4] = 0.0;
      Tmp2[5] = 0.0;
      Tmp2[6] = 0.0;
      Tmp2[7] = 0.0;

      for (int i = 0; i < this->NbrBands; ++i)
	{
	  Tmp1[0] += LocalBasis[band][i] * Conj(LocalBasisIncX[band][i]);
	  Tmp1[1] += LocalBasis[band][i] * Conj(LocalBasisDecX[band][i]);
	  Tmp1[2] += LocalBasis[band][i] * Conj(LocalBasisIncY[band][i]);
	  Tmp1[3] += LocalBasis[band][i] * Conj(LocalBasisDecY[band][i]);

	  Tmp2[0] += Conj(LocalBasisIncX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[1] += Conj(LocalBasisDecX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[2] += Conj(LocalBasisIncX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[3] += Conj(LocalBasisDecX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[4] += Conj(LocalBasisIncY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[5] += Conj(LocalBasisDecY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[6] += Conj(LocalBasisIncY[band][i]) * LocalBasisDecX[band][i];
	  Tmp2[7] += Conj(LocalBasisDecY[band][i]) * LocalBasisDecX[band][i];
	}

      TmpChernNumber += (Tmp1[2] * Conj(Tmp1[0]) * Tmp2[0]);
      TmpChernNumber -= (Tmp1[2] * Conj(Tmp1[1]) * Tmp2[1]);
      TmpChernNumber -= (Tmp1[3] * Conj(Tmp1[0]) * Tmp2[2]);
      TmpChernNumber += (Tmp1[3] * Conj(Tmp1[1]) * Tmp2[3]);
	  
      TmpChernNumber -= (Tmp1[0] * Conj(Tmp1[2]) * Tmp2[4]);
      TmpChernNumber += (Tmp1[0] * Conj(Tmp1[3]) * Tmp2[5]);
      TmpChernNumber += (Tmp1[1] * Conj(Tmp1[2]) * Tmp2[6]);
      TmpChernNumber -= (Tmp1[1] * Conj(Tmp1[3]) * Tmp2[7]);

    }
  TmpChernNumber /= 8.0 * M_PI;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
  cout << "Chern number computed in  " << Dt << "s" << endl;
  return TmpChernNumber.Im;
}

// compute the Chern number of several bands
//
// bands = band indices
// nbrBands = number of bands that have to be taken into account
// return value = Chern number

double Abstract2DTightBindingModel::ComputeChernNumber(int* bands, int nbrBands)
{
  cout << "warning, this code is not working !" << endl;
  return 0.0;
  if (this->HaveOneBodyBasis() == false)
    {
      cout << "error, the tight binding model does not provide the one body basis" << endl;
      return 0.0;
    }
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  Complex TmpChernNumber = 0.0;
  Complex Tmp1[4];
  Complex Tmp2[8];
  for (long LinearizedMomentumIndex = 0l; LinearizedMomentumIndex < this->NbrStatePerBand; ++LinearizedMomentumIndex)
    {
      int Kx;
      int Ky;
      this->GetLinearizedMomentumIndex(LinearizedMomentumIndex, Kx, Ky);
      int LinearizedMomentumIndexIncX = this->GetLinearizedMomentumIndex((Kx + 1) % this->NbrSiteX, Ky);
      int LinearizedMomentumIndexDecX;
      if (Kx > 0)
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex((Kx - 1) % this->NbrSiteX, Ky);
      else
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex(this->NbrSiteX - 1, Ky);
      int LinearizedMomentumIndexIncY = this->GetLinearizedMomentumIndex(Kx, (Ky + 1) % this->NbrSiteY);
      int LinearizedMomentumIndexDecY;
      if (Ky > 0)
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, (Ky - 1) % this->NbrSiteY);
      else
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, this->NbrSiteY - 1);

      ComplexMatrix& LocalBasis = this->OneBodyBasis[LinearizedMomentumIndex];
      ComplexMatrix& LocalBasisIncX = this->OneBodyBasis[LinearizedMomentumIndexIncX];
      ComplexMatrix& LocalBasisDecX = this->OneBodyBasis[LinearizedMomentumIndexDecX];
      ComplexMatrix& LocalBasisIncY = this->OneBodyBasis[LinearizedMomentumIndexIncY];
      ComplexMatrix& LocalBasisDecY = this->OneBodyBasis[LinearizedMomentumIndexDecY];  
      Tmp1[0] = 0.0;
      Tmp1[1] = 0.0;
      Tmp1[2] = 0.0;
      Tmp1[3] = 0.0;

      Tmp2[0] = 0.0;
      Tmp2[1] = 0.0;
      Tmp2[2] = 0.0;
      Tmp2[3] = 0.0;
      Tmp2[4] = 0.0;
      Tmp2[5] = 0.0;
      Tmp2[6] = 0.0;
      Tmp2[7] = 0.0;

      for (int i = 0; i < this->NbrBands; ++i)
	{
	  Tmp1[0] += LocalBasis[bands[0]][i] * Conj(LocalBasisIncX[bands[0]][i]);
	  Tmp1[1] += LocalBasis[bands[0]][i] * Conj(LocalBasisDecX[bands[0]][i]);
	  Tmp1[2] += LocalBasis[bands[0]][i] * Conj(LocalBasisIncY[bands[0]][i]);
	  Tmp1[3] += LocalBasis[bands[0]][i] * Conj(LocalBasisDecY[bands[0]][i]);

	  Tmp2[0] += Conj(LocalBasisIncX[bands[0]][i]) * LocalBasisIncY[bands[0]][i];
	  Tmp2[1] += Conj(LocalBasisDecX[bands[0]][i]) * LocalBasisIncY[bands[0]][i];
	  Tmp2[2] += Conj(LocalBasisIncX[bands[0]][i]) * LocalBasisDecY[bands[0]][i];
	  Tmp2[3] += Conj(LocalBasisDecX[bands[0]][i]) * LocalBasisDecY[bands[0]][i];
	  Tmp2[4] += Conj(LocalBasisIncY[bands[0]][i]) * LocalBasisIncX[bands[0]][i];
	  Tmp2[5] += Conj(LocalBasisDecY[bands[0]][i]) * LocalBasisIncX[bands[0]][i];
	  Tmp2[6] += Conj(LocalBasisIncY[bands[0]][i]) * LocalBasisDecX[bands[0]][i];
	  Tmp2[7] += Conj(LocalBasisDecY[bands[0]][i]) * LocalBasisDecX[bands[0]][i];
	}

      TmpChernNumber += (Tmp1[2] * Conj(Tmp1[0]) * Tmp2[0]);
      TmpChernNumber -= (Tmp1[2] * Conj(Tmp1[1]) * Tmp2[1]);
      TmpChernNumber -= (Tmp1[3] * Conj(Tmp1[0]) * Tmp2[2]);
      TmpChernNumber += (Tmp1[3] * Conj(Tmp1[1]) * Tmp2[3]);
	  
      TmpChernNumber -= (Tmp1[0] * Conj(Tmp1[2]) * Tmp2[4]);
      TmpChernNumber += (Tmp1[0] * Conj(Tmp1[3]) * Tmp2[5]);
      TmpChernNumber += (Tmp1[1] * Conj(Tmp1[2]) * Tmp2[6]);
      TmpChernNumber -= (Tmp1[1] * Conj(Tmp1[3]) * Tmp2[7]);

    }
  TmpChernNumber /= 8.0 * M_PI;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
  cout << "Chern number computed in  " << Dt << "s" << endl;
  return TmpChernNumber.Im;
}

// compute the Berry curvature  of a given band
//
// band = band index
// fileName = name of the output file 
// return value = Chern number

double Abstract2DTightBindingModel::ComputeBerryCurvature(int band, char* fileName)
{
  if (this->HaveOneBodyBasis() == false)
    {
      cout << "error, the tight binding model does not provide the one body basis" << endl;
      return 0.0;
    }
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  Complex TmpChernNumber = 0.0;
  Complex Tmp1[4];
  Complex Tmp2[8];
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky    Berry_curvature\n";
  double Fluctations = 0.0;
  for (long LinearizedMomentumIndex = 0l; LinearizedMomentumIndex < this->NbrStatePerBand; ++LinearizedMomentumIndex)
    {
      int Kx;
      int Ky;
      this->GetLinearizedMomentumIndex(LinearizedMomentumIndex, Kx, Ky);
      int LinearizedMomentumIndexIncX = this->GetLinearizedMomentumIndex((Kx + 1) % this->NbrSiteX, Ky);
      int LinearizedMomentumIndexDecX;
      if (Kx > 0)
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex((Kx - 1) % this->NbrSiteX, Ky);
      else
	LinearizedMomentumIndexDecX = this->GetLinearizedMomentumIndex(this->NbrSiteX - 1, Ky);
      int LinearizedMomentumIndexIncY = this->GetLinearizedMomentumIndex(Kx, (Ky + 1) % this->NbrSiteY);
      int LinearizedMomentumIndexDecY;
      if (Ky > 0)
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, (Ky - 1) % this->NbrSiteY);
      else
	LinearizedMomentumIndexDecY = this->GetLinearizedMomentumIndex(Kx, this->NbrSiteY - 1);

      ComplexMatrix& LocalBasis = this->OneBodyBasis[LinearizedMomentumIndex];
      ComplexMatrix& LocalBasisIncX = this->OneBodyBasis[LinearizedMomentumIndexIncX];
      ComplexMatrix& LocalBasisDecX = this->OneBodyBasis[LinearizedMomentumIndexDecX];
      ComplexMatrix& LocalBasisIncY = this->OneBodyBasis[LinearizedMomentumIndexIncY];
      ComplexMatrix& LocalBasisDecY = this->OneBodyBasis[LinearizedMomentumIndexDecY];  
      Tmp1[0] = 0.0;
      Tmp1[1] = 0.0;
      Tmp1[2] = 0.0;
      Tmp1[3] = 0.0;

      Tmp2[0] = 0.0;
      Tmp2[1] = 0.0;
      Tmp2[2] = 0.0;
      Tmp2[3] = 0.0;
      Tmp2[4] = 0.0;
      Tmp2[5] = 0.0;
      Tmp2[6] = 0.0;
      Tmp2[7] = 0.0;

      for (int i = 0; i < this->NbrBands; ++i)
	{
	  Tmp1[0] += LocalBasis[band][i] * Conj(LocalBasisIncX[band][i]);
	  Tmp1[1] += LocalBasis[band][i] * Conj(LocalBasisDecX[band][i]);
	  Tmp1[2] += LocalBasis[band][i] * Conj(LocalBasisIncY[band][i]);
	  Tmp1[3] += LocalBasis[band][i] * Conj(LocalBasisDecY[band][i]);

	  Tmp2[0] += Conj(LocalBasisIncX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[1] += Conj(LocalBasisDecX[band][i]) * LocalBasisIncY[band][i];
	  Tmp2[2] += Conj(LocalBasisIncX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[3] += Conj(LocalBasisDecX[band][i]) * LocalBasisDecY[band][i];
	  Tmp2[4] += Conj(LocalBasisIncY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[5] += Conj(LocalBasisDecY[band][i]) * LocalBasisIncX[band][i];
	  Tmp2[6] += Conj(LocalBasisIncY[band][i]) * LocalBasisDecX[band][i];
	  Tmp2[7] += Conj(LocalBasisDecY[band][i]) * LocalBasisDecX[band][i];
	}

      Complex TmpCurvature = 0.0;
      TmpCurvature += (Tmp1[2] * Conj(Tmp1[0]) * Tmp2[0]);
      TmpCurvature -= (Tmp1[2] * Conj(Tmp1[1]) * Tmp2[1]);
      TmpCurvature -= (Tmp1[3] * Conj(Tmp1[0]) * Tmp2[2]);
      TmpCurvature += (Tmp1[3] * Conj(Tmp1[1]) * Tmp2[3]);
	  
      TmpCurvature -= (Tmp1[0] * Conj(Tmp1[2]) * Tmp2[4]);
      TmpCurvature += (Tmp1[0] * Conj(Tmp1[3]) * Tmp2[5]);
      TmpCurvature += (Tmp1[1] * Conj(Tmp1[2]) * Tmp2[6]);
      TmpCurvature -= (Tmp1[1] * Conj(Tmp1[3]) * Tmp2[7]);

      TmpCurvature *= 0.25;

      Fluctations += (TmpCurvature.Im - (2.0 * M_PI / ((double) this->NbrStatePerBand))) * (TmpCurvature.Im - (2.0 * M_PI / ((double) this->NbrStatePerBand)));

      File << Kx << " " << Ky << " " << TmpCurvature.Im << endl;

      TmpChernNumber += TmpCurvature;
    }
//  Fluctations *= ((double) this->NbrStatePerBand) ;
  cout << "Berry curvature fluctuations = " << Fluctations << " " <<  sqrt(Fluctations) << " " << (TmpChernNumber.Im * TmpChernNumber.Im) << " " << sqrt(Fluctations - (TmpChernNumber.Im * TmpChernNumber.Im)) 
       << "( " << (sqrt(Fluctations - (TmpChernNumber.Im * TmpChernNumber.Im)) / (2.0 * M_PI) )<< ") in 2 pi units" << endl;
  TmpChernNumber /= 2.0 * M_PI;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
  cout << "Chern number computed in  " << Dt << "s" << endl;
//  cout << "Berry curvature fluctuations = " << sqrt ()<< endl;

  File.close();

  return TmpChernNumber.Im;
}


// compute the complex eigenvalues of the D(ky) matrix (in order to compute the Z2 invariant)
//
// bandIndex = band index (corresponds to two bands that are related by time reversal symmetry)
// nbrOccupiedBands = dimension of the D matrix
// DMatrixEigenvalues = array of eigenvalues of the D Matrix, for all values of ky
// kyMin = minimal value of ky for which the D matrix has to be diagonalized
// kyMax = maximal value of ky for which the D matrix has to be diagonalized
// nbrKy = number of ky values for which the D matrix has to be diagonalized
// return value = array of eigenvalues of the D matrix

Complex** Abstract2DTightBindingModel::ComputeDMatrixEigenvalues(int nbrOccupiedBands, int kyMin, int kyMax, int nbrKy)
{
  Complex** DMatrixEigenvalues;
  DMatrixEigenvalues = new Complex*[this->NbrSiteY];
  for (int i = 0; i < this->NbrSiteY; ++i)
    DMatrixEigenvalues[i] = new Complex[nbrOccupiedBands];
  ComplexMatrix TmpDMatrix(nbrOccupiedBands, nbrOccupiedBands, true);
  ComplexMatrix FMatrix(nbrOccupiedBands, nbrOccupiedBands, true);
  
//   ComplexMatrix Rotation(this->NbrBands, this->NbrBands, true);
//   Rotation.SetMatrixElement(0, 0, M_SQRT1_2);
//   Rotation.SetMatrixElement(0, 1, Complex(0.0, M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 0, Complex(0.0, -1.0*M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 1, -1.0*M_SQRT1_2);
//   Rotation.SetMatrixElement(2, 2, 1.0);
//   Rotation.SetMatrixElement(3, 3, 1.0);
  
//   ComplexMatrix Rotation1(this->NbrBands, this->NbrBands, true);
//   Rotation.SetMatrixElement(0, 0, M_SQRT1_2);
//   Rotation.SetMatrixElement(0, 1, Complex(0.0, M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 0, Complex(0.0, -1.0*M_SQRT1_2));
//   Rotation.SetMatrixElement(1, 1, -1.0*M_SQRT1_2);
//   Rotation1.SetMatrixElement(2, 2, 1.0);
//   Rotation1.SetMatrixElement(3, 3, 1.0);
  
  
  for (int ky = kyMin; ky <= kyMax; ++ky)
  {
    TmpDMatrix.SetToIdentity();
    for (int i = 0; i < this->NbrSiteX; ++i)
    {
      double KX = (double) i *2.0 * M_PI / ((double) this->NbrSiteX);
      double KX1 = (double) (i + 1) *2.0 * M_PI / ((double) this->NbrSiteX);
//       cout << cos(KX) << " " << sin(KX) << endl;
//       Rotation.SetMatrixElement(0, 0, Complex(M_SQRT1_2*cos(KX), M_SQRT1_2*sin(KX)));
//       Rotation.SetMatrixElement(0, 1, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 0, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 1, Complex(-1.0*M_SQRT1_2*cos(KX), M_SQRT1_2*sin(KX)));
// 
//       
//       Rotation1.SetMatrixElement(0, 0, Complex(M_SQRT1_2*cos(KX1), M_SQRT1_2*sin(KX1)));
//       Rotation1.SetMatrixElement(0, 1, M_SQRT1_2);
//       Rotation1.SetMatrixElement(1, 0, M_SQRT1_2);
//       Rotation1.SetMatrixElement(1, 1, Complex(-1.0*M_SQRT1_2*cos(KX1), M_SQRT1_2*sin(KX1)));
//       Rotation.SetMatrixElement(0, 0, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 1, -1.0*M_SQRT1_2);
//       Rotation.SetMatrixElement(0, 1, M_SQRT1_2);
//       Rotation.SetMatrixElement(1, 0, M_SQRT1_2);
 
	int LinearizedMomentumIndex1 = this->GetLinearizedMomentumIndex(i, ky);
	int LinearizedMomentumIndex2 = this->GetLinearizedMomentumIndex((i + 1) % this->NbrSiteX, ky);
// 	cout << i << " " << LinearizedMomentumIndex1 << " " << LinearizedMomentumIndex2 << " " << endl;
	ComplexMatrix& LocalBasis = this->OneBodyBasis[LinearizedMomentumIndex1];
	ComplexMatrix& LocalBasisIncX = this->OneBodyBasis[LinearizedMomentumIndex2];
// 	ComplexMatrix LocalBasis(this->NbrBands, this->NbrBands, true);
// 	ComplexMatrix LocalBasisIncX(this->NbrBands, this->NbrBands, true);
	
// 	LocalBasis = TmpLocalBasis*Rotation;
// 	LocalBasisIncX = TmpLocalBasisIncX*Rotation1;
	
// 	LocalBasis = TmpLocalBasis;
// 	LocalBasisIncX = TmpLocalBasisIncX;
	
	for (int n = 0; n < nbrOccupiedBands; ++n)
	  {
	    for (int m = 0; m < nbrOccupiedBands; ++m)
	      {
// 		Complex Tmp = 0.0;
// 		for (int alpha = 0; alpha < this->NbrBands; ++alpha)
// 		{
// 		  Tmp += Conj(LocalBasis[n][alpha]) * LocalBasisIncX[m][alpha];
// 		}
// 		FMatrix.SetMatrixElement(n, m, Tmp);
		FMatrix.SetMatrixElement(n, m, LocalBasis[n] * LocalBasisIncX[m]);
	      }
	  }
// 	  cout << i << endl;
// 	  cout << FMatrix << endl;
// 	  ComplexMatrix TmpMatrix = TmpDMatrix;
// 	  TmpDMatrix = TmpMatrix*FMatrix;
	  TmpDMatrix.Multiply(FMatrix);
      }
    
// 	  if (ky == 0)
// 	  {
// 	   cout << "ky = 0" << endl;
// 	   cout << TmpDMatrix << endl; 
// 	  }
    ComplexDiagonalMatrix TmpDiag(nbrOccupiedBands);
#ifdef __LAPACK__
    TmpDMatrix.LapackDiagonalize(TmpDiag);
#else
    TmpDMatrix.Diagonalize(TmpDiag);
#endif
//     cout << TmpDiag << endl;
    
    for (int j = 0; j < nbrOccupiedBands; ++j)
    {
      DMatrixEigenvalues[ky][j] = TmpDiag[j] ;
//       cout << DMatrixEigenvalues[ky][j] << endl;
    }
  }
  
  return DMatrixEigenvalues;
}

// compute the Z2 topological invariant for a system with time reversal symmetry
//
// nbrOccupiedBands = number of occupied bands
// return value = Z2 invariant

int Abstract2DTightBindingModel::ComputeZ2Invariant(int nbrOccupiedBands)
{
  int z2Invariant = 0;
  double referenceLine = 0.9267;
  
  double distancePlus;
  double distanceMoins;
  double distanceMod2PiPlus;
  double distanceMod2PiMoins;
  
  int ModPiPlus = 0;
  int ModPiMoins = 0;
  
  double theta1;
  double theta2;
  
  Complex** Lambda = this->ComputeDMatrixEigenvalues(nbrOccupiedBands, 0, this->NbrSiteY - 1, this->NbrSiteY); 
  double** Theta = new double*[this->NbrSiteY];
  for (int ky = 0; ky < this->NbrSiteY; ++ky)
  {
   Theta[ky] = new double[2];
   for (int i = 0; i < 2; ++i)
   {
     theta1 = atan2(Lambda[ky][nbrOccupiedBands - 2].Im,Lambda[ky][nbrOccupiedBands - 2].Re);
     theta2 = atan2(Lambda[ky][nbrOccupiedBands - 1].Im,Lambda[ky][nbrOccupiedBands - 1].Re);
     Theta[ky][0] = max(theta1, theta2);
     Theta[ky][1] = min(theta1, theta2); 
   }
  }
//   double referenceLine = 0.9267 + theta[0][0];
  for (int ky = 0; ky < this->NbrSiteY  - 1; ++ ky)
  {
    distancePlus = fabs(Theta[ky][0] - Theta[ky + 1][0]);
    distanceMod2PiPlus = fabs(Theta[ky][0] - Theta[ky + 1][1] - 2*M_PI);
    distanceMoins = fabs(Theta[ky][1] - Theta[ky + 1][1]);
    distanceMod2PiMoins = fabs(Theta[ky][1] - Theta[ky + 1][0] + 2*M_PI);
    
    if (distanceMod2PiPlus < distancePlus)
    {
     ModPiPlus += 1;
     double Tmp = Theta[ky + 1][0];
     Theta[ky + 1][0] = Theta[ky + 1][1] + 2*M_PI;
     Theta[ky + 1][1] = Tmp;
    }
    
    if (distanceMod2PiMoins < distanceMoins)
    {
     ModPiMoins += 1;
     double Tmp = Theta[ky + 1][1];
     Theta[ky + 1][1] = Theta[ky + 1][0] - 2*M_PI;
     Theta[ky + 1][0] = Tmp;
    }
  }
//   cout << ModPiPlus << " " << ModPiMoins << endl;
  for (int ky = 0 ; ky < this->NbrSiteY/2 ; ++ky)
  {
//     cout << ky << " " << Theta[ky][0] << " " << Theta[ky][1] << endl;
    for (int i = 0; i <= ModPiPlus; ++i)
    {
//       cout << Theta[ky + 1][0] - (referenceLine + i*2*M_PI) << " " << Theta[ky][0] - (referenceLine + i*2*M_PI) << " " << (Theta[ky + 1][1] - (referenceLine - i*2*M_PI)) << " " << (Theta[ky][1] - (referenceLine - i*2*M_PI)) << endl;
      if ((Theta[ky + 1][0] - (referenceLine + i*2*M_PI)) * (Theta[ky][0] - (referenceLine + i*2*M_PI)) < 0)
      {
	z2Invariant += 1;
      }
    }
    for (int i = 0; i <= ModPiMoins; ++i)
    {
      if ((Theta[ky + 1][1] - (referenceLine - i*2*M_PI)) * (Theta[ky][1] - (referenceLine - i*2*M_PI)) < 0)
      {
	z2Invariant += 1;
      }
    }
  }
  return (z2Invariant % 2); 
}

//computes all the values of the projected momentum and stores them in a double array
//

void Abstract2DTightBindingModel::ComputeAllProjectedMomenta()
{
  this->ProjectedMomenta = new double* [this->NbrStatePerBand];
  for (int i = 0; i < this->NbrStatePerBand; ++i)
    this->ProjectedMomenta[i] = new double [2];
  double projectedMomentum1;
  double projectedMomentum2;
  if (this->NbrConnectedOrbitals != 0)
    {
      for (int kx = 0; kx < this->NbrSiteX; ++kx)
	{
	  for (int ky = 0; ky < this->NbrSiteY; ++ky)
	    {
	      double kx_trans = kx + this->Offset * ky;
	      double ky_trans = ky;
	      projectedMomentum1 = 2.0 * M_PI * ((double) kx_trans * (double) this->Ny2 - (double) ky_trans * (double) this->Ny1) / ((double) (this->NbrSiteX * this->NbrSiteY));
	      projectedMomentum2 = 2.0 * M_PI * ((double) kx_trans * (double) (-this->Nx2) + (double) ky_trans * (double)this->Nx1) / ((double) (this->NbrSiteX * this->NbrSiteY));
	      this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky)][0] = projectedMomentum1;
	      this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky)][1] = projectedMomentum2;
	    }
	}
    }
  else
    {
      for (int kx = 0; kx < this->NbrSiteX; ++kx)
	{
	  for (int ky = 0; ky < this->NbrSiteY; ++ky)
	    {
	      double kx_trans = kx + this->Offset * ky + this->GammaX;
	      double ky_trans = ky + this->GammaY;
	      projectedMomentum1 = 2.0 * M_PI * ((double) kx_trans * (double) this->Ny2 - (double) ky_trans * (double) this->Ny1) / ((double) (this->NbrSiteX * this->NbrSiteY));
	      projectedMomentum2 = 2.0 * M_PI * ((double) kx_trans * (double) (-this->Nx2) + (double) ky_trans * (double)this->Nx1) / ((double) (this->NbrSiteX * this->NbrSiteY));
	      this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky)][0] = projectedMomentum1;
	      this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky)][1] = projectedMomentum2;
	    }
	}
    }
}

//set the value of the embedding vectors from an external ascii file
//
//embeddingFileName = name of the ascii file that defines the embedding
//
bool Abstract2DTightBindingModel::SetEmbeddingFromAsciiFile(char* embeddingFileName)
{
  MultiColumnASCIIFile Parser;
  Parser.Parse(embeddingFileName);
//   cout << Parser.GetNbrColumns() <<  "  " << Parser.GetNbrLines() << endl;
  if (Parser.GetNbrColumns() != 2 || Parser.GetNbrLines() != this->NbrBands)
  {
    cout << "Embedding file should have " << this->NbrBands << " lines and 2 columns" << endl;
    return false;
  }
  else
  {
    double** coefficients = new double*[2];
    for (int i = 0; i < 2; ++i)
      coefficients[i] = Parser.GetAsDoubleArray(i);

    for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EmbeddingX[i] = coefficients[0][i];
//       cout << this->EmbeddingX[i] << endl;
      this->EmbeddingY[i] = coefficients[1][i];
//       cout << this->EmbeddingY[i] << endl;
    }
    delete[] coefficients;
    return true;
  }
 
}

// build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
//
// nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
// orbitalIndices = array that gives the orbital indices of the connected orbitals
// spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
// hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
// return value = tight binding hamiltonian in real space 

HermitianMatrix Abstract2DTightBindingModel::BuildTightBindingHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes)
{
  HermitianMatrix TmpHamiltonian(this->NbrBands * this->NbrSiteX * this->NbrSiteY, true);
  for (int i = 0; i < this->NbrSiteX; ++i)
    {
      for (int j = 0; j < this->NbrSiteY; ++j)
	{
	  for (int k = 0; k < this->NbrBands; ++k)
	    {
	      int Index2 = this->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, k);
	      for (int l = 0; l < nbrConnectedOrbitals[k]; ++l)
		{
		  int Index1 = this->GetRealSpaceTightBindingLinearizedIndexSafe(spatialIndices[k][l << 1] + i, spatialIndices[k][(l << 1) + 1] + j, orbitalIndices[k][l]);
		  if(Index1 >= Index2)
		    TmpHamiltonian.AddToMatrixElement(Index1, Index2, hoppingAmplitudes[k][l]);
		}
	    }
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

ComplexMatrix Abstract2DTightBindingModel::BuildTightBindingNonHermitianHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes)
{
  ComplexMatrix TmpHamiltonian(this->NbrBands * this->NbrSiteX * this->NbrSiteY, this->NbrBands * this->NbrSiteX * this->NbrSiteY, true);
  for (int i = 0; i < this->NbrSiteX; ++i)
    {
      for (int j = 0; j < this->NbrSiteY; ++j)
	{
	  for (int k = 0; k < this->NbrBands; ++k)
	    {
	      int Index2 = this->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, k);
	      for (int l = 0; l < nbrConnectedOrbitals[k]; ++l)
		{
		  int Index1 = this->GetRealSpaceTightBindingLinearizedIndexSafe(spatialIndices[k][l << 1] + i, spatialIndices[k][(l << 1) + 1] + j, orbitalIndices[k][l]);
		  TmpHamiltonian.AddToMatrixElement(Index1, Index2, hoppingAmplitudes[k][l]);
		}
	    }
	}      
    }
  return TmpHamiltonian;
}

// build the tight binding hamiltonian in recirpocal space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
//
// kx = momentum along the x direction (in 2pi /N_x unit) for which the hamiltonian in recirpocal space has to be computed
// ky = momentum along the y direction (in 2pi /N_y unit) for which the hamiltonian in recirpocal space has to be computed
// nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
// orbitalIndices = array that gives the orbital indices of the connected orbitals
// spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
// hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
// return value = tight binding hamiltonian in real space 

HermitianMatrix Abstract2DTightBindingModel::BuildTightBindingHamiltonianReciprocalSpace(int kx, int ky, int* nbrConnectedOrbitals, int** orbitalIndices, 
											 int** spatialIndices, Complex** hoppingAmplitudes)
{
  HermitianMatrix TmpHamiltonian(this->NbrBands, true);
  double TmpKx = this->GetProjectedMomentum(kx, ky, 0);
  double TmpKy = this->GetProjectedMomentum(kx, ky, 1);
  int p;
  int q;
  
  for (int k = 0; k < this->NbrBands; ++k)
    {
      for (int l = 0; l < nbrConnectedOrbitals[k]; ++l)
	{
	  this->GetRealSpaceIndex(spatialIndices[k][l << 1], spatialIndices[k][(l << 1) + 1], p, q);
	  double TmpPhase = ((TmpKx * ((double) p)) 
			      + (TmpKy * ((double) q)));
	  if (k >= orbitalIndices[k][l])
	    {
	      TmpHamiltonian.AddToMatrixElement(k, orbitalIndices[k][l], Conj(hoppingAmplitudes[k][l]) * Phase(TmpPhase));
	    }
	}
    }
  return TmpHamiltonian; 
}


// generate a tight-binding Density-Density interaction in real space for the current Tight-Binding Model
// nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsPotentials = intensity of each density-density term 
RealSymmetricMatrix Abstract2DTightBindingModel::GenerateDensityDensityInteraction(int *NbrInteractingOrbitals, int **InteractingOrbitalsOrbitalIndices, int **InteractingOrbitalsSpatialIndices,  double **InteractingOrbitalsPotentials)
{
  RealSymmetricMatrix interaction(this->GetNbrBands() * this->GetNbrStatePerBand(), true);
  for (int x = 0; x < this->NbrSiteX; ++x)
    {
      for (int y = 0; y < this->NbrSiteY; ++y)
	{
	  for (int OrbitalIndex = 0; OrbitalIndex < this->GetNbrBands(); ++OrbitalIndex)
	    {
	      for (int k = 0; k < NbrInteractingOrbitals[OrbitalIndex]; ++k)
		{
		  interaction.AddToMatrixElement(this->GetRealSpaceTightBindingLinearizedIndexSafe(x, y, OrbitalIndex), 
							       this->GetRealSpaceTightBindingLinearizedIndexSafe(x + InteractingOrbitalsSpatialIndices[OrbitalIndex][2 * k], 
														 y + InteractingOrbitalsSpatialIndices[OrbitalIndex][(2 * k) + 1], 
														 InteractingOrbitalsOrbitalIndices[OrbitalIndex][k]), 
							       InteractingOrbitalsPotentials[OrbitalIndex][k]);
				  
		}
	    }
	}
    }
  return interaction;
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void Abstract2DTightBindingModel::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  this->FindConnectedOrbitals();
  if (this->NbrConnectedOrbitals != 0)
    {
      if (nbrStates == 0l)
	nbrStates = this->NbrStatePerBand;
      long MaxStateIndex = minStateIndex + nbrStates;
      double KX;
      double KY;
      for (int kx = 0; kx < this->NbrSiteX; ++kx)
	{
	  for (int ky = 0; ky < this->NbrSiteY; ++ky)
	    {
	      int Index = this->GetLinearizedMomentumIndex(kx, ky);
	      if ((Index >= minStateIndex) && (Index < MaxStateIndex))
		{
		  HermitianMatrix TmpOneBodyHamiltonian = this->BuildTightBindingHamiltonianReciprocalSpace(kx, ky, this->NbrConnectedOrbitals, this->ConnectedOrbitalIndices,
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
		      this->OneBodyBasis[Index] = TmpMatrix;
		      for (int i = 0; i < this->NbrBands; ++i)
		      {
			this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
// 			cout << i << " " << Index << " " << TmpDiag(i,i) << endl;
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
			this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		    }
		}
	    }
	}
    }
}

// compute the band structure at a single point of the Brillouin zone
//
// kx = momentum along the x axis
// ky = momentum along the x axis
// energies = array where the energies will be stored

void Abstract2DTightBindingModel::ComputeBandStructureSinglePoint(double kx, double ky, double* energies)
{
  HermitianMatrix TmpOneBodyHamiltonian = this->ComputeBlochHamiltonian(kx, ky);
  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
  TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif
  for (int i = 0; i < this->NbrBands; ++i)
    energies[i] = TmpDiag(i, i);
}

// compute the Bloch hamiltonian at a point of the Brillouin zone
//
// kx = momentum along the x axis
// ky = momentum along the x axis
// return value = Bloch hamiltonian

HermitianMatrix Abstract2DTightBindingModel::ComputeBlochHamiltonian(double kx, double ky)
{
  HermitianMatrix TmpOneBodyHamiltonian;
  return TmpOneBodyHamiltonian;
}

// get the high symmetry points 
//
// pointNames = name of each high symmetry point
// pointCoordinates = coordinates in the fist Brillouin zone of the high symmetry points
// return value = number of high symmetry points

int Abstract2DTightBindingModel::GetHighSymmetryPoints(char**& pointNames, double**& pointCoordinates)
{
  pointNames = 0;
  pointCoordinates = 0;
  return 0;
}

// compute the distance between two points in the first Brillouin zone, changing the coordinates the second one by a reciprocal lattice vector if needed
//
// kx1 = momentum of the first point along the x axis
// ky1 = momentum of the first point along the y axis
// kx2 = reference on the momentum of the second point along the x axis
// ky2 = reference on the momentum of the second point along the y axis
// return value = distance between the two points

double Abstract2DTightBindingModel::GetDistanceReciprocalSpace(double kx1, double ky1, double& kx2, double& ky2)
{
  double MinDistance = sqrt (((kx1 - kx2) * (kx1 - kx2)) + ((ky1 - ky2) * (ky1 - ky2)));
  double MinKx2 = kx2;
  double MinKy2 = ky2;
  for (int i = -1; i <= 1; ++i)
    {
      double TmpKx2  = kx2 + (2.0 * ((double) i) * M_PI);
      for (int j = -1; j <= 1; ++j)
	{
	  double TmpKy2  = ky2 + (2.0 * ((double) j) * M_PI);	  
	  double TmpDistance = sqrt (((kx1 - TmpKx2) * (kx1 - TmpKx2)) + ((ky1 - TmpKy2) * (ky1 - TmpKy2)));
	  if (TmpDistance < MinDistance)
	    {
	      MinDistance = TmpDistance;
	      MinKx2 = TmpKx2;
	      MinKy2 = TmpKy2;
	    }
	}
    }
  kx2 = MinKx2;
  ky2 = MinKy2;
  return MinDistance;
}

// evaluate the two point correlation function 
//
// x = linearized position index of the first point
// y = linearized position index of the second point
// occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
// nbrOccupiedMomenta = number of occupied momenta
// bandIndex = index of the band to consider
// return value = value of the two point correlation function 

Complex Abstract2DTightBindingModel::EvaluateTwoPointCorrelationFunction(int x, int y, int* occupiedMomenta, int nbrOccupiedMomenta, int bandIndex)
{
  int TmpXx;
  int TmpXy;
  int TmpXOrbital;
  this->GetRealSpaceTightBindingLinearizedIndex(x, TmpXx, TmpXy, TmpXOrbital);
  int TmpYx;
  int TmpYy;
  int TmpYOrbital;
  this->GetRealSpaceTightBindingLinearizedIndex(y, TmpYx, TmpYy, TmpYOrbital);
  Complex Tmp = 0.0;
  int TmpMomentumX;
  int TmpMomentumY;
  for (int i = 0; i < nbrOccupiedMomenta; ++i)
    {
      this->GetLinearizedMomentumIndex(occupiedMomenta[i], TmpMomentumX, TmpMomentumY);
      Tmp += (Phase((this->KxFactor * ((double) (TmpMomentumX * (TmpXx - TmpYx))))
		    + (this->KyFactor * ((double) (TmpMomentumY * (TmpXy - TmpYy)))))
	      * Conj(this->OneBodyBasis[occupiedMomenta[i]][bandIndex][TmpXOrbital]) * this->OneBodyBasis[occupiedMomenta[i]][bandIndex][TmpYOrbital]);
    }
  Tmp /= ((double) (this->NbrSiteX * this->NbrSiteY));
  return Tmp;
}

// evaluate the two point correlation function in a given region
//
// maxX = x coordinate of the region upper right corner 
// maxY = y coordinate of the region upper right corner 
// occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
// nbrOccupiedMomenta = number of occupied momenta
// bandIndex = index of the band to consider
// return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)

HermitianMatrix Abstract2DTightBindingModel::EvaluateFullTwoPointCorrelationFunction(int maxX, int maxY, int* occupiedMomenta, int nbrOccupiedMomenta, int bandIndex)
{
  int TotalNbrSites = maxX * maxY * this->NbrBands;
  int TmpMomentumX;
  int TmpMomentumY;
  HermitianMatrix EntanglementHamiltonian(TotalNbrSites, true);

  int TmpMaxX2 = 2 * maxX + 1;
  Complex** TmpPhaseFactorX = new Complex*[TmpMaxX2];
  for (int i = 0; i < TmpMaxX2; ++i)
    {
      TmpPhaseFactorX[i] = new Complex[nbrOccupiedMomenta];
      for (int j = 0; j < nbrOccupiedMomenta; ++j)
	{
	  this->GetLinearizedMomentumIndex(occupiedMomenta[j], TmpMomentumX, TmpMomentumY);
	  TmpPhaseFactorX[i][j] = Phase(this->KxFactor * ((double) (TmpMomentumX * (maxX - i))));
	}
    }
  int TmpMaxY2 = 2 * maxY + 1;
  Complex** TmpPhaseFactorY = new Complex*[TmpMaxY2];
  for (int i = 0; i < TmpMaxY2; ++i)
    {
      TmpPhaseFactorY[i] = new Complex[nbrOccupiedMomenta];
      for (int j = 0; j < nbrOccupiedMomenta; ++j)
	{
	  this->GetLinearizedMomentumIndex(occupiedMomenta[j], TmpMomentumX, TmpMomentumY);
	  TmpPhaseFactorY[i][j] = Phase(this->KyFactor * ((double) (TmpMomentumY * (maxY - i))));
	}
    }
  
  Complex*** TmpFormFactors = new Complex** [nbrOccupiedMomenta];
  for (int i = 0; i < nbrOccupiedMomenta; ++i)
    {
      TmpFormFactors[i] = new Complex*[this->NbrBands];
      for (int j = 0; j < this->NbrBands; ++j)
	{
	  TmpFormFactors[i][j] = new Complex[this->NbrBands];
	  for (int k = 0; k < this->NbrBands; ++k)
	    {	      
	      TmpFormFactors[i][j][k] = (Conj(this->OneBodyBasis[occupiedMomenta[i]][bandIndex][j])
					 * this->OneBodyBasis[occupiedMomenta[i]][bandIndex][k]);
	    }
	}
    }

  for (int TmpX1 = 0; TmpX1 < maxX; ++TmpX1)
    {
      for (int TmpY1 = 0; TmpY1 < maxY; ++TmpY1)
	{	  
	  for (int TmpOrbital1 = 0; TmpOrbital1 < this->NbrBands; ++TmpOrbital1)
	    {
	      int TmpLinearizedIndex1 = this->GetRealSpaceTightBindingLinearizedIndex(TmpX1, TmpY1, TmpOrbital1);
	      int TmpReducedLinearizedIndex1 = TmpOrbital1 + ((TmpY1  + TmpX1 * maxY) * this->NbrBands);
	      for (int TmpX2 = 0; TmpX2 < maxX; ++TmpX2)
		{
		  for (int TmpY2 = 0; TmpY2 < maxY; ++TmpY2)
		    {	  
		      for (int TmpOrbital2 = 0; TmpOrbital2 < this->NbrBands; ++TmpOrbital2)
			{
			  int TmpLinearizedIndex2 = this->GetRealSpaceTightBindingLinearizedIndex(TmpX2, TmpY2, TmpOrbital2);
			  int TmpReducedLinearizedIndex2 = TmpOrbital2 + ((TmpY2  + TmpX2 * maxY) * this->NbrBands);
			  if (TmpReducedLinearizedIndex1 <= TmpReducedLinearizedIndex2)		  
			    {
			      Complex Tmp = 0.0;
			      Complex* DiffX = TmpPhaseFactorX[TmpX2 - TmpX1 + maxX];
			      Complex* DiffY = TmpPhaseFactorY[TmpY2 - TmpY1 + maxY];
			      for (int i = 0; i < nbrOccupiedMomenta; ++i)
				{
				  Tmp += (DiffX[i] * DiffY[i] * TmpFormFactors[i][TmpOrbital1][TmpOrbital2]);
				}
			      EntanglementHamiltonian.SetMatrixElement(TmpReducedLinearizedIndex1, TmpReducedLinearizedIndex2, Tmp);
			    }
			}
		    }
		}
	    }
	}
    }

  for (int i = 0; i < TmpMaxX2; ++i)
    delete[] TmpPhaseFactorX[i];
  delete[] TmpPhaseFactorX;
  for (int i = 0; i < TmpMaxY2; ++i)
    delete[] TmpPhaseFactorY[i];
  delete[] TmpPhaseFactorY;
  for (int i = 0; i < nbrOccupiedMomenta; ++i)
    {
      for (int j = 0; j < this->NbrBands; ++j)
	delete[] TmpFormFactors[i][j];
      delete[] TmpFormFactors[i];
    }
  delete[] TmpFormFactors;
  EntanglementHamiltonian /= ((double) (this->NbrSiteX * this->NbrSiteY));
  return EntanglementHamiltonian;
}

// evaluate the two point correlation function in a given region
//
// maxX = x coordinate of the region upper right corner 
// maxY = y coordinate of the region upper right corner 
// occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
// bandIndices = indices of the band corresponding ot each occupied state
// nbrOccupiedStates = number of occupied states
// bandIndex = index of the band to consider
// return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)

HermitianMatrix Abstract2DTightBindingModel::EvaluateFullTwoPointCorrelationFunction(int maxX, int maxY, int* occupiedMomenta, int* bandIndices, int nbrOccupiedStates)
{
  int TotalNbrSites = maxX * maxY * this->NbrBands;
  int TmpMomentumX;
  int TmpMomentumY;
  HermitianMatrix EntanglementHamiltonian(TotalNbrSites, true);

  int TmpMaxX2 = 2 * maxX + 1;
  Complex** TmpPhaseFactorX = new Complex*[TmpMaxX2];
  for (int i = 0; i < TmpMaxX2; ++i)
    {
      TmpPhaseFactorX[i] = new Complex[nbrOccupiedStates];
      for (int j = 0; j < nbrOccupiedStates; ++j)
	{
	  this->GetLinearizedMomentumIndex(occupiedMomenta[j], TmpMomentumX, TmpMomentumY);
	  TmpPhaseFactorX[i][j] = Phase(this->KxFactor * ((double) (TmpMomentumX * (maxX - i))));
	}
    }
  int TmpMaxY2 = 2 * maxY + 1;
  Complex** TmpPhaseFactorY = new Complex*[TmpMaxY2];
  for (int i = 0; i < TmpMaxY2; ++i)
    {
      TmpPhaseFactorY[i] = new Complex[nbrOccupiedStates];
      for (int j = 0; j < nbrOccupiedStates; ++j)
	{
	  this->GetLinearizedMomentumIndex(occupiedMomenta[j], TmpMomentumX, TmpMomentumY);
	  TmpPhaseFactorY[i][j] = Phase(this->KyFactor * ((double) (TmpMomentumY * (maxY - i))));
	}
    }
  
  Complex*** TmpFormFactors = new Complex** [nbrOccupiedStates];
  for (int i = 0; i < nbrOccupiedStates; ++i)
    {
      TmpFormFactors[i] = new Complex*[this->NbrBands];
      for (int j = 0; j < this->NbrBands; ++j)
	{
	  TmpFormFactors[i][j] = new Complex[this->NbrBands];
	  for (int k = 0; k < this->NbrBands; ++k)
	    {	      
	      TmpFormFactors[i][j][k] = (Conj(this->OneBodyBasis[occupiedMomenta[i]][bandIndices[i]][j])
					 * this->OneBodyBasis[occupiedMomenta[i]][bandIndices[i]][k]);
	    }
	}
    }

  for (int TmpX1 = 0; TmpX1 < maxX; ++TmpX1)
    {
      for (int TmpY1 = 0; TmpY1 < maxY; ++TmpY1)
	{	  
	  for (int TmpOrbital1 = 0; TmpOrbital1 < this->NbrBands; ++TmpOrbital1)
	    {
	      int TmpLinearizedIndex1 = this->GetRealSpaceTightBindingLinearizedIndex(TmpX1, TmpY1, TmpOrbital1);
	      int TmpReducedLinearizedIndex1 = TmpOrbital1 + ((TmpY1  + TmpX1 * maxY) * this->NbrBands);
	      for (int TmpX2 = 0; TmpX2 < maxX; ++TmpX2)
		{
		  for (int TmpY2 = 0; TmpY2 < maxY; ++TmpY2)
		    {	  
		      for (int TmpOrbital2 = 0; TmpOrbital2 < this->NbrBands; ++TmpOrbital2)
			{
			  int TmpLinearizedIndex2 = this->GetRealSpaceTightBindingLinearizedIndex(TmpX2, TmpY2, TmpOrbital2);
			  int TmpReducedLinearizedIndex2 = TmpOrbital2 + ((TmpY2  + TmpX2 * maxY) * this->NbrBands);
			  if (TmpReducedLinearizedIndex1 <= TmpReducedLinearizedIndex2)		  
			    {
			      Complex Tmp = 0.0;
			      Complex* DiffX = TmpPhaseFactorX[TmpX2 - TmpX1 + maxX];
			      Complex* DiffY = TmpPhaseFactorY[TmpY2 - TmpY1 + maxY];
			      for (int i = 0; i < nbrOccupiedStates; ++i)
				{
				  Tmp += (DiffX[i] * DiffY[i] * TmpFormFactors[i][TmpOrbital1][TmpOrbital2]);
				}
			      EntanglementHamiltonian.SetMatrixElement(TmpReducedLinearizedIndex1, TmpReducedLinearizedIndex2, Tmp);
			    }
			}
		    }
		}
	    }
	}
    }

  for (int i = 0; i < TmpMaxX2; ++i)
    delete[] TmpPhaseFactorX[i];
  delete[] TmpPhaseFactorX;
  for (int i = 0; i < TmpMaxY2; ++i)
    delete[] TmpPhaseFactorY[i];
  delete[] TmpPhaseFactorY;
  for (int i = 0; i < nbrOccupiedStates; ++i)
    {
      for (int j = 0; j < this->NbrBands; ++j)
	delete[] TmpFormFactors[i][j];
      delete[] TmpFormFactors[i];
    }
  delete[] TmpFormFactors;
  EntanglementHamiltonian /= ((double) (this->NbrSiteX * this->NbrSiteY));
  return EntanglementHamiltonian;
}

// evaluate the mixed two point correlation function in a given region, assuming translation invariance along one direction
//
// maxX = length along the borken translation direction of the region 
// ky = momentum along the translation invariant direction
// occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
// bandIndices = array that gives the band index of each occupied state
// nbrOccupiedMomenta = number of occupied momenta
// return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)

HermitianMatrix Abstract2DTightBindingModel::EvaluateFullMixedTwoPointCorrelationFunctionWithK(int maxX, int ky, int* occupiedMomenta, int* bandIndices, int nbrOccupiedMomenta)
{
  int TotalNbrSites = maxX * this->NbrBands;
  int TmpMomentumX;
  int TmpMomentumY;
  HermitianMatrix EntanglementHamiltonian(TotalNbrSites, true);

  int TmpNbrKeptStates = 0;
  for (int j = 0; j < nbrOccupiedMomenta; ++j)
    {
      this->GetLinearizedMomentumIndex(occupiedMomenta[j], TmpMomentumX, TmpMomentumY);
      if (TmpMomentumY == ky)
	{
	  ++TmpNbrKeptStates;
	}
    }
  if (TmpNbrKeptStates == 0)
    {
      return EntanglementHamiltonian;
    }
  int* TmpKeptStates = new int[TmpNbrKeptStates];
  TmpNbrKeptStates = 0;
  for (int j = 0; j < nbrOccupiedMomenta; ++j)
   {
     this->GetLinearizedMomentumIndex(occupiedMomenta[j], TmpMomentumX, TmpMomentumY);
     if (TmpMomentumY == ky)
	{
	  TmpKeptStates[TmpNbrKeptStates] = j;
	  ++TmpNbrKeptStates;
	}
   }

  int TmpMaxX2 = 2 * maxX + 1;
  Complex** TmpPhaseFactorX = new Complex*[TmpMaxX2];
  for (int i = 0; i < TmpMaxX2; ++i)
    {
      TmpPhaseFactorX[i] = new Complex[TmpNbrKeptStates];
      for (int j = 0; j < TmpNbrKeptStates; ++j)
	{
	  this->GetLinearizedMomentumIndex(occupiedMomenta[TmpKeptStates[j]], TmpMomentumX, TmpMomentumY);
	  TmpPhaseFactorX[i][j] = Phase(this->KxFactor * ((double) (TmpMomentumX * (maxX - i))));
	}
    }
  
  
  Complex*** TmpFormFactors = new Complex** [TmpNbrKeptStates];
  for (int i = 0; i < TmpNbrKeptStates; ++i)
    {
      TmpFormFactors[i] = new Complex*[this->NbrBands];
      for (int j = 0; j < this->NbrBands; ++j)
	{
	  TmpFormFactors[i][j] = new Complex[this->NbrBands];
	  for (int k = 0; k < this->NbrBands; ++k)
	    {	      
	      TmpFormFactors[i][j][k] = (Conj(this->OneBodyBasis[occupiedMomenta[TmpKeptStates[i]]][bandIndices[TmpKeptStates[i]]][j])
					 * this->OneBodyBasis[occupiedMomenta[TmpKeptStates[i]]][bandIndices[TmpKeptStates[i]]][k]);
	    }
	}
    }


  for (int TmpX1 = 0; TmpX1 < maxX; ++TmpX1)
    {
      for (int TmpOrbital1 = 0; TmpOrbital1 < this->NbrBands; ++TmpOrbital1)
	{
	  int TmpReducedLinearizedIndex1 = TmpOrbital1 + (TmpX1 * this->NbrBands);
	  for (int TmpX2 = 0; TmpX2 < maxX; ++TmpX2)
	    {
	      for (int TmpOrbital2 = 0; TmpOrbital2 < this->NbrBands; ++TmpOrbital2)
		{
		  int TmpReducedLinearizedIndex2 = TmpOrbital2 + (TmpX2 * this->NbrBands);
		  if (TmpReducedLinearizedIndex1 <= TmpReducedLinearizedIndex2)		  
		    {
		      Complex Tmp = 0.0;
		      Complex* DiffX = TmpPhaseFactorX[TmpX2 - TmpX1 + maxX];
		      for (int i = 0; i < TmpNbrKeptStates; ++i)
			{
			  Tmp += (DiffX[i] * TmpFormFactors[i][TmpOrbital1][TmpOrbital2]);
			}
		      EntanglementHamiltonian.SetMatrixElement(TmpReducedLinearizedIndex1, TmpReducedLinearizedIndex2, Tmp);
		    }
		}
	    }
	}
    }

  for (int i = 0; i < TmpMaxX2; ++i)
    delete[] TmpPhaseFactorX[i];
  delete[] TmpPhaseFactorX;
  for (int i = 0; i < TmpNbrKeptStates; ++i)
    {
      for (int j = 0; j < this->NbrBands; ++j)
	delete[] TmpFormFactors[i][j];
      delete[] TmpFormFactors[i];
    }
  delete[] TmpFormFactors;
  EntanglementHamiltonian /= ((double) (this->NbrSiteX));
  return EntanglementHamiltonian;
}


// compute the form factor for the density operator 
// 
// kx = momentum along x of annihilation operator
// ky = momentum along y of creation operator
// qx = momentum transfer along x direction
// qy = momentum transfer along y direction
// valleyIndex = valley index of density operator

Complex Abstract2DTightBindingModel::ComputeDensityFormFactor(int kx, int ky, int qx, int qy, int valleyIndex)
{
    cout << "Warning: using dummy method Abstract2DTightBindingModel::ComputeDensityFormFactor" << endl;
    return 0.0;
}
