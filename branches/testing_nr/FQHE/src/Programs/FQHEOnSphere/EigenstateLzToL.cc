#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Tools/FQHESpectrum/QHEOnSphereLzSortedSpectrum.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


// subdivide a basis into different L orthonormal eigenstate
//
// vectors = reference on the matrix whose columns are the vectors that span the basis
// oper = pointer to operator whivh allow to evaluate total L matrix elements
// totalMaxLz = maximum L value that can be reach by the system (-1/2 if L is half integer)
// subspaceSize = reference on the array that will be filled  with the dimension of each fixed L subspace (index equal to L if L is integer, L-1/2 if L is half integer)
// subspacePositions = reference on the array that will be filled  with the position of the first occurence of a vector with a given fixed L subspace 
//                     (index equal to L if L is integer, L-1/2 if L is half integer)
// lError = allowed error on L determination, if relative error between <L> and int(<L>) is greater than lError, no  L orthonormal eigenstates van be extracted
// vectorNames = array that contains file name of each vectors (used to identify vector in error message, 0 if vector index has to be used instead)
void LSortBasis(RealMatrix& vectors, AbstractOperator* oper, int totalMaxLz, int* subspaceSize, int* subspacePositions, double lError, char** vectorNames = 0);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("EigenstateLzToL" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* DataGroup = new OptionGroup ("data options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += DataGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz", "twice the total lz value of the system", 0);
  (*SystemGroup) += new SingleStringOption  ('s', "statistics", "particle statistics (boson or fermion, try to guess it from file name if not defined)");

  (*DataGroup) += new SingleStringOption  ('f', "file-prefix", "prefix for all data files (vectors have to be of the form prefix_x.y.vec and Lz spectrum prefix.dat)");
  (*DataGroup) += new SingleDoubleOption  ('\n', "l-error", "error above which a vector is no more considerated as an eigenvector of L^2", 1e-12);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHENBodyQuasiHoleOverlap -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  double LError = ((SingleDoubleOption*) Manager["l-error"])->GetDouble();

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["lz"])->GetInteger();
  bool FermionFlag = false;
  if (((SingleStringOption*) Manager["statistics"])->GetString() == 0)
    FermionFlag = true;
  if (FQHEOnSphereFindSystemInfoFromFileName(((SingleStringOption*) Manager["file-prefix"])->GetString(), NbrParticles, LzMax, FermionFlag) == false)
    {
      return -1;
    }
  if (((SingleStringOption*) Manager["statistics"])->GetString() != 0)
    if ((strcmp ("fermions", ((SingleStringOption*) Manager["statistics"])->GetString()) == 0))
      {
	FermionFlag = true;
      }
    else
      if ((strcmp ("fermions", ((SingleStringOption*) Manager["statistics"])->GetString()) == 0))
	{
	  FermionFlag = false;
	}
      else
	{
	  cout << ((SingleStringOption*) Manager["statistics"])->GetString() << " is an undefined statistics" << endl;
	}  
  int Parity = Lz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

  char* SpectrumFilename = AddExtensionToFileName(((SingleStringOption*) Manager["file-prefix"])->GetString(), "dat");
  cout << SpectrumFilename << endl;
  if (IsFile(SpectrumFilename) == false)
    {
      cout << "Spectrum " << SpectrumFilename << " does not exist or can't be opened" << endl;
      return -1;           
    }
  QHEOnSphereLzSortedSpectrum Spectrum (SpectrumFilename, 1e-12);
  if (Spectrum.IsSpectrumValid() == false)
    {
      cout << "Spectrum " << SpectrumFilename << " is not valid" << endl;
      return -1;           
    }
  
  char* VectorPrefix = new char [strlen(((SingleStringOption*) Manager["file-prefix"])->GetString()) + 24];
  sprintf (VectorPrefix , "%s_%d.", ((SingleStringOption*) Manager["file-prefix"])->GetString(), Lz);
  char** VectorFiles;
  int NbrVectorFiles = GetAllFilesDirectories(VectorPrefix, VectorFiles, ".vec");
  if (NbrVectorFiles == 0)
    {
      cout << "no available vector" << endl;
      return -1;
    }

  long MemorySpace = 9l << 20;
  ParticleOnSphere* Space;
  if (FermionFlag == true)
    {
#ifdef __64_BITS__
      if (LzMax <= 63)
	{
	  Space = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  Space = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#else
      if (LzMax <= 31)
	{
	  Space = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  Space = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#endif
    }
  else
    {
      Space = new BosonOnSphere(NbrParticles, Lz, LzMax);
    }

  int TotalMaxLz = LzMax * NbrParticles;
  if (FermionFlag == true)
    {
      TotalMaxLz = (LzMax - NbrParticles + 1) * NbrParticles;
    }

   ParticleOnSphereSquareTotalMomentumOperator oper(Space, LzMax);

  int CurrentVector = 0;
  int* NbrSortedVectors = new int [TotalMaxLz + 1];
  int* VectorPosition = new int [TotalMaxLz + 1];
  int* GlobalPosition = new int [TotalMaxLz + 1];
  for (int i = 0; i <= TotalMaxLz; ++i)
    GlobalPosition[i] = 0;

  char* OutputVectorPrefix = new char [strlen(((SingleStringOption*) Manager["file-prefix"])->GetString()) + 64];
  while (CurrentVector < NbrVectorFiles)
    {
      int TmpDegeneracy = Spectrum.GetDegeneracy(Lz,CurrentVector );
      RealVector* TmpVectors = new RealVector[TmpDegeneracy]; 
      for (int j = 0; j < TmpDegeneracy; ++j)
	{
	  if (TmpVectors[j].ReadVector(VectorFiles[j + CurrentVector]) == false)
	    {
	      cout << "error while reading " << VectorFiles[j + CurrentVector] << endl;
	      return -1;
	    }
	}
      RealMatrix DiagonalBasis (TmpVectors, TmpDegeneracy);      
      LSortBasis(DiagonalBasis, &oper, TotalMaxLz, NbrSortedVectors, VectorPosition, LError, VectorFiles + CurrentVector);   
      int Pos = 0;
      for (int j = 0; j <= TotalMaxLz; ++j)
	if (NbrSortedVectors[j] > 0)
	  {
	    Pos = VectorPosition[j];
	    for  (int k = 0; k < NbrSortedVectors[j]; ++k)
	      {
		sprintf (OutputVectorPrefix , "%s_%d_l_%d.%d.vec", ((SingleStringOption*) Manager["file-prefix"])->GetString(), Lz, ((2 * j) + Parity), GlobalPosition[j]);
		cout << "write " << OutputVectorPrefix << endl;
		TmpVectors[Pos].WriteVector(OutputVectorPrefix);
		++GlobalPosition[j];
		++Pos;
	      }
	  }
      CurrentVector += TmpDegeneracy;
    }

  for (int i = 0; i < NbrVectorFiles; ++i)
    delete[] VectorFiles[i];
  delete[] VectorFiles;
  delete[] NbrSortedVectors;
  delete[] VectorPosition;
  delete[] GlobalPosition;
  return 0;
}

// subdivide a basis into different L orthonormal eigenstate
//
// vectors = reference on the matrix whose columns are the vectors that span the basis
// oper = pointer to operator whivh allow to evaluate total L matrix elements
// totalMaxLz = maximum L value that can be reach by the system (-1/2 if L is half integer)
// subspaceSize = reference on the array that will be filled  with the dimension of each fixed L subspace (index equal to L if L is integer, L-1/2 if L is half integer)
// subspacePositions = reference on the array that will be filled  with the position of the first occurence of a vector with a given fixed L subspace 
//                     (index equal to L if L is integer, L-1/2 if L is half integer)
// lError = allowed error on L determination, if relative error between <L> and int(<L>) is greater than lError, no  L orthonormal eigenstates van be extracted
// vectorNames = array that contains file name of each vectors (used to identify vector in error message, 0 if vector index has to be used instead)

void LSortBasis(RealMatrix& vectors, AbstractOperator* oper, int totalMaxLz, int* subspaceSize, int* subspacePositions, double lError, char** vectorNames)
{
  for (int j = 0; j <= totalMaxLz; ++j)    
    { 
      subspacePositions[j] = totalMaxLz + 1;
      subspaceSize[j] = 0;
    }
  if (vectors.GetNbrColumn() > 1)
    {
      RealSymmetricMatrix HRep (vectors.GetNbrColumn(), vectors.GetNbrColumn());
      for (int k = 0; k < vectors.GetNbrColumn(); ++k)
	for (int l = k; l < vectors.GetNbrColumn(); ++l)
	  HRep.SetMatrixElement(l, k,  oper->MatrixElement(vectors[k], vectors[l]).Re);
      RealMatrix TmpEigenvector (vectors.GetNbrColumn(), vectors.GetNbrColumn(), true);
      for (int l = 0; l < vectors.GetNbrColumn(); ++l)
	TmpEigenvector(l, l) = 1.0;
      RealTriDiagonalSymmetricMatrix TmpTriDiag (vectors.GetNbrColumn());
      HRep.Householder(TmpTriDiag, 1e-7, TmpEigenvector);
      TmpTriDiag.Diagonalize(TmpEigenvector);
      TmpTriDiag.SortMatrixUpOrder(TmpEigenvector);
      int TmpAngularMomentum;
      double RawTmpAngularMomentum;
      for (int k = 0; k < vectors.GetNbrColumn(); ++k)	    
	{
	  RawTmpAngularMomentum = (sqrt ((4.0 * TmpTriDiag.DiagonalElement(k)) + 1.0) - 1.0);
	  if (RawTmpAngularMomentum < 0.0)
	    RawTmpAngularMomentum = 0.0;
	  TmpAngularMomentum = ((int) round(RawTmpAngularMomentum));
	  if (((TmpAngularMomentum == 0) && (fabs(RawTmpAngularMomentum) > lError)) ||
	      (fabs(RawTmpAngularMomentum - round(RawTmpAngularMomentum)) > (lError * fabs(RawTmpAngularMomentum))))
	    {
	      cout << "error while processing vectors ";
	      cout << vectorNames[0] << " ";
	      for (int i = 1; i < vectors.GetNbrColumn(); ++i)
		cout << "," << vectorNames[i] << " ";
	      cout << endl << "<2L> = " << RawTmpAngularMomentum << endl;
	    }
	  if ((TmpAngularMomentum & 1) == 0)
	    {
	      TmpAngularMomentum >>= 1;
	    }
	  else
	    {		  
	      TmpAngularMomentum -= 1;
	      TmpAngularMomentum >>= 1;		  
	    }
	  subspaceSize[TmpAngularMomentum]++;
	  if (subspacePositions[TmpAngularMomentum] > k)
	    subspacePositions[TmpAngularMomentum] = k;
	}
      vectors.Multiply(TmpEigenvector);
    }
  else
    {
      double RawTmpAngularMomentum = (sqrt ((4.0 * oper->MatrixElement(vectors[0], vectors[0]).Re) + 1.0) - 1.0);
      if (RawTmpAngularMomentum < 0.0)
	RawTmpAngularMomentum = 0.0;
      int TmpAngularMomentum = ((int) round(RawTmpAngularMomentum));
      if (((TmpAngularMomentum == 0) && (fabs(RawTmpAngularMomentum) > lError)) ||
	  (fabs(RawTmpAngularMomentum - round(RawTmpAngularMomentum)) > (lError * fabs(RawTmpAngularMomentum))))
	{
	  cout << "error while processing vectors " << vectorNames[0] << endl << "<2L> = " << RawTmpAngularMomentum << endl;
	}
      if ((TmpAngularMomentum & 1) == 0)
	{
	  TmpAngularMomentum >>= 1;
	}
      else
	{		  
	  TmpAngularMomentum -= 1;
	  TmpAngularMomentum >>= 1;		  
	}
      subspaceSize[TmpAngularMomentum]++;	    
      subspacePositions[TmpAngularMomentum] = 0;
    }
}

