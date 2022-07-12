#include "config.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

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

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/ArrayTools.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "HilbertSpace/BosonOnSphere.h"

#include "Hamiltonian/ParticleOnSphereDeltaHamiltonian.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHENBodyQuasiHoleOverlap" , "0.01");
  OptionGroup* MainGroup = new OptionGroup ("main options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += MainGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*MainGroup) += new SingleStringOption  ('\n', "input-file", "name of the file which contains definition of overlaps to evaluate");
  (*MainGroup) += new SingleStringOption  ('\n', "output-file", "name of the file which contains definition of overlaps to evaluate");
  (*MainGroup) += new BooleanOption  ('\n', "lsort", "diagonalize in each fixed total l sector");
  (*MainGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
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

  bool LSortFlag = ((BooleanOption*) Manager["lsort"])->GetBoolean();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;

  ConfigurationParser OverlapDefinition;
  if (OverlapDefinition.Parse(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      OverlapDefinition.DumpErrors(cout) << endl;
      return -1;
    }

  int LzMax = 0;
  int NbrParticles = 0;
  if ((OverlapDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax < 0))
    {
      cout << "LzMax is not defined or as a wrong value" << endl;
      return -1;
    }
  if ((OverlapDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
    {
      cout << "NbrParticles is not defined or as a wrong value" << endl;
      return -1;
    }

  if (OverlapDefinition["Degeneracy"] == 0)
    {
      cout << "no Degeneracy defined in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }

  int MaxNbrLz;
  int* Degeneracy;
  if (OverlapDefinition.GetAsIntegerArray("Degeneracy", ' ', Degeneracy, MaxNbrLz) == false)
    {
      cout << "error while parsing Degeneracy in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }
  
  
  char* VectorBasisName = new char [strlen(OverlapDefinition["Vectors"]) + 24];
  strcpy (VectorBasisName, OverlapDefinition["Vectors"]);
  char* VectorBasisName2 = VectorBasisName + strlen(VectorBasisName);

  int TotalMaxLz = LzMax * NbrParticles;
  bool Parity = true;
  if ((TotalMaxLz & 1) != 0)
    Parity = false;

  ofstream File;
  if (((SingleStringOption*) Manager["output-file"])->GetString() != 0)
    {
       File.open(((SingleStringOption*) Manager["output-file"])->GetString(), ios::binary | ios::out);
       File.precision(14);     
    }

  if (LSortFlag == true)
    {
      int* NbrSortedTestVectors = new int [2 * TotalMaxLz];
      int* TestVectorPosition = new int [2 * TotalMaxLz];

      for (int i = 0; i < MaxNbrLz; ++i)
	{
	  int Lz = i << 1;
	  if (Parity == false)
	    {
	      ++Lz;
	      cout << "Lz = " << Lz << "/2" << endl;
	    }
	  else
	    {
	      cout << "Lz = " << i << endl;
	    }
	  RealVector* TestVectors = new RealVector[Degeneracy[i]]; 
	  for (int j = 0; j < TotalMaxLz; ++j)  
	    {  
	      TestVectorPosition[j] = TotalMaxLz;
	      NbrSortedTestVectors[j] = 0;
	    }
	  
	  for (int j = 0; j < Degeneracy[i]; ++j)
	    {
	      sprintf (VectorBasisName2, "%d.%d.vec", Lz, j);	      
	      if (TestVectors[j].ReadVector(VectorBasisName) == false)
		{
		  cout << "error while reading " << VectorBasisName << endl;
		  return -1;
		}
	    }
	  
	  RealMatrix DiagonalBasis (TestVectors, Degeneracy[i]);
	  BosonOnSphere Space (NbrParticles, Lz, LzMax);
	  
	  if (Degeneracy[i] > 1)
	    {
	      RealSymmetricMatrix HRep (Degeneracy[i], Degeneracy[i]);
	      ParticleOnSphereSquareTotalMomentumOperator oper(&Space, LzMax);
	      for (int k = 0; k < Degeneracy[i]; ++k)
		{
		  for (int l = k; l < Degeneracy[i]; ++l)
		    {	
		      HRep.SetMatrixElement(l, k,  oper.MatrixElement(TestVectors[k], TestVectors[l]).Re);
		    }
		}
	      cout << endl << endl;
	      
	      RealMatrix TmpEigenvector (Degeneracy[i], Degeneracy[i], true);
	      for (int l = 0; l < Degeneracy[i]; ++l)
		TmpEigenvector(l, l) = 1.0;
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Degeneracy[i]);
	      HRep.Householder(TmpTriDiag, 1e-7, TmpEigenvector);
	      TmpTriDiag.Diagonalize(TmpEigenvector);
	      TmpTriDiag.SortMatrixUpOrder(TmpEigenvector);
	      cout << "angular momentum of test eigenvectors = ";
	      int TmpAngularMomentum;
	      for (int k = 0; k < Degeneracy[i]; ++k)	    
		{
		  TmpAngularMomentum = ((int) round((sqrt ((4.0 * TmpTriDiag.DiagonalElement(k)) + 1.0) - 1.0)));
		  if (Parity == false)
		    cout << TmpAngularMomentum << "/2 ";	      
		  else
		    cout << (TmpAngularMomentum >> 1) << " ";	      
		  NbrSortedTestVectors[TmpAngularMomentum]++;
		  if (TestVectorPosition[TmpAngularMomentum] > k)
		    TestVectorPosition[TmpAngularMomentum] = k;
		}
	      
	      cout << endl;
	      DiagonalBasis.Multiply(TmpEigenvector);
	    }
	  else
	    if (Degeneracy[i] == 1)
	      {
		ParticleOnSphereSquareTotalMomentumOperator oper(&Space, Lz, LzMax);
		int TmpAngularMomentum = ((int) round((sqrt ((4.0 * oper.MatrixElement(TestVectors[0], TestVectors[0]).Re) + 1.0) - 1.0)));
		cout << "angular momentum of test eigenvectors = ";
		if (Parity == false)
		  cout << TmpAngularMomentum  << "/2 " << endl;
		else
		  cout << (TmpAngularMomentum >> 1)  << endl;
		NbrSortedTestVectors[TmpAngularMomentum]++;
		TestVectorPosition[TmpAngularMomentum] = 0;
	      }
	  
	  Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());
	  AbstractQHEOnSphereHamiltonian* Hamiltonian = new ParticleOnSphereDeltaHamiltonian(&Space, NbrParticles, LzMax, Architecture.GetArchitecture(), Memory);
	  for (int j = 0; j < MaxNbrLz; ++j)
	    {
	      int Pos = j;
	      if (Parity == false)
		Pos = (j << 1) + 1;
	      else
		Pos = (j << 1);
	      if (NbrSortedTestVectors[Pos] > 1)
		{
		  int Shift = TestVectorPosition[Pos];
		  RealSymmetricMatrix HRep (NbrSortedTestVectors[Pos], NbrSortedTestVectors[Pos]);
		  RealVector TmpVector (DiagonalBasis[Shift].GetVectorDimension());
		  for (int k = 0; k < NbrSortedTestVectors[Pos]; ++k)
		    {
		      Hamiltonian->Multiply(DiagonalBasis[k + Shift], TmpVector);
		      for (int l = k; l < NbrSortedTestVectors[Pos]; ++l)
			{	
			  HRep.SetMatrixElement(l, k, (TmpVector * DiagonalBasis[l + Shift]));
			}
		    }
		  RealTriDiagonalSymmetricMatrix TmpTriDiag (NbrSortedTestVectors[Pos]);
		  HRep.Householder(TmpTriDiag, 1e-7);
		  TmpTriDiag.Diagonalize();
		  TmpTriDiag.SortMatrixUpOrder();
		  for (int k = 0; k < NbrSortedTestVectors[Pos]; ++k)
		    cout << j << " " << TmpTriDiag.DiagonalElement(k) << endl;
		  if ((i == 0) && (((SingleStringOption*) Manager["output-file"])->GetString() != 0))
		    for (int k = 0; k < NbrSortedTestVectors[Pos]; ++k)
		      File << j << " " << TmpTriDiag.DiagonalElement(k) << endl;		    
		}
	      else
		if (NbrSortedTestVectors[Pos] == 1)
		  {
		    int Shift = TestVectorPosition[Pos];
		    RealVector TmpVector (DiagonalBasis[Shift].GetVectorDimension());
		    Hamiltonian->Multiply(DiagonalBasis[Shift], TmpVector);
		    double Scalar = (TmpVector * DiagonalBasis[Shift]);
		    cout << j << " " << Scalar << endl;
		    if ((i == 0) && (((SingleStringOption*) Manager["output-file"])->GetString() != 0))
		      File << j << " " << Scalar << endl;
		  }
	    }
	  delete Hamiltonian;
	}
      
      delete[] Degeneracy;
      delete[] VectorBasisName;
      delete[] NbrSortedTestVectors;
      delete[] TestVectorPosition;
    }
  else
    {
      for (int i = 0; i < MaxNbrLz; ++i)
	{
	  int Lz = i << 1;
	  if (Parity == false)
	    {
	      ++Lz;
	      cout << "Lz = " << Lz << "/2" << endl;
	    }
	  else
	    {
	      cout << "Lz = " << i << endl;
	    }
	  BosonOnSphere Space (NbrParticles, Lz, LzMax);
	  Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());
	  AbstractQHEOnSphereHamiltonian* Hamiltonian = new ParticleOnSphereDeltaHamiltonian(&Space, NbrParticles, LzMax, Architecture.GetArchitecture(), Memory);
	  if (Degeneracy[i] > 1)
	    {
	      RealSymmetricMatrix HRep (Degeneracy[i]);
	      RealVector TmpVector;
	      for (int k = 0; k < Degeneracy[i]; ++k)
		{
		  RealVector TestVector;
		  sprintf (VectorBasisName2, "%d.%d.vec", Lz, k);	      
		  if (TestVector.ReadVector(VectorBasisName) == false)
		    {
		      cout << "error while reading " << VectorBasisName << endl;
		      return -1;
		    }
		  RealVector TmpVector (TestVector.GetVectorDimension());
		  Hamiltonian->Multiply(TestVector, TmpVector);
		  HRep.SetMatrixElement(k, k, (TmpVector * TestVector));
		  for (int l = k + 1; l < Degeneracy[i]; ++l)
		    {	
		      sprintf (VectorBasisName2, "%d.%d.vec", Lz, l);	      
		      if (TestVector.ReadVector(VectorBasisName) == false)
			{
			  cout << "error while reading " << VectorBasisName << endl;
			  return -1;
			}
		      HRep.SetMatrixElement(l, k, (TmpVector * TestVector));
		    }
		}
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Degeneracy[i]);
	      HRep.Householder(TmpTriDiag, 1e-7);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      for (int k = 0; k < Degeneracy[i]; ++k)
		cout << i << " " << TmpTriDiag.DiagonalElement(k) << endl;
	      if (((SingleStringOption*) Manager["output-file"])->GetString() != 0)
		for (int k = 0; k < Degeneracy[i]; ++k)
		  File << i << " " << TmpTriDiag.DiagonalElement(k) << endl;
	    }
	  else
	    {
	      RealVector TestVector;
	      sprintf (VectorBasisName2, "%d.%d.vec", Lz, 0);	      
	      if (TestVector.ReadVector(VectorBasisName) == false)
		{
		  cout << "error while reading " << VectorBasisName << endl;
		  return -1;
		}
	      RealVector TmpVector (TestVector.GetVectorDimension());
	      Hamiltonian->Multiply(TestVector, TmpVector);
	      double Scalar = (TmpVector * TestVector);
	      cout << i << " " << Scalar << endl;
	      if (((SingleStringOption*) Manager["output-file"])->GetString() != 0)
		File << i << " " << Scalar << endl;
	    }
	  delete Hamiltonian;
	}
    }
  if (((SingleStringOption*) Manager["output-file"])->GetString() != 0)
    {
      File.close();
    }
  return 0;
}


