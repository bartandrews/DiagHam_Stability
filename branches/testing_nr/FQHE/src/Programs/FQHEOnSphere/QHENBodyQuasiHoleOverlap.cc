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

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/ArrayTools.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Tools/FQHESpectrum/QHEOnSphereLzSortedSpectrum.h"

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
void LSortBasis(RealMatrix& vectors, AbstractOperator* oper, int totalMaxLz, int* subspaceSize, int* subspacePositions);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHENBodyQuasiHoleOverlap" , "0.01");
  OptionGroup* MainGroup = new OptionGroup ("main options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += MiscGroup;
  Manager += MainGroup;
 
  (*MainGroup) += new SingleStringOption  ('\n', "input-file", "name of the file which contains definition of overlaps to evaluate");
  (*MainGroup) += new SingleDoubleOption  ('\n', "ortho-error", "scalar product value below which two states are considered as orthogonal", 1e-12);
  (*MainGroup) += new SingleIntegerOption  ('\n', "output-precision", "numerical display precision", 14, true, 2, true, 14);
  (*MainGroup) += new BooleanOption  ('\n', "latex-output", "display final results as a latex table formatted string");
  (*MainGroup) += new BooleanOption  ('\n', "global-overlap", "display the global overlap for the whole quasihole subspace");
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

  double OrthogonalityError = ((SingleDoubleOption*) Manager["ortho-error"])->GetDouble();

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
  int* LzDegeneracy;
  if (OverlapDefinition.GetAsIntegerArray("Degeneracy", ' ', LzDegeneracy, MaxNbrLz) == false)
    {
      cout << "error while parsing Degeneracy in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }
  
  if (OverlapDefinition["Spectrum"] == 0)
    {
      cout << "no Spectrum defined in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }
  QHEOnSphereLzSortedSpectrum Spectrum (OverlapDefinition["Spectrum"]);
  if (Spectrum.IsSpectrumValid() == false)
    {
      cout << "Spectrum " << OverlapDefinition["Spectrum"] << " is unreadable or is not valid" << endl;
      return -1;           
    }
  
  char* InputVectors = new char [strlen(OverlapDefinition["ExactStates"]) + 24];
  char* OutputVectors = new char [strlen(OverlapDefinition["QuasiholeStates"]) + 24];
  strcpy (InputVectors, OverlapDefinition["ExactStates"]);
  strcpy (OutputVectors, OverlapDefinition["QuasiholeStates"]);
  char* InputVectors2 = InputVectors + strlen(InputVectors);
  char* OutputVectors2 = OutputVectors + strlen(OutputVectors);
  bool FermionFlag = false;
  if ((OverlapDefinition["Statistics"] != 0) && (strcmp ("fermions", OverlapDefinition["Statistics"]) == 0))
    {
      FermionFlag = true;
    }
  long MemorySpace = 9l << 20;
  int* LDegeneracy = new int [MaxNbrLz];
  LDegeneracy[MaxNbrLz - 1] = LzDegeneracy[MaxNbrLz - 1];
  for (int i = MaxNbrLz - 2; i >= 0; --i)
    LDegeneracy[i] = LzDegeneracy[i] - LzDegeneracy[i + 1];

  int TotalMaxLz = LzMax * NbrParticles;
  if (FermionFlag == true)
    {
      TotalMaxLz = (LzMax - NbrParticles + 1) * NbrParticles;
    }
  RealVector** SortedTestVectors = new RealVector*[MaxNbrLz];
  int* NbrSortedTestVectors = new int [TotalMaxLz + 1];
  int* ReferenceNbrSortedVectors = new int [TotalMaxLz + 1];
  bool Parity = true;
  if ((TotalMaxLz & 1) != 0)
    Parity = false;
  double* BestOverlaps = new double [MaxNbrLz + 1];
  int* TestVectorPosition = new int [TotalMaxLz + 1];
  int* ReferenceVectorPosition = new int [TotalMaxLz + 1];

  for (int i = MaxNbrLz - 1; i >= 0; --i)
    if (LDegeneracy[i] > 0)
      {
	cout << "----------------------------------------------------" << endl;
	int Lz = i << 1;
	if (Parity == false)
	  {
	    ++Lz;
	    cout << "L = " << Lz << "/2" << endl;
	  }
	else
	  {
	    cout << "L = " << i << endl;
	  }
	int* NbrVectorWithMomentum = new int [TotalMaxLz]; 
	RealVector* TestVectors = new RealVector[LzDegeneracy[i]]; 
	for (int j = 0; j < TotalMaxLz; ++j)    
	  { 
	    NbrVectorWithMomentum[j] = 0;
	  }
	for (int j = 0; j < LzDegeneracy[i]; ++j)
	  {
	    sprintf (OutputVectors2, "%d.%d.vec", Lz, j);	      
	    if (TestVectors[j].ReadVector(OutputVectors) == false)
	      {
		cout << "error while reading " << OutputVectors << endl;
		return -1;
	      }
	  }
	RealMatrix DiagonalBasis (TestVectors, LzDegeneracy[i]);

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
	ParticleOnSphereSquareTotalMomentumOperator oper(Space, LzMax);
	LSortBasis(DiagonalBasis, &oper, TotalMaxLz, NbrSortedTestVectors, TestVectorPosition);
	cout << endl;

	bool OpenFlag = true;
	int NbrReferenceStates = 0;
	int VectorPosition = 0;
	RealVector* ReferenceVectors = new RealVector[200];	
	while ((OpenFlag == true) && (NbrReferenceStates < LDegeneracy[i]))
	  {
	    int ReferenceDegeneracy = Spectrum.GetDegeneracy(Lz, VectorPosition);
	    if (ReferenceDegeneracy > 0)
	      {
		RealVector* TmpReferenceVectors = new RealVector[ReferenceDegeneracy];	
		int k = 0;
		for (; (k < ReferenceDegeneracy) && (OpenFlag == true); ++k)  
		  {
		    sprintf (InputVectors2, "%d.%d.vec", Lz, VectorPosition + k);
		    if (TmpReferenceVectors[k].ReadVector(InputVectors) == false)
		      {
			OpenFlag = false;
		      }
		  }
		if (OpenFlag == true)
		  {
		    RealMatrix ReferenceDiagonalBasis (TmpReferenceVectors, ReferenceDegeneracy);
		    LSortBasis(ReferenceDiagonalBasis, &oper, TotalMaxLz, ReferenceNbrSortedVectors, ReferenceVectorPosition);
		    for (k = 0; k < ReferenceNbrSortedVectors[i]; ++k)  
		      {
			ReferenceVectors[NbrReferenceStates] = ReferenceDiagonalBasis[ReferenceVectorPosition[i] + k];
			++NbrReferenceStates;
		      }
		    if ( ReferenceNbrSortedVectors[i] > 0)
		      for (k = 0; k < ReferenceDegeneracy; ++k)  
			{
			  sprintf (InputVectors2, "%d.%d.vec", Lz, VectorPosition + k);
			  cout << InputVectors << endl;
			}
		  }
		VectorPosition += ReferenceDegeneracy;
	      }
	    else
	      OpenFlag = false;
	  }
	if (NbrReferenceStates == 0)
	  {
	    cout << "no possible overlap calculation" << endl;
	  }
	else
	  if (NbrReferenceStates == NbrSortedTestVectors[i])
	    {
	      double Scalar = 0.0;
	      double TmpScalar;
	      int Shift = TestVectorPosition[i];
	      cout << "scalar products: ";
	      for (int k = 0; k < NbrReferenceStates; ++k)
		for (int l = 0; l < NbrReferenceStates; ++l)		
		  {
		    TmpScalar = ReferenceVectors[k] * DiagonalBasis[l + Shift];
		    TmpScalar *= TmpScalar;
		    cout << TmpScalar << " ";
		    Scalar += TmpScalar;
		  }
	      cout << endl;
	      Scalar /= (double) NbrReferenceStates;
	      BestOverlaps[i] = Scalar;
	      cout << "overlap = " << Scalar << endl;
	    }	  
	  else
	    if (NbrSortedTestVectors[i] == 1)
	      {
		int Shift = TestVectorPosition[i];
		BestOverlaps[i] = 0.0;
		for (int k = 0; k < ReferenceNbrSortedVectors[i]; ++k)
		  {
		    double Scalar = ReferenceVectors[ReferenceVectorPosition[i] + k] * DiagonalBasis[Shift];
		    if (fabs(Scalar) > OrthogonalityError)
		      cout << "    " << OutputVectors << "  = " << Scalar << endl;
		    BestOverlaps[i] += Scalar * Scalar;
		  }
		cout << "overlap = " << BestOverlaps[i]  << endl;
	      }
	    else
	      {
		cout << "can't compare vectors if both subspaces don't have the same dimension " << endl;		  
		BestOverlaps[i] = -1.0;
	      }
	delete Space;
      }

  
  cout << "----------------------------------------------------" << endl;
  cout << "final results" << endl;
  cout << "----------------------------------------------------" << endl << endl;
  cout.precision(((SingleIntegerOption*) Manager["output-precision"])->GetInteger());
  double GlobalOverlap = 0.0;
  if (((BooleanOption*) Manager["global-overlap"])->GetBoolean() == true)
    {
      cout << endl << "global overlap = ";
      int TotalDegeneracy = 0;
      int TmpLValue = 1;
      if (Parity == false)
	{
	  TmpLValue = 2;
	}
      for (int i = 0; i < MaxNbrLz; ++i)
	{
	  if (LDegeneracy[i] > 0)
	    {
	      GlobalOverlap +=  ((double) (LDegeneracy[i] * TmpLValue)) * BestOverlaps[i];
	      TotalDegeneracy += LDegeneracy[i] * TmpLValue;
	    }
	  TmpLValue += 2;
	}
      GlobalOverlap /= ((double) TotalDegeneracy);
      cout << "total overlap =" << GlobalOverlap << endl;
    }
  for (int i = 0; i < MaxNbrLz; ++i)
    {
      if (LDegeneracy[i] > 0)
	{
	  int Lz = i << 1;
	  if (Parity == false)
	    {
	      ++Lz;
	      cout << "L = " << Lz << "/2 : ";
	    }
	  else
	    {
	      cout << "L = " << i <<" : ";
	    }
	  cout << BestOverlaps[i] << endl;
	}
    }
  if (((BooleanOption*) Manager["latex-output"])->GetBoolean() == true)
    {
      cout << endl << "latex output:" << endl;
      if (((BooleanOption*) Manager["global-overlap"])->GetBoolean() == true)
	{
	  cout << "$" << GlobalOverlap << "$ & " ;
	}
      for (int i = 0; i < (MaxNbrLz - 1); ++i)
	{
	  if (LDegeneracy[i] > 0)
	    {
	      cout << "$" << BestOverlaps[i] << "$" ;
	    }
	  cout << " & ";
	}
      if (LDegeneracy[MaxNbrLz - 1] > 0)
	{
	  cout << "$" << BestOverlaps[MaxNbrLz - 1] << "$" ;
	}
      cout << " \\\\" << endl;
    }


  delete[] BestOverlaps;
  delete[] LzDegeneracy;
  delete[] InputVectors;
  delete[] OutputVectors;
  delete[] SortedTestVectors;
  delete[] NbrSortedTestVectors;
  delete[] TestVectorPosition;
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

void LSortBasis(RealMatrix& vectors, AbstractOperator* oper, int totalMaxLz, int* subspaceSize, int* subspacePositions)
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
      for (int k = 0; k < vectors.GetNbrColumn(); ++k)	    
	{
	  TmpAngularMomentum = ((int) round((sqrt ((4.0 * TmpTriDiag.DiagonalElement(k)) + 1.0) - 1.0)));
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
      int TmpAngularMomentum = ((int) round((sqrt ((4.0 * oper->MatrixElement(vectors[0], vectors[0]).Re) + 1.0) - 1.0)));
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

