#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithSpin.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorusWithSpin.h"

#include "Operator/ParticleOnTorusKxOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


// invert states from the (kx,ky) basis to the ky basis
//
// manager = pointer to the option manager
// return value = 0 if no error occured
int InvertU1KxKyStates(OptionManager* manager);

// invert states from the (kx,ky) basis to the ky basis with SU(2) degree of freedom
//
// manager = pointer to the option manager
// return value = 0 if no error occured
int InvertSU2KxKyStates(OptionManager* manager);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHETorusComputeKx" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption('i', "input-state", "name of the file containing the state whose Kx momentum has to be computed");
  (*SystemGroup) += new SingleStringOption('\n', "degenerated-states", "name of the file containing a list of states (override input-state)");
  (*SystemGroup) += new SingleDoubleOption   ('r', "ratio", "ratio between lengths along the x and y directions", 1);
  (*SystemGroup) += new BooleanOption ('\n',  "compute-eigenstate", "compute the eigenstates of the Kx operator in the given basis");
  (*SystemGroup) += new SingleStringOption ('\n',  "interaction-name", "name that should be inserted in the output file names", "dummy");
  (*SystemGroup) += new BooleanOption ('\n',  "no-convertion", "do not convert the final vectors to the (Kx,Ky) n-body basis");
  (*SystemGroup) += new BooleanOption ('\n',  "invert", "assume the input states are in the (Kx,Ky) and express them in the Ky basis");
  (*SystemGroup) += new BooleanOption ('\n',  "invert-real", "when using the invert option, assume that the final vector is real");
  (*SystemGroup) += new BooleanOption ('\n',  "export-transformation", "export the transformation matrix in a ascii file (one per momentum sector)");
  (*SystemGroup) += new BooleanOption ('\n',  "export-bintransformation", "export the transformation matrix in a binary file (one per momentum sector)");
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusComputeKx -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int MaxMomentum = 0;
  int YMomentum = -1;
  int TotalSz = 0;
  bool Statistics = true;
  bool SU2Flag = false;
  bool SU3Flag = false;
  bool SU4Flag = false;
  if ((Manager.GetBoolean("invert")) || (Manager.GetBoolean("invert-real")))
    {
      return InvertU1KxKyStates(&Manager);
    }
  

  RealVector* InputStates = 0;
  ComplexVector* ComplexInputStates = 0;
  int NbrInputStates = 0;
  bool ComplexVectorFlag = false;
  if (Manager.GetString("degenerated-states") == 0)
    {
      if (Manager.GetString("input-state") == 0)
	{
	  cout << "error, either input-state or degenerated-states has to be provided" << endl;
	  return -1;
	}
      RealVector TestVector;
      if (TestVector.ReadVectorTest(Manager.GetString("input-state")) == false)
	{
	  ComplexVectorFlag = true;
	}
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						      NbrParticles, MaxMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("ground-state") << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Ky=" << YMomentum << " ";
      NbrInputStates = 1;
      if (ComplexVectorFlag == false)
	{
	  InputStates = new RealVector [NbrInputStates];
	  if (InputStates[0].ReadVector(Manager.GetString("input-state")) == false)
	    {
	      cout << "error while reading " << Manager.GetString("input-state") << endl;
	      return -1;
	    }
	}
      else
	{
	  ComplexInputStates = new ComplexVector [NbrInputStates];
	  if (ComplexInputStates[0].ReadVector(Manager.GetString("input-state")) == false)
	    {
	      cout << "error while reading " << Manager.GetString("input-state") << endl;
	      return -1;
	    }
	}
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-states")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	} 
      NbrInputStates = DegeneratedFile.GetNbrLines();
      if (NbrInputStates < 1)
	{
	  cout << "no state found in " << Manager.GetString("degenerated-states") << endl;
	}
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, 0),
						      NbrParticles, MaxMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, 0) << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Ky=" << YMomentum << " ";
      RealVector TestVector;
      if (TestVector.ReadVectorTest(DegeneratedFile(0, 0)) == false)
	ComplexVectorFlag = true;
     
      if (ComplexVectorFlag == false)
	{
	  InputStates = new RealVector [NbrInputStates];
	  if (InputStates[0].ReadVector(DegeneratedFile(0, 0)) == false)
	    {
	      cout << "error while reading " << DegeneratedFile(0, 0) << endl;
	      return -1;
	    }
	}
       else
	{
	  ComplexInputStates = new ComplexVector [NbrInputStates];
	  if (ComplexInputStates[0].ReadVector(DegeneratedFile(0, 0)) == false)
	    {
	      cout << "error while reading " << DegeneratedFile(0, 0) << endl;
	      return -1;
	    }
	}
      for (int i = 1; i < NbrInputStates; ++i)
	{
	  int TmpNbrParticles = 0;
	  int TmpMaxMomentum = 0;
	  int TmpYMomentum = -1;
	  bool TmpStatistics = true;
	  if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, i),
							  TmpNbrParticles, TmpMaxMomentum, TmpYMomentum, TmpStatistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, i) << endl;
	      return -1;
	    }
	  if ((TmpNbrParticles != NbrParticles) || (TmpMaxMomentum != MaxMomentum) || 
	      (TmpYMomentum != YMomentum) || (Statistics != TmpStatistics))
	    {
	      cout << "error, " << DegeneratedFile(0, i) << " has different system parameters than " << DegeneratedFile(0, 0) 
		   << ", N=" << TmpNbrParticles << "(" << NbrParticles << "), N_phi=" << TmpMaxMomentum << "(" << MaxMomentum 
		   << "), Ky=" << TmpYMomentum << "(" << YMomentum << ")" << endl;
	    }
	  if (ComplexVectorFlag == false)
	    {
	      if (InputStates[i].ReadVector(DegeneratedFile(0, i)) == false)
		{
		  cout << "error while reading " << DegeneratedFile(0, i) << endl;
		  return -1;
		}
	      if (InputStates[i].GetVectorDimension() != InputStates[0].GetVectorDimension())
		{
		  cout << "error, " << DegeneratedFile(0, i) << " has different dimension than " << DegeneratedFile(0, 0) << endl;	      
		}
	    }
	  else
	    {
	      if (ComplexInputStates[i].ReadVector(DegeneratedFile(0, i)) == false)
		{
		  cout << "error while reading " << DegeneratedFile(0, i) << endl;
		  return -1;
		}
	      if (ComplexInputStates[i].GetVectorDimension() != ComplexInputStates[0].GetVectorDimension())
		{
		  cout << "error, " << DegeneratedFile(0, i) << " has different dimension than " << DegeneratedFile(0, 0) << endl;	      
		}
	    }
	}
    }
  int MomentumModulo = FindGCD(NbrParticles, MaxMomentum);
  ParticleOnTorus* TotalSpace = 0;
  ParticleOnTorusWithMagneticTranslations** TargetSpaces = 0;
  char* OutputNamePrefix = new char [256 + strlen(Manager.GetString("interaction-name"))];

  if (Statistics == false)
    {
#ifdef  __64_BITS__
      if ((MaxMomentum + NbrParticles - 1) < 63)
#else
	if ((MaxMomentum + NbrParticles - 1) < 31)	
#endif
	  {
	    TotalSpace = new BosonOnTorusShort(NbrParticles, MaxMomentum, YMomentum);
	  }
	else
	  {
	    TotalSpace = new BosonOnTorus(NbrParticles, MaxMomentum, YMomentum);
	  }
      sprintf (OutputNamePrefix, "bosons_torus_%s_n_%d_2s_%d_ratio_%.6f", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum, Manager.GetDouble("ratio"));
      if (Manager.GetBoolean("compute-eigenstate") == true)
	{
	  TargetSpaces = new ParticleOnTorusWithMagneticTranslations*[MomentumModulo];
	  for (int i = 0; i < MomentumModulo; ++i)
	    TargetSpaces[i] = new BosonOnTorusWithMagneticTranslationsShort(NbrParticles, MaxMomentum, i, YMomentum);
	}
    }
  else
    {
      TotalSpace = new FermionOnTorus (NbrParticles, MaxMomentum, YMomentum);
      sprintf (OutputNamePrefix, "fermions_torus_%s_n_%d_2s_%d_ratio_%.6f", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum, Manager.GetDouble("ratio"));
      if (Manager.GetBoolean("compute-eigenstate") == true)
	{
	  TargetSpaces = new ParticleOnTorusWithMagneticTranslations*[MomentumModulo];
	  for (int i = 0; i < MomentumModulo; ++i)
	    TargetSpaces[i] = new FermionOnTorusWithMagneticTranslations(NbrParticles, MaxMomentum, i, YMomentum);
	}
    }
  if (((ComplexVectorFlag == false) && (InputStates[0].GetVectorDimension() != TotalSpace->GetHilbertSpaceDimension()))
      || ((ComplexVectorFlag == true) && (ComplexInputStates[0].GetVectorDimension() != TotalSpace->GetHilbertSpaceDimension())))
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputStates[0].GetVectorDimension() 
	   << " "<< TotalSpace->GetHilbertSpaceDimension() << endl;
      return -1;
    }

  Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());

  ComplexMatrix KxRep (NbrInputStates, NbrInputStates, true);
  ParticleOnTorusKxOperator KxOperator(TotalSpace);


  if (ComplexVectorFlag == false)
    {
      RealVector TmpVector(TotalSpace->GetHilbertSpaceDimension());
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  VectorOperatorMultiplyOperation Operation(&KxOperator, &(InputStates[i]), &TmpVector);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  for (int j = 0; j < NbrInputStates; ++j)
	    {
	      KxRep.SetMatrixElement(j, i, (InputStates[j] * TmpVector));
	    }
	}
    }
  else
    {
      ComplexVector TmpVector(TotalSpace->GetHilbertSpaceDimension());
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  VectorOperatorMultiplyOperation Operation(&KxOperator, &(ComplexInputStates[i]), &TmpVector);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  for (int j = 0; j < NbrInputStates; ++j)
	    {
	      KxRep.SetMatrixElement(j, i, (TmpVector * ComplexInputStates[j]));
	    }
	}
    }
  char* OutputName = new char [256 + strlen(OutputNamePrefix)];
  sprintf (OutputName, "%s_ky_%d.dat", OutputNamePrefix, YMomentum);
  ofstream File;
  File.open(OutputName, ios::binary | ios::out);
  File.precision(14);
  File << "# eigenvalue Norm Arg (GCD(N,N_phi)*Arg/2pi) Kx round(Kx)" << endl;

  if (Manager.GetBoolean("compute-eigenstate") == false)
    {
      ComplexDiagonalMatrix Eigenvalues(NbrInputStates, true);
#ifdef __LAPACK__
      KxRep.LapackDiagonalize(Eigenvalues);
#else
  cout << "error, lapack is required" << endl;
#endif
      double Factor = ((double) MomentumModulo) / (2.0 * M_PI);
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  Complex TmpValue = Eigenvalues[i];
	  File << TmpValue << " " << Norm(TmpValue) << " " << Arg(TmpValue);
	  double TmpValue2 = Arg(TmpValue);
	  if (fabs(TmpValue2) < MACHINE_PRECISION)
	    TmpValue2 = 0.0;
	  if (TmpValue2 < 0.0)
	    TmpValue2 += 2.0 * M_PI;
	  if (TmpValue2 >= (2.0 * M_PI))
	    TmpValue2 -= 2.0 * M_PI;
	  TmpValue2 *= Factor;
	  int TmpKx = round(TmpValue2);
	  if (TmpKx == MomentumModulo)
	    TmpKx = 0;	    
	  File << " " << TmpValue2 << " " << TmpKx << endl;
	}      
    }
  else
    {
      ComplexDiagonalMatrix Eigenvalues(NbrInputStates, true);
      ComplexMatrix Eigenstates(NbrInputStates, NbrInputStates);
#ifdef __LAPACK__
      KxRep.LapackDiagonalize(Eigenvalues, Eigenstates);
#endif
      double Factor = ((double) MomentumModulo) / (2.0 * M_PI);
      int* NbrStatePerKxSector = new int [MomentumModulo];
      int* KxValues = new int [NbrInputStates] ;
      for (int i = 0; i < MomentumModulo; ++i)
	NbrStatePerKxSector[i] = 0;
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  Complex TmpValue = Eigenvalues[i];
	  File << TmpValue << " " << Norm(TmpValue) << " " << Arg(TmpValue);
	  double TmpValue2 = Arg(TmpValue);
	  if (fabs(TmpValue2) < MACHINE_PRECISION)
	    TmpValue2 = 0.0;
	  if (TmpValue2 < 0.0)
	    TmpValue2 += 2.0 * M_PI;
	  if (TmpValue2 >= (2.0 * M_PI))
	    TmpValue2 -= 2.0 * M_PI;
	  TmpValue2 *= Factor;
	  int TmpKx = round(TmpValue2);
	  if (TmpKx == MomentumModulo)
	    TmpKx = 0;	    
	  File << " " << TmpValue2 << " " << TmpKx << endl;
	  KxValues[i] = TmpKx;
	  NbrStatePerKxSector[TmpKx]++;
	}
      ComplexVector** TmpVectors = new ComplexVector*[MomentumModulo];
      for (int i = 0; i < MomentumModulo; ++i)
	{
	  TmpVectors[i]  = new ComplexVector[NbrStatePerKxSector[i]];
	  NbrStatePerKxSector[i] = 0;
	}
      ComplexVector TmpVector (TotalSpace->GetHilbertSpaceDimension());
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  Complex TmpValue = Eigenvalues[i];
	  if (ComplexVectorFlag == false)
	    {
	      for (int j = 0; j < TotalSpace->GetHilbertSpaceDimension(); ++j)
		{
		  Complex Tmp = 0.0;
		  for (int k = 0; k < NbrInputStates; ++k)
		    {
		      Tmp += Conj(Eigenstates[i][k]) * InputStates[k][j];
		    }
		  TmpVector[j] = Tmp;		  
		}
	    }
	  else
	    {
	      for (int j = 0; j < TotalSpace->GetHilbertSpaceDimension(); ++j)
		{
		  Complex Tmp = 0.0;
		  for (int k = 0; k < NbrInputStates; ++k)
		    {
		      Tmp += Conj(Eigenstates[i][k]) * ComplexInputStates[k][j];
		    }
		  TmpVector[j] = Tmp;		  
		}
	    }
	  if (Manager.GetBoolean("no-convertion") == false)
	    TmpVectors[KxValues[i]][NbrStatePerKxSector[KxValues[i]]] = (TargetSpaces[KxValues[i]]->ConvertToKxKyBasis(TmpVector, TotalSpace));
	  else
	    TmpVectors[KxValues[i]][NbrStatePerKxSector[KxValues[i]]] = ComplexVector(TmpVector, true);
	  NbrStatePerKxSector[KxValues[i]]++;	      
	}	  
      
      ComplexMatrix* TransformationMatrices = new ComplexMatrix [MomentumModulo]; 
      for (int i = 0; i < MomentumModulo; ++i)
	{
	  ComplexMatrix TmpMatrix (TmpVectors[i], NbrStatePerKxSector[i]);
	  TmpMatrix.OrthoNormalizeColumns(TransformationMatrices[i]);
	  for (int j = 0; j < NbrStatePerKxSector[i]; ++j)
	    {
	      char* VectorOutputName = new char [256 + strlen(OutputNamePrefix)];
	      sprintf (VectorOutputName, "%s_kx_%d_ky_%d.%d.vec", OutputNamePrefix, i, YMomentum, j);
	      if (TmpMatrix[j].WriteVector(VectorOutputName) == false)
		{
		  cout << "error, can't write vector " << VectorOutputName << endl;
		}
	      delete[] VectorOutputName;	  	      
	    }
	  //	      delete[] TmpVectors[i];
	}	  
      if ((Manager.GetBoolean("export-transformation") == true) || 
	  (Manager.GetBoolean("export-bintransformation") == true))
	{
	  for (int i = 0; i < MomentumModulo; ++i)
	    {
	      if (NbrStatePerKxSector[i] > 0)
		{
		  ComplexMatrix TmpMatrix(NbrInputStates, NbrStatePerKxSector[i]);
		  int Tmp = 0;
		  for (int j = 0; j < NbrInputStates; ++j)
		    {
		      if (KxValues[j] == i)
			{
			  TmpMatrix[Tmp] = Eigenstates[j];
			  ++Tmp;
			}
		    }
		  TmpMatrix.ComplexConjugate();
		  TmpMatrix.Multiply(TransformationMatrices[i]);
		  char* MatrixOutputName =  new char [256 + strlen(OutputNamePrefix)];
		  if (Manager.GetBoolean("export-transformation") == true)
		    {
		      sprintf (MatrixOutputName, "%s_kx_%d_ky_%d.mat.txt", OutputNamePrefix, i, YMomentum);
		      TmpMatrix.WriteAsciiMatrix(MatrixOutputName);
		    }
		  else
		    {
		      sprintf (MatrixOutputName, "%s_kx_%d_ky_%d.mat", OutputNamePrefix, i, YMomentum);
		      TmpMatrix.WriteMatrix(MatrixOutputName);
		    }
		  delete[] MatrixOutputName;
		}
	    }
	}
      delete[] TmpVectors;
    }
  File.close();
  return 0;

}



// invert states from the (kx,ky) basis to the ky basis
//
// manager = pointer to the option manager
// return value = 0 if no error occured

int InvertU1KxKyStates(OptionManager* manager)
{
  int NbrParticles = 0;
  int MaxMomentum = 0;
  int YMomentum = -1;
  bool Statistics = true;
  bool SU2Flag = false;
  bool SU3Flag = false;
  bool SU4Flag = false;
  int XMomentum = -1;
  ComplexVector* InputStates = 0;
  int NbrInputStates = 0;
  char* OutputNamePrefix = new char [256 + strlen(manager->GetString("interaction-name"))];
  char** OutputFileNames = 0;
  if (manager->GetString("degenerated-states") == 0)
    {
      if (manager->GetString("input-state") == 0)
	{
	  cout << "error, either input-state or degenerated-states has to be provided" << endl;
	  return -1;
	}
      FQHEOnSphereFindInternalSymmetryGroupFromFileName(manager->GetString("input-state"), SU2Flag, SU3Flag, SU4Flag);
      if (SU2Flag == true)
	{
	  return InvertSU2KxKyStates(manager);
	}
      
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(manager->GetString("input-state"),
						      NbrParticles, MaxMomentum, XMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << manager->GetString("ground-state") << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Kx=" << XMomentum << " Ky=" << YMomentum << " ";
      if (Statistics == false)
	{
	  sprintf (OutputNamePrefix, "bosons_torus_kysym_%s_n_%d_2s_%d_ratio_%.6f", manager->GetString("interaction-name"), NbrParticles, MaxMomentum, manager->GetDouble("ratio"));
	}
      else
	{
	  sprintf (OutputNamePrefix, "fermions_torus_kysym_%s_n_%d_2s_%d_ratio_%.6f", manager->GetString("interaction-name"), NbrParticles, MaxMomentum, manager->GetDouble("ratio"));
	}
      NbrInputStates = 1;
      InputStates = new ComplexVector [NbrInputStates];
      if (InputStates[0].ReadVector(manager->GetString("input-state")) == false)
	{
	  cout << "error while reading " << manager->GetString("input-state") << endl;
	  return -1;
	}
      OutputFileNames = new char* [1];
      OutputFileNames[0] = ReplaceString(manager->GetString("input-state"), "_torus_", "_torus_kysym_");
      if (OutputFileNames[0] == 0)
	{
	  OutputFileNames[0] = new char [256 + strlen(OutputNamePrefix)];
	  sprintf (OutputFileNames[0], "%s_kx_%d_ky_%d.0.vec", OutputNamePrefix, XMomentum, YMomentum);
	  cout << "warning, can't deduce output file name from " << manager->GetString("input-state") << " will use " << OutputFileNames[0] << " instead" << endl;
	}
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(manager->GetString("degenerated-states")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	} 
      NbrInputStates = DegeneratedFile.GetNbrLines();
      if (NbrInputStates < 1)
	{
	  cout << "no state found in " << manager->GetString("degenerated-states") << endl;
	}
      InputStates = new ComplexVector [NbrInputStates];
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, 0),
						      NbrParticles, MaxMomentum, XMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, 0) << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Ky=" << YMomentum << " ";
      if (InputStates[0].ReadVector(DegeneratedFile(0, 0)) == false)
	{
	  cout << "error while reading " << DegeneratedFile(0, 0) << endl;
	  return -1;
	}
      if (Statistics == false)
	{
	  sprintf (OutputNamePrefix, "bosons_torus_kysym_%s_n_%d_2s_%d_ratio_%.6f", manager->GetString("interaction-name"), NbrParticles, MaxMomentum, manager->GetDouble("ratio"));
	}
      else
	{
	  sprintf (OutputNamePrefix, "fermions_torus_kysym_%s_n_%d_2s_%d_ratio_%.6f", manager->GetString("interaction-name"), NbrParticles, MaxMomentum, manager->GetDouble("ratio"));
	}
      OutputFileNames = new char* [NbrInputStates];
      OutputFileNames[0] = ReplaceString(manager->GetString("input-state"), "_torus_", "_torus_kysym_");
      if (OutputFileNames[0] == 0)
	{
	  OutputFileNames[0] = new char [256 + strlen(OutputNamePrefix)];
	  sprintf (OutputFileNames[0], "%s_kx_%d_ky_%d.0.vec", OutputNamePrefix, XMomentum, YMomentum);
	  cout << "warning, can't deduce output file name from " << manager->GetString("input-state") << " will use " << OutputFileNames[0] << " instead" << endl;
	}
      for (int i = 1; i < NbrInputStates; ++i)
	{
	  int TmpNbrParticles = 0;
	  int TmpMaxMomentum = 0;
	  int TmpYMomentum = -1;
	  int TmpXMomentum = -1;
	  bool TmpStatistics = true;
	  if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, i),
							  TmpNbrParticles, TmpMaxMomentum, TmpXMomentum, TmpYMomentum, TmpStatistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, i) << endl;
	      return -1;
	    }
	  if ((TmpNbrParticles != NbrParticles) || (TmpMaxMomentum != MaxMomentum) || 
	      (TmpYMomentum != YMomentum) || (TmpXMomentum != XMomentum) || (Statistics != TmpStatistics))
	    {
	      cout << "error, " << DegeneratedFile(0, i) << " has different system parameters than " << DegeneratedFile(0, 0) 
		   << ", N=" << TmpNbrParticles << "(" << NbrParticles << "), N_phi=" << TmpMaxMomentum << "(" << MaxMomentum 
		   << "), Kx=" << TmpXMomentum << "(" << XMomentum << "), Ky=" << TmpYMomentum << "(" << YMomentum << ")" << endl;
	    }
	  if (InputStates[i].ReadVector(DegeneratedFile(0, i)) == false)
	    {
	      cout << "error while reading " << DegeneratedFile(0, i) << endl;
	      return -1;
	    }
	  if (InputStates[i].GetVectorDimension() != InputStates[0].GetVectorDimension())
	    {
	      cout << "error, " << DegeneratedFile(0, i) << " has different dimension than " << DegeneratedFile(0, 0) << endl;	      
	    }
	  OutputFileNames[i] = ReplaceString(manager->GetString("input-state"), "_torus_", "_torus_kysym_");
	  if (OutputFileNames[i] == 0)
	    {
	      OutputFileNames[i] = new char [256 + strlen(OutputNamePrefix)];
	      sprintf (OutputFileNames[i], "%s_kx_%d_ky_%d.%d.vec", OutputNamePrefix, XMomentum, YMomentum, i);
	      cout << "warning, can't deduce output file name from " << manager->GetString("input-state") << " will use " << OutputFileNames[0] << " instead" << endl;
	    }
	}
    }
  ParticleOnTorus* TotalSpace = 0;
  ParticleOnTorusWithMagneticTranslations* TargetSpace = 0;
  if (Statistics == false)
    {
      TargetSpace = new BosonOnTorusWithMagneticTranslationsShort(NbrParticles, MaxMomentum, XMomentum, YMomentum);
      TotalSpace = new BosonOnTorusShort(NbrParticles, MaxMomentum, YMomentum);
    }
  else
    {
      TotalSpace = new FermionOnTorus (NbrParticles, MaxMomentum, YMomentum);
      TargetSpace = new FermionOnTorusWithMagneticTranslations(NbrParticles, MaxMomentum, XMomentum, YMomentum);
    }
  
  if (InputStates[0].GetLargeVectorDimension() != TargetSpace->GetLargeHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputStates[0].GetLargeVectorDimension() 
	   << " "<< TargetSpace->GetLargeHilbertSpaceDimension() << endl;
      return -1;
    }
  for (int i = 0; i < NbrInputStates; ++i)
    {
      ComplexVector TmpVector = TargetSpace->ConvertFromKxKyBasis(InputStates[i], TotalSpace);
      if (manager->GetBoolean("invert-real") == false)
	{
	  if (TmpVector.WriteVector(OutputFileNames[i]) == false)
	    {
	      cout << "error, can't write vector " << OutputFileNames[i] << endl;
	    }
	}
      else
	{
	  Complex TmpPhase = TmpVector.GlobalPhase();
	  cout << "global phase = " << TmpPhase << endl;
	  TmpVector /= TmpPhase;
	  RealVector TmpVector2(TmpVector);
	  TmpVector2.Normalize();
	  if (TmpVector2.WriteVector(OutputFileNames[i]) == false)
	    {
	      cout << "error, can't write vector " << OutputFileNames[i] << endl;
	    }	
	}
    }
  return 0;
}

// invert states from the (kx,ky) basis to the ky basis with SU(2) degree of freedom
//
// manager = pointer to the option manager
// return value = 0 if no error occured

int InvertSU2KxKyStates(OptionManager* manager)
{
  int NbrParticles = 0;
  int MaxMomentum = 0;
  int YMomentum = -1;
  int TotalSz = 0;
  bool Statistics = true;
  bool SU2Flag = false;
  bool SU3Flag = false;
  bool SU4Flag = false;
  int XMomentum = -1;
  ComplexVector* InputStates = 0;
  int NbrInputStates = 0;
  char* OutputNamePrefix = new char [256 + strlen(manager->GetString("interaction-name"))];
  if (manager->GetString("degenerated-states") == 0)
    {
      if (manager->GetString("input-state") == 0)
	{
	  cout << "error, either input-state or degenerated-states has to be provided" << endl;
	  return -1;
	}
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(manager->GetString("input-state"),
						      NbrParticles, MaxMomentum, XMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << manager->GetString("ground-state") << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Kx=" << XMomentum << " Ky=" << YMomentum << " ";
      NbrInputStates = 1;
      InputStates = new ComplexVector [NbrInputStates];
      if (InputStates[0].ReadVector(manager->GetString("input-state")) == false)
	{
	  cout << "error while reading " << manager->GetString("input-state") << endl;
	  return -1;
	}
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(manager->GetString("degenerated-states")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	} 
      NbrInputStates = DegeneratedFile.GetNbrLines();
      if (NbrInputStates < 1)
	{
	  cout << "no state found in " << manager->GetString("degenerated-states") << endl;
	}
      InputStates = new ComplexVector [NbrInputStates];
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, 0),
						      NbrParticles, MaxMomentum, XMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, 0) << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Ky=" << YMomentum << " ";
      if (InputStates[0].ReadVector(DegeneratedFile(0, 0)) == false)
	{
	  cout << "error while reading " << DegeneratedFile(0, 0) << endl;
	  return -1;
	}
      for (int i = 1; i < NbrInputStates; ++i)
	{
	  int TmpNbrParticles = 0;
	  int TmpMaxMomentum = 0;
	  int TmpYMomentum = -1;
	  int TmpXMomentum = -1;
	  bool TmpStatistics = true;
	  if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, i),
							  TmpNbrParticles, TmpMaxMomentum, TmpXMomentum, TmpYMomentum, TmpStatistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, i) << endl;
	      return -1;
	    }
	  if ((TmpNbrParticles != NbrParticles) || (TmpMaxMomentum != MaxMomentum) || 
	      (TmpYMomentum != YMomentum) || (TmpXMomentum != XMomentum) || (Statistics != TmpStatistics))
	    {
	      cout << "error, " << DegeneratedFile(0, i) << " has different system parameters than " << DegeneratedFile(0, 0) 
		   << ", N=" << TmpNbrParticles << "(" << NbrParticles << "), N_phi=" << TmpMaxMomentum << "(" << MaxMomentum 
		   << "), Kx=" << TmpXMomentum << "(" << XMomentum << "), Ky=" << TmpYMomentum << "(" << YMomentum << ")" << endl;
	    }
	  if (InputStates[i].ReadVector(DegeneratedFile(0, i)) == false)
	    {
	      cout << "error while reading " << DegeneratedFile(0, i) << endl;
	      return -1;
	    }
	  if (InputStates[i].GetVectorDimension() != InputStates[0].GetVectorDimension())
	    {
	      cout << "error, " << DegeneratedFile(0, i) << " has different dimension than " << DegeneratedFile(0, 0) << endl;	      
	    }
	}
    }
  ParticleOnSphereWithSpin* TotalSpace = 0;
  ParticleOnTorusWithSpinAndMagneticTranslations* TargetSpace = 0;
  if (Statistics == false)
    {
      TargetSpace = new BosonOnTorusWithSpinAndMagneticTranslations(NbrParticles, TotalSz, MaxMomentum, XMomentum, YMomentum);
      TotalSpace = new BosonOnTorusWithSpin(NbrParticles, MaxMomentum, TotalSz, YMomentum);
      sprintf (OutputNamePrefix, "bosons_torus_kysym_%s_n_%d_2s_%d_ratio_%.6f_sz_%d", manager->GetString("interaction-name"), NbrParticles, MaxMomentum, manager->GetDouble("ratio"), TotalSz);
    }
  else
    {
      TotalSpace = new FermionOnTorusWithSpin (NbrParticles, MaxMomentum, TotalSz, YMomentum);
      TargetSpace = new FermionOnTorusWithSpinAndMagneticTranslations(NbrParticles, TotalSz, MaxMomentum, XMomentum, YMomentum);
      sprintf (OutputNamePrefix, "fermions_torus_kysym_%s_n_%d_2s_%d_ratio_%.6f_sz_%d", manager->GetString("interaction-name"), NbrParticles, MaxMomentum, manager->GetDouble("ratio"), TotalSz);
    }
  
  if (InputStates[0].GetVectorDimension() != TargetSpace->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputStates[0].GetVectorDimension() 
	   << " "<< TargetSpace->GetHilbertSpaceDimension() << endl;
      return -1;
    }
  for (int i = 0; i < NbrInputStates; ++i)
    {
      char* VectorOutputName = new char [256 + strlen(OutputNamePrefix)];
      sprintf (VectorOutputName, "%s_kx_%d_ky_%d.%d.vec", OutputNamePrefix, XMomentum, YMomentum, i);
      ComplexVector TmpVector = TargetSpace->ConvertFromKxKyBasis(InputStates[i], TotalSpace);
      if (manager->GetBoolean("invert-real") == false)
	{
	  if (TmpVector.WriteVector(VectorOutputName) == false)
	    {
	      cout << "error, can't write vector " << VectorOutputName << endl;
	    }
	}
      else
	{
	  Complex TmpPhase = TmpVector.GlobalPhase();
	  TmpVector /= TmpPhase;
	  RealVector TmpVector2(TmpVector);
	  TmpVector2.Normalize();
	  if (TmpVector2.WriteVector(VectorOutputName) == false)
	    {
	      cout << "error, can't write vector " << VectorOutputName << endl;
	    }	      
	}
      delete[] VectorOutputName;
    }
  return 0;
}
