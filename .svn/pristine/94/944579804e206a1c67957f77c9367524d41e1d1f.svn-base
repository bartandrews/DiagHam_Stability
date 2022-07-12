#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereFullShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/BosonOnSphereLong.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisLong.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereBosonsFullEntanglementEntropy" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-la", "minimum size of the subsystem whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-la", "maximum size of the subsystem whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices to a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "na-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed number of particles", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*OutputGroup) += new SingleStringOption ('\n', "full-densitymatrix", "store full density matrices to a given file");
  (*PrecalculationGroup) += new BooleanOption  ('\n', "load-hilbert", "load Hilbert space description from a file (only available for the Haldane basis)");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsFullEntanglementEntropy -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("ground-file") == 0)
    {
      cout << "error, a ground state decription should be provided. See man page for option syntax or type FQHESphereBosonsFullEntanglementEntropy -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("ground-file")) == false)
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }


  int NbrParticles = 0;
  int LzMax = 0;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  int FilterNa = Manager.GetInteger("na-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  int* TotalLz = 0;
  bool Statistics = true;
  int NbrSpaces = 1;
  BosonOnSphereShort** Spaces = 0;
  RealVector* GroundStates = 0;
  char** GroundStateFiles = 0;

  MultiColumnASCIIFile DegeneratedFile;
  if (DegeneratedFile.Parse(Manager.GetString("ground-file")) == false)
    {
      DegeneratedFile.DumpErrors(cout);
      return -1;
    }
  NbrSpaces = DegeneratedFile.GetNbrLines();
  GroundStateFiles = new char* [NbrSpaces];
  TotalLz = new int[NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
      strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalLz[i] = 0;
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(GroundStateFiles[i],
						       NbrParticles, LzMax, TotalLz[i], Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
      if (Statistics == true)
	{
	  cout << GroundStateFiles[i] << " is not a bosonic state" << endl;
	  return -1;
	}
      if (((NbrParticles * LzMax) & 1) != (TotalLz[i] & 1))
	{
	  cout << "incompatible values for nbr-particles, nbr-flux and total-lz for ground state file " << GroundStateFiles[i] << endl;
	  return -1;
	}
    }

  double* RealGroundStateCoefficients = DegeneratedFile.GetAsDoubleArray(1);
  GroundStates = new RealVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
      {
	cout << "can't open vector file " << GroundStateFiles[i] << endl;
	return -1;      
      }

  Spaces = new BosonOnSphereShort* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (Manager.GetBoolean("haldane") == false)
	{
	  Spaces[i] = new BosonOnSphereShort (NbrParticles, TotalLz[i], LzMax);
	}
      else
	{
	  int* ReferenceState = 0;
	  if (Manager.GetString("reference-file") == 0)
	    {
	      cout << "error, a reference file is needed" << endl;
	      return 0;
	    }
	  ConfigurationParser ReferenceStateDefinition;
	  if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
	    {
	      ReferenceStateDefinition.DumpErrors(cout) << endl;
	      return 0;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
	    {
	      cout << "NbrParticles is not defined or as a wrong value" << endl;
	      return 0;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax < 0))
	    {
	      cout << "LzMax is not defined or as a wrong value" << endl;
	      return 0;
	    }
	  int MaxNbrLz;
	  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	    {
	      cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
	      return 0;     
	    }
	  if (MaxNbrLz != (LzMax + 1))
	    {
	      cout << "wrong LzMax value in ReferenceState" << endl;
	      return 0;     
	    }
	  if (Manager.GetBoolean("load-hilbert") == true)
	    {
	      Spaces[i] = new BosonOnSphereHaldaneBasisShort(Manager.GetString("load-hilbert"));
	    }
	  else
	    {
	      Spaces[i] = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz[i], LzMax, ReferenceState);	  
	    }
	}
      if (Spaces[i]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }


  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "# l_a    N    <Lz>    lambda" << endl;
      DensityMatrixFile.close();
    }

  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "ent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);
  int MeanSubsystemSize = LzMax >> 1;
  if ((LzMax & 1) != 0)
    ++MeanSubsystemSize;
  if (Manager.GetInteger("max-la") > 0)
    {
      MeanSubsystemSize = Manager.GetInteger("max-la");
      if (MeanSubsystemSize > LzMax)
	MeanSubsystemSize = LzMax;
    }
  int SubsystemSize = Manager.GetInteger("min-la");
  if (SubsystemSize < 1)
    SubsystemSize = 1;
  BinomialCoefficients Coefs(MeanSubsystemSize + NbrParticles - 1);
  for (; SubsystemSize <= MeanSubsystemSize; ++SubsystemSize)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      int MaxSubsystemNbrParticles = NbrParticles;
      if (MaxSubsystemNbrParticles > SubsystemSize)
	MaxSubsystemNbrParticles = SubsystemSize;
      long MaximumSize = 0;
      for (int i = 0; i <= NbrParticles; ++i)
	MaximumSize += Coefs(SubsystemSize + i - 1, i);
      double* TmpDensityMatrixEigenvalues = new double [MaximumSize];
      long TmpDensityMatrixEigenvaluePosition = 0;
      for (int SubsystemNbrParticles = 0; SubsystemNbrParticles <= NbrParticles; ++SubsystemNbrParticles)
	{
	  cout << "processing subsystem size=" << SubsystemSize << "  subsystem nbr of particles=" << SubsystemNbrParticles << endl;
	  BosonOnSphereFullShort TmpSpace (SubsystemNbrParticles, SubsystemSize - 1);
	  RealSymmetricMatrix PartialDensityMatrix(TmpSpace.GetHilbertSpaceDimension(), true);
	  RealMatrix TmpPartialDensityMatrix(TmpSpace.GetHilbertSpaceDimension(), true);
	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      TmpSpace.EvaluatePartialDensityMatrix(TmpPartialDensityMatrix, Spaces[i], GroundStates[i], Spaces[i], GroundStates[i]);
	      TmpPartialDensityMatrix *= (0.5 * RealGroundStateCoefficients[i] * RealGroundStateCoefficients[i]);
	      PartialDensityMatrix += TmpPartialDensityMatrix;
	      for (int j = i + 1; j < NbrSpaces; ++j)
		{
		  TmpSpace.EvaluatePartialDensityMatrix(TmpPartialDensityMatrix, Spaces[i], GroundStates[i], Spaces[j], GroundStates[j]);
		  TmpPartialDensityMatrix *= (RealGroundStateCoefficients[i] * RealGroundStateCoefficients[j]);
		  PartialDensityMatrix += TmpPartialDensityMatrix;
		}
	    }
	  if ((Manager.GetString("full-densitymatrix") != 0) && (FilterNa == SubsystemNbrParticles))
	    {
	      ofstream FullDensityMatrixFile;
	      FullDensityMatrixFile.open(Manager.GetString("full-densitymatrix"), ios::binary | ios::out | ios::app); 
	      FullDensityMatrixFile << PartialDensityMatrix;
	      FullDensityMatrixFile.close();
	    }
	  if (PartialDensityMatrix.GetNbrRow() > 1)
	    {
	      RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
	      if (LapackFlag == true)
		{
// 		  if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles))
// 		    {
// 		      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
// 							PartialDensityMatrix.GetNbrRow(), true);
// 		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
// 			TmpEigenstates[i][i] = 1.0;
// 		      PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
// 		      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
// 		      char* TmpEigenstateName = new char[512];
// 		      int MaxNbrEigenstates = NbrEigenstates;
// 		      if (NbrEigenstates == 0)
// 			MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
// 		      for (int i = 0; i < MaxNbrEigenstates; ++i)
// 			{
// 			  if (TmpDiag[i] > 1e-14)
// 			    {
// 			      sprintf (TmpEigenstateName,
// 				       "bosons_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
// 				       NbrParticles, LzMax, TotalLz[0], SubsystemSize,
// 				       SubsystemNbrParticles, SubsystemTotalLz, i);
// 			      TmpEigenstates[i].WriteVector(TmpEigenstateName);
// 			    }
// 			}
// 		      delete[] TmpEigenstateName;
// 		    }
// 		  else
//		    {
		      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		      //		    }
		}
	      else
		{
// 		  if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles))
// 		    {
// 		      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
// 						PartialDensityMatrix.GetNbrRow(), true);
// 		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
// 			TmpEigenstates[i][i] = 1.0;
// 		      PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
// 		      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
// 		      char* TmpEigenstateName = new char[512];
// 		      int MaxNbrEigenstates = NbrEigenstates;
// 		      if (NbrEigenstates == 0)
// 			MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
// 		      for (int i = 0; i < MaxNbrEigenstates; ++i)
// 			{
// 			  if (TmpDiag[i] > 1e-14)
// 			    {
// 			      sprintf (TmpEigenstateName,
// 				       "bosons_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
// 				       NbrParticles, LzMax, TotalLz[0], SubsystemSize,
// 				       SubsystemNbrParticles, SubsystemTotalLz, i);
// 			      TmpEigenstates[i].WriteVector(TmpEigenstateName);
// 			    }
// 			}
// 		      delete[] TmpEigenstateName;
// 		    }
// 		  else
// 		    {
		      PartialDensityMatrix.Diagonalize(TmpDiag);
		      //		    }
		}
#else
	      PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
	      TmpDiag.SortMatrixDownOrder();
	      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		TmpDensityMatrixEigenvalues[TmpDensityMatrixEigenvaluePosition++] = TmpDiag[i];
	      if (DensityMatrixFileName != 0)
		{
		  ofstream DensityMatrixFile;
		  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		  DensityMatrixFile.precision(14);
		  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		    DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << " " << TmpDiag[i] << endl;
		  DensityMatrixFile.close();
		}
	    }
	  else
	    if (PartialDensityMatrix.GetNbrRow() == 1)
	      {
		double TmpValue = PartialDensityMatrix(0,0);
		TmpDensityMatrixEigenvalues[TmpDensityMatrixEigenvaluePosition++] = TmpValue;
		if (DensityMatrixFileName != 0)
		  {
		    ofstream DensityMatrixFile;
		    DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		    DensityMatrixFile.precision(14);
		    DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << " " << TmpValue << endl;
		    DensityMatrixFile.close();
		  }		  
	      }
	}
      EntanglementEntropy = 0.0;
      DensitySum = 0.0;
      cout << "sorting density matrix eigenvalues and computing entanglement entropy" << endl;
      SortArrayDownOrdering(TmpDensityMatrixEigenvalues, TmpDensityMatrixEigenvaluePosition);
      unsigned TmpPos = 0;
      for (; (TmpPos < TmpDensityMatrixEigenvaluePosition) && (DensitySum < 1.0); ++TmpPos)
	{
	  if (TmpDensityMatrixEigenvalues[TmpPos] > 1e-14)
	    {
	      EntanglementEntropy += TmpDensityMatrixEigenvalues[TmpPos] * log(TmpDensityMatrixEigenvalues[TmpPos]);
	      DensitySum += TmpDensityMatrixEigenvalues[TmpPos];
	    }
	}
      double DensitySumError = 0.0;
      for (; TmpPos < TmpDensityMatrixEigenvaluePosition; ++TmpPos)
	DensitySumError += TmpDensityMatrixEigenvalues[TmpPos];
      delete[] TmpDensityMatrixEigenvalues;
      File << SubsystemSize << " " << (-EntanglementEntropy) << " " << DensitySum << " " << DensitySumError << endl;
    }
  File.close();
}

