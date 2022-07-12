#include "Options/Options.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceFixedParity.h"

#include "Hamiltonian/ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/ArrayTools.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"


#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("HubbardFloquetHeisenbergModel" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sites", "number of sites", 4);
  (*SystemGroup) += new  SingleDoubleOption('\n', "tau", "evolution time", 1.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "multiple-taus", "compute the specrtum of the evolution operator at several times in a single shot");
  (*SystemGroup) += new  SingleDoubleOption ('\n', "j-xx", "coupling constant between neighboring sites along the x axis in the spin language", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "j-yy", "coupling constant between neighboring sites along the y axis in the spin language", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "j-zz", "coupling constant between neighboring sites along the z axis in the spin language", 1.0);
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hz-value", "amplitude of the Zeeman term along the z axis", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-hzvalue", "amplitude of the random Zeeman term on each site", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-gaussianhzvalue", "amplitude of the random Zeeman term on each site, using a gaussian disrtibution with zero mean value and a given standard deviation", 0.0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "run-id", "add an additional run id to the file name when using the --random-hzvalue option", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "fullhz-values", "name of the file that contains the Zeeman term amplitudes for each site");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleDoubleOption ('\n', "testhermitian-error", "error threshold when testing hermiticy (0 for machine accuracy)", 0.0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardFloquetHeisenbergModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSites = Manager.GetInteger("nbr-sites"); 
 
  double Tau = Manager.GetDouble("tau");
 
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [64];
  sprintf (StatisticPrefix, "fermions_floquet_heisenberg"); 

  char* BoundaryName = new char [16];
  if (Manager.GetBoolean("use-periodic") == false)
    sprintf (BoundaryName, "open");
  else
    sprintf (BoundaryName, "closed");
  
  char* FilePrefix = new char [256];
  sprintf (FilePrefix, "%s_%s_ns_%d", StatisticPrefix, BoundaryName, NbrSites);
  
  char* TmpFileParameterString = new char [256];
  char* FileParameterString = 0;
  if ((Manager.GetDouble("hz-value") == 0.0) && (Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
    {
      sprintf (TmpFileParameterString, "jxx_%.6f_jyy_%.6f_jzz_%.6f", Manager.GetDouble("j-xx"), Manager.GetDouble("j-yy"), Manager.GetDouble("j-zz"));
    }
  else
    {
      if ((Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	{
	  sprintf (TmpFileParameterString, "jxx_%.6f_jyy_%.6f_jzz_%.6f_hz_%.6f", 
		   Manager.GetDouble("j-xx"), Manager.GetDouble("j-yy"), Manager.GetDouble("j-zz"), Manager.GetDouble("hz-value"));
	}
      else
	{
	  if (Manager.GetDouble("random-gaussianhzvalue") == 0.0)
	    {
	      sprintf (TmpFileParameterString, "jxx_%.6f_jyy_%.6f_jzz_%.6f_hz_%.6f_randomhz_%.6f", 
		       Manager.GetDouble("j-xx"), Manager.GetDouble("j-yy"), Manager.GetDouble("j-zz"), 
		       Manager.GetDouble("hz-value"), Manager.GetDouble("random-hzvalue"));
	    }
	  else
	    {
	      sprintf (TmpFileParameterString, "jxx_%.6f_jyy_%.6f_jzz_%.6f_hz_%.6f_gaussianrandomhz_%.6f", 
		       Manager.GetDouble("j-xx"), Manager.GetDouble("j-yy"), Manager.GetDouble("j-zz"), 
		       Manager.GetDouble("hz-value"), Manager.GetDouble("random-gaussianhzvalue"));
	    }
	}
    }
  if (Tau == 0.0)
    {
      FileParameterString = new char[strlen(TmpFileParameterString) + 64];
      if ((Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	{
	  strcpy(FileParameterString, TmpFileParameterString);
	}
      else
	{
	  sprintf (FileParameterString, "%s_runid_%ld", TmpFileParameterString, Manager.GetInteger("run-id"));
	}  
    }
  else
    {
      FileParameterString = new char[strlen(TmpFileParameterString) + 64];
      if ((Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	{
	  sprintf (FileParameterString, "%s_tau_%.6f", TmpFileParameterString, Tau);
	}
      else
	{
	  sprintf (FileParameterString, "%s_tau_%.6f_runid_%ld", TmpFileParameterString, Tau, Manager.GetInteger("run-id"));
	}
    }
  //%.6f_tK_%.6f_tKphi_%.6f_j1_%.6f_j2_%.6f", Manager.GetDouble("isotropic-t"), Manager.GetDouble("anisotropic-t"), Manager.GetDouble("anisotropic-phase"), Manager.GetDouble("j1"), Manager.GetDouble("j2"));

  char* CommentLine = new char [256];
  sprintf (CommentLine, "Z2 ");

  char* EigenvalueOutputFile = new char [512];
  sprintf(EigenvalueOutputFile, "%s_%s.dat", FilePrefix, FileParameterString);

  bool FirstRunFlag = true;

  double* HzValues = 0;
  if (Manager.GetString("fullhz-values") != 0)
    {
      MultiColumnASCIIFile HFieldFile;
      if (HFieldFile.Parse(Manager.GetString("fullhz-values")) == false)
	{
	  HFieldFile.DumpErrors(cout);
	  return -1;
	}
      if (HFieldFile.GetNbrLines() == NbrSites)
	{
	  HzValues = HFieldFile.GetAsDoubleArray(0);
	}
      else
	{
	  if (HFieldFile.GetNbrLines() > NbrSites)
	    {
	      cout << "warning, " << Manager.GetString("fullhz-values") << " has more hz values than the number of sites" << endl;
	      HzValues = HFieldFile.GetAsDoubleArray(0);
	    }
	  else
	    {
	      cout << "error, " << Manager.GetString("fullhz-values") << " has less hz values than the number of sites" << endl;
	      return 0;
	    }	  
	}
    }
  else
    {
      if ((Manager.GetDouble("hz-value") != 0.0) || (Manager.GetDouble("random-hzvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	{
	  HzValues = new double [NbrSites];
	  HzValues[0] = Manager.GetDouble("hz-value");
	  for (int i = 1; i < NbrSites; ++i)
	    HzValues[i] = HzValues[0];
	  if ((Manager.GetDouble("random-hzvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	    {
	      AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	      RandomNumber->UseTimeSeed();
	      char* TmpExtention = new char [16];
	      sprintf (TmpExtention, ".hzvalues");
	      char* HzOutputFileName = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	      ofstream File;
	      File.open(HzOutputFileName, ios::binary | ios::out); 
	      File.precision(14); 
	      for (int i = 0; i < NbrSites; ++i)
		{
		  double Tmp;
		  if (Manager.GetDouble("random-hzvalue") != 0.0)
		    {
		      Tmp = Manager.GetDouble("random-hzvalue") * (2.0 * RandomNumber->GetRealRandomNumber() - 1.0);
		    }
		  else
		    {
		      Tmp = RandomNumber->GetGaussianRandomNumber(0.0, Manager.GetDouble("random-gaussianhzvalue"));
		    }
		  HzValues[i] += Tmp;
		  File << Tmp << endl;
		}
	      File.close();	      
	    }
	}
    }


  double Jxx = Manager.GetDouble("j-xx");
  double Jyy = Manager.GetDouble("j-yy");
  double Jzz = Manager.GetDouble("j-zz");
      
  if ((Tau == 0.0) && (Manager.GetString("multiple-taus") == 0))
    {
      int ParitySector = 0;
      int MaxParitySector = 1;
      for (; ParitySector <= MaxParitySector; ++ParitySector)
	{
	  ParticleOnSphere* Space = 0;
	  AbstractHamiltonian* Hamiltonian = 0;
	  Space = new FermionOnLatticeRealSpaceFixedParity (NbrSites, ParitySector);
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  Hamiltonian = new ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian(Space, ParitySector, Jxx, Jyy, Jzz, HzValues, Manager.GetBoolean("use-periodic"),
											      Architecture.GetArchitecture(), Memory);
	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d", ParitySector);
	  char* EigenstateOutputFile;
	  char* TmpExtention = new char [512];
	  sprintf (TmpExtention, "_");
	  EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	  
	  GenericRealMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	  FirstRunFlag = false;
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  cout << "------------------------------------" << endl;
	  delete Hamiltonian;
	  delete Space;
	  delete[] EigenstateOutputFile;
	  delete[] ContentPrefix;
	}  
    }
  else
    {
      int NbrTauValues  = 1;
      double* TauValues = 0;
      char** TauOutputFileNames = 0;
      if (Manager.GetString("multiple-taus") != 0)
	{
	  MultiColumnASCIIFile TauFile;
	  if (TauFile.Parse(Manager.GetString("multiple-taus")) == false)
	    {
	      TauFile.DumpErrors(cout);
	      return false;
	    }
	  NbrTauValues = TauFile.GetNbrLines();
	  TauValues = TauFile.GetAsDoubleArray(0);
	  TauOutputFileNames = new char* [NbrTauValues];
	  for (int i = 0; i < NbrTauValues; ++i)
	    {
	      FileParameterString = new char[strlen(TmpFileParameterString) + 64];
	      if ((Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
		{
		  sprintf (FileParameterString, "%s_tau_%.6f", TmpFileParameterString, TauValues[i]);
		}
	      else
		{
		  sprintf (FileParameterString, "%s_tau_%.6f_runid_%ld", TmpFileParameterString, TauValues[i], Manager.GetInteger("run-id"));
		}
	      TauOutputFileNames[i] = new char [512];
	      sprintf(TauOutputFileNames[i], "%s_%s.dat", FilePrefix, FileParameterString);
	      delete[] FileParameterString;
	    }
	}
      else
	{
	  TauValues = new double[1];
	  TauValues[0] = Tau;
	  TauOutputFileNames = new char* [1];
	  TauOutputFileNames[0] = new char [512];
	  sprintf(TauOutputFileNames[0], "%s_%s.dat", FilePrefix, FileParameterString);
	}

      for (int i = 0; i < NbrTauValues; ++i)
	{
	  ofstream File;
	  File.open(TauOutputFileNames[i], ios::out);
	  File.precision(14);  
	  File << "# Z2 lambda norm(lambda) phase(lambda)" << endl;
	  File.close();
	}
      ComplexDiagonalMatrix TmpDiagUnitaryEvolution;
      int ParitySector = 0;
      int MaxParitySector = 1;
      for (; ParitySector <= MaxParitySector; ++ParitySector)
	{
	  ParticleOnSphere* Space = new FermionOnLatticeRealSpaceFixedParity (NbrSites, ParitySector);
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  double Dt;
	  cout << "building H1 Hamiltonian" <<  endl;
	  gettimeofday (&(TotalStartingTime), 0);
	  int StartTimeSecond = TotalStartingTime.tv_sec;
	  ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian* Hamiltonian1 = new ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian(Space, ParitySector, 0.0, 0.0, Jzz, HzValues, Manager.GetBoolean("use-periodic"),
																				  Architecture.GetArchitecture(), Memory);
	  RealSymmetricMatrix HRep1 (Hamiltonian1->GetHilbertSpaceDimension(), true);
	  RealDiagonalMatrix TmpDiag1 (Hamiltonian1->GetHilbertSpaceDimension());
	  Hamiltonian1->GetHamiltonian(HRep1);
	  for (int i = 0; i < HRep1.GetNbrRow(); ++i)
	    {
	      HRep1.GetMatrixElement(i, i, TmpDiag1[i]);
	    }
	  delete Hamiltonian1;
	  gettimeofday (&(TotalEndingTime), 0);
	  Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			 ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "done in " << Dt << "sec" << endl;
	  
	  cout << "building H2 Hamiltonian" <<  endl;
	  gettimeofday (&(TotalStartingTime), 0);
	  ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian* Hamiltonian2 = new ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian(Space, ParitySector, 0.0, Jyy, 0.0, 0, Manager.GetBoolean("use-periodic"),
																				  Architecture.GetArchitecture(), Memory);
	  RealMatrix Basis2(Hamiltonian2->GetHilbertSpaceDimension(), Hamiltonian2->GetHilbertSpaceDimension());
	  RealSymmetricMatrix HRep2 (Hamiltonian2->GetHilbertSpaceDimension(), true);
	  RealDiagonalMatrix TmpDiag2 (Hamiltonian2->GetHilbertSpaceDimension());
	  Hamiltonian2->GetHamiltonian(HRep2);
#ifdef __LAPACK__     
	  HRep2.LapackDiagonalize(TmpDiag2, Basis2);
#else
	  HRep2.Diagonalize(TmpDiag2, Basis2);
#endif	  
	  delete Hamiltonian2;
	  RealMatrix TransposedBasis2;
	  TransposedBasis2.Copy(Basis2);
	  TransposedBasis2.Transpose();
	  gettimeofday (&(TotalEndingTime), 0);
	  Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			 ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "done in " << Dt << "sec" << endl;
	  
	  cout << "building H3 Hamiltonian" <<  endl;
	  gettimeofday (&(TotalStartingTime), 0);
	  ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian*Hamiltonian3 = new ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian(Space, ParitySector, Jxx, 0.0, 0.0, 0, Manager.GetBoolean("use-periodic"),
											       Architecture.GetArchitecture(), Memory);
	  RealSymmetricMatrix HRep3 (Hamiltonian3->GetHilbertSpaceDimension(), true);
	  RealMatrix Basis3(Hamiltonian3->GetHilbertSpaceDimension(), Hamiltonian3->GetHilbertSpaceDimension());
	  RealDiagonalMatrix TmpDiag3(Hamiltonian3->GetHilbertSpaceDimension());
	  Hamiltonian3->GetHamiltonian(HRep3);
#ifdef __LAPACK__
	  HRep3.LapackDiagonalize(TmpDiag3, Basis3);
#else
	  HRep3.Diagonalize(TmpDiag3, Basis3);
#endif
	  delete Hamiltonian3;
	  RealMatrix TransposedBasis3;
	  TransposedBasis3.Copy(Basis3);
	  TransposedBasis3.Transpose();
	  gettimeofday (&(TotalEndingTime), 0);
	  Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			 ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "done in " << Dt << "sec" << endl;
	  
	  for (int j = 0; j < NbrTauValues; ++j)
	    {
	      cout << "building unitary evolution operator at time tau=" <<  TauValues[j] << endl;
	      double TauValueHx = - TauValues[j] / 3.0;
	      double TauValueHy = - TauValues[j] / 3.0;
	      double TauValueHz = - TauValues[j] / 3.0;
	      gettimeofday (&(TotalStartingTime), 0);
	      ComplexMatrix UnitaryEvolution(Basis2);
	      for (int i = 0; i < UnitaryEvolution.GetNbrColumn(); ++i)
		{
		  UnitaryEvolution[i] *= Phase(TauValueHx * TmpDiag2[i]);
		}
 	      UnitaryEvolution.Multiply(TransposedBasis2);
	      for (int i = 0; i < UnitaryEvolution.GetNbrColumn(); ++i)
		{
		  UnitaryEvolution[i] *= Phase(TauValueHz * TmpDiag1[i]);
		}
	      gettimeofday (&(TotalEndingTime), 0);
	      
	      UnitaryEvolution.Multiply(Basis3);
	      for (int i = 0; i < UnitaryEvolution.GetNbrColumn(); ++i)
		{
		  UnitaryEvolution[i] *= Phase(TauValueHy * TmpDiag3[i]);
		}      
 	      UnitaryEvolution.Multiply(TransposedBasis3);
	      
	      Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			     ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "done in " << Dt << "sec" << endl;
	      
	      cout << "diagonalizing unitary evolution operator" <<  endl;
	      gettimeofday (&(TotalStartingTime), 0);
	      TmpDiagUnitaryEvolution = ComplexDiagonalMatrix (UnitaryEvolution.GetNbrColumn());
#ifdef __LAPACK__     
	      UnitaryEvolution.LapackDiagonalize(TmpDiagUnitaryEvolution);
#else
	      UnitaryEvolution.Diagonalize(TmpDiagUnitaryEvolution);
#endif
	      gettimeofday (&(TotalEndingTime), 0);
	      Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			     ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "done in " << Dt << "sec" << endl;
	      int Lim = TmpDiagUnitaryEvolution.GetNbrColumn();
	      double* TmpPhases = new double[Lim];
	      for (int i = 0; i < Lim; ++i)
		{
		  TmpPhases[i] = Arg(TmpDiagUnitaryEvolution[i]);
		}
	      SortArrayUpOrdering<Complex>(TmpPhases, TmpDiagUnitaryEvolution.GetDiagonalElements(), Lim);
	      ofstream File;
	      File.open(TauOutputFileNames[j], ios::out | ios::app);
	      File.precision(14);  
	      for (int i = 0; i < Lim; ++i)
		{
		  File << ParitySector << " " << i << " " << TmpDiagUnitaryEvolution[i] << " " << Norm(TmpDiagUnitaryEvolution[i]) << " " << TmpPhases[i] << endl;
		}	  
	      File.close();
	      delete[] TmpPhases;
	    }    
	  delete Space;	  
	}
   }
  return 0;
}
