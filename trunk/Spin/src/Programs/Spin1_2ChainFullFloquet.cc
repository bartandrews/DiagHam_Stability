#include "Hamiltonian/SpinChainHamiltonian.h"
#include "Hamiltonian/SpinChainFullHamiltonian.h"
#include "Hamiltonian/SpinChainRealFullHamiltonian.h"

#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainFixedParity.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/GenericRealMainTask.h"

#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
 
#include "Options/Options.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("Spin1_2ChainFullFloquet" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 8);
  (*SystemGroup) += new  SingleDoubleOption('\n', "tau", "coupling constant along the x axis", 1.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "multiple-taus", "compute the specrtum of the evolution operator at several times in a single shot");
  (*SystemGroup) += new  SingleDoubleOption ('x', "jx-value", "coupling constant in the x direction between neighboring sites", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('y', "jy-value", "coupling constant in the x direction between neighboring sites", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('z', "jz-value", "coupling constant in the x direction between neighboring sites", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hx-value", "amplitude of the Zeeman term along the x axis", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hy-value", "amplitude of the Zeeman term along the y axis", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hz-value", "amplitude of the Zeeman term along the z axis", 0.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "step-description", "an ASCII column formatted text file that describes each step of the time evolution instead of the default model");
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-hxvalue", "amplitude of the random Zeeman term along the x direction on each site", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-gaussianhxvalue", "amplitude of the random Zeeman term along the x direction on each site, using a gaussian disrtibution with zero mean value and a given standard deviation", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-hyvalue", "amplitude of the random Zeeman term along the y direction on each site", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-gaussianhyvalue", "amplitude of the random Zeeman term along the y direction on each site, using a gaussian disrtibution with zero mean value and a given standard deviation", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-hzvalue", "amplitude of the random Zeeman term along the z direction on each site", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-gaussianhzvalue", "amplitude of the random Zeeman term along the z direction on each site, using a gaussian disrtibution with zero mean value and a given standard deviation", 0.0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "run-id", "add an additional run id to the file name when using the --random-hzvalue option", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "fullh-values", "name of the file that contains the Zeeman term amplitudes for each site");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*OutputGroup) += new  SingleStringOption ('\n', "output-suffix", "apprend and extra suffix to the string describing the system in the output file name");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type Spin1_2ChainFullFloquet -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSpins = Manager.GetInteger("nbr-spin");
  double HxValue = Manager.GetDouble("hx-value");
  double HyValue = Manager.GetDouble("hy-value");
  double HzValue = Manager.GetDouble("hz-value");
  double JxValue = Manager.GetDouble("jx-value");
  double JyValue = Manager.GetDouble("jy-value");
  double JzValue = Manager.GetDouble("jz-value");


  char* BoundaryName = new char [16];
  if (Manager.GetBoolean("use-periodic") == false)
    sprintf (BoundaryName, "open");
  else
    sprintf (BoundaryName, "closed");
  
  char* FilePrefix = new char [256];
  sprintf (FilePrefix, "spin_1_2_floquet_%s_n_%d", BoundaryName, NbrSpins);

  char* OutputParameterFileName = new char [256];
  sprintf (OutputParameterFileName, "jx_%.6f_jy_%.6f_jz_%.6f", Manager.GetDouble("jx-value"), Manager.GetDouble("jy-value"), Manager.GetDouble("jz-value"));    
  if ((Manager.GetDouble("hx-value") != 0.0) || (Manager.GetDouble("hy-value") != 0.0) || (Manager.GetDouble("hz-value") != 0.0))
    {
      sprintf (OutputParameterFileName + strlen(OutputParameterFileName), "_hx_%.6f_hy_%.6f_hz_%.6f", Manager.GetDouble("hx-value"), 
	       Manager.GetDouble("hy-value"), Manager.GetDouble("hz-value"));
    }
  if ((Manager.GetDouble("random-hxvalue") != 0.0) || (Manager.GetDouble("random-hyvalue") != 0.0) || (Manager.GetDouble("random-hzvalue") != 0.0))
    {
      sprintf (OutputParameterFileName + strlen(OutputParameterFileName), "_randomhx_%.6f_randomhy_%.6f_randomhz_%.6f", Manager.GetDouble("random-hxvalue"), 
	       Manager.GetDouble("random-hyvalue"), Manager.GetDouble("random-hzvalue"));
    }
  else
    { 
      if ((Manager.GetDouble("random-gaussianhxvalue") != 0.0) || (Manager.GetDouble("random-gaussianhyvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	{
	  sprintf (OutputParameterFileName + strlen(OutputParameterFileName), "_grandomhx_%.6f_grandomhy_%.6f_grandomhz_%.6f", Manager.GetDouble("random-gaussianhxvalue"), 
		   Manager.GetDouble("random-gaussianhyvalue"), Manager.GetDouble("random-gaussianhzvalue"));
	}
    }
  


  char* TmpFileParameterString = new char [strlen(FilePrefix) + strlen(OutputParameterFileName) + 64];
  sprintf (TmpFileParameterString, "%s_%s", FilePrefix, OutputParameterFileName);
  double* JxValues = new double [NbrSpins];
  JxValues[0] = Manager.GetDouble("jx-value");
  for (int i = 1; i < NbrSpins; ++i)
    JxValues[i] = JxValues[0];
  double* JyValues = new double [NbrSpins];
  JyValues[0] = Manager.GetDouble("jy-value");
  for (int i = 1; i < NbrSpins; ++i)
    JyValues[i] = JyValues[0];
  double* JzValues = new double [NbrSpins];
  JzValues[0] = Manager.GetDouble("jz-value");
  for (int i = 1; i < NbrSpins; ++i)
    JzValues[i] = JzValues[0];
  double* HxValues = 0;
  double* HyValues = 0;
  double* HzValues = 0;
  if (Manager.GetString("fullh-values") != 0)
    {
      MultiColumnASCIIFile HFieldFile;
      if (HFieldFile.Parse(Manager.GetString("fullh-values")) == false)
	{
	  HFieldFile.DumpErrors(cout);
	  return -1;
	}
      if (HFieldFile.GetNbrLines() == NbrSpins)
	{
	  HxValues = HFieldFile.GetAsDoubleArray(0);
	  HyValues = HFieldFile.GetAsDoubleArray(1);
	  HzValues = HFieldFile.GetAsDoubleArray(2);
	}
      else
	{
	  if (HFieldFile.GetNbrLines() > NbrSpins)
	    {
	      cout << "warning, " << Manager.GetString("fullh-values") << " has more hz values than the number of sites" << endl;
	      HzValues = HFieldFile.GetAsDoubleArray(0);
	    }
	  else
	    {
	      cout << "error, " << Manager.GetString("fullh-values") << " has less hz values than the number of sites" << endl;
	      return 0;
	    }	  
	}
    }
  else
    {
      HxValues = new double [NbrSpins];
      HxValues[0] = Manager.GetDouble("hx-value");
      for (int i = 1; i < NbrSpins; ++i)
	HxValues[i] = HxValues[0];
      if ((Manager.GetDouble("random-hxvalue") != 0.0) || (Manager.GetDouble("random-gaussianhxvalue") != 0.0))
	{
	  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	  RandomNumber->UseTimeSeed();
	  for (int i = 0; i < NbrSpins; ++i)
	    {
	      double Tmp;
	      if (Manager.GetDouble("random-hxvalue") != 0.0)
		{
		  Tmp = Manager.GetDouble("random-hxvalue") * (2.0 * RandomNumber->GetRealRandomNumber() - 1.0);
		}
	      else
		{
		  Tmp = RandomNumber->GetGaussianRandomNumber(0.0, Manager.GetDouble("random-gaussianhxvalue"));
		}
	      HxValues[i] += Tmp;
	    }
	}
      HyValues = new double [NbrSpins];
      HyValues[0] = Manager.GetDouble("hy-value");
      for (int i = 1; i < NbrSpins; ++i)
	HyValues[i] = HyValues[0];
      if ((Manager.GetDouble("random-hyvalue") != 0.0) || (Manager.GetDouble("random-gaussianhyvalue") != 0.0))
	{
	  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	  RandomNumber->UseTimeSeed();
	  for (int i = 0; i < NbrSpins; ++i)
	    {
	      double Tmp;
	      if (Manager.GetDouble("random-hyvalue") != 0.0)
		{
		  Tmp = Manager.GetDouble("random-hyvalue") * (2.0 * RandomNumber->GetRealRandomNumber() - 1.0);
		}
	      else
		{
		  Tmp = RandomNumber->GetGaussianRandomNumber(0.0, Manager.GetDouble("random-gaussianhyvalue"));
		}
	      HyValues[i] += Tmp;
	    }
	}
      HzValues = new double [NbrSpins];
      HzValues[0] = Manager.GetDouble("hz-value");
      for (int i = 1; i < NbrSpins; ++i)
	HzValues[i] = HzValues[0];
      if ((Manager.GetDouble("random-hzvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	{
	  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	  RandomNumber->UseTimeSeed();
	  for (int i = 0; i < NbrSpins; ++i)
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
	    }
	}
      if ((Manager.GetDouble("random-hxvalue") != 0.0) || (Manager.GetDouble("random-hyvalue") != 0.0) || (Manager.GetDouble("random-hzvalue") != 0.0) ||
	  (Manager.GetDouble("random-gaussianhxvalue") != 0.0) || (Manager.GetDouble("random-gaussianhyvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	{
	  char* HOutputFileName = new char [strlen(TmpFileParameterString) + 64];
	  sprintf (HOutputFileName, "%s_runid_%ld.hvalues", TmpFileParameterString, Manager.GetInteger("run-id"));
	  ofstream File;
	  File.open(HOutputFileName, ios::binary | ios::out); 
	  File.precision(14); 
	  for (int i = 0; i < NbrSpins; ++i)
	    {
	      File << HxValues[i] << " " << HyValues[i] << " " << HzValues[i] << endl;
	    }
	  File.close();	      
	}
    }

  double* NegativeHxValues = new double [NbrSpins];
  double* NegativeHyValues = new double [NbrSpins];
  double* NullValues = new double [NbrSpins];
  for (int i = 0; i < NbrSpins; ++i)
    {
      NegativeHxValues[i] = -HxValues[i];
      NegativeHyValues[i] = -HyValues[i];
      NullValues[i] = 0.0;
    }


  int NbrTauValues  = 1;
  double* TauValues = 0;
  char** TauOutputFileNames = 0;
  char* FileParameterString = 0;
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
    }
  else
    {
      TauValues = new double[1];
      TauValues[0] = Manager.GetDouble("tau");
    }
  TauOutputFileNames = new char* [NbrTauValues];
  for (int i = 0; i < NbrTauValues; ++i)
    {
      TauOutputFileNames[i] = new char[strlen(TmpFileParameterString) + 64];
      if ((Manager.GetDouble("random-hxvalue") == 0.0) && (Manager.GetDouble("random-gaussianhxvalue") == 0.0)
	  && (Manager.GetDouble("random-hyvalue") == 0.0) && (Manager.GetDouble("random-gaussianhyvalue") == 0.0)
	  && (Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	{
	  sprintf (TauOutputFileNames[i], "%s_tau_%.6f.dat", TmpFileParameterString, TauValues[i]);
	}
      else
	{
	  sprintf (TauOutputFileNames[i], "%s_tau_%.6f_runid_%ld.dat", TmpFileParameterString, TauValues[i], Manager.GetInteger("run-id"));
	}
    }
  
  for (int i = 0; i < NbrTauValues; ++i)
    {
      ofstream File;
      File.open(TauOutputFileNames[i], ios::out);
      File.precision(14);  
      File << "# lambda norm(lambda) phase(lambda)" << endl;
      File.close();
    }


  AbstractSpinChain* Space = new Spin1_2ChainFull (NbrSpins);
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  double Dt;

  int NbrSteps = 0;
  double* TauRatios;
  double* JxFactors;
  double* JyFactors;
  double* JzFactors;
  double* HxFactors;
  double* HyFactors;
  double* HzFactors;
  if (Manager.GetString("step-description") != 0)
    {
      MultiColumnASCIIFile StepFile;
      if (StepFile.Parse(Manager.GetString("step-description")) == false)
	{
	  StepFile.DumpErrors(cout);
	  return -1;
	}
      NbrSteps = StepFile.GetNbrLines();
      TauRatios = StepFile.GetAsDoubleArray(0);
      JxFactors = StepFile.GetAsDoubleArray(1);
      JyFactors = StepFile.GetAsDoubleArray(2);
      JzFactors = StepFile.GetAsDoubleArray(3);
      HxFactors = StepFile.GetAsDoubleArray(4);
      HyFactors = StepFile.GetAsDoubleArray(5);
      HzFactors = StepFile.GetAsDoubleArray(6);
    }
  else
    {
      NbrSteps = 3;

      TauRatios = new double[NbrSteps];
      JxFactors = new double[NbrSteps];
      JyFactors = new double[NbrSteps];
      JzFactors = new double[NbrSteps];
      HxFactors = new double[NbrSteps];
      HyFactors = new double[NbrSteps];
      HzFactors = new double[NbrSteps];

      TauRatios[0] = 1.0 / 3.0;
      JxFactors[0] = 1.0;
      JyFactors[0] = 0.0;
      JzFactors[0] = 0.0;
      HxFactors[0] = 1.0;
      HyFactors[0] = 1.0;
      HzFactors[0] = 0.0;
      
      TauRatios[1] = 1.0 / 3.0;
      JxFactors[1] = 0.0;
      JyFactors[1] = 0.0;
      JzFactors[1] = 1.0;
      HxFactors[1] = 0.0;
      HyFactors[1] = 0.0;
      HzFactors[1] = 1.0;
      
      TauRatios[2] = 1.0 / 3.0;
      JxFactors[2] = 0.0;
      JyFactors[2] = 1.0;
      JzFactors[2] = 0.0;
      HxFactors[2] = -1.0;
      HyFactors[2] = 0.0;
      HzFactors[2] = 0.0;
    }

  RealDiagonalMatrix* HamiltonianEigenvalues = new RealDiagonalMatrix[NbrSteps]; 

  ComplexMatrix* StepMatrices = new ComplexMatrix[NbrSteps + 1];
  for (int j = 0; j < NbrSteps; ++j)
    {
      cout << "building H" << j << " Hamiltonian" <<  endl;
      gettimeofday (&(TotalStartingTime), 0);
      double* LocalJxValues = new double[NbrSpins];
      double* LocalJyValues = new double[NbrSpins];
      double* LocalJzValues = new double[NbrSpins];
      double* LocalHxValues = new double[NbrSpins];
      double* LocalHyValues = new double[NbrSpins];
      double* LocalHzValues = new double[NbrSpins];
      for (int i = 0; i < NbrSpins; ++i)
	{
	  LocalJxValues[i] = JxValues[i] * JxFactors[j];
	  LocalJyValues[i] = JyValues[i] * JyFactors[j];
	  LocalJzValues[i] = JzValues[i] * JzFactors[j];
	  LocalHxValues[i] = HxValues[i] * HxFactors[j];
	  LocalHyValues[i] = HyValues[i] * HyFactors[j];
	  LocalHzValues[i] = HzValues[i] * HzFactors[j];
	}
      if (HyFactors[j] == 0.0)
	{
	  if ((HxFactors[j] == 0.0) && (JxFactors[j] == 0.0) && (JyFactors[j] == 0.0)) 
	    {
	      SpinChainHamiltonian* Hamiltonian = new SpinChainHamiltonian(Space, NbrSpins, NullValues, LocalJzValues, 
									   LocalHzValues, Manager.GetBoolean("use-periodic"));
	      RealSymmetricMatrix HRep (Hamiltonian->GetHilbertSpaceDimension(), true);
	      HamiltonianEigenvalues[j] = RealDiagonalMatrix (Hamiltonian->GetHilbertSpaceDimension());
	      Hamiltonian->GetHamiltonian(HRep);
	      for (int i = 0; i < HRep.GetNbrRow(); ++i)
		{
		  HRep.GetMatrixElement(i, i, HamiltonianEigenvalues[j][i]);
		}
	      delete Hamiltonian;
	      if (j == 0)
		{
		  StepMatrices[0] = ComplexMatrix(Space->GetHilbertSpaceDimension(), Space->GetHilbertSpaceDimension());
		  StepMatrices[0].SetToIdentity();
		}
	    }
	  else
	    {
	      SpinChainRealFullHamiltonian* Hamiltonian = new SpinChainRealFullHamiltonian(Space, NbrSpins, LocalJxValues, LocalJyValues, LocalJzValues, 
											   LocalHxValues, LocalHzValues, Manager.GetBoolean("use-periodic"));
	      RealSymmetricMatrix HRep (Hamiltonian->GetHilbertSpaceDimension(), true);
	      RealMatrix TmpBasis (Hamiltonian->GetHilbertSpaceDimension(), Hamiltonian->GetHilbertSpaceDimension());
	      HamiltonianEigenvalues[j] = RealDiagonalMatrix(Hamiltonian->GetHilbertSpaceDimension());
	      Hamiltonian->GetHamiltonian(HRep);
	      cout << "done in " << Dt << "sec" << endl;      
	      cout << "diagonalizing H" << j << " Hamiltonian" <<  endl;
	      gettimeofday (&(TotalStartingTime), 0);
	      HRep.LapackDiagonalize(HamiltonianEigenvalues[j], TmpBasis);
	      delete Hamiltonian;
	      if (StepMatrices[j].GetNbrRow() == 0)
		{
		  StepMatrices[j] = ComplexMatrix(TmpBasis);
		}
	      else
		{
		  StepMatrices[j].Multiply(TmpBasis);
		}
	      TmpBasis.Transpose();
	      StepMatrices[j + 1] = ComplexMatrix(TmpBasis);
	    }
	}
      else
	{
	  SpinChainFullHamiltonian* Hamiltonian = new SpinChainFullHamiltonian(Space, NbrSpins, LocalJxValues, LocalJyValues, LocalJzValues,
									       LocalHxValues, LocalHyValues, LocalHzValues, Manager.GetBoolean("use-periodic"));
	  StepMatrices[j + 1] = ComplexMatrix (Hamiltonian->GetHilbertSpaceDimension(), Hamiltonian->GetHilbertSpaceDimension());
	  HermitianMatrix HRep (Hamiltonian->GetHilbertSpaceDimension(), true);
	  HamiltonianEigenvalues[j] = RealDiagonalMatrix(Hamiltonian->GetHilbertSpaceDimension());
	  Hamiltonian->GetHamiltonian(HRep);
	  gettimeofday (&(TotalEndingTime), 0);
	  Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			 ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "done in " << Dt << "sec" << endl;      
	  cout << "diagonalizing H" << j << " Hamiltonian" <<  endl;
	  gettimeofday (&(TotalStartingTime), 0);
	  HRep.LapackDiagonalize(HamiltonianEigenvalues[j], StepMatrices[j + 1]);
	  delete Hamiltonian;
	  if (StepMatrices[j].GetNbrRow() == 0)
	    {
	      StepMatrices[j].Copy(StepMatrices[j + 1]);
	    }
	  else
	    {
	      StepMatrices[j].Multiply(StepMatrices[j + 1]);
	    }
	  StepMatrices[j + 1].HermitianTranspose();
	}
      delete[] LocalJxValues;
      delete[] LocalJyValues;
      delete[] LocalJzValues;
      delete[] LocalHxValues;
      delete[] LocalHyValues;
      delete[] LocalHzValues;
      gettimeofday (&(TotalEndingTime), 0);
      Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
		     ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
      cout << "done in " << Dt << "sec" << endl;      
    }

  ComplexDiagonalMatrix TmpDiagUnitaryEvolution;
  for (int j = 0; j < NbrTauValues; ++j)
    {
      cout << "building unitary evolution operator at time tau=" <<  TauValues[j] << endl;
      double* TauSteps = new double [NbrSteps];
      for (int i = 0; i < NbrSteps; ++i)
	{
	  TauSteps[i] = -TauValues[j] * TauRatios[i];
	} 
      gettimeofday (&(TotalStartingTime), 0);
      ComplexMatrix UnitaryEvolution;
      UnitaryEvolution.Copy(StepMatrices[0]);
      for (int i = 0; i < NbrSteps; ++i)
	{
 	  for (int k = 0; k < UnitaryEvolution.GetNbrColumn(); ++k)
 	    {
 	      UnitaryEvolution[k] *= Phase(TauSteps[i] * HamiltonianEigenvalues[i][k]);
 	    }
	  if (StepMatrices[i + 1].GetNbrRow() != 0)
	    {
	      UnitaryEvolution.Multiply(StepMatrices[i + 1]);
	    }
	}
      gettimeofday (&(TotalEndingTime), 0);
      Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
		     ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
      cout << "done in " << Dt << "sec" << endl;
      
      cout << "diagonalizing unitary evolution operator" <<  endl;
      gettimeofday (&(TotalStartingTime), 0);
      TmpDiagUnitaryEvolution = ComplexDiagonalMatrix (UnitaryEvolution.GetNbrColumn());
      ComplexMatrix TmpUnitaryEvolutionEigenstates;
      if (((Manager["all-eigenstates"] != 0) && (Manager.GetBoolean("all-eigenstates") == true)) ||
	  ((Manager["eigenstate"] != 0) && (Manager.GetBoolean("eigenstate") == true)))
	{
	  TmpUnitaryEvolutionEigenstates = ComplexMatrix(UnitaryEvolution.GetNbrColumn(), UnitaryEvolution.GetNbrColumn());
	  UnitaryEvolution.LapackDiagonalize(TmpDiagUnitaryEvolution, TmpUnitaryEvolutionEigenstates);
	}
      else
	{
	  UnitaryEvolution.LapackDiagonalize(TmpDiagUnitaryEvolution);
	}
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
      if (TmpUnitaryEvolutionEigenstates.GetNbrRow() == 0)
	{
	  SortArrayUpOrdering<Complex>(TmpPhases, TmpDiagUnitaryEvolution.GetDiagonalElements(), Lim);
	}
      else
	{
	  int* TmpArray = new int[Lim];
	  for  (int i = 0; i < Lim; ++i)
	    {
	      TmpArray[i] = i;
	    }
	  SortArrayUpOrdering<int>(TmpPhases, TmpArray, Lim);
	  ComplexVector* TmpVectors = new ComplexVector[Lim];
	  for  (int i = 0; i < Lim; ++i)
	    {
	      TmpVectors[i] = TmpUnitaryEvolutionEigenstates[TmpArray[i]];		  
	    }
	  ComplexDiagonalMatrix TmpDiagUnitaryEvolution2(Lim);
	  for  (int i = 0; i < Lim; ++i)
	    {
	      TmpDiagUnitaryEvolution2[i] = TmpDiagUnitaryEvolution[TmpArray[i]];
	    }
	  TmpDiagUnitaryEvolution = TmpDiagUnitaryEvolution2;
	  ComplexMatrix TmpUnitaryEvolutionEigenstates2(TmpVectors, Lim);
	  if ((Manager["all-eigenstates"] != 0) && (Manager.GetBoolean("all-eigenstates") == true))
	    {
	      char* TmpVectorName = ReplaceExtensionToFileName (TauOutputFileNames[j], "dat", "eigenvec.mat");
	      TmpUnitaryEvolutionEigenstates2.WriteMatrix(TmpVectorName );
	      delete[] TmpVectorName;
	    }
	  else
	    {
	      int FirstEigenstateIndex = 0;
	      if (Manager["first-eigenstate"] != 0)
		FirstEigenstateIndex = Manager.GetInteger("first-eigenstate");
	      int NbrEigenvalues = Manager.GetInteger("nbr-eigen");
	      int LastEigenstateIndex = FirstEigenstateIndex + NbrEigenvalues;
	      for (; FirstEigenstateIndex < LastEigenstateIndex; ++FirstEigenstateIndex)
		{
		  char* TmpExtension = new char [16];
		  sprintf (TmpExtension, "%d.vec", FirstEigenstateIndex);
		  char* TmpVectorName = ReplaceExtensionToFileName (TauOutputFileNames[j], "dat", TmpExtension);
		  TmpUnitaryEvolutionEigenstates2[FirstEigenstateIndex].WriteVector(TmpVectorName);
		  delete[] TmpExtension;
		  delete[] TmpVectorName;
		}
	    }
	}

      ofstream File;
      File.open(TauOutputFileNames[j], ios::out | ios::app);
      File.precision(14);  
      for (int i = 0; i < Lim; ++i)
	{
	  File << i << " " << TmpDiagUnitaryEvolution[i] << " " << Norm(TmpDiagUnitaryEvolution[i]) << " " << TmpPhases[i] << endl;
	}	  
      File.close();
      
      delete[] TmpPhases;
    }    
  delete Space;	  

  return 0;
}
