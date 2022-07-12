#include "Hamiltonian/SpinChainXYZHamiltonian.h"
#include "Hamiltonian/SpinChainPureHFieldHamiltonian.h"

#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainFixedParity.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/GenericRealMainTask.h"

#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("Spin1_2ChainFloquet" , "0.01");
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
  (*SystemGroup) += new  SingleDoubleOption('\n', "jz-value", "coupling strength along the z axis between two neighboring spins", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('\n', "hx-value", "Zeeman strength along the x axis", 0.0);
  (*SystemGroup) += new  SingleDoubleOption('\n', "hy-value", "Zeeman strength along the y axis", 0.0);
  (*SystemGroup) += new  SingleDoubleOption('\n', "hz-value", "Zeeman strength along the z axis", 0.0);
  (*SystemGroup) += new  BooleanOption  ('\n', "use-parity", "take the parity into account when computing the spectrum");
  (*SystemGroup) += new  BooleanOption  ('\n', "read-spectrum", "read the unitary matrix spectrum from the disk instead of recomputing it");
  (*SystemGroup) += new  SingleDoubleOption('\n', "hy-tau", "fraction of the time tau where the hy hamiltonian is turned on (if negative, it is set to 1/3)", -1.0);
  (*SystemGroup) += new  SingleDoubleOption('\n', "hz-tau", "fraction of the time tau where the hz hamiltonian is turned on (if negative, it is set to 1/2 if hy-value=0, 1/3 otherwise)", -1.0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*OutputGroup) += new  SingleStringOption ('\n', "output-suffix", "apprend and extra suffix to the string describing the system in the output file name");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type Spin1_2ChainFloquet -h" << endl;
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
  double JzValue = Manager.GetDouble("jz-value");
  double TauValueHx = -0.5 * Manager.GetDouble("tau");
  double TauValueHy = 0.0;
  double TauValueHz = -0.5 * Manager.GetDouble("tau");
  if (HyValue != 0.0)
    {
      TauValueHx = - Manager.GetDouble("tau") / 3.0;
      TauValueHy = - Manager.GetDouble("tau") / 3.0;
      TauValueHz = - Manager.GetDouble("tau") / 3.0;
      if ((Manager.GetDouble("hz-tau") >= 0.0) && (Manager.GetDouble("hy-tau") >= 0.0))
	{
	  TauValueHx = - Manager.GetDouble("tau") * (1.0 - Manager.GetDouble("hz-tau") - Manager.GetDouble("hy-tau"));
	  TauValueHy = - Manager.GetDouble("tau") * Manager.GetDouble("hy-tau");
	  TauValueHz = - Manager.GetDouble("tau") * Manager.GetDouble("hz-tau");
	  
	}
    }
  else
    {
      if (Manager.GetDouble("hz-tau") >= 0.0)
	{
	  TauValueHx = - Manager.GetDouble("tau") * (1.0 - Manager.GetDouble("hz-tau"));
	  TauValueHz = - Manager.GetDouble("tau") * Manager.GetDouble("hz-tau");
	  
	}
    }
  char* OutputFileName = new char[1024];
  if (Manager.GetBoolean("use-parity") == false)
    {
      if ((Manager.GetDouble("hz-tau") >= 0.0) || (Manager.GetDouble("hy-tau") >= 0.0))
	{
	  sprintf (OutputFileName, "spin_1_2_floquet_n_%d_jz_%.6f_hx_%.6f_hy_%.6f_hz_%.6f_tau_%.6f_tauhy_%.6f_tauhz_%.6f.dat", 
		   NbrSpins, JzValue, HxValue, HyValue, HzValue, Manager.GetDouble("tau"), Manager.GetDouble("hy-tau"), Manager.GetDouble("hz-tau"));     
	}
      else
	{
	  sprintf (OutputFileName, "spin_1_2_floquet_n_%d_jz_%.6f_hx_%.6f_hy_%.6f_hz_%.6f_tau_%.6f.dat", NbrSpins, JzValue, HxValue, HyValue, HzValue, Manager.GetDouble("tau"));     
	}
    }
  else
    {
      if ((Manager.GetDouble("hz-tau") >= 0.0) || (Manager.GetDouble("hy-tau") >= 0.0))
	{
	  sprintf (OutputFileName, "spin_1_2_floquet_parity_n_%d_jz_%.6f_hx_%.6f_hy_%.6f_hz_%.6f_tau_%.6f_tauhy_%.6f_tauhz_%.6f.dat", 
		   NbrSpins, JzValue, HxValue, HyValue, HzValue, Manager.GetDouble("tau"), Manager.GetDouble("hy-tau"), Manager.GetDouble("hz-tau"));
    	}
      else
	{
	  sprintf (OutputFileName, "spin_1_2_floquet_parity_n_%d_jz_%.6f_hx_%.6f_hy_%.6f_hz_%.6f_tau_%.6f.dat", NbrSpins, JzValue, HxValue, HyValue, HzValue, Manager.GetDouble("tau"));     
	}
    }
  ofstream File;
  if (Manager.GetBoolean("read-spectrum") == false)
    {
      File.open(OutputFileName, ios::out);
      File.precision(14);  
    }


  Spin1_2Chain* Chain;
  int ParitySector = 0;
  int MaxParitySector = 1;
  if (Manager.GetBoolean("use-parity") == false)
    {
      MaxParitySector = 0;
    }
  double Min = 0.0;
  double Max = 0.0;
  while (ParitySector <= MaxParitySector)
    {
      ComplexDiagonalMatrix TmpDiagUnitaryEvolution;
      if (Manager.GetBoolean("read-spectrum") == false)
	{
	  if (Manager.GetBoolean("use-parity") == false)
	    {
	      Chain = new Spin1_2ChainFull (NbrSpins);
	    }
	  else
	    {
	      cout << "computing parity sector " << ParitySector << endl;
	      Chain = new Spin1_2ChainFixedParity (NbrSpins, ParitySector);
	    }
	  
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  double Dt;
	  cout << "building H1 Hamiltonian" <<  endl;
	  gettimeofday (&(TotalStartingTime), 0);
	  int StartTimeSecond = TotalStartingTime.tv_sec;
	  SpinChainXYZHamiltonian* Hamiltonian1 = new SpinChainXYZHamiltonian (Chain, NbrSpins, 0.0, 0.0, JzValue, 2.0 * HzValue);
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
	  SpinChainPureHFieldHamiltonian* Hamiltonian2 = new SpinChainPureHFieldHamiltonian (Chain, NbrSpins, HxValue, 0.0, 0.0);
	  RealMatrix Basis2(Hamiltonian2->GetHilbertSpaceDimension(), Hamiltonian2->GetHilbertSpaceDimension());
	  RealSymmetricMatrix HRep2 (Hamiltonian2->GetHilbertSpaceDimension(), true);
	  RealDiagonalMatrix TmpDiag2 (Hamiltonian2->GetHilbertSpaceDimension());
	  Hamiltonian2->GetHamiltonian(HRep2);
	  HRep2.LapackDiagonalize(TmpDiag2, Basis2);
	  delete Hamiltonian2;
	  gettimeofday (&(TotalEndingTime), 0);
	  Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			 ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "done in " << Dt << "sec" << endl;
	  
	  ComplexMatrix Basis3;
	  RealDiagonalMatrix TmpDiag3;
	  if (HyValue != 0.0)
	    {
	      cout << "building H3 Hamiltonian" <<  endl;
	      gettimeofday (&(TotalStartingTime), 0);
	      SpinChainPureHFieldHamiltonian* Hamiltonian3 = new SpinChainPureHFieldHamiltonian (Chain, NbrSpins, 0.0, HyValue, 0.0);
	      Basis3 = ComplexMatrix(Hamiltonian3->GetHilbertSpaceDimension(), Hamiltonian3->GetHilbertSpaceDimension());
	      HermitianMatrix HRep3 (Hamiltonian3->GetHilbertSpaceDimension(), true);
	      TmpDiag3 = RealDiagonalMatrix(Hamiltonian3->GetHilbertSpaceDimension());
	      Hamiltonian3->GetHamiltonian(HRep3);
	      HRep3.LapackDiagonalize(TmpDiag3, Basis3);
	      delete Hamiltonian3;
	      gettimeofday (&(TotalEndingTime), 0);
	      Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			     ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "done in " << Dt << "sec" << endl;
	    }
	  
	  cout << "building unitary evolution operator" <<  endl;
	  gettimeofday (&(TotalStartingTime), 0);
	  ComplexMatrix UnitaryEvolution(Basis2);
	  for (int i = 0; i < UnitaryEvolution.GetNbrColumn(); ++i)
	    {
	      UnitaryEvolution[i] *= Phase(TauValueHx * TmpDiag2[i]);
	    }
	  Basis2.Transpose();
	  UnitaryEvolution.Multiply(Basis2);
	  for (int i = 0; i < UnitaryEvolution.GetNbrColumn(); ++i)
	    {
	      UnitaryEvolution[i] *= Phase(TauValueHz * TmpDiag1[i]);
	    }
	  gettimeofday (&(TotalEndingTime), 0);
	  
	  if (HyValue != 0.0)
	    {
	      UnitaryEvolution.Multiply(Basis3);
	      for (int i = 0; i < UnitaryEvolution.GetNbrColumn(); ++i)
		{
		  UnitaryEvolution[i] *= Phase(TauValueHy * TmpDiag3[i]);
		}      
	      Basis3.HermitianTranspose();
	      UnitaryEvolution.Multiply(Basis3);
	    }
	  Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			 ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "done in " << Dt << "sec" << endl;
      
	  cout << "diagonalizing unitary evolution operator" <<  endl;
	  gettimeofday (&(TotalStartingTime), 0);
	  TmpDiagUnitaryEvolution = ComplexDiagonalMatrix (UnitaryEvolution.GetNbrColumn());
	  UnitaryEvolution.LapackDiagonalize(TmpDiagUnitaryEvolution);
	  gettimeofday (&(TotalEndingTime), 0);
	  Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			 ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "done in " << Dt << "sec" << endl;
	  delete Chain;	  
	}
      else
	{
	  MultiColumnASCIIFile SpectrumFile;
	  if (SpectrumFile.Parse(OutputFileName) == false)
	    {
	      SpectrumFile.DumpErrors(cout);
	      return -1;
	    }
	  if (Manager.GetBoolean("use-parity") == false)
	    {
	      int TmpSize = SpectrumFile.GetNbrLines();
	      Complex* TmpSpectrum = SpectrumFile.GetAsComplexArray(1);
	      TmpDiagUnitaryEvolution = ComplexDiagonalMatrix (TmpSpectrum, TmpSize);
	    }
	  else
	    {
	      int TmpSize = SpectrumFile.GetNbrLines();
	      int* TmpParity = SpectrumFile.GetAsIntegerArray(0);
	      Complex* TmpSpectrum = SpectrumFile.GetAsComplexArray(2);
	      int TmpActualSize = 0;
	      for (int i = 0; i < TmpSize; ++i)
		{
		  if (TmpParity[i] == ParitySector)
		    ++TmpActualSize;
		}
	      Complex* TmpSpectrum2 = new Complex [TmpActualSize];
	      TmpActualSize = 0;
	      for (int i = 0; i < TmpSize; ++i)
		{
		  if (TmpParity[i] == ParitySector)
		    {
		      TmpSpectrum2[TmpActualSize] = TmpSpectrum[i];
		      ++TmpActualSize;
		    }
		}
	      TmpDiagUnitaryEvolution = ComplexDiagonalMatrix (TmpSpectrum2, TmpActualSize);
	      delete[] TmpSpectrum;
	      delete[] TmpParity;
	    }
	}
            
      int Lim = TmpDiagUnitaryEvolution.GetNbrColumn();
      double* TmpPhases = new double[Lim];
      for (int i = 0; i < Lim; ++i)
	{
	  TmpPhases[i] = Arg(TmpDiagUnitaryEvolution[i]);
	}
      if (Manager.GetBoolean("read-spectrum") == false)
	{
	  SortArrayUpOrdering<Complex>(TmpPhases, TmpDiagUnitaryEvolution.GetDiagonalElements(), TmpDiagUnitaryEvolution.GetNbrColumn());
	  if (Manager.GetBoolean("use-parity") == false)
	    {
	      for (int i = 0; i < Lim; ++i)
		{
		  File << i << " " << TmpDiagUnitaryEvolution[i] << " " << Norm(TmpDiagUnitaryEvolution[i]) << " " << TmpPhases[i] << endl;
		}
	    }
	  else
	    {
	      for (int i = 0; i < Lim; ++i)
		{
		  File << ParitySector << " " << i << " " << TmpDiagUnitaryEvolution[i] << " " << Norm(TmpDiagUnitaryEvolution[i]) << " " << TmpPhases[i] << endl;
		}
	    }
	}
      --Lim;
      double TmpInfDiff;
      double TmpSupDiff;
      for (int i = 1; i < Lim; ++i)
	{
	  TmpInfDiff = TmpPhases[i] - TmpPhases[i - 1];
	  TmpSupDiff = TmpPhases[i + 1] - TmpPhases[i];
	  if (TmpInfDiff > TmpSupDiff)
	    {
	      Min += TmpSupDiff;
	      Max += TmpInfDiff;
	    }
	  else
	    {
	      Max += TmpSupDiff;
	      Min += TmpInfDiff;
	    }
	}
      delete[] TmpPhases;
      ++ParitySector;
    }
  cout << "r=" << (Min / Max) << endl;
  if (Manager.GetBoolean("read-spectrum") == false)
    {
      File.close();
    }
  return 0;
}
