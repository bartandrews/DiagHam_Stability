#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "Hamiltonian/ParticleOnSphereCoulombHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "MathTools/ClebschGordanCoefficients.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <stdlib.h>
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
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 40);
  SingleIntegerOption LzMaxOption ('l', "lzmax", "twice the maximum momentum for a single particle", 14);
  SingleIntegerOption NbrFermionOption ('p', "nbr-particles", "number of particles", 8);
  SingleDoubleOption FrequencyShiftOption ('f', "frequency-shift", "frequency shift from FQHE", 0.0);
  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &NbrFermionOption;
  OptionList += &LzMaxOption;
  OptionList += &FrequencyShiftOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type ExplicitMatrixExample -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }
  bool SMPFlag = SMPOption.GetBoolean();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  int NbrFermions = NbrFermionOption.GetInteger();
  int LzMax = LzMaxOption.GetInteger();

//  int L = 0;
  double GroundStateEnergy = 0.0;
  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "fermions_coulomb_n_%d_2s_%d_lz.dat", NbrFermions, LzMax);

  int Shift = 0;
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  int Max = ((LzMax - NbrFermions + 1) * NbrFermions);
  int TotalSize = 0;
  double** Eigenvalues = new double* [Max + 1];
  int* Dimensions = new int [Max + 1];
  int  L = 0;
  if ((abs(Max) & 1) != 0)
     L = 1;
  for (; L <= Max; L += 2)
    {
      cout << "----------------------------------------------------------------" << endl;
      cout << " LzTotal = " << L << endl;
      FermionOnSphere Space (NbrFermions, L, LzMax);
      cout << " Hilbert space dimension = " << Space.GetHilbertSpaceDimension() << endl;
      TotalSize += Space.GetHilbertSpaceDimension();
/*     for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
	{
	  cout << i << " = ";
	  Space.PrintState(cout, i) << endl;
	}*/
      ParticleOnSphereCoulombHamiltonian Hamiltonian(&Space, NbrFermions, LzMax, 0);
      if (Hamiltonian.GetHilbertSpaceDimension() < 200)
	{
	  Dimensions[L >> 1] = Hamiltonian.GetHilbertSpaceDimension();
	  Eigenvalues[L >> 1] = new double [Hamiltonian.GetHilbertSpaceDimension()];
	  RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
	  Hamiltonian.GetHamiltonian(HRep);
//	  cout << HRep << endl;
	  if (Hamiltonian.GetHilbertSpaceDimension() > 1)
	    {
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian.GetHilbertSpaceDimension());
	      HRep.Householder(TmpTriDiag, 1e-7);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      if (L == 0)
		GroundStateEnergy = TmpTriDiag.DiagonalElement(0);
	      //	  cout << "eigenvalues : " << endl;
	      for (int j = 0; j < Hamiltonian.GetHilbertSpaceDimension() ; j++)
		{
		  Eigenvalues[L >> 1][j] = TmpTriDiag.DiagonalElement(j);
		  cout << TmpTriDiag.DiagonalElement(j) << " ";
		  File << (L / 2) << " " << TmpTriDiag.DiagonalElement(j) << endl;
// (TmpTriDiag.DiagonalElement(j) - GroundStateEnergy) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
	      Eigenvalues[L >> 1][0] = HRep(0, 0);
	      //	      GroundStateEnergy = HRep(0, 0);
	      File << (L / 2) << " " << HRep(0, 0) << endl;// - GroundStateEnergy) / (4 * M_PI)) << endl;
	    }
	}
      else
	{
	  AbstractArchitecture* Architecture = 0;
	  if (SMPFlag == false)
	    Architecture = new MonoProcessorArchitecture;
	  else
	    Architecture = new SMPArchitecture(2);
	  FullReorthogonalizedLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = NbrEigenvalue + 3;
	  Lanczos.SetHamiltonian(&Hamiltonian);
	  Lanczos.InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  double Dt;
	  gettimeofday (&(TotalStartingTime), 0);
	  Lanczos.RunLanczosAlgorithm(NbrEigenvalue + 2);
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while ((Precision > 1e-14) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
	    {
	      Lanczos.RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos.GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue);//Lanczos.GetGroundStateEnergy();
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << " "<< endl;
	    }
	  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
	    {
	      cout << "too much Lanczos iterations" << endl;
	      File << "too much Lanczos iterations" << endl;
	      File.close();
	      exit(0);
	    }
	  GroundStateEnergy = Lowest;
	  cout << endl;
	  cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	       << CurrentNbrIterLanczos << endl;
	  for (int i = 0; i <= NbrEigenvalue; ++i)
	    {
	      cout << TmpMatrix.DiagonalElement(i) << " ";
	      File << (L / 2) << " " << TmpMatrix.DiagonalElement(i) << endl;
	    }
	  cout << endl;
//	  File << (L / 2) << " " << Lowest << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	}
      cout << "----------------------------------------------------------------" << endl;
      cout << " Total Hilbert space dimension = " << TotalSize << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrFermions) << endl;
    }
  File.close();

  return 0;
}
