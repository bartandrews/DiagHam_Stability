#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/BosonOnDisk.h"
#include "Hamiltonian/ParticleOnDiskDeltaHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ios;

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEBosonsDiskDelta" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += MiscGroup;

 (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
 (*SystemGroup) += new SingleIntegerOption  ('l', "maximum-momentum", "maximum single particle momentum to study", 10, true, 1);
 (*SystemGroup) += new SingleIntegerOption  ('\n', "minimum-momentum", "minimum single particle momentum to study", 1, true, 1);

  (*LanczosGroup) += new SingleIntegerOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup) += new SingleIntegerOption ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup) += new SingleIntegerOption ('\n', "full-diag", 
					      "maximum Hilbert space dimension for which full diagonalization is applied", 
					      500, true, 100);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsDiskDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int MaxNbrIterLanczos = ((SingleIntegerOption*) Manager["iter-max"])->GetInteger();
  int NbrEigenvalue = ((SingleIntegerOption*) Manager["nbr-eigen"])->GetInteger();
  int NbrBosons = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int MMin = ((SingleIntegerOption*) Manager["minimum-momentum"])->GetInteger();
  int MMax = ((SingleIntegerOption*) Manager["maximum-momentum"])->GetInteger();
  if (MMax < MMin)
    MMax = MMin;
  int FullDiagonalizationLimit = ((SingleIntegerOption*) Manager["full-diag"])->GetInteger();

  char* OutputName = new char [1024];
  sprintf (OutputName, "bosons_disk_delta_n_%d_l_%d.dat", NbrBosons, MMax);
  ofstream File;
  File.open(OutputName, ios::binary | ios::out);

  for (int  L = MMin; L <= MMax; ++L)
    {
      BosonOnDisk Space (NbrBosons, L);
      cout << "Nbr bosons = " << NbrBosons << "    L = " << L << "    Dimension = " << Space.GetHilbertSpaceDimension() << endl;
      Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());
      ParticleOnDiskDeltaHamiltonian* Hamiltonian = new ParticleOnDiskDeltaHamiltonian(&Space, L);
      if (Hamiltonian->GetHilbertSpaceDimension() < FullDiagonalizationLimit)
	{
	  RealSymmetricMatrix HRep (Hamiltonian->GetHilbertSpaceDimension());
	  Hamiltonian->GetHamiltonian(HRep);
	  if (Hamiltonian->GetHilbertSpaceDimension() > 1)
	    {
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian->GetHilbertSpaceDimension());
	      HRep.Householder(TmpTriDiag, 1e-7);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      for (int j = 0; j < Hamiltonian->GetHilbertSpaceDimension(); j++)
		{
		  File << L << " " << TmpTriDiag.DiagonalElement(j) << endl;
		}
	    }
	  else
	    {
	      double TmpVal = HRep(0, 0);
	      File << L << " " << TmpVal << endl;
	    }
	}
      else
	{
	  AbstractLanczosAlgorithm* Lanczos;
	  if (NbrEigenvalue == 1)
	    {
	       Lanczos = new BasicLanczosAlgorithm(Architecture.GetArchitecture(), NbrEigenvalue, MaxNbrIterLanczos);
	    }
	  else
	    {
	      Lanczos = new FullReorthogonalizedLanczosAlgorithm (Architecture.GetArchitecture(), NbrEigenvalue, MaxNbrIterLanczos);
	    }
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = 0;
	  Lanczos->SetHamiltonian(Hamiltonian);
	  Lanczos->InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  double Dt;
	  gettimeofday (&(TotalStartingTime), 0);
	  Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	  CurrentNbrIterLanczos = NbrEigenvalue + 3;
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while ((Lanczos->TestConvergence() == false) &&  (CurrentNbrIterLanczos < MaxNbrIterLanczos))
	    {
	      ++CurrentNbrIterLanczos;
	      Lanczos->RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
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
	  for (int i = 0; i <= NbrEigenvalue; ++i)
	    {
	      cout << TmpMatrix.DiagonalElement(i) << " ";
	      File << (int) (L / 2) << " " << (TmpMatrix.DiagonalElement(i)) << endl;
	    }
	  cout << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	  delete Lanczos;
	}
      cout << endl;
      cout << "//////////////////////////////////////////////////////" << endl;
    }
  File.close();
  return 0;
}

