#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Hamiltonian/ParticleOnTorusCoulombHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusDeltaHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "QuantumNumber/AbstractQuantumNumber.h"

#include "GeneralTools/ListIterator.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("QHEBosonsTorus" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption ('l', "max-momentum", "maximum momentum for a single particle", 15);
  (*SystemGroup) += new SingleIntegerOption ('m', "momentum", "constraint on the total momentum modulo the maximum momentum (negative if none)", -1);
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the two torus lengths", 0.57735026919);
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");

  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag", 
						"maximum Hilbert space dimension for which full diagonalization is applied", 
						500, true, 100);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);  
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "limit-time", "use limit in time instead of a number of lanczos iteration (0 if none, time in seconds)", 0);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false); 
  (*LanczosGroup) += new SingleStringOption ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
  (*LanczosGroup) += new BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsNBodyHardCore -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  bool GroundFlag = ((BooleanOption*) Manager["ground"])->GetBoolean();
  int NbrBosons = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int MaxMomentum = ((SingleIntegerOption*) Manager["max-momentum"])->GetInteger();
  int Momentum = ((SingleIntegerOption*) Manager["momentum"])->GetInteger();
  double XRatio = ((SingleDoubleOption*) Manager["ratio"])->GetDouble();

  int MaxNbrIterLanczos = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int NbrEigenvalue = ((SingleIntegerOption*) Manager["nbr-eigen"])->GetInteger();
  int FullDiagonalizationLimit = ((SingleIntegerOption*) Manager["full-diag"])->GetInteger();

  int L = 0;
  double GroundStateEnergy = 0.0;

  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "bosons_torus_delta_n_%d_2s_%d_ratio_%f.dat", NbrBosons, MaxMomentum, XRatio);
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);

  int Max = (MaxMomentum - 1);
  if (Momentum < 0)
    Momentum = 0;
  else
    Max = Momentum;
  
  for (; Momentum <= Max; ++Momentum)
    {     
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
      BosonOnTorus Space (NbrBosons, MaxMomentum, Momentum);
      cout << " Hilbert space dimension = " << Space.GetHilbertSpaceDimension() << endl;

      Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());
      AbstractHamiltonian* Hamiltonian = new ParticleOnTorusDeltaHamiltonian (&Space, NbrBosons, MaxMomentum, XRatio, Architecture.GetArchitecture());

      if (Hamiltonian->GetHilbertSpaceDimension() < FullDiagonalizationLimit)
	{
	  RealSymmetricMatrix HRep (Hamiltonian->GetHilbertSpaceDimension());
	  Hamiltonian->GetHamiltonian(HRep);
	  cout << HRep << endl;
	  if (Hamiltonian->GetHilbertSpaceDimension() > 1)
	    {
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian->GetHilbertSpaceDimension());
	      HRep.Householder(TmpTriDiag, 1e-7);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      if (L == 0)
		GroundStateEnergy = TmpTriDiag.DiagonalElement(0);
	      for (int j = 0; j < Hamiltonian->GetHilbertSpaceDimension() ; j++)
		{
		  File << Momentum << " " << TmpTriDiag.DiagonalElement(j) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
	      File << Momentum << " " << HRep(0, 0) << endl;
	    }
	}
      else
	{
	  AbstractLanczosAlgorithm* Lanczos = 0;
	  if (NbrEigenvalue > 1)
	    Lanczos = new FullReorthogonalizedLanczosAlgorithm (Architecture.GetArchitecture(), NbrEigenvalue, MaxNbrIterLanczos);
	  else
	    Lanczos = new BasicLanczosAlgorithm (Architecture.GetArchitecture(), NbrEigenvalue, MaxNbrIterLanczos);
	    
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = NbrEigenvalue + 3;
	  Lanczos->SetHamiltonian(Hamiltonian);
	  Lanczos->InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  double Dt;
	  gettimeofday (&(TotalStartingTime), 0);
	  Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while ((Precision > 1e-14) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
	    {
	      Lanczos->RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue);
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
	      File << Momentum << " " << TmpMatrix.DiagonalElement(i) << endl;
	    }
	  cout << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	  delete Lanczos;
	}
      cout << "----------------------------------------------------------------" << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrBosons) << endl;
      delete Hamiltonian;
    }
  File.close();

  return 0;
}
