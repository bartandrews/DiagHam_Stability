#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/ExplicitHamiltonian.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "Hamiltonian/ParticleOnSphereDeltaHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"

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
  BooleanOption LanczosOption ('l', "lanczos", "enable lanczos diagonalization algorithm", false);
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption VerboseOption ('v', "verbose", "verbose mode");
  BooleanOption EigenstateOption ('e', "eigenstate", "evaluate eigenstates");
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 40);
  List<AbstractOption*> OptionList;
  OptionList += &LanczosOption;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &VerboseOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &EigenstateOption;
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
  bool LanczosFlag = LanczosOption.GetBoolean();
  bool SMPFlag = SMPOption.GetBoolean();
  bool VerboseFlag = VerboseOption.GetBoolean();
  bool EigenstateFlag = EigenstateOption.GetBoolean();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();


  // initialize matrix associated to the hamiltonian

  // delete this part (just an example)
  BosonOnSphere Space (6, 0, 6);
  ParticleOnSphereDeltaHamiltonian Hamiltonian2(&Space, 6, 6, 0);
  RealSymmetricMatrix HRep (Hamiltonian2.GetHilbertSpaceDimension());
  Hamiltonian2.GetHamiltonian(HRep);
 // stop delete

  // enter here hilbert space dimension
  int Dimension = Hamiltonian2.GetHilbertSpaceDimension();

  RealSymmetricMatrix HamiltonianRepresentation (Dimension);


  cout << "dim = " << (Hamiltonian2.GetHilbertSpaceDimension()) << endl;
  for (int i = 0; i < Dimension; ++i)
    for (int j = 0; j <= i; ++j)
      {
	// insert here your algorithm to fill the hamiltonian representation
	HamiltonianRepresentation(j, i) = HRep(j, i);
      }

  // store matrix in an hamiltonian class and create corresponding hilbert space
  UndescribedHilbertSpace HilbertSpace(Dimension);
  ExplicitHamiltonian Hamiltonian(&HilbertSpace, &HamiltonianRepresentation);
  RealVector* Eigenstates = 0;
  double* Eigenvalues = 0;

  // find the eigenvalues (and eigenvectors if needed)
  if (LanczosFlag == false)
    {
      if (Hamiltonian.GetHilbertSpaceDimension() > 1)
	{
	  // diagonalize the hamiltonian
	  RealTriDiagonalSymmetricMatrix TmpTriDiag (Dimension);

	  if (EigenstateFlag == false)
	    {
	      ((RealSymmetricMatrix*) Hamiltonian.GetHamiltonian())->Householder(TmpTriDiag, MACHINE_PRECISION);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	    }
	  else
	    {
	      RealMatrix TmpEigenvectors (Dimension, Dimension);
	      ((RealSymmetricMatrix*) Hamiltonian.GetHamiltonian())->Householder(TmpTriDiag, MACHINE_PRECISION, TmpEigenvectors);
	      TmpTriDiag.Diagonalize(TmpEigenvectors);
	      TmpTriDiag.SortMatrixUpOrder(TmpEigenvectors);
	      Eigenstates = new RealVector [NbrEigenvalue];
	      for (int i = 0; i < NbrEigenvalue; ++i)
		{
		  Eigenstates[i] = TmpEigenvectors[i];
		}
	    }

	  // store eigenvalues
	  int Max = Hamiltonian.GetHilbertSpaceDimension();
	  if (Max > NbrEigenvalue)
	    Max = NbrEigenvalue;
	  Eigenvalues = new double [NbrEigenvalue];
	  for (int j = 0; j < Max ; j++)
	    {
	      Eigenvalues[j]= TmpTriDiag.DiagonalElement(j);
	    }
	}
      else
	{
	  cout << (*(Hamiltonian.GetHamiltonian()))(0, 0) << endl;
	}
    }
  else
    {
      
      // architecture type (i.e. 1 CPU or multi CPU)
      AbstractArchitecture* Architecture;
      if (SMPFlag == true)
	Architecture = new SMPArchitecture(2);
      else
	Architecture = new MonoProcessorArchitecture;


      // type of lanczos algorithm (with or without reorthogonalization)
//	  BasicLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
      FullReorthogonalizedLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
	
      // initialization of lanczos algorithm
      double Precision = 1.0;
      double PreviousLowest = 1e50;
      double Lowest = PreviousLowest;
      int CurrentNbrIterLanczos = NbrEigenvalue + 3;
      Lanczos.SetHamiltonian(&Hamiltonian);
      Lanczos.InitializeLanczosAlgorithm();
      Lanczos.RunLanczosAlgorithm(NbrEigenvalue + 2);
      RealTriDiagonalSymmetricMatrix TmpMatrix;

      // run Lancos algorithm up to desired precision on the n-th eigenvalues
      while ((Precision > 1e-14) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
	{
	  Lanczos.RunLanczosAlgorithm(1);
	  TmpMatrix.Copy(Lanczos.GetDiagonalizedMatrix());
	  TmpMatrix.SortMatrixUpOrder();
	  Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
	  if (VerboseFlag == true)
	    cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << endl;
	  Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	  PreviousLowest = Lowest; 
	}      
      if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
	{
	  cout << "too much Lanczos iterations" << endl;
	  exit(0);
	}
      
      // store eigenvalues 
      Eigenvalues = new double [NbrEigenvalue];
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  Eigenvalues[i] = TmpMatrix.DiagonalElement(i);
	}
      cout << endl;
      
      //compute eigenstates
      if (EigenstateFlag == true)
	Eigenstates = (RealVector*) Lanczos.GetEigenstates(NbrEigenvalue);
    }


  // insert here your code using the eigenvalues and the eigenvectors
  for (int i = 0; i < NbrEigenvalue; ++i)
    {
      cout << Eigenvalues[i] << " ";
    }
  cout << endl;
  if ((EigenstateFlag == true) && (Eigenstates != 0))
    {
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  RealVector TmpVector (Dimension);
	  TmpVector.Multiply(HRep, Eigenstates[i]);	  
	  cout << (TmpVector * (RealVector&) Eigenstates[i]) << " ";
	}
      cout << endl;
    }

  return 0;
}
