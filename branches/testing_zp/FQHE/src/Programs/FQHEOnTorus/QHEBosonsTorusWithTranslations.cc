#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "Hamiltonian/ExplicitHamiltonian.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"

#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusState.h"
#include "Hamiltonian/ParticleOnTorusCoulombHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusDeltaWithMagneticTranslationsHamiltonian.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

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
#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption DiskOption ('\n', "disk", "enable disk resume capabilities", false);
  BooleanOption ResumeOption ('\n', "resume", "resume from disk datas", false);
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  SingleIntegerOption IterationOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrIterationOption ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 40);
  BooleanOption GroundOption ('g', "ground", "restrict to the largest subspace");
  SingleIntegerOption NbrBosonOption ('p', "nbr-particles", "number of particles", 3);
  SingleIntegerOption MaxMomentumOption ('l', "max-momentum", "maximum momentum for a single particle", 6);
  SingleIntegerOption MemoryOption ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  SingleIntegerOption MaxFullDiagonalizationOption ('f', "max-full", "maximum hilbert space size allowed to use full diagonalization", 300);
  SingleDoubleOption RatioOption ('r', "ratio", "ratio between the two torus lengths", 0.57735026919);
  SingleIntegerOption VectorMemoryOption ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  SingleStringOption SavePrecalculationOption ('\n', "save-precalculation", "save precalculation in a file",0);
  SingleStringOption LoadPrecalculationOption ('\n', "load-precalculation", "load precalculation from a file",0);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &GroundOption;
  OptionList += &SMPNbrProcessorOption;
  OptionList += &IterationOption;
  OptionList += &NbrIterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &NbrBosonOption;
  OptionList += &MaxMomentumOption;
  OptionList += &MemoryOption;
  OptionList += &MaxFullDiagonalizationOption;
  OptionList += &RatioOption;
  OptionList += &VectorMemoryOption;
  OptionList += &DiskOption;
  OptionList += &ResumeOption;
  OptionList += &LoadPrecalculationOption;
  OptionList += &SavePrecalculationOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsTorus -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }

/*  int Dim = 20;
  HermitianMatrix Mat (Dim, true);
  
  for (int i = 0; i < Dim; ++i)
    {
      Mat.SetMatrixElement(i, i, drand48());
      for (int j = i + 1; j < Dim; ++j)
        Mat.SetMatrixElement(i, j, Complex(drand48(), drand48()));
    }

  UndescribedHilbertSpace Sp (Dim); 
  ExplicitHamiltonian Hamil (&Sp, &Mat);
//  cout << Mat << endl;
  for (int i = 0; i < Dim; ++i)
    {
      cout << " i = " << i << endl;
      ComplexVector V1 (Dim, true);
      ComplexVector V2 (Dim, true);
      ComplexVector V3 (Dim, true);
      V1.Re(i) = drand48();
      V1.Im(i) = drand48();
      Hamil.LowLevelMultiply(V1, V2, 0, 3);
      Hamil.LowLevelAddMultiply(V1, V2, 3, 2);
      Hamil.LowLevelAddMultiply(V1, V2, 5, 5);
      Hamil.LowLevelAddMultiply(V1, V2, 10, 10);
      V3.Multiply (Mat, V1);
      for (int j = 0; j < Dim; ++j)
	{
	  if ((V2.Re(j) != V3.Re(j)) || (V2.Im(j) != V3.Im(j)))
	    cout << "error at " << j << " " << V2.Re(j) << " " << V3.Re(j) << " " << V2.Im(j) << " " << V3.Im(j) << endl;
	}
//      cout << V2 << endl << endl << V3 << endl << endl;
   }
  return 0;*/


  bool GroundFlag = GroundOption.GetBoolean();
  bool SMPFlag = SMPOption.GetBoolean();
  int NbrProcessor = SMPNbrProcessorOption.GetInteger();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrIterLanczos = NbrIterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  int NbrBosons = NbrBosonOption.GetInteger();
  int MaxMomentum = MaxMomentumOption.GetInteger();
  int MaxFullDiagonalization = MaxFullDiagonalizationOption.GetInteger();
  double XRatio = RatioOption.GetDouble();
  long Memory = -1;
  if (MemoryOption.GetInteger() > 0)
    Memory = MemoryOption.GetInteger()<< 20;
  bool ResumeFlag = ResumeOption.GetBoolean();
  bool DiskFlag = DiskOption.GetBoolean();
  int VectorMemory = VectorMemoryOption.GetInteger();
  char* LoadPrecalculationFileName = LoadPrecalculationOption.GetString();
  char* SavePrecalculationFileName = SavePrecalculationOption.GetString();

  int L = 0;
  double GroundStateEnergy = 0.0;

/*  int NbrState = 9;
  int ReducedNbrState = NbrState >> 2;
  int NbrStateRemainder = NbrState - (ReducedNbrState << 2);
  if (NbrStateRemainder == 0)
    {
      NbrStateRemainder = 4;
      --ReducedNbrState;
    }
  for (int k = 0; k < NbrState; ++k)
    {
      BosonOnTorusState State (ReducedNbrState + 1);
      for (int i = 0; i < NbrState; ++i)
	State.SetOccupation(i, 65 | (i << 1));
      BosonOnTorusState TmpState(State, ReducedNbrState + 1);
      State.LeftShiftState(ReducedNbrState, NbrStateRemainder, k);
      for (int i = 0; i < NbrState; ++i)
	if (((i >=  k) && (State.GetOccupation(i - k) != TmpState.GetOccupation(i)))
	    || ((i <  k) && (State.GetOccupation(NbrState + i -  k) != TmpState.GetOccupation(i))))
	  cout << "error " << i << endl;
      State.PrintState(cout, ReducedNbrState, NbrStateRemainder) << endl;
    }
  return 0;*/

  char* OutputNameLz = new char [512];
  sprintf (OutputNameLz, "bosons_torus_delta_n_%d_2s_%d_ratio_%f.dat", NbrBosons, MaxMomentum, XRatio);
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);


  
  AbstractArchitecture* Architecture = 0;
  if (SMPFlag == false)
    Architecture = new MonoProcessorArchitecture;
  else
    Architecture = new SMPArchitecture(NbrProcessor);

  int MomentumModulo = FindGCD(NbrBosons, MaxMomentum);
  //  MomentumModulo = 1;
  for (int x = 0; x < MomentumModulo; ++x)
    for (int y = 0; y < MomentumModulo; ++y)
    {     
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
//      BosonOnTorus TotalSpace (NbrBosons, MaxMomentum, y);
      BosonOnTorusWithMagneticTranslations TotalSpace2 (NbrBosons, MaxMomentum, x, y);
      cout << " Total Hilbert space dimension = " << TotalSpace2.GetHilbertSpaceDimension() << endl;
      cout << "momentum = (" << x << "," << y << ")" << endl;
//       for (int i = 0; i < TotalSpace2.GetHilbertSpaceDimension(); ++i)
// 	{
// 	  cout << i << " = ";
// 	  TotalSpace2.PrintState(cout, i) << endl;
// 	}
      cout << endl << endl;
      }
  return 0;
  for (int x = 0; x < MomentumModulo; ++x)
    for (int y = 0; y < MomentumModulo; ++y)
    {     
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
//      BosonOnTorus TotalSpace (NbrBosons, MaxMomentum, y);
      BosonOnTorusWithMagneticTranslations TotalSpace (NbrBosons, MaxMomentum, x, y);
      cout << " Total Hilbert space dimension = " << TotalSpace.GetHilbertSpaceDimension() << endl;
/*      cout << "momentum = (" << x << "," << y << ")" << endl;
      for (int i = 0; i < TotalSpace.GetHilbertSpaceDimension(); ++i)
	{
	  cout << i << " = ";
	  TotalSpace.PrintState(cout, i) << endl;
	}
      cout << endl << endl;*/
/*      for (int i = 0; i < TotalSpace.GetHilbertSpaceDimension(); ++i)
	{
	  cout << "---------------------------------------------" << endl;
	  cout << i << " = " << endl;;
	  for (int m1 = 0; m1 < MaxMomentum; ++m1)
	    for (int m2 = 0; m2 < m1; ++m2)
	      for (int m3 = 0; m3 < MaxMomentum; ++m3)
		{
		  int m4 = m1 + m2 - m3;
		  if (m4 < 0)
		    m4 += MaxMomentum;
		  else
		    if (m4 >= MaxMomentum)
		      m4 -= MaxMomentum;
		  if (m3 > m4)
		    {
		      double Coefficient = 0.0;
		      int NbrTranslations = 0;
		      TotalSpace.AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslations);
		    }
		}
	}*/
      AbstractArchitecture* Architecture = 0;
      if (SMPFlag == false)
	Architecture = new MonoProcessorArchitecture;
      else
	Architecture = new SMPArchitecture(NbrProcessor);
      AbstractHamiltonian* Hamiltonian = new ParticleOnTorusDeltaWithMagneticTranslationsHamiltonian (&TotalSpace, NbrBosons, MaxMomentum, x,
												      XRatio, Architecture, Memory);
      if (Hamiltonian->GetHilbertSpaceDimension() < MaxFullDiagonalization)
	{
	  HermitianMatrix HRep2 (Hamiltonian->GetHilbertSpaceDimension());
	  Hamiltonian->GetHamiltonian(HRep2);
//	  cout << HRep2 << endl;
	  RealSymmetricMatrix HRep (HRep2.ConvertToSymmetricMatrix());
	  if (Hamiltonian->GetHilbertSpaceDimension() > 1)
	    {
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian->GetHilbertSpaceDimension(), true);
	      HRep.Householder(TmpTriDiag, 1e-7);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      if (L == 0)
		GroundStateEnergy = TmpTriDiag.DiagonalElement(0);
	      for (int j = 0; j < Hamiltonian->GetHilbertSpaceDimension() ; j++)
		{
		  File << x << " " << y << " " << TmpTriDiag.DiagonalElement(2 * j) << endl;
		  cout << x << " " << y << " " << TmpTriDiag.DiagonalElement(2 * j) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
	      File << x << " " << y << " " << HRep(0, 0) << endl;
	    }
	}
      else
	{
	  int MaxNbrIterLanczos = 4000;
	  AbstractLanczosAlgorithm* Lanczos;
	  if (NbrEigenvalue == 1)
	    {
	      if (DiskFlag == false)
		Lanczos = new ComplexBasicLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	      else
		Lanczos = new ComplexBasicLanczosAlgorithmWithDiskStorage(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	    }
	  else
	    {
	      if (DiskFlag == false)
		Lanczos = new FullReorthogonalizedComplexLanczosAlgorithm (Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	      else
		Lanczos = new FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage (Architecture, NbrEigenvalue, VectorMemory, MaxNbrIterLanczos);
	    }
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = NbrEigenvalue + 3;
	  Lanczos->SetHamiltonian(Hamiltonian);
	  if ((DiskFlag == true) && (ResumeFlag == true))
	    Lanczos->ResumeLanczosAlgorithm();
	  else
	    Lanczos->InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  double Dt;
	  gettimeofday (&(TotalStartingTime), 0);
	  if (ResumeFlag == false)
	    {
	      Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	      CurrentNbrIterLanczos = NbrEigenvalue + 3;
	    }
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while ((Lanczos->TestConvergence() == false) && (((DiskFlag == true) && (CurrentNbrIterLanczos < NbrIterLanczos)) ||
							   ((DiskFlag == false) && (CurrentNbrIterLanczos < MaxNbrIterLanczos))))
	    {
	      Lanczos->RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << " "<< endl;
	      ++CurrentNbrIterLanczos;
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
	      File << x << " " << y << " " << TmpMatrix.DiagonalElement(i) << endl;
	    }
	  cout << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	}
      cout << "----------------------------------------------------------------" << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrBosons) << endl;
      delete Hamiltonian;

    }
  File.close();

  return 0;
}
