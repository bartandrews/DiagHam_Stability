#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnTorusWithSpinNew.h"
#include "Hamiltonian/ParticleOnTorusCoulombHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithSpinHamiltonian.h"
//#include "Hamiltonian/ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/Options.h"

#include "MainTask/GenericRealMainTask.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);
    
  // some running options and help
  OptionManager Manager ("FQHETorusFermionsWithSpinAndTranslations" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "max-momentum", "maximum momentum for a single particle", 18);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-spin", "total spin of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "constraint on the total momentum in the y direction (negative if none)", -1);
  (*SystemGroup) += new SingleDoubleOption   ('R', "ratio", 
					      "ratio between lengths along the x and y directions (-1 if has to be taken equal to nbr-particles/4)", -1.0);
  (*SystemGroup) += new SingleDoubleOption   ('d', "layerSeparation", 
					      "for bilayer simulations: layer separation in magnetic lengths", 0.0);
  (*SystemGroup) += new  BooleanOption  ('\n', "redundantYMomenta", "Calculate all subspaces up to YMomentum = MaxMomentum-1", false);

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusFermionsWithSpinAndTranslations -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int TotalSpin = Manager.GetInteger("total-spin");
  int NbrFermions = Manager.GetInteger("nbr-particles");
int MaxMomentum = Manager.GetInteger("max-momentum");
  int YMomentum = Manager.GetInteger("y-momentum");
  double XRatio = NbrFermions / 4.0;
  if (Manager.GetDouble("ratio") > 0)
    {
       XRatio = Manager.GetDouble("ratio");
    }
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  double LayerSeparation = Manager.GetDouble("layerSeparation");

  char* OutputFileName = new char [512];
  if (LayerSeparation==0.0)
    sprintf (OutputFileName, "fermions_torus_su2_coulomb_n_%d_2s_%d_Sz_%d_ratio_%f.dat", NbrFermions, MaxMomentum, TotalSpin, XRatio);
  else
    sprintf (OutputFileName, "fermions_torus_d_%f_coulomb_n_%d_2s_%d_Sz_%d_ratio_%f.dat", LayerSeparation, 
NbrFermions,MaxMomentum, TotalSpin, XRatio);
  ofstream File;
  File.open(OutputFileName, ios::binary | ios::out);
  File.precision(14);


  int MomentumModulo = FindGCD(NbrFermions, MaxMomentum);
  int YMaxMomentum;
  if (Manager.GetBoolean("redundantYMomenta"))
    YMaxMomentum = (MaxMomentum - 1);
  else
    YMaxMomentum = (MomentumModulo - 1);
  if (YMomentum < 0)
    YMomentum = 0;
  else
    YMaxMomentum = YMomentum; 

  char *SubspaceStr = new char[50];  
  char *SubspaceLegend= new char[50];
  sprintf(SubspaceLegend,"Ky");
  bool FirstRun=true;
  for (int YMomentum2 = YMomentum; YMomentum2 <= YMaxMomentum; ++YMomentum2)
    {
      sprintf(SubspaceStr,"%d", YMomentum2);
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
      FermionOnTorusWithSpinNew TotalSpace (NbrFermions, TotalSpin, MaxMomentum, YMomentum2);	
      cout << " Total Hilbert space dimension = " << TotalSpace.GetHilbertSpaceDimension() << endl;
      //      cout << "momentum = " << Momentum << endl;
      cout << "momentum Ky = " << YMomentum2 << endl;
/*	for (int i = 0; i < TotalSpace.GetHilbertSpaceDimension(); ++i)
	  {
	    cout << i << " = ";
	    TotalSpace.PrintState(cout, i) << endl;
	  }
	cout << endl << endl;
	exit(0);*/
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
	
	Architecture.GetArchitecture()->SetDimension(TotalSpace.GetHilbertSpaceDimension());
	
	AbstractHamiltonian* Hamiltonian = // new ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian
	  new ParticleOnTorusCoulombWithSpinHamiltonian(&TotalSpace, NbrFermions, MaxMomentum, XRatio,
							/*Zeeman*/ 0.0, LayerSeparation, Memory);


	char* EigenvectorName = 0;
	if ( Manager.GetBoolean("eigenstate") == true)	
	  {
	    EigenvectorName = new char [100];
	    if (LayerSeparation==0.0)
	      sprintf (OutputFileName, "fermions_torus_su2_coulomb_n_%d_2s_%d_Sz_%d_ratio_%f_k_%d", NbrFermions, MaxMomentum,
			TotalSpin, XRatio,YMomentum2);
	    else
	      sprintf (OutputFileName, "fermions_torus_d_%f_coulomb_n_%d_2s_%d_Sz_%d_ratio_%f_k_%d", LayerSeparation, NbrFermions,
		       MaxMomentum, TotalSpin, XRatio,YMomentum2);
	  }
	
	GenericRealMainTask Task (&Manager, &TotalSpace, &Lanczos, Hamiltonian, SubspaceStr, SubspaceLegend,
				  /*Shift*/ 0.0, OutputFileName, FirstRun, EigenvectorName);
	MainTaskOperation TaskOperation (&Task);
	TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	if (EigenvectorName != 0)
	  {
	    delete[] EigenvectorName;
	  }
	
	if (FirstRun == true)
	  FirstRun = false;
    
	delete Hamiltonian;
    }
  File.close();
  delete[] OutputFileName;
  return 0;
}
