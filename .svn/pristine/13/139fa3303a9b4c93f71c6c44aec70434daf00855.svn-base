#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinPartialPolarization.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinPartialPolarization.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Operator/ParticleOnSphereWithSpinDensityOperator.h"
#include "Operator/ParticleOnSphereWithSpinDensityDensityOperator.h"

#include "FunctionBasis/ParticleOnCylinderFunctionBasis.h"
#include "Hamiltonian/ParticleOnCylinderStructureFactor.h"
#include "Hamiltonian/ParticleOnCylinderPairAmplitude.h"
#include "Hamiltonian/ParticleOnCylinderDensityDensity.h"

#include "Tools/FQHEFiles/FQHEOnCylinderFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "Options/OptionManager.h"
#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

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

  // some running options and help
  OptionManager Manager ("FQHECylinderWithSU2SpinDensity" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('\n', "input-state", "name of the file containing the state to convert");
  (*SystemGroup) += new SingleStringOption('\n', "degenerate-states", "name of the file containing a list of states (override input-state)");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state of the output state");
  (*SystemGroup) += new BooleanOption  ('\n', "all-sz", "assume a hilbert space including all sz values");
  (*SystemGroup) += new BooleanOption  ('\n', "use-alt", "use alternative Hilbert space for  bosonic states");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-points", "number of points along the cylinder axis", 100);
  (*OutputGroup) += new SingleDoubleOption  ('\n', "offset", "additional length along the cylinder axis on each side of the [-Lx/2,Lx/2] region where the density should be computed", 5.0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderWithSU2SpinDensity -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int KyMax = 0;
  int TotalKy = 0;
  int TotalSz = 0;
  bool FermionFlag = true;  
  bool SzSymmetrizedBasis = false;
  bool SzMinusParity = false;
  bool LzSymmetrizedBasis = false;
  bool LzMinusParity = false;
  int NbrPoints = Manager.GetInteger("nbr-points");
  double Ratio = 0.0;
  double Perimeter = 0.0;
  int LandauLevel = 0;
  int NbrPolarizedOrbitals = 0;
  int NbrInputStates = 0;
  char** InputStateNames = 0;

  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0) )
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type FQHECylinderWithSU2SpinDensity -h" << endl;
      return -1;
    }

  if (Manager.GetString("input-state") != 0)
    {
      if (FQHEOnCylinderWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, KyMax, TotalKy, TotalSz, 
								 FermionFlag, Ratio, Perimeter) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state")  << endl;
	  return -1;
	}
      if (FQHEOnSphereWithSpinFindSystemPolarizationFromFileName(Manager.GetString("input-state"), NbrPolarizedOrbitals) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state")  << endl;
	}
      NbrInputStates = 1;
      InputStateNames = new char*[NbrInputStates];
      InputStateNames[0] = new char [strlen(Manager.GetString("input-state")) + 1];
      strcpy (InputStateNames[0], Manager.GetString("input-state"));
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      NbrInputStates = DegenerateFile.GetNbrLines();
      InputStateNames = new char*[NbrInputStates];
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  InputStateNames[i] = new char [strlen(DegenerateFile(0, i)) + 1];
	  strcpy (InputStateNames[i], DegenerateFile(0, i));
	}
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  if (FQHEOnCylinderWithSpinFindSystemInfoFromVectorFileName(InputStateNames[i], NbrParticles, KyMax, TotalKy, TotalSz, FermionFlag, Ratio, Perimeter) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << InputStateNames[i] << endl;
	      return -1;
	    }
	  if (FQHEOnSphereWithSpinFindSystemPolarizationFromFileName(InputStateNames[i], NbrPolarizedOrbitals) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state")  << endl;
	    }
	}
    }

  if ((Perimeter == 0.0) & (Ratio == 0.0))
    {
      cout << "error, neither the aspect ratio or the cylinder perimeter is defined" << endl;
      return -1;
    }
  if (Perimeter == 0.0)
    {
      Perimeter = sqrt(2.0 * M_PI * (KyMax + 1) * Ratio);
    }
  else
    {
      Ratio = (Perimeter * Perimeter) / (2.0 * M_PI * (KyMax + 1));
    }
  double CylinderLength = Perimeter / Ratio;

  ParticleOnSphereWithSpin* Space;
  if (FermionFlag == true)
    {
      if (Manager.GetBoolean("all-sz"))
	{
	  if (LzSymmetrizedBasis == false)
	    Space = new FermionOnSphereWithSpinAllSz (NbrParticles, TotalKy, KyMax);
	  else
	    Space = new FermionOnSphereWithSpinAllSzLzSymmetry (NbrParticles, KyMax, LzMinusParity);
	}
      else if ((SzSymmetrizedBasis == false) && (LzSymmetrizedBasis == false))
	  {
#ifdef __64_BITS__
	  if (KyMax <= 31)
#else
	    if (KyMax <= 15)
#endif
	      {
		if (NbrPolarizedOrbitals == 0)
		  {
		    Space = new FermionOnSphereWithSpin(NbrParticles, TotalKy, KyMax, TotalSz);
		  }
		else
		  {
		    Space = new FermionOnSphereWithSpinPartialPolarization(NbrParticles, TotalKy, KyMax, TotalSz, NbrPolarizedOrbitals);
		  }
	      }
	    else
	      {
#ifdef __128_BIT_LONGLONG__
		if (KyMax <= 63)
#else
		  if (KyMax <= 31)
#endif
		    {
		      Space = new FermionOnSphereWithSpinLong(NbrParticles, TotalKy, KyMax, TotalSz);
		    }
		  else
		    {
		      cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		      return -1;
		    }	
	      }
	  }
      else
	{
#ifdef __128_BIT_LONGLONG__
	  if (KyMax >= 61)
#else
	    if (KyMax >= 29)
#endif
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
	  if (SzSymmetrizedBasis == true) 
	    if (LzSymmetrizedBasis == false)
	      {
#ifdef __64_BITS__
		if (KyMax <= 28)
#else
		  if (KyMax <= 13)
#endif
		    {
		      if (Manager.GetString("load-hilbert") == 0)
			Space = new FermionOnSphereWithSpinSzSymmetry(NbrParticles, TotalKy, KyMax, SzMinusParity);
		      else
			Space = new FermionOnSphereWithSpinSzSymmetry(Manager.GetString("load-hilbert"));
		    }
		  else
		    {
		      if (Manager.GetString("load-hilbert") == 0)
			Space = new FermionOnSphereWithSpinSzSymmetryLong(NbrParticles, TotalKy, KyMax, SzMinusParity);
		      else
			Space = new FermionOnSphereWithSpinSzSymmetryLong(Manager.GetString("load-hilbert"));
		    }
		  }
	    else
#ifdef __64_BITS__
	      if (KyMax <= 28)
#else
		if (KyMax <= 13)
#endif
		  {
		    if (Manager.GetString("load-hilbert") == 0)
		      {
			Space = new FermionOnSphereWithSpinLzSzSymmetry(NbrParticles, KyMax, SzMinusParity,
									LzMinusParity);
		      }
		    else
		      Space = new FermionOnSphereWithSpinLzSzSymmetry(Manager.GetString("load-hilbert"));
		  }
		else
		  {
		    if (Manager.GetString("load-hilbert") == 0)
		      {
			Space = new FermionOnSphereWithSpinLzSzSymmetryLong(NbrParticles, KyMax, SzMinusParity,
									    LzMinusParity);
		      }
		    else
		      Space = new FermionOnSphereWithSpinLzSzSymmetryLong(Manager.GetString("load-hilbert"));
		    
		  }
	      else
#ifdef __64_BITS__
		if (KyMax <= 28)
#else
		  if (KyMax <= 13)
#endif
		    {
		      if (Manager.GetString("load-hilbert") == 0)
			Space = new FermionOnSphereWithSpinLzSymmetry(NbrParticles, KyMax, TotalSz, LzMinusParity);
		      else
			Space = new FermionOnSphereWithSpinLzSymmetry(Manager.GetString("load-hilbert"));	      
		    }
		  else
		    {
		      if (Manager.GetString("load-hilbert") == 0)
			Space = new FermionOnSphereWithSpinLzSymmetryLong(NbrParticles, KyMax, TotalSz, LzMinusParity);
		      else
			Space = new FermionOnSphereWithSpinLzSymmetryLong(Manager.GetString("load-hilbert"));	      
		    }
	}
    }
  else
    {
      if (Manager.GetBoolean("all-sz"))
	{
	  Space = new BosonOnSphereWithSpinAllSz (NbrParticles, TotalKy, KyMax);
	}
      else
	{
	  if (NbrPolarizedOrbitals == 0)
	    {
	      if (Manager.GetBoolean("use-alt") == false)
		{
		  Space = new BosonOnSphereWithSpin(NbrParticles, TotalKy, KyMax, TotalSz);
		}
	      else
		{
		  Space = new BosonOnSphereWithSU2Spin(NbrParticles, TotalKy, KyMax, TotalSz);
		}
	    }
	  else
	    {
	      Space = new BosonOnSphereWithSU2SpinPartialPolarization(NbrParticles, TotalKy, KyMax, TotalSz, NbrPolarizedOrbitals);
	    }
	}
    }
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

  
  RealVector* InputStates = new RealVector[NbrInputStates];
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (InputStates[i].ReadVector (InputStateNames[i]) == false)
	{
	  cout << "can't open vector file " << InputStateNames[i] << endl;
	  return -1;      
		}	  
      if (Space->GetHilbertSpaceDimension() != InputStates[i].GetVectorDimension())
	{
	  cout << "error, " << InputStateNames[i] << " does not have the correct dimension (is "
	       << InputStates[i].GetVectorDimension() << ", should be " << Space->GetHilbertSpaceDimension() << ")" << endl;
	  return -1;
	}
    }

  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = ReplaceExtensionToFileName(InputStateNames[0], "vec", "rho.dat");
    }


  ofstream File;
  File.precision(14);
  File.open(OutputFileName, ios::binary | ios::out);
  File << "# Length = " << CylinderLength << ", Perimeter = " << Perimeter  << ", Aspect ratio = " << Ratio << endl;
  HermitianMatrix* OccupationMatrixElementsUp = new HermitianMatrix[KyMax + 1];
  HermitianMatrix* OccupationMatrixElementsDown = new HermitianMatrix[KyMax + 1];
  
      
  for (int i = 0; i <= KyMax; ++i)
    {
      OccupationMatrixElementsUp[i] = HermitianMatrix(NbrInputStates, true);
      OccupationMatrixElementsDown[i] = HermitianMatrix(NbrInputStates, true);
      ParticleOnSphereWithSpinDensityOperator TmpOperatorUp (Space, i, 0, i, 0);
      ParticleOnSphereWithSpinDensityOperator TmpOperatorDown (Space, i, 1, i, 1);
      RealVector TmpStateUp(Space->GetHilbertSpaceDimension(), true);
      RealVector TmpStateDown(Space->GetHilbertSpaceDimension(), true);
      for (int j = 0; j < NbrInputStates; ++j)
	{
	  VectorOperatorMultiplyOperation OperationUp (&TmpOperatorUp, &InputStates[j], &TmpStateUp);
	  OperationUp.ApplyOperation(Architecture.GetArchitecture());
	  VectorOperatorMultiplyOperation OperationDown (&TmpOperatorDown, &InputStates[j], &TmpStateDown);
	  OperationDown.ApplyOperation(Architecture.GetArchitecture());
	  for (int k = j; k <  NbrInputStates; ++k)
	    {
	      Complex TmpUp = InputStates[k] * TmpStateUp;
	      OccupationMatrixElementsUp[i].SetMatrixElement(k, j, TmpUp);
	      Complex TmpDown = InputStates[k] * TmpStateDown;
	      OccupationMatrixElementsDown[i].SetMatrixElement(k, j, TmpDown);
	      File << "# " << i << " " << j << " " << k << " " << TmpUp << " " << TmpDown << endl;
	    }
	}
    }
  File << "# x rho_up integral(rho_up) rho_down integral(rho_down)" << endl; 
  ParticleOnCylinderFunctionBasis Basis (KyMax, LandauLevel, Ratio);
  HermitianMatrix TmpMatrixUp (NbrInputStates, true);
  HermitianMatrix TmpMatrixUp2 (NbrInputStates, true);
  HermitianMatrix TmpMatrixDown (NbrInputStates, true);
  HermitianMatrix TmpMatrixDown2 (NbrInputStates, true);
  RealDiagonalMatrix TmpEigenvaluesUp(NbrInputStates);
  RealDiagonalMatrix TmpEigenvaluesUp2(NbrInputStates);
  RealDiagonalMatrix TmpEigenvaluesDown(NbrInputStates);
  RealDiagonalMatrix TmpEigenvaluesDown2(NbrInputStates);
  double Offset = Manager.GetDouble("offset");
  double XPosition = -(0.5 * CylinderLength) - Offset;
  double XStep = -2.0 * XPosition / ((double) (NbrPoints + 1));
  double TmpPrefactor1 =  1.0 / sqrt(M_PI);
  double TmpPrefactor2 =  2.0 * M_PI / Perimeter;
  double TmpShift = 0.5 * ((double) KyMax) * TmpPrefactor2;
  double* TotalChargesUp = new double[NbrInputStates];
  double* TotalChargesDown = new double[NbrInputStates];
  for (int i = 0; i < NbrInputStates; ++i)
    {   
      TotalChargesUp[i] = 0.0;
      TotalChargesDown[i] = 0.0;
    }
  for (int TmpX = 0; TmpX <= NbrPoints; ++TmpX)
    { 
      TmpMatrixUp.ClearMatrix();
      TmpMatrixUp2.ClearMatrix();
      TmpMatrixDown.ClearMatrix();
      TmpMatrixDown2.ClearMatrix();
      for (int i = 0; i <= KyMax; ++i)
	{
	  double TmpFactor = TmpPrefactor1 * exp(-(XPosition + TmpShift - (TmpPrefactor2 * ((double) i))) * (XPosition + TmpShift - (TmpPrefactor2 * ((double) i))));
	  TmpMatrixUp.AddLinearCombination(TmpFactor, OccupationMatrixElementsUp[i]);
	  TmpMatrixDown.AddLinearCombination(TmpFactor, OccupationMatrixElementsDown[i]);
	  TmpFactor = 0.5 * (1.0 + erf(XPosition + TmpShift - (TmpPrefactor2 * ((double) i))));
	  TmpMatrixUp2.AddLinearCombination(TmpFactor, OccupationMatrixElementsUp[i]);
	  TmpMatrixDown2.AddLinearCombination(TmpFactor, OccupationMatrixElementsDown[i]);
	}
      TmpMatrixUp.LapackDiagonalize(TmpEigenvaluesUp);
      TmpMatrixUp2.LapackDiagonalize(TmpEigenvaluesUp2);
      TmpMatrixDown.LapackDiagonalize(TmpEigenvaluesDown);
      TmpMatrixDown2.LapackDiagonalize(TmpEigenvaluesDown2);
      File << XPosition;
      for (int i = 0; i < NbrInputStates; ++i)
	{   
	  File << " " << TmpEigenvaluesUp[i] << " " << TmpEigenvaluesUp2[i] << " " << TmpEigenvaluesDown[i] << " " << TmpEigenvaluesDown2[i];
	  TotalChargesUp[i] +=  TmpEigenvaluesUp[i] * XStep;
	  TotalChargesDown[i] +=  TmpEigenvaluesDown[i] * XStep;
	}      
      File << endl;
      XPosition += XStep;
    }

  cout << "total charge up : " << endl;
  for (int i = 0; i < NbrInputStates; ++i)
    {   
      cout << TotalChargesUp[i] << " ";
    }
  cout << endl;
  cout << "total charge down : " << endl;
  for (int i = 0; i < NbrInputStates; ++i)
    {   
      cout << TotalChargesDown[i] << " ";
    }
  cout << endl;
  File.close();

  return 0;
}

