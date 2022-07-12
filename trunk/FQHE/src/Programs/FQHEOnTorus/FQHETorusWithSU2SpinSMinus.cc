#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithSpin.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorusWithSpin.h"

#include "Operator/ParticleOnSphereWithSpinSMinusOperator.h"
#include "Operator/ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

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
  OptionManager Manager ("FQHETorusWithSU2SpinSMinus" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-file", "name of the file corresponding to the input state");
  (*SystemGroup) += new  BooleanOption ('\n', "fully-unpolarized", "apply the S- minus as many times as need to reach the Sz=0 (or Sz=1/2) sector");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (changing the Sz value)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusWithSU2SpinSMinus -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-file") == 0)
    {
      cout << "error, an input state file should be provided. See man page for option syntax or type FQHETorusWithSU2SpinSMinus -h" << endl;
      return -1;
    }
  if ((Manager.GetString("input-file") != 0) && 
      (IsFile(Manager.GetString("input-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("input-file") << endl;
      return -1;
    }

  int NbrParticles = 0; 
  int KyMax = 0; 
  int TotalKx = 0;
  int TotalKy = 0;
  int TotalSz = 0;
  bool ComplexFlag = false;
  bool Statistics = true;
  bool KxFlag = false;
  bool MagneticTranslationFlag = false;

  char* InputStateFile = new char [strlen(Manager.GetString("input-file")) + 1];
  strcpy(InputStateFile, Manager.GetString("input-file"));

  if (FQHEOnTorusWithSpinFindSystemInfoFromVectorFileName(InputStateFile, NbrParticles, KyMax,  TotalKx, TotalKy, TotalSz, Statistics) == false)
    {
      if (FQHEOnTorusWithSpinFindSystemInfoFromVectorFileName(InputStateFile, NbrParticles, KyMax, TotalKy, TotalSz, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << InputStateFile << endl;
	  return -1;
	}
    }
  else
    {
      MagneticTranslationFlag = true;
    }
  
  int MinSz = TotalSz - 2;
  int CurrentSz = MinSz;
  if (Manager.GetBoolean("fully-unpolarized"))
    {
      if ((NbrParticles & 1) == 0)
	MinSz = 0;
      else
	MinSz = 1;
    }
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") == 0)
    {
      char* TmpString = strstr(InputStateFile, "_sz_");
      if (TmpString == 0)
	{
	  cout << InputStateFile << " file name does not contain _sz_, can't create output file name" << endl;
	  return 0;
	}
      int SizeString = 4;
      while ((TmpString[SizeString] != '\0') && (TmpString[SizeString] != '_') && (TmpString[SizeString] != '.') &&
	     (((TmpString[SizeString] >= '0') && (TmpString[SizeString] <= '9')) || (TmpString[SizeString] != '-')))
	++SizeString;
       if (((TmpString[SizeString] == '_') || (TmpString[SizeString] == '.')) && (SizeString != 4))
	 {
	   char TmpChar = TmpString[4];
	   TmpString[4] = '\0';
	   OutputFileName = new char [strlen(InputStateFile) + 16 + strlen(TmpString + SizeString)];
	   sprintf (OutputFileName, "%s%d%s", InputStateFile, MinSz, TmpString + SizeString);
	   TmpString[4] = TmpChar;
	 }
       else
	 {
	   cout << "error while replacing sz value in " << InputStateFile << endl;
	   return 0;
	 }     
    }
  else
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }

  RealVector InputState;
  ComplexVector ComplexInputState;
  ComplexFlag = ComplexInputState.ReadVectorTest(InputStateFile);

  if(ComplexFlag == false)
    {
      if (InputState.ReadVector (InputStateFile) == false)
	{
	  cout << "can't open vector file " << InputStateFile << endl;
	  return -1;      
	}
    }
  else
    {
      if (ComplexInputState.ReadVector (InputStateFile) == false)
	{
	  cout << "can't open vector file " << InputStateFile << endl;
	  return -1;      
	}
    }

  if (MagneticTranslationFlag == false)
    {
      ParticleOnSphereWithSpin* InputSpace;
      if (Statistics == true)
	{
	  InputSpace = new FermionOnTorusWithSpin (NbrParticles, KyMax, TotalSz, TotalKy);
	}
      else
	{
	  InputSpace = new BosonOnTorusWithSpin (NbrParticles, KyMax, TotalSz, TotalKy);
	}
      
      if (ComplexFlag == false)
	{
	  if (InputSpace->GetLargeHilbertSpaceDimension() != InputState.GetLargeVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and ground state" << endl;
	      return 0;
	    }
	}
      else
	{
	  if (InputSpace->GetLargeHilbertSpaceDimension() != ComplexInputState.GetLargeVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and ground state" << endl;
	      return 0;
	    }
	}
      
      
      while (CurrentSz >= MinSz)
	{
	  ParticleOnSphereWithSpin* OutputSpace;
	  if (Statistics == true)
	    {
	      InputSpace = new FermionOnTorusWithSpin (NbrParticles, KyMax, CurrentSz, TotalKy);
	    }
	  else
	    {
	      OutputSpace = new BosonOnTorusWithSpin (NbrParticles, KyMax, CurrentSz, TotalKy);
	    }
	  InputSpace->SetTargetSpace(OutputSpace);
	  ParticleOnSphereWithSpinSMinusOperator SMinusOperator(InputSpace);
	  if (ComplexFlag == false)
	    {
	      RealVector TmpVector (OutputSpace->GetHilbertSpaceDimension());
	      VectorOperatorMultiplyOperation Operation(&SMinusOperator, &InputState, &TmpVector);
	      Operation.ApplyOperation(Architecture.GetArchitecture());
	      InputState = TmpVector;
	      InputState /= InputState.Norm();
	    }
	  else
	    {
	      ComplexVector TmpVector (OutputSpace->GetHilbertSpaceDimension());
	      VectorOperatorMultiplyOperation Operation(&SMinusOperator, &ComplexInputState, &TmpVector);
	      Operation.ApplyOperation(Architecture.GetArchitecture());
	      ComplexInputState = TmpVector;
	      ComplexInputState /= ComplexInputState.Norm();
	    }
	  delete InputSpace;
	  InputSpace = OutputSpace;
	  CurrentSz -= 2;
	}
      if (ComplexFlag == false)
	{
	  if (InputState.WriteVector(OutputFileName) == false)
	    {
	      cout << "can't write " << OutputFileName << endl;
	    }
	}
      else
	{
	  if (ComplexInputState.WriteVector(OutputFileName) == false)
	    {
	      cout << "can't write " << OutputFileName << endl;
	    }
	}
      delete InputSpace;
    }
  else
    {
      ParticleOnTorusWithSpinAndMagneticTranslations* InputSpace;
      cout << NbrParticles << " " << KyMax << " " << TotalSz << " " <<  TotalKx << " " <<  TotalKy << endl;
      if (Statistics == true)
	{
	  cout << "error : fermions are not yet supported" << endl;
	  return 0;
	}
      else
	{
	  InputSpace = new BosonOnTorusWithSpinAndMagneticTranslations (NbrParticles, TotalSz, KyMax, TotalKx, TotalKy);
	}
      
      if (InputSpace->GetLargeHilbertSpaceDimension() != ComplexInputState.GetLargeVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space and ground state (" << InputSpace->GetLargeHilbertSpaceDimension() << " and " << ComplexInputState.GetLargeVectorDimension() << ")" <<endl;
	  return 0;
	}
      
      
      while (CurrentSz >= MinSz)
	{
	  ParticleOnTorusWithSpinAndMagneticTranslations* OutputSpace;
	  if (Statistics == true)
	    {
	      cout << "error : fermions are not yet supported" << endl;
	      return 0;
	    }
	  else
	    {
	      OutputSpace = new BosonOnTorusWithSpinAndMagneticTranslations (NbrParticles, CurrentSz, KyMax, TotalKx, TotalKy);
	    }
	  InputSpace->SetTargetSpace(OutputSpace);
	  ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator SMinusOperator(InputSpace, KyMax, TotalKx);
	  ComplexVector TmpVector (OutputSpace->GetHilbertSpaceDimension());
	  VectorOperatorMultiplyOperation Operation(&SMinusOperator, &ComplexInputState, &TmpVector);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  ComplexInputState = TmpVector;
	  ComplexInputState /= ComplexInputState.Norm();
	  delete InputSpace;
	  InputSpace = OutputSpace;
	  CurrentSz -= 2;
	}
      if (ComplexInputState.WriteVector(OutputFileName) == false)
	{
	  cout << "can't write " << OutputFileName << endl;
	}
      delete InputSpace;
    }
  return 0;
}
