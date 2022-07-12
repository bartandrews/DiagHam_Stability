#include "config.h"

#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"


#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;



int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereFermionsWithSpinRotation" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the file that contains the state");
  (*SystemGroup) += new SingleDoubleOption  ('a', "angle", "Rotation angle", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "sz", "Sz sector of the initial state", 1000000);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lz", "Lz sector of the initial state", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "u1", "if the original state belongs to a u1 space", false);
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsWithSpinRotation -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(Manager.GetString("state") == 0)
    {
      cout << "no input state " << endl << "see man page for option syntax or type FQHESphereFermionsWithSpinRotation -h" << endl;
      return -1;
    }

  int NbrParticles = 0;
  int LzMax = 0;
  double theta = Manager.GetDouble("angle");
  bool IsU1 = Manager.GetBoolean("u1");
  int TotalSz = Manager.GetInteger("sz");
  bool IsSU2 = (TotalSz != 1000000 );
  int TotalLz = Manager.GetInteger("lz");

  bool LzSymmetrizedBasis = false;
  bool LzMinusParity = false;
  bool FermionFlag = true;

  char* StateFileName = Manager.GetString("state");

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(StateFileName, NbrParticles, LzMax, TotalLz, FermionFlag) == false)
    {
      cout << "error while retrieving system informations from file name " << Manager.GetString("state") << endl;
      return -1;
    }

  cout << "N=" << NbrParticles << "  LzMax=" << LzMax << "  TotalLz=" << TotalLz << endl;

  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }

  RealVector State;
  if (State.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }


  long MemorySpace = 9l << 20;
    
  char* OldExtension = new char [512];
  char* NewExtension = new char [512];
  if (IsSU2)
  {
    sprintf (OldExtension, "n_%d_2s_%d_sz_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalSz, TotalLz);  
    sprintf (NewExtension, "theta_%g_n_%d_2s_%d_lz_%d.0.vec", theta, NbrParticles, LzMax, TotalLz);   
  }
  else
  {
    sprintf (OldExtension, "n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);  
    sprintf (NewExtension, "theta_%g_n_%d_2s_%d_lz_%d.0.vec", theta, NbrParticles, LzMax, TotalLz);   
  }
  char* OutputName = ReplaceExtensionToFileName(Manager.GetString("state"), OldExtension, NewExtension);
   
//  printf("%s\n",OutputName);  
    
  
  ComplexVector InitialState ;
  InitialState = State;

  FermionOnSphereWithSpinAllSz* Space = 0;
  Space = new FermionOnSphereWithSpinAllSz(NbrParticles, TotalLz, LzMax, MemorySpace);
  ComplexVector SpinRotatedState(Space->GetHilbertSpaceDimension());

  if (IsU1)
  {
    FermionOnSphere* U1Space = 0;
    U1Space = new FermionOnSphere(NbrParticles, TotalLz, LzMax);
    if(InitialState.GetVectorDimension() != U1Space->GetHilbertSpaceDimension())
    {
      cout << "Error: Given U1 state has the wrong dimension. It is " <<  InitialState.GetVectorDimension() << ", while it should be " << U1Space->GetHilbertSpaceDimension() <<  "." << endl; 
      exit(1);
    }
    InitialState = Space->U1ToSU2AllSz(InitialState, *U1Space);
    delete U1Space;
  }
  else 
  {
    if (IsSU2)
    {
      FermionOnSphereWithSpin* SU2Space = 0;
      SU2Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz);
      if(InitialState.GetVectorDimension() != SU2Space->GetHilbertSpaceDimension())
	{
	  cout << "Error: Given SU2 state has the wrong dimension. It is " <<  InitialState.GetVectorDimension() << ", while it should be " << SU2Space->GetHilbertSpaceDimension() <<  "." << endl; 
	  exit(1);
	}
      InitialState = Space->SU2ToSU2AllSz(InitialState, *SU2Space);
      delete SU2Space;
    }
    else
    {
      if ( InitialState.GetVectorDimension() != Space->GetHilbertSpaceDimension() )
      {
	cout << "Error: Given state has the wrong dimension. It is " <<  InitialState.GetVectorDimension() << ", while it should be " << Space->GetHilbertSpaceDimension() <<  "." << endl; 
	exit(1);
      }
    }
  }

  ComplexMatrix* RotationMatrices;
  RotationMatrices = new ComplexMatrix [LzMax+1];

  // Create rotation matrix
  ComplexMatrix RotationMatrix(2,2);
  RotationMatrix.SetMatrixElement(0,0,cos(theta/2));
  RotationMatrix.SetMatrixElement(0,1,-sin(theta/2));
  RotationMatrix.SetMatrixElement(1,0,sin(theta/2));
  RotationMatrix.SetMatrixElement(1,1,cos(theta/2));

  for( int k=0; k<=LzMax ; ++k) *(RotationMatrices+k) = RotationMatrix;

//  cout << *(RotationMatrices+LzMax);

  Space->TransformOneBodyBasis(InitialState, SpinRotatedState, RotationMatrices);
      
  RealVector OutState;
  OutState = SpinRotatedState;
  OutState.WriteVector(OutputName);

  delete[] RotationMatrices;
  delete Space;
  delete[] OldExtension;
  delete[] NewExtension;
  delete[] OutputName;

  return 0;
}

