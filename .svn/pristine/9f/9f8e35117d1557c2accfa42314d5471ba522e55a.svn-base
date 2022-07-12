#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Hamiltonian/ParticleOnTorusCoulombHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusGenericHamiltonian.h"

#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/FQHEOnTorusMainTask.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "QuantumNumber/AbstractQuantumNumber.h"

#include "GeneralTools/ListIterator.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

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

  OptionManager Manager ("FQHETorusBosonsTwoBodyGeneric" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption ('l', "max-momentum", "maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption ('y', "ky-momentum", "constraint on the total momentum modulo the maximum momentum (negative if none)", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  BooleanOption  ('\n', "redundant-kymomenta", "Calculate all subspaces up to Ky  = MaxMomentum-1", false);
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the two torus lengths", 1.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");


  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusBosonsTwoBodyGeneric -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrParticles = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  int Momentum = Manager.GetInteger("ky-momentum");
  double XRatio = Manager.GetDouble("ratio");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  bool FirstRun = true;
  
  double* PseudoPotentials;
  int NbrPseudoPotentials = 0;
  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      if (FQHETorusGetPseudopotentials(Manager.GetString("interaction-file"), NbrPseudoPotentials, PseudoPotentials) == false)
	return -1;
    }
  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "bosons_torus_kysym_%s_n_%d_2s_%d_ratio_%f.dat", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum, XRatio);
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

      ParticleOnTorus* Space = 0;
#ifdef  __64_BITS__
      if ((MaxMomentum + NbrParticles - 1) < 63)
#else
	if ((MaxMomentum + NbrParticles - 1) < 31)	
#endif
	  {
	    Space = new BosonOnTorusShort(NbrParticles, MaxMomentum, Momentum);	    
	  }
	else
	  {
	    Space = new BosonOnTorus(NbrParticles, MaxMomentum, Momentum);
	  }
      cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;

      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();

      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnTorusGenericHamiltonian (Space, NbrParticles, MaxMomentum, XRatio, NbrPseudoPotentials, PseudoPotentials,
										   Architecture.GetArchitecture(), Memory);

      if (Manager.GetString("energy-expectation") != 0 )
	{
          cout<<"Warning: assumes input vector is complex, while Hamiltonian is real."<<endl;

	  char* StateFileName = Manager.GetString("energy-expectation");
	  if (IsFile(StateFileName) == false)
	    {
	      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
	      return -1;           
	    }
	  ComplexVector InputState;
 	  if (InputState.ReadVector(StateFileName) == false)
	    {
	      cout << "error while reading " << StateFileName << endl;
	      return -1;
	    }
          RealVector InputStateRe (InputState.GetVectorDimension(), true);
          RealVector InputStateIm (InputState.GetVectorDimension(), true);
          for (int i = 0; i < InputState.GetVectorDimension(); i++)
            {
                Complex Tmp = InputState[i];
                InputStateRe[i] = Tmp.Re;
                InputStateIm[i] = Tmp.Im;
            }

	  if (InputState.GetVectorDimension()!=Space->GetHilbertSpaceDimension())
	    {
	      cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
	      return -1;
	    }
	  RealVector TmpStateRe(Space->GetHilbertSpaceDimension(), true);
	  RealVector TmpStateIm(Space->GetHilbertSpaceDimension(), true);

	  VectorHamiltonianMultiplyOperation OperationRe (Hamiltonian, &InputStateRe, &TmpStateRe);
	  OperationRe.ApplyOperation(Architecture.GetArchitecture());
	  VectorHamiltonianMultiplyOperation OperationIm (Hamiltonian, &InputStateIm, &TmpStateIm);
	  OperationIm.ApplyOperation(Architecture.GetArchitecture());

	  Complex EnergyValue = InputStateRe * TmpStateRe + InputStateIm * TmpStateIm;
          cout << "<Energy>= " << EnergyValue.Re << " " << EnergyValue.Im << endl;
	  return 0;
	}


      double Shift = -10.0;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [256];
	  sprintf (EigenvectorName, "bosons_torus_kysym_%s_n_%d_2s_%d_ratio_%f_ky_%d", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum, XRatio, Momentum);
	}
      FQHEOnTorusMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, Momentum, Shift, OutputNameLz, FirstRun, EigenvectorName);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      if (FirstRun == true)
	FirstRun = false;
      delete Hamiltonian;
      delete Space;
    }
  File.close();

  return 0;
}
