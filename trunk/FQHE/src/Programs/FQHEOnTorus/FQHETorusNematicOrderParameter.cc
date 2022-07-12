#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/IntegerAlgebraTools.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Hamiltonian/ParticleOnTorusNematicParameterWithMagneticTranslationsHamiltonian.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHETorusNematicOrderParameter" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption('i', "input-state", "name of the file containing the state whose Kx momentum has to be computed");
  (*SystemGroup) += new SingleStringOption('\n', "degenerated-states", "name of the file containing a list of states (override input-state)");
  (*SystemGroup) += new SingleStringOption ('\n',  "interaction-name", "name that should be inserted in the output file names", "nematic_parameter");
  (*SystemGroup) += new SingleDoubleOption ('\n',  "length", "characteristic length (in magnetic lentgh unit) used in the order parameter", 1.0);
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusNematicOrderParameter -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int MaxMomentum = 0;
  int YMomentum = -1;
  int MomentumModulo = -1;
  int* YMomenta = 0;
  bool Statistics = true;
  int XMomentum = -1;
  double Ratio = 0.0;
  ComplexVector* InputStates = 0;
  int NbrInputStates = 0;
  char* OutputNamePrefix = new char [256 + strlen(Manager.GetString("interaction-name"))];
  if (Manager.GetString("degenerated-states") == 0)
    {
      if (Manager.GetString("input-state") == 0)
	{
	  cout << "error, either input-state or degenerated-states has to be provided" << endl;
	  return -1;
	}
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						      NbrParticles, MaxMomentum, XMomentum, YMomentum, Ratio, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("ground-state") << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Kx=" << XMomentum << " " << " Ky=" << YMomentum << endl;
      NbrInputStates = 1;
      InputStates = new ComplexVector [NbrInputStates];
      YMomenta = new int [NbrInputStates];
      YMomenta[0] = YMomentum;
      MomentumModulo = FindGCD(NbrParticles, MaxMomentum);
      if (InputStates[0].ReadVector(Manager.GetString("input-state")) == false)
	{
	  cout << "error while reading " << Manager.GetString("input-state") << endl;
	  return -1;
	}
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-states")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	} 
      NbrInputStates = DegeneratedFile.GetNbrLines();
      if (NbrInputStates < 1)
	{
	  cout << "no state found in " << Manager.GetString("degenerated-states") << endl;
	}
      InputStates = new ComplexVector [NbrInputStates];
      YMomenta = new int [NbrInputStates];
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, 0),
						      NbrParticles, MaxMomentum, XMomentum, YMomentum, Ratio, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, 0) << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Kx=" << XMomentum << " " << " Ky=" << YMomentum << endl;
      MomentumModulo = FindGCD(NbrParticles, MaxMomentum);
      if (InputStates[0].ReadVector(DegeneratedFile(0, 0)) == false)
	{
	  cout << "error while reading " << DegeneratedFile(0, 0) << endl;
	  return -1;
	}
      YMomenta[0] = YMomentum;
      for (int i = 1; i < NbrInputStates; ++i)
	{
	  int TmpNbrParticles = 0;
	  int TmpMaxMomentum = 0;
	  YMomenta[i] = -1;
	  int TmpXMomentum = -1;
	  bool TmpStatistics = true;
	  double TmpRatio = 0.0;
	  if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, i),
							  TmpNbrParticles, TmpMaxMomentum, TmpXMomentum, YMomenta[i], TmpRatio, TmpStatistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, i) << endl;
	      return -1;
	    }
	  if ((TmpNbrParticles != NbrParticles) || (TmpMaxMomentum != MaxMomentum) || 
	      (YMomenta[i] != YMomentum) || (TmpXMomentum != XMomentum) || (Statistics != TmpStatistics))
	    {
	      cout << "error, " << DegeneratedFile(0, i) << " has different system parameters than " << DegeneratedFile(0, 0) 
		   << ", N=" << TmpNbrParticles << "(" << NbrParticles << "), N_phi=" << TmpMaxMomentum << "(" << MaxMomentum 
		   << "), Kx=" << TmpXMomentum << "(" << XMomentum << "), Ky=" << YMomenta[i] << " (" << YMomentum << ")" << endl;
	    }
	  if (InputStates[i].ReadVector(DegeneratedFile(0, i)) == false)
	    {
	      cout << "error while reading " << DegeneratedFile(0, i) << endl;
	      return -1;
	    }
	  if (InputStates[i].GetVectorDimension() != InputStates[0].GetVectorDimension())
	    {
	      cout << "error, " << DegeneratedFile(0, i) << " has different dimension than " << DegeneratedFile(0, 0) << endl;	      
	    }
	}
    }
  int MomentumStep = MaxMomentum / MomentumModulo;
  
  ParticleOnTorusWithMagneticTranslations* Space = 0;
  if (Statistics == false)
    {
      Space = new BosonOnTorusWithMagneticTranslationsShort(NbrParticles, MaxMomentum, XMomentum, YMomentum);
      sprintf (OutputNamePrefix, "bosons_torus_%s_l_%.6f_n_%d_2s_%d", Manager.GetString("interaction-name"), Manager.GetDouble("length"), NbrParticles, MaxMomentum);
    }
  else
    {
      Space = new FermionOnTorusWithMagneticTranslations(NbrParticles, MaxMomentum, XMomentum, YMomentum);
      sprintf (OutputNamePrefix, "fermions_torus_%s_l_%.6f_n_%d_2s_%d", Manager.GetString("interaction-name"), Manager.GetDouble("length"), NbrParticles, MaxMomentum);
    }

  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (InputStates[i].GetVectorDimension() != Space->GetHilbertSpaceDimension())
	{
	  cout << "error: vector and Hilbert-space have unequal dimensions " << InputStates[i].GetVectorDimension() 
	       << " "<< Space->GetHilbertSpaceDimension() << endl;
	  return -1;
	}
    }


  char* OutputName = new char [256 + strlen(OutputNamePrefix)];
  sprintf (OutputName, "%s_ratio_%.6f_kx_%d_ky_%d.dat", OutputNamePrefix, Ratio, XMomentum, YMomentum);

  ParticleOnTorusNematicParameterWithMagneticTranslationsHamiltonian NematicParameter(Space, NbrParticles, MaxMomentum, XMomentum, Ratio, Manager.GetDouble("length"), Architecture.GetArchitecture(), 0l);
//   HermitianMatrix NematicParameterRep (Space->GetHilbertSpaceDimension(),true);
//   NematicParameter.GetHamiltonian(NematicParameterRep);
  HermitianMatrix NematicParameterRep (NbrInputStates,true);
  ComplexVector TmpVector(Space->GetHilbertSpaceDimension());
  for (int i = 0; i < NbrInputStates; ++i)
    {
      VectorHamiltonianMultiplyOperation Operation (&NematicParameter, &(InputStates[i]), &TmpVector);
      Operation.ApplyOperation(Architecture.GetArchitecture());	  
      for (int j = 0; j < NbrInputStates; ++j)
	{
	  NematicParameterRep.SetMatrixElement(j, i, (TmpVector * InputStates[j]));
	}
    }
  ofstream File;
  File.open(OutputName, ios::binary | ios::out);
  File.precision(14);
  File << "# eigenvalue Norm Arg (2 * Arg/pi) round(C4)" << endl;

  RealDiagonalMatrix Eigenvalues(NematicParameterRep.GetNbrColumn(), true);
#ifdef __LAPACK__
  NematicParameterRep.LapackDiagonalize(Eigenvalues);
#else
  cout << "error, lapack is required" << endl;
#endif
  for (int i = 0; i < NematicParameterRep.GetNbrColumn(); ++i)
    {
      double TmpValue = Eigenvalues[i];
      File << TmpValue << endl;
    }
  File.close();
  return 0;

}
