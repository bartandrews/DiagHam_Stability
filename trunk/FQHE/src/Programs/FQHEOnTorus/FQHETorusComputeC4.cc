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
  OptionManager Manager ("FQHETorusComputeC4" , "0.01");
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
  (*SystemGroup) += new SingleStringOption ('\n',  "interaction-name", "name that should be inserted in the output file names", "dummy");
  (*SystemGroup) += new BooleanOption ('\n',  "apply-c4", "apply the C4 rotation to the input states instead of computing the C4 eigenvalue(s)");
  (*SystemGroup) += new SingleIntegerOption ('\n',  "target-ky", "force the ky sector of the rotated state, if negative use -kx mod GCD(Ne, Nphi)", -1);
  (*SystemGroup) += new BooleanOption ('\n',  "clockwise", "apply the C4 rotation clockwise");
  (*SystemGroup) += new BooleanOption ('\n',  "export-transformation", "export the transformation matrix in a ascii file (one per momentum sector)");
  (*SystemGroup) += new BooleanOption ('\n',  "export-bintransformation", "export the transformation matrix in a binary file (one per momentum sector)");
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusComputeC4 -h" << endl;
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
						      NbrParticles, MaxMomentum, XMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("ground-state") << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Kx=" << XMomentum << " " << " Ky=" << YMomentum << endl;
      if ((Manager.GetBoolean("apply-c4") == false) && ((XMomentum != YMomentum) || (!((XMomentum == 0) || (((NbrParticles & 1) == 0) && (XMomentum == (NbrParticles / 2)))))))
	{
	  cout << "C4 symmetry can only be computed in the (0,0) or (pi,pi) sectors" << endl;
	  return -1;
	}
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
						      NbrParticles, MaxMomentum, XMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, 0) << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Kx=" << XMomentum << " " << " Ky=" << YMomentum << endl;
      MomentumModulo = FindGCD(NbrParticles, MaxMomentum);
      if ((Manager.GetBoolean("apply-c4") == false) && ((XMomentum != (YMomentum % MomentumModulo)) || (!((XMomentum == 0) || (((NbrParticles & 1) == 0) && (XMomentum == (NbrParticles / 2)))))))
	{
	  cout << "C4 symmetry can only be computed in the (0,0) or (pi,pi) sectors" << endl;
	  return -1;
	}
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
	  if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, i),
							  TmpNbrParticles, TmpMaxMomentum, TmpXMomentum, YMomenta[i], TmpStatistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, i) << endl;
	      return -1;
	    }
	  if ((TmpNbrParticles != NbrParticles) || (TmpMaxMomentum != MaxMomentum) || 
	      ((YMomenta[i] % MomentumModulo) != (YMomentum % MomentumModulo)) || (TmpXMomentum != XMomentum) || (Statistics != TmpStatistics))
	    {
	      cout << "error, " << DegeneratedFile(0, i) << " has different system parameters than " << DegeneratedFile(0, 0) 
		   << ", N=" << TmpNbrParticles << "(" << NbrParticles << "), N_phi=" << TmpMaxMomentum << "(" << MaxMomentum 
		   << "), Kx=" << TmpXMomentum << "(" << XMomentum << "), Ky=" << YMomenta[i] << " mod " << MomentumModulo << " (" << YMomentum << " mod " << MomentumModulo << ")" << endl;
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
  int* UsedYMomenta = new int [MaxMomentum];
  for (int i = 0; i < MaxMomentum; ++i)
    UsedYMomenta[i] = 0;
  for (int i = 0; i < NbrInputStates; ++i)
    UsedYMomenta[YMomenta[i]]++;
  int MomentumStep = MaxMomentum / MomentumModulo;
  for (int i = 0; i < MomentumModulo; ++i)
    {
      for (int j = 1; j < MomentumStep; ++j)
	{
	  if (UsedYMomenta[i + j * MomentumModulo] != UsedYMomenta[i])
	    {
	      cout << "Warning, not providing the same number of states in sector " << i << " and " << (i + j * MomentumModulo) << ", C4 eigenvalue calculation might be wrong" << endl;
	    }
	}
    }
  bool PiPiSectorFlag = false;
  if (XMomentum != 0)
    PiPiSectorFlag = true;
  ParticleOnTorusWithMagneticTranslations** Spaces = new ParticleOnTorusWithMagneticTranslations*[MaxMomentum];
  ParticleOnTorusWithMagneticTranslations** OutputSpaces = new ParticleOnTorusWithMagneticTranslations*[MaxMomentum];
  int TargetXMomentum = (MaxMomentum - YMomentum) % MomentumModulo;
  int TargetYMomentum  = XMomentum;
  if (Manager.GetBoolean("clockwise"))
    {
      TargetXMomentum = YMomentum & MomentumModulo;
      TargetYMomentum = (MaxMomentum - XMomentum) % MomentumModulo;
    }
  if (Manager.GetInteger("target-ky") >= 0)
    {
      TargetYMomentum = Manager.GetInteger("target-ky") % MaxMomentum;
    }
  for (int i = 0; i < MaxMomentum; ++i)
    {
      int LocalTargetYMomentum = TargetYMomentum + (i / MomentumModulo) * MomentumModulo;
      if (UsedYMomenta[i] > 0)
	{
	  if (Statistics == false)
	    {
	      Spaces[i] = new BosonOnTorusWithMagneticTranslationsShort(NbrParticles, MaxMomentum, XMomentum, i);
	      OutputSpaces[LocalTargetYMomentum] = new BosonOnTorusWithMagneticTranslationsShort(NbrParticles, MaxMomentum, TargetXMomentum, LocalTargetYMomentum);	
	      sprintf (OutputNamePrefix, "bosons_torus_%s_n_%d_2s_%d", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum);
	    }
	  else
	    {
	      Spaces[i] = new FermionOnTorusWithMagneticTranslations(NbrParticles, MaxMomentum, XMomentum, i);
	      OutputSpaces[LocalTargetYMomentum] = new FermionOnTorusWithMagneticTranslations(NbrParticles, MaxMomentum, TargetXMomentum, LocalTargetYMomentum);
	      sprintf (OutputNamePrefix, "fermions_torus_%s_n_%d_2s_%d", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum);
	    }
	  Architecture.GetArchitecture()->SetDimension(Spaces[i]->GetHilbertSpaceDimension());
	}
    }

  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (InputStates[i].GetVectorDimension() != Spaces[YMomenta[i]]->GetHilbertSpaceDimension())
	{
	  cout << "error: vector and Hilbert-space have unequal dimensions " << InputStates[i].GetVectorDimension() 
	       << " "<< Spaces[YMomenta[i]]->GetHilbertSpaceDimension() << endl;
	  return -1;
	}
    }


  char* OutputName = new char [256 + strlen(OutputNamePrefix)];
  sprintf (OutputName, "%s_c4_kx_%d_ky_%d.dat", OutputNamePrefix, YMomentum, XMomentum);
  ComplexMatrix C4Rep (NbrInputStates, NbrInputStates, true);
  if (Manager.GetBoolean("apply-c4"))
    {
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  ComplexVector TmpVector = OutputSpaces[TargetYMomentum]->C4Rotation(InputStates[i], Spaces[YMomenta[i]], Manager.GetBoolean("clockwise"), Architecture.GetArchitecture());
	  double TmpNorm = TmpVector.Norm();
	  cout << "Norm of rotated state " << i << " : " << TmpNorm << endl;
	  TmpVector /= TmpNorm;
	  char* VectorOutputName = new char [256 + strlen(OutputNamePrefix)];
	  sprintf (VectorOutputName, "%s_c4_kx_%d_ky_%d.%d.vec", OutputNamePrefix, TargetXMomentum, TargetYMomentum, i);
	  if (TmpVector.WriteVector(VectorOutputName) == false)
	    {
	      cout << "error, can't write vector " << VectorOutputName << endl;	      
	    }
	  delete[] VectorOutputName;
	}
      return 0;
    }
  else
    {
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  for (int j = 0; j < NbrInputStates; ++j)
	    {
	      ComplexVector TmpVector = OutputSpaces[YMomenta[j]]->C4Rotation(InputStates[i], Spaces[YMomenta[i]], false, Architecture.GetArchitecture());
	      double TmpNorm = TmpVector.Norm();
	      C4Rep.SetMatrixElement(j, i, (TmpVector * InputStates[j]));
	    }
	}
    }

  ofstream File;
  File.open(OutputName, ios::binary | ios::out);
  File.precision(14);
  File << "# eigenvalue Norm Arg (2 * Arg/pi) round(C4)" << endl;

  ComplexDiagonalMatrix Eigenvalues(NbrInputStates, true);
#ifdef __LAPACK__
  C4Rep.LapackDiagonalize(Eigenvalues);
#endif
  double Factor = 2.0 /  M_PI;
  for (int i = 0; i < NbrInputStates; ++i)
    {
      Complex TmpValue = Eigenvalues[i];
      File << TmpValue << " " << Norm(TmpValue) << " " << Arg(TmpValue);
      double TmpValue2 = Arg(TmpValue);
      if (fabs(TmpValue2) < MACHINE_PRECISION)
	TmpValue2 = 0.0;
      if (TmpValue2 < 0.0)
	TmpValue2 += 2.0 * M_PI;
      if (TmpValue2 >= (2.0 * M_PI))
	TmpValue2 -= 2.0 * M_PI;
      TmpValue2 *= Factor;
      File << " " << TmpValue2 << " " << round(TmpValue2) << endl;
      cout << "C4 = " << round(TmpValue2) << " (" << TmpValue2 << ")" << endl;
    }
  File.close();
  return 0;

}
