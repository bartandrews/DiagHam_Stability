#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"

#include "Operator/ParticleOnTorusWithMagneticTranslationsDensityDensityOperator.h"
#include "Operator/ParticleOnTorusWithMagneticTranslationsDensityOperator.h"
#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundState.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithEigenstates.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundStateFastDisk.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

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
  OptionManager Manager ("FQHETorusStructureFactor" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  Manager += ToolsGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "vector file that corresponds to the input state");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusStructureFactor -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles;
  int NbrFluxQuanta;
  int KxMomentum;
  int KyMomentum;
  double Ratio = 1.0;
  bool Statistics = true;
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						  NbrParticles, NbrFluxQuanta, KxMomentum, KyMomentum, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
      return -1;
    }
  cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << NbrFluxQuanta << " Ky=" << KyMomentum << " ";
  if (Statistics == false)
    {
      cout << "bosons";
    }
  else
    {
      cout << "fermions";
    }
  cout << endl;

 
  ParticleOnTorusWithMagneticTranslations* SpaceWithTranslations = 0;
  if (Statistics == false)
    {
       SpaceWithTranslations = new BosonOnTorusWithMagneticTranslationsShort(NbrParticles, NbrFluxQuanta, KxMomentum, KyMomentum);
    }
  else
    {
      SpaceWithTranslations = new FermionOnTorusWithMagneticTranslations (NbrParticles, NbrFluxQuanta, KxMomentum, KyMomentum);
    }
//   ParticleOnTorus* Space = 0;
//   if (Statistics == false)
//     {
//       Space = new BosonOnTorusShort(NbrParticles, NbrFluxQuanta, KyMomentum);
//     }
//   else
//     {
//       Space = new FermionOnTorus (NbrParticles, NbrFluxQuanta, KyMomentum);
//     }
  
  Architecture.GetArchitecture()->SetDimension(SpaceWithTranslations->GetHilbertSpaceDimension());

  char* StateFileName = Manager.GetString("input-state");
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }

  ComplexVector InputStateWithTranslations;
  if (InputStateWithTranslations.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }
  if (InputStateWithTranslations.GetVectorDimension() != SpaceWithTranslations->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputStateWithTranslations.GetVectorDimension() 
	   << " " << SpaceWithTranslations->GetHilbertSpaceDimension() << endl;
      return -1;
    }

  //  ComplexVector InputState = SpaceWithTranslations->ConvertFromKxKyBasis(InputStateWithTranslations, Space);  
 
  Complex**** DensityDensityCoefficients = new Complex***[NbrFluxQuanta];
  double Sign = 1.0;
  if (Statistics == true)
    {
      Sign = -1.0;
    }
  for (int m1 = 0; m1 < NbrFluxQuanta; ++m1)
    {
      DensityDensityCoefficients[m1] = new Complex**[NbrFluxQuanta];
      for (int m2 = 0; m2 < NbrFluxQuanta; ++m2)
	{
	  DensityDensityCoefficients[m1][m2] = new Complex*[NbrFluxQuanta];
	  for (int n1 = 0; n1 < NbrFluxQuanta; ++n1)
	    {
	      DensityDensityCoefficients[m1][m2][n1] = new Complex[NbrFluxQuanta];
	    }
	}
    }

  for (int m1 = 0; m1 < NbrFluxQuanta; ++m1)
    {
      for (int m2 = m1; m2 < NbrFluxQuanta; ++m2)
	{
	  for (int n1 = 0; n1 < NbrFluxQuanta; ++n1)
	    {
	      for (int n2 = n1; n2 < NbrFluxQuanta; ++n2)
		{
		  if (((n2 + n1) % NbrFluxQuanta) == ((m2 + m1) % NbrFluxQuanta))
		    {
		      if ((Statistics == true) && ((n1 == n2) || (m1 == m2)))
			{
			  DensityDensityCoefficients[m1][m2][n1][n2] = 0.0;
			  DensityDensityCoefficients[m2][m1][n1][n2] = 0.0;
			  DensityDensityCoefficients[m1][m2][n2][n1] = 0.0;
			  DensityDensityCoefficients[m2][m1][n2][n1] = 0.0;
			}
		      else
			{
			  //			  ParticleOnSphereDensityDensityOperator TmpOperator(Space, m1, m2, n1, n2);
			  ParticleOnTorusWithMagneticTranslationsDensityDensityOperator TmpOperator(SpaceWithTranslations, m1, m2, n1, n2);
			  OperatorMatrixElementOperation Operation (&TmpOperator, InputStateWithTranslations, InputStateWithTranslations);
			  Operation.ApplyOperation(Architecture.GetArchitecture());
			  DensityDensityCoefficients[m1][m2][n1][n2] = Operation.GetScalar();
			  DensityDensityCoefficients[m2][m1][n1][n2] = Sign * Operation.GetScalar();
			  DensityDensityCoefficients[m1][m2][n2][n1] = Sign * Operation.GetScalar();
			  DensityDensityCoefficients[m2][m1][n2][n1] = Operation.GetScalar();
			}
		    }
		}
	    }
	}
    }
  
  Complex* DensityCoefficients = new Complex[NbrFluxQuanta];
  for (int m1 = 0; m1 < NbrFluxQuanta; ++m1)
    {
//      ParticleOnSphereDensityOperator TmpOperator(Space, m1);
//       OperatorMatrixElementOperation Operation (&TmpOperator, InputState, InputState);
      ParticleOnTorusWithMagneticTranslationsDensityOperator TmpOperator(SpaceWithTranslations, m1, m1);
      OperatorMatrixElementOperation Operation (&TmpOperator, InputStateWithTranslations, InputStateWithTranslations);
      Operation.ApplyOperation(Architecture.GetArchitecture());
      DensityCoefficients[m1] = Operation.GetScalar();
    }

  double InvRatio = 1.0 / Ratio;
  double Lx = sqrt(2.0 * M_PI * ((double) NbrFluxQuanta) * Ratio);
  double Ly = sqrt(2.0 * M_PI * ((double) NbrFluxQuanta) * InvRatio);
  double InvNPhi = 2.0 * M_PI / ((double) NbrFluxQuanta);

  char* OutputName = ReplaceExtensionToFileName(Manager.GetString("input-state"), "vec", "sf");
  ofstream File;
  File.open(OutputName, ios::binary | ios::out);
  File.precision(14);
  File << "# qx qy qx_full qy_full q^2 SF.Re SF.Im" << endl;

  Complex Sum = 0.0;
  int MaxQ = NbrFluxQuanta / 2;
  for (int Qx = -MaxQ; Qx <= MaxQ; ++Qx)
    {
      double QxFull = 2.0 * M_PI * ((double) Qx) / Lx;
      Complex AverageDensitySquare = 0.0;
      for (int ky1 = 0; ky1 < NbrFluxQuanta; ++ky1)
 	{
 	  AverageDensitySquare += DensityCoefficients[ky1] * Phase (-InvNPhi *  ((double) (ky1 * Qx)));
 	}      
      AverageDensitySquare = SqrNorm(AverageDensitySquare);
      for (int Qy = -MaxQ; Qy <= MaxQ; ++Qy)
	{
	  
	  double QyFull = 2.0 * M_PI * ((double) Qy) / Ly;
	  double Q2 = (QxFull * QxFull) + (QyFull * QyFull);
	  double GaussianFactor = exp(-0.5 * Q2);
	  Complex StructureFactor = 0.0;
	  for (int ky1 = 0; ky1 < NbrFluxQuanta; ++ky1)
	    {
	      int m1 = ky1 + Qy;
	      while (m1 < 0)
		m1 += NbrFluxQuanta;
	      while (m1 >= NbrFluxQuanta)
		m1 -= NbrFluxQuanta;
	      for (int ky2 = 0; ky2 < NbrFluxQuanta; ++ky2)
		{
		  int m2 = ky2 - Qy;
		  while (m2 < 0)
		    m2 += NbrFluxQuanta;
		  while (m2 >= NbrFluxQuanta)
		    m2 -= NbrFluxQuanta;
		  StructureFactor += DensityDensityCoefficients[m1][m2][ky1][ky2] * Phase (-InvNPhi *  ((double) ((ky1 - ky2 + Qy) * Qx)));
		}
	    }
	  StructureFactor -= (double) NbrParticles;
	  //	  StructureFactor *= GaussianFactor;
	  Sum += StructureFactor * Q2 * GaussianFactor;
	  if ((Qy == 0) && ((abs(Qx % SpaceWithTranslations->GetMaxXMomentum())) == 0))
	    {
	      StructureFactor += AverageDensitySquare;
	    }
	  StructureFactor *= -1.0 / ((double) NbrFluxQuanta);
	  File << Qx << " " << Qy << " " << QxFull << " " << QyFull << " " << Q2 << " " << StructureFactor.Re << " " << StructureFactor.Im << endl;
	}
      File << endl;
    }
  File.close();
  cout << "checksum = " << Sum << endl;

  for (int m1 = 0; m1 < NbrFluxQuanta; ++m1)
    {
       for (int m2 = m1; m2 < NbrFluxQuanta; ++m2)
	{
	  for (int n1 = 0; n1 < NbrFluxQuanta; ++n1)
	    {
	      delete[] DensityDensityCoefficients[m1][m2][n1];
	    }
	  delete[] DensityDensityCoefficients[m1][m2];
	}
       delete[] DensityDensityCoefficients[m1];
    }
  delete[] DensityDensityCoefficients;

  return 0;
}
