#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Operator/ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator.h"


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


// find the gauge transformation to get a smooth gauge choice between two set of vectors
//
// leftVectors = vectors for which tha gauge transformation has to be computed
// rightVectors = reference vectors that we have to smoothly evolve from
// rightMatrix = additional transformation that has to be applied to the reference vectors
// return value = gauge transformation
ComplexMatrix FindGaugeTransformation (ComplexVector* leftVectors, ComplexVector* rightVectors, int nbrVectors, ComplexMatrix& rightMatrix);


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FCIRealSpaceManyBodyChernNumber" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "name of the file at zero twist");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-gammax", "number of discretization step in the x direction", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-gammay", "number of discretization step in the y direction", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-states", "number of states corresponding to the low energy state manifold", 1);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIRealSpaceManyBodyChernNumber -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, a state file should be provided. See man page for option syntax or type FCIRealSpaceManyBodyChernNumber -h" << endl;
      return -1;
    }

  int NbrGammaX = Manager.GetInteger("nbr-gammax");
  int NbrGammaY = Manager.GetInteger("nbr-gammay");
  double GammaXStep = 1.0 / ((double) NbrGammaX);
  double GammaYStep = 1.0 / ((double) NbrGammaY);
  int NbrStates = Manager.GetInteger("nbr-states");

  double GammaX = 0.0;
  double GammaY = 0.0;
  char* TmpGammaFileName = new char[512];
  char* TmpGammaFileName0 = new char[512];
  char* TmpIndexFileName = new char[512];
  char* TmpIndexFileName0 = new char[512];
  sprintf (TmpGammaFileName0, "gx_%f_gy_%f", GammaX, GammaY);
  sprintf (TmpIndexFileName0, ".0.vec");
  cout << "checking the " << ((NbrGammaX + 1) * (NbrGammaY + 1) * NbrStates) << " vector files" << endl;
  for (int GammaXIndex = 0; GammaXIndex <= NbrGammaX; ++GammaXIndex)
    {
      GammaY = 0.0;
      for (int GammaYIndex = 0; GammaYIndex <= NbrGammaY; ++GammaYIndex)
	{
	  sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	  char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	  for (int i = 0; i < NbrStates; ++i)
	    {
	      sprintf (TmpIndexFileName, ".%d.vec", i);
	      char* TmpFileName2 = ReplaceString(TmpFileName, TmpIndexFileName0, TmpIndexFileName);
	      if (IsFile(TmpFileName2) == false)
		{
		  cout << "can't find file " << TmpFileName2 << endl;
		  return 0;		  
		}
	      delete[] TmpFileName2;
	    }
	  delete[] TmpFileName;
	  GammaY  += GammaYStep;
	} 
      GammaX  += GammaXStep;
    } 
  cout << "check done" << endl;
  int NbrParticles = 0;
  int NbrSites = 0;
  bool Statistics = true;
  bool GutzwillerFlag = false;
  int MomentumFlag = false;
  int KxMomentum = 0;
  int XPeriodicity = 0;
  int KyMomentum = 0;
  int YPeriodicity = 0;
  bool FixedSzFlag = false;
  int TotalSz = 0;

  cout << "Smoothing gauge" << endl;
  if (NbrStates == 0)
    {
      ComplexMatrix GaugeFixingPhases (NbrGammaX + 1 , NbrGammaY + 1);
      GaugeFixingPhases.SetMatrixElement(0, 0, 1.0);
      ComplexVector TmpVector1;
      ComplexVector TmpVector2;
      if (TmpVector1.ReadVector(Manager.GetString("input-state")) == false)
	{
	  cout << "can't read " << Manager.GetString("input-state") << endl;
	  return 0;
	}
      GammaX = 0.0;
      GammaY = GammaYStep;
      for (int GammaYIndex = 1; GammaYIndex <= NbrGammaY; ++GammaYIndex)
	{
	  sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	  char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	  if (TmpVector2.ReadVector(TmpFileName) == false)
	    {
	      cout << "can't read " << TmpFileName << endl;
	      return 0;
	    }
	  Complex Tmp = TmpVector2 * TmpVector1;
	  Tmp /= Norm(Tmp);
	  Complex Tmp2;
	  GaugeFixingPhases.GetMatrixElement(0, GammaYIndex - 1, Tmp2);	  
	  Tmp *= Tmp2;
	  GaugeFixingPhases.SetMatrixElement(0, GammaYIndex, Tmp);	  
	  ComplexVector TmpVector3 = TmpVector1;
	  TmpVector1 = TmpVector2;
	  TmpVector2 = TmpVector3;
	  delete[] TmpFileName;
	  GammaY  += GammaYStep;
	} 
      for (int GammaXIndex = 1; GammaXIndex <= NbrGammaX; ++GammaXIndex)
	{
	  GammaY = 0.0;
	  sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	  char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	  if (TmpVector1.ReadVector(TmpFileName) == false)
	    {
	      cout << "can't read " << TmpFileName << endl;
	      return 0;
	    }
	  GammaX  += GammaXStep;
	  sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	  TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	  if (TmpVector2.ReadVector(TmpFileName) == false)
	    {
	      cout << "can't read " << TmpFileName << endl;
	      return 0;
	    }
	  Complex Tmp = TmpVector2 * TmpVector1;
	  Tmp /= Norm(Tmp);
	  Complex Tmp2;
	  GaugeFixingPhases.GetMatrixElement(GammaXIndex - 1, 0, Tmp2);	  
	  Tmp *= Tmp2;
	  GaugeFixingPhases.SetMatrixElement(GammaXIndex, 0, Tmp);	  
	  GammaY = GammaYStep;
	  ComplexVector TmpVector3 = TmpVector1;
	  TmpVector1 = TmpVector2;
	  TmpVector2 = TmpVector3;
	  for (int GammaYIndex = 1; GammaYIndex <= NbrGammaY; ++GammaYIndex)
	    {
	      sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	      TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	      if (TmpVector2.ReadVector(TmpFileName) == false)
		{
		  cout << "can't read " << TmpFileName << endl;
		  return 0;
		}
	      Tmp = TmpVector2 * TmpVector1;
	      Tmp /= Norm(Tmp);
	      Complex Tmp2;
	      GaugeFixingPhases.GetMatrixElement(GammaXIndex, GammaYIndex - 1, Tmp2);	  
	      Tmp *= Tmp2;
	      GaugeFixingPhases.SetMatrixElement(GammaXIndex, GammaYIndex, Tmp);	  
	      TmpVector3 = TmpVector1;
	      TmpVector1 = TmpVector2;
	      TmpVector2 = TmpVector3;
	      delete[] TmpFileName;
	      GammaY  += GammaYStep;
	    } 
	} 

      cout << "computing many-body Chern number" << endl;
      double ChernNumber = 0.0;
      ComplexVector* VectorGammaX1 = new ComplexVector[NbrGammaY + 1];
      ComplexVector* VectorGammaX2 = new ComplexVector[NbrGammaY + 1];
      ComplexVector* VectorGammaX3 = new ComplexVector[NbrGammaY + 1];
      GammaX = 0.0;
      GammaY = 0.0;
      for (int GammaYIndex = 0; GammaYIndex <= NbrGammaY; ++GammaYIndex)
	{
	  sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	  char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	  if (VectorGammaX1[GammaYIndex].ReadVector(TmpFileName) == false)
	    {
	      cout << "can't read " << TmpFileName << endl;
	      return 0;
	    }
	  delete[] TmpFileName;
	  GammaY  += GammaYStep;
	}
      GammaX  += GammaXStep;
      GammaY = 0.0;
      for (int GammaYIndex = 0; GammaYIndex <= NbrGammaY; ++GammaYIndex)
	{
	  sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	  char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	  if (VectorGammaX2[GammaYIndex].ReadVector(TmpFileName) == false)
	    {
	      cout << "can't read " << TmpFileName << endl;
	      return 0;
	    }
	  delete[] TmpFileName;
	  GammaY  += GammaYStep;
	}
      GammaX  += GammaXStep;
      Complex TmpXDec;
      Complex TmpXInc;
      Complex TmpYDec;	      
      Complex TmpYInc;	      
      Complex Tmp;

      GaugeFixingPhases.GetMatrixElement(0, 0, TmpXDec);
      GaugeFixingPhases.GetMatrixElement(1, 0, TmpXInc);
      GaugeFixingPhases.GetMatrixElement(0, 0, TmpYDec);
      GaugeFixingPhases.GetMatrixElement(0, 1, TmpYInc);
      Tmp  = (VectorGammaX2[0] * VectorGammaX1[1]) * TmpYInc * Conj(TmpXInc);
      Tmp += (VectorGammaX1[0] * VectorGammaX1[0]) * TmpYDec * Conj(TmpXDec);
      Tmp -= (VectorGammaX1[0] * VectorGammaX1[1]) * TmpYInc * Conj(TmpXDec);
      Tmp -= (VectorGammaX2[0] * VectorGammaX1[0]) * TmpYDec * Conj(TmpXInc);
      ChernNumber += Tmp.Im;
      for (int GammaYIndex = 1; GammaYIndex < NbrGammaY; ++GammaYIndex)
	{
	  GaugeFixingPhases.GetMatrixElement(0, GammaYIndex, TmpXDec);
	  GaugeFixingPhases.GetMatrixElement(1, GammaYIndex, TmpXInc);
	  GaugeFixingPhases.GetMatrixElement(0, GammaYIndex - 1, TmpYDec);
	  GaugeFixingPhases.GetMatrixElement(0, GammaYIndex + 1, TmpYInc);
	  Tmp  = (VectorGammaX2[GammaYIndex] * VectorGammaX1[GammaYIndex + 1]) * TmpYInc * Conj(TmpXInc);
	  Tmp += (VectorGammaX1[GammaYIndex] * VectorGammaX1[GammaYIndex - 1]) * TmpYDec * Conj(TmpXDec);
	  Tmp -= (VectorGammaX2[GammaYIndex] * VectorGammaX1[GammaYIndex - 1]) * TmpYDec * Conj(TmpXInc);
	  Tmp -= (VectorGammaX1[GammaYIndex] * VectorGammaX1[GammaYIndex + 1]) * TmpYInc * Conj(TmpXDec);
	  ChernNumber += 0.5 * Tmp.Im;
	}
      GaugeFixingPhases.GetMatrixElement(0, NbrGammaY, TmpXDec);
      GaugeFixingPhases.GetMatrixElement(1, NbrGammaY, TmpXInc);
      GaugeFixingPhases.GetMatrixElement(0, NbrGammaY - 1 , TmpYDec);
      GaugeFixingPhases.GetMatrixElement(0, NbrGammaY , TmpYInc);
      Tmp  = (VectorGammaX2[NbrGammaY] * VectorGammaX1[NbrGammaY]) * TmpYInc * Conj(TmpXInc);
      Tmp += (VectorGammaX1[NbrGammaY] * VectorGammaX1[NbrGammaY - 1]) * TmpYDec * Conj(TmpXDec);
      Tmp -= (VectorGammaX2[NbrGammaY] * VectorGammaX1[NbrGammaY - 1]) * TmpYDec * Conj(TmpXInc);
      Tmp -= (VectorGammaX1[NbrGammaY] * VectorGammaX1[NbrGammaY]) * TmpYInc * Conj(TmpXDec);
      ChernNumber += Tmp.Im;


      for (int GammaXIndex = 1; GammaXIndex < NbrGammaX; ++GammaXIndex)
	{
	  GammaY = 0.0;
	  for (int GammaYIndex = 0; GammaYIndex <= NbrGammaY; ++GammaYIndex)
	    {
	      sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	      char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	      if (VectorGammaX3[GammaYIndex].ReadVector(TmpFileName) == false)
		{
		  cout << "can't read " << TmpFileName << endl;
		  return 0;
		}
	      delete[] TmpFileName;
	      GammaY  += GammaYStep;
	    }
	  GaugeFixingPhases.GetMatrixElement(GammaXIndex - 1, 0, TmpXDec);
	  GaugeFixingPhases.GetMatrixElement(GammaXIndex + 1, 0, TmpXInc);
	  GaugeFixingPhases.GetMatrixElement(GammaXIndex, 0, TmpYDec);
	  GaugeFixingPhases.GetMatrixElement(GammaXIndex, 1, TmpYInc);
	  Tmp  = (VectorGammaX1[0] * VectorGammaX2[0]) * TmpYDec * Conj(TmpXDec);
	  Tmp += (VectorGammaX3[0] * VectorGammaX2[1]) * TmpYInc * Conj(TmpXInc);
	  Tmp -= (VectorGammaX1[0] * VectorGammaX2[1]) * TmpYInc * Conj(TmpXDec);
	  Tmp -= (VectorGammaX3[0] * VectorGammaX2[0]) * TmpYDec * Conj(TmpXInc);
	  ChernNumber += 0.5 * Tmp.Im;
	  for (int GammaYIndex = 1; GammaYIndex < NbrGammaY; ++GammaYIndex)
	    {
	      GaugeFixingPhases.GetMatrixElement(GammaXIndex - 1, GammaYIndex, TmpXDec);
	      GaugeFixingPhases.GetMatrixElement(GammaXIndex + 1, GammaYIndex, TmpXInc);
	      GaugeFixingPhases.GetMatrixElement(GammaXIndex, GammaYIndex - 1, TmpYDec);
	      GaugeFixingPhases.GetMatrixElement(GammaXIndex, GammaYIndex + 1, TmpYInc);
	      Tmp  = (VectorGammaX1[GammaYIndex] * VectorGammaX2[GammaYIndex - 1]) * TmpYDec * Conj(TmpXDec);
	      Tmp += (VectorGammaX3[GammaYIndex] * VectorGammaX2[GammaYIndex + 1]) * TmpYInc * Conj(TmpXInc);
	      Tmp -= (VectorGammaX1[GammaYIndex] * VectorGammaX2[GammaYIndex + 1]) * TmpYInc * Conj(TmpXDec);
	      Tmp -= (VectorGammaX3[GammaYIndex] * VectorGammaX2[GammaYIndex - 1]) * TmpYDec * Conj(TmpXInc);
	      ChernNumber += 0.25 * Tmp.Im;
	    }
	  GaugeFixingPhases.GetMatrixElement(GammaXIndex - 1, NbrGammaY, TmpXDec);
	  GaugeFixingPhases.GetMatrixElement(GammaXIndex + 1, NbrGammaY, TmpXInc);
	  GaugeFixingPhases.GetMatrixElement(GammaXIndex,NbrGammaY - 1 , TmpYDec);
	  GaugeFixingPhases.GetMatrixElement(GammaXIndex,NbrGammaY , TmpYInc);
	  Tmp  = (VectorGammaX1[NbrGammaY] * VectorGammaX2[NbrGammaY - 1]) * TmpYDec * Conj(TmpXDec);
	  Tmp += (VectorGammaX3[NbrGammaY] * VectorGammaX2[NbrGammaY]) * TmpYInc * Conj(TmpXInc);
	  Tmp -= (VectorGammaX1[NbrGammaY] * VectorGammaX2[NbrGammaY]) * TmpYInc * Conj(TmpXDec);
	  Tmp -= (VectorGammaX3[NbrGammaY] * VectorGammaX2[NbrGammaY- 1 ]) * TmpYDec * Conj(TmpXInc);
	  ChernNumber += 0.5 * Tmp.Im;
	  ComplexVector* VectorGammaX4 = VectorGammaX1;
	  VectorGammaX1 = VectorGammaX2;
	  VectorGammaX2 = VectorGammaX3;
	  VectorGammaX3 = VectorGammaX4;
	  GammaX  += GammaXStep;      
	}
      
      GaugeFixingPhases.GetMatrixElement(NbrGammaX - 1, 0, TmpXDec);
      GaugeFixingPhases.GetMatrixElement(NbrGammaX, 0, TmpXInc);
      GaugeFixingPhases.GetMatrixElement(NbrGammaX, 0, TmpYDec);
      GaugeFixingPhases.GetMatrixElement(NbrGammaX, 1, TmpYInc);
      Tmp  = (VectorGammaX1[0] * VectorGammaX2[1]) * TmpYInc * Conj(TmpXDec);
      Tmp += (VectorGammaX2[0] * VectorGammaX2[0]) * TmpYDec * Conj(TmpXInc);
      Tmp -= (VectorGammaX1[0] * VectorGammaX2[0]) * TmpYDec * Conj(TmpXDec);
      Tmp -= (VectorGammaX2[0] * VectorGammaX2[1]) * TmpYInc * Conj(TmpXInc);
      ChernNumber += Tmp.Im;
      for (int GammaYIndex = 1; GammaYIndex < NbrGammaY; ++GammaYIndex)
	{
	  GaugeFixingPhases.GetMatrixElement(NbrGammaX - 1, GammaYIndex, TmpXDec);
	  GaugeFixingPhases.GetMatrixElement(NbrGammaX, GammaYIndex, TmpXInc);
	  GaugeFixingPhases.GetMatrixElement(NbrGammaX, GammaYIndex - 1, TmpYDec);
	  GaugeFixingPhases.GetMatrixElement(NbrGammaX, GammaYIndex + 1, TmpYInc);
	  Tmp  = (VectorGammaX1[GammaYIndex] * VectorGammaX2[GammaYIndex + 1]) * TmpYInc * Conj(TmpXDec);
	  Tmp += (VectorGammaX2[GammaYIndex] * VectorGammaX2[GammaYIndex - 1]) * TmpYDec * Conj(TmpXInc);
	  Tmp -= (VectorGammaX1[GammaYIndex] * VectorGammaX2[GammaYIndex - 1]) * TmpYDec * Conj(TmpXDec);
	  Tmp -= (VectorGammaX2[GammaYIndex] * VectorGammaX2[GammaYIndex + 1]) * TmpYInc * Conj(TmpXInc);
	  ChernNumber += 0.5 * Tmp.Im;
	}
      GaugeFixingPhases.GetMatrixElement(NbrGammaX - 1, NbrGammaY, TmpXDec);
      GaugeFixingPhases.GetMatrixElement(NbrGammaX, NbrGammaY, TmpXInc);
      GaugeFixingPhases.GetMatrixElement(NbrGammaX, NbrGammaY - 1 , TmpYDec);
      GaugeFixingPhases.GetMatrixElement(NbrGammaX, NbrGammaY , TmpYInc);
      Tmp  = (VectorGammaX1[NbrGammaY] * VectorGammaX2[NbrGammaY]) * TmpYInc * Conj(TmpXDec);
      Tmp += (VectorGammaX2[NbrGammaY] * VectorGammaX2[NbrGammaY - 1]) * TmpYDec * Conj(TmpXInc);
      Tmp -= (VectorGammaX1[NbrGammaY] * VectorGammaX2[NbrGammaY - 1]) * TmpYDec * Conj(TmpXDec);
      Tmp -= (VectorGammaX2[NbrGammaY] * VectorGammaX2[NbrGammaY]) * TmpYInc * Conj(TmpXInc);
      ChernNumber += Tmp.Im;
            
      ChernNumber /= 2.0 * M_PI;
      ChernNumber *= 2.0;
      cout << "C=" << ChernNumber << endl;
    }
  else
    {
      ComplexMatrix** GaugeFixingMatrices = new ComplexMatrix* [NbrGammaX + 1];
      for (int i = 0; i <= NbrGammaX; ++i)
	GaugeFixingMatrices[i] = new ComplexMatrix[NbrGammaY + 1];
      GaugeFixingMatrices[0][0] = ComplexMatrix(NbrStates, NbrStates);
      GaugeFixingMatrices[0][0].SetToIdentity();
      ComplexVector* TmpVectors1 = new ComplexVector[NbrStates];
      ComplexVector* TmpVectors2 = new ComplexVector[NbrStates];
      for (int i = 0; i < NbrStates ; ++i)
	{
	  sprintf (TmpIndexFileName, ".%d.vec", i);
	  char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpIndexFileName0, TmpIndexFileName);
	  if (TmpVectors1[i].ReadVector(TmpFileName) == false)
	    {
	      cout << "can't read " << TmpFileName << endl;
	      return 0;
	    }
	  delete[] TmpFileName;
	}
      GammaX = 0.0;
      GammaY = GammaYStep;
      for (int GammaYIndex = 1; GammaYIndex <= NbrGammaY; ++GammaYIndex)
	{
	  sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	  char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	  for (int i = 0; i < NbrStates ; ++i)
	    {
	      sprintf (TmpIndexFileName, ".%d.vec", i);
	      char* TmpFileName2 = ReplaceString(TmpFileName, TmpIndexFileName0, TmpIndexFileName);
	      if (TmpVectors2[i].ReadVector(TmpFileName2) == false)
		{
		  cout << "can't read " << TmpFileName2 << endl;
		  return 0;
		}
	      delete[] TmpFileName2;
	    }
	  GaugeFixingMatrices[0][GammaYIndex] = FindGaugeTransformation(TmpVectors2, TmpVectors1, NbrStates, GaugeFixingMatrices[0][GammaYIndex - 1]);
	  ComplexVector* TmpVectors3 = TmpVectors1;
	  TmpVectors1 = TmpVectors2;
	  TmpVectors2 = TmpVectors3;
	  delete[] TmpFileName;
	  GammaY  += GammaYStep;
	} 
      for (int GammaXIndex = 1; GammaXIndex <= NbrGammaX; ++GammaXIndex)
	{
	  GammaY = 0.0;
	  sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	  char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	  for (int i = 0; i < NbrStates ; ++i)
	    {
	      sprintf (TmpIndexFileName, ".%d.vec", i);
	      char* TmpFileName2 = ReplaceString(TmpFileName, TmpIndexFileName0, TmpIndexFileName);
	      if (TmpVectors1[i].ReadVector(TmpFileName2) == false)
		{
		  cout << "can't read " << TmpFileName2 << endl;
		  return 0;
		}
	      delete[] TmpFileName2;
	    }
	  GammaX  += GammaXStep;
	  sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	  TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	  for (int i = 0; i < NbrStates ; ++i)
	    {
	      sprintf (TmpIndexFileName, ".%d.vec", i);
	      char* TmpFileName2 = ReplaceString(TmpFileName, TmpIndexFileName0, TmpIndexFileName);
	      if (TmpVectors2[i].ReadVector(TmpFileName2) == false)
		{
		  cout << "can't read " << TmpFileName2 << endl;
		  return 0;
		}
	      delete[] TmpFileName2;
	    }
	  GaugeFixingMatrices[GammaXIndex][0] = FindGaugeTransformation(TmpVectors2, TmpVectors1, NbrStates, GaugeFixingMatrices[GammaXIndex - 1][0]);
	  ComplexVector* TmpVectors3 = TmpVectors1;
	  TmpVectors1 = TmpVectors2;
	  TmpVectors2 = TmpVectors3;
	  GammaY = GammaYStep;
	  for (int GammaYIndex = 1; GammaYIndex <= NbrGammaY; ++GammaYIndex)
	    {
	      sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	      TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	      for (int i = 0; i < NbrStates ; ++i)
		{
		  sprintf (TmpIndexFileName, ".%d.vec", i);
		  char* TmpFileName2 = ReplaceString(TmpFileName, TmpIndexFileName0, TmpIndexFileName);
		  if (TmpVectors2[i].ReadVector(TmpFileName2) == false)
		    {
		      cout << "can't read " << TmpFileName2 << endl;
		      return 0;
		    }
		  delete[] TmpFileName2;
		}
	      GaugeFixingMatrices[GammaXIndex][GammaYIndex] = FindGaugeTransformation(TmpVectors2, TmpVectors1, NbrStates, GaugeFixingMatrices[GammaXIndex][GammaYIndex - 1]);
	      TmpVectors3 = TmpVectors1;
	      TmpVectors1 = TmpVectors2;
	      TmpVectors2 = TmpVectors3;
	      delete[] TmpFileName;
	      GammaY  += GammaYStep;
	    } 
	} 


      cout << "computing many-body Chern number" << endl;
      double ChernNumber = 0.0;
      ComplexMatrix* VectorGammaX1 = new ComplexMatrix[NbrGammaY + 1];
      ComplexMatrix* VectorGammaX2 = new ComplexMatrix[NbrGammaY + 1];
      ComplexMatrix* VectorGammaX3 = new ComplexMatrix[NbrGammaY + 1];

      GammaX = 0.0;
      GammaY = 0.0;
      for (int GammaYIndex = 0; GammaYIndex <= NbrGammaY; ++GammaYIndex)
	{
	  sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	  char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	  ComplexVector* TmpVectors = new ComplexVector[NbrStates];
	  for (int i = 0; i < NbrStates ; ++i)
	    {
	      sprintf (TmpIndexFileName, ".%d.vec", i);
	      char* TmpFileName2 = ReplaceString(TmpFileName, TmpIndexFileName0, TmpIndexFileName);
	      if (TmpVectors[i].ReadVector(TmpFileName2) == false)
		{
		  cout << "can't read " << TmpFileName2 << endl;
		  return 0;
		}
	      delete[] TmpFileName2;
	    }
	  VectorGammaX1[GammaYIndex] = ComplexMatrix(TmpVectors, NbrStates);
	  VectorGammaX1[GammaYIndex].Multiply(GaugeFixingMatrices[0][GammaYIndex]);
	  delete[] TmpFileName;
	  GammaY  += GammaYStep;
	}
      GammaX += GammaXStep;
      GammaY = 0.0;
      for (int GammaYIndex = 0; GammaYIndex <= NbrGammaY; ++GammaYIndex)
	{
	  sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	  char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	  ComplexVector* TmpVectors = new ComplexVector[NbrStates];
	  for (int i = 0; i < NbrStates ; ++i)
	    {
	      sprintf (TmpIndexFileName, ".%d.vec", i);
	      char* TmpFileName2 = ReplaceString(TmpFileName, TmpIndexFileName0, TmpIndexFileName);
	      if (TmpVectors[i].ReadVector(TmpFileName2) == false)
		{
		  cout << "can't read " << TmpFileName2 << endl;
		  return 0;
		}
	      delete[] TmpFileName2;
	    }
	  VectorGammaX2[GammaYIndex] = ComplexMatrix(TmpVectors, NbrStates);
	  VectorGammaX2[GammaYIndex].Multiply(GaugeFixingMatrices[1][GammaYIndex]);
	  VectorGammaX3[GammaYIndex] = ComplexMatrix(TmpVectors[0].GetVectorDimension(), NbrStates);
	  delete[] TmpFileName;
	  GammaY  += GammaYStep;
	}
      GammaX  += GammaXStep;
      Complex TmpXDec;
      Complex TmpXInc;
      Complex TmpYDec;	      
      Complex TmpYInc;	      
      Complex Tmp;

       for (int i = 0; i < NbrStates ; ++i)
 	{
 	  Tmp  = (VectorGammaX2[0][i] * VectorGammaX1[1][i]);
 	  Tmp += (VectorGammaX1[0][i] * VectorGammaX1[0][i]);
 	  Tmp -= (VectorGammaX1[0][i] * VectorGammaX1[1][i]);
 	  Tmp -= (VectorGammaX2[0][i] * VectorGammaX1[0][i]);
	  ChernNumber += Tmp.Im;
 	}
      for (int GammaYIndex = 1; GammaYIndex < NbrGammaY; ++GammaYIndex)
	{
	  for (int i = 0; i < NbrStates ; ++i)
	    {
	      Tmp  = (VectorGammaX2[GammaYIndex][i] * VectorGammaX1[GammaYIndex + 1][i]);
	      Tmp += (VectorGammaX1[GammaYIndex][i] * VectorGammaX1[GammaYIndex - 1][i]);
	      Tmp -= (VectorGammaX2[GammaYIndex][i] * VectorGammaX1[GammaYIndex - 1][i]);
	      Tmp -= (VectorGammaX1[GammaYIndex][i] * VectorGammaX1[GammaYIndex + 1][i]);
	      ChernNumber += 0.5 * Tmp.Im;
	    }
	}
      for (int i = 0; i < NbrStates ; ++i)
 	{
	  Tmp  = (VectorGammaX2[NbrGammaY][i] * VectorGammaX1[NbrGammaY][i]);
	  Tmp += (VectorGammaX1[NbrGammaY][i] * VectorGammaX1[NbrGammaY - 1][i]);
	  Tmp -= (VectorGammaX2[NbrGammaY][i] * VectorGammaX1[NbrGammaY - 1][i]);
	  Tmp -= (VectorGammaX1[NbrGammaY][i] * VectorGammaX1[NbrGammaY][i]);
	  ChernNumber += Tmp.Im;
	}


      for (int GammaXIndex = 1; GammaXIndex < NbrGammaX; ++GammaXIndex)
	{
	  GammaY = 0.0;
	  for (int GammaYIndex = 0; GammaYIndex <= NbrGammaY; ++GammaYIndex)
	    {
	      sprintf (TmpGammaFileName, "gx_%f_gy_%f", GammaX, GammaY);
	      char* TmpFileName = ReplaceString(Manager.GetString("input-state"), TmpGammaFileName0, TmpGammaFileName);
	      for (int i = 0; i < NbrStates ; ++i)
		{
		  sprintf (TmpIndexFileName, ".%d.vec", i);
		  char* TmpFileName2 = ReplaceString(TmpFileName, TmpIndexFileName0, TmpIndexFileName);
		  if (VectorGammaX3[GammaYIndex][i].ReadVector(TmpFileName2) == false)
		    {
		      cout << "can't read " << TmpFileName2 << endl;
		      return 0;
		    }
		  delete[] TmpFileName2;
		}
	      VectorGammaX3[GammaYIndex].Multiply(GaugeFixingMatrices[GammaXIndex + 1][GammaYIndex]);
	      delete[] TmpFileName;
	      GammaY  += GammaYStep;
	    }
	  for (int i = 0; i < NbrStates ; ++i)
	    {
	      Tmp  = (VectorGammaX1[0][i] * VectorGammaX2[0][i]);
	      Tmp += (VectorGammaX3[0][i] * VectorGammaX2[1][i]);
	      Tmp -= (VectorGammaX1[0][i] * VectorGammaX2[1][i]);
	      Tmp -= (VectorGammaX3[0][i] * VectorGammaX2[0][i]);
	      ChernNumber += 0.5 * Tmp.Im;
	    }
	  for (int GammaYIndex = 1; GammaYIndex < NbrGammaY; ++GammaYIndex)
	    {
	      for (int i = 0; i < NbrStates ; ++i)
		{
		  Tmp  = (VectorGammaX1[GammaYIndex][i] * VectorGammaX2[GammaYIndex - 1][i]);
		  Tmp += (VectorGammaX3[GammaYIndex][i] * VectorGammaX2[GammaYIndex + 1][i]);
		  Tmp -= (VectorGammaX1[GammaYIndex][i] * VectorGammaX2[GammaYIndex + 1][i]);
		  Tmp -= (VectorGammaX3[GammaYIndex][i] * VectorGammaX2[GammaYIndex - 1][i]);
		  ChernNumber += 0.25 * Tmp.Im;
		}
// 	      ComplexMatrix Ax = HermitianMultiply(VectorGammaX2[GammaYIndex], VectorGammaX3[GammaYIndex]) -  HermitianMultiply(VectorGammaX2[GammaYIndex], VectorGammaX1[GammaYIndex]);
// 	      ComplexMatrix Ay = HermitianMultiply(VectorGammaX2[GammaYIndex], VectorGammaX2[GammaYIndex + 1]) -  HermitianMultiply(VectorGammaX2[GammaYIndex], VectorGammaX2[GammaYIndex - 1]);
// 	      ComplexMatrix TmpA2 = (Ax*  Ax) + (Ay * Ay);
// 	      cout << TmpA2.ComplexTr() << endl;
// 	      ChernNumber += 0.25 * (TmpA2.ComplexTr()).Im;
	      
	    }
	  for (int i = 0; i < NbrStates ; ++i)
	    {
	      Tmp  = (VectorGammaX1[NbrGammaY][i] * VectorGammaX2[NbrGammaY - 1][i]);
	      Tmp += (VectorGammaX3[NbrGammaY][i] * VectorGammaX2[NbrGammaY][i]);
	      Tmp -= (VectorGammaX1[NbrGammaY][i] * VectorGammaX2[NbrGammaY][i]);
	      Tmp -= (VectorGammaX3[NbrGammaY][i] * VectorGammaX2[NbrGammaY- 1][i]);
	      ChernNumber += 0.5 * Tmp.Im;
	    }
	  ComplexMatrix* VectorGammaX4 = VectorGammaX1;
	  VectorGammaX1 = VectorGammaX2;
	  VectorGammaX2 = VectorGammaX3;
	  VectorGammaX3 = VectorGammaX4;
	  GammaX  += GammaXStep;      
	}
      
      for (int i = 0; i < NbrStates ; ++i)
	{
	  Tmp  = (VectorGammaX1[0][i] * VectorGammaX2[1][i]);
	  Tmp += (VectorGammaX2[0][i] * VectorGammaX2[0][i]);
	  Tmp -= (VectorGammaX1[0][i] * VectorGammaX2[0][i]);
	  Tmp -= (VectorGammaX2[0][i] * VectorGammaX2[1][i]);
	  ChernNumber += Tmp.Im;
	}
      for (int GammaYIndex = 1; GammaYIndex < NbrGammaY; ++GammaYIndex)
	{
	  for (int i = 0; i < NbrStates ; ++i)
	    {
	      Tmp  = (VectorGammaX1[GammaYIndex][i] * VectorGammaX2[GammaYIndex + 1][i]);
	      Tmp += (VectorGammaX2[GammaYIndex][i] * VectorGammaX2[GammaYIndex - 1][i]);
	      Tmp -= (VectorGammaX1[GammaYIndex][i] * VectorGammaX2[GammaYIndex - 1][i]);
	      Tmp -= (VectorGammaX2[GammaYIndex][i] * VectorGammaX2[GammaYIndex + 1][i]);
	      ChernNumber += 0.5 * Tmp.Im;
	    }
	}
      for (int i = 0; i < NbrStates ; ++i)
 	{
	  Tmp  = (VectorGammaX1[NbrGammaY][i] * VectorGammaX2[NbrGammaY][i]);
	  Tmp += (VectorGammaX2[NbrGammaY][i] * VectorGammaX2[NbrGammaY - 1][i]);
	  Tmp -= (VectorGammaX1[NbrGammaY][i] * VectorGammaX2[NbrGammaY - 1][i]);
	  Tmp -= (VectorGammaX2[NbrGammaY][i] * VectorGammaX2[NbrGammaY][i]);
	  ChernNumber += Tmp.Im;
	}
            
      ChernNumber /= 2.0 * M_PI;
      ChernNumber *= 2.0;
      ChernNumber /= (double) NbrStates;
      cout << "C=" << ChernNumber << endl;
    }
  return 0;
}

// find the gauge transformation to get a smooth gauge choice between two set of vectors
//
// leftVectors = vectors for which tha gauge transformation has to be computed
// rightVectors = reference vectors that we have to smoothly evolve from
// rightMatrix = additional transformation that has to be applied to the reference vectors
// return value = gauge transformation

ComplexMatrix FindGaugeTransformation (ComplexVector* leftVectors, ComplexVector* rightVectors, int nbrVectors, ComplexMatrix& rightMatrix)
{
  ComplexMatrix TmpMatrix (nbrVectors, nbrVectors);
  for (int i = 0; i < nbrVectors; ++i)
    {
      for (int j = 0; j < nbrVectors; ++j)
	{
	  TmpMatrix.SetMatrixElement(i, j, leftVectors[i] * rightVectors[j]);
	}      
    }
  ComplexMatrix TmpMatrix2;
  TmpMatrix2.Copy(rightMatrix);
  TmpMatrix2.HermitianTranspose();
  TmpMatrix.Multiply(TmpMatrix2);
  TmpMatrix.OrthoNormalizeColumns();
  return TmpMatrix;
}
