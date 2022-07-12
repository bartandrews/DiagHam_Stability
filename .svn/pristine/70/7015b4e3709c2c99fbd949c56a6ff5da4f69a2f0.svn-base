#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string.h>


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// compute all the quasihole states for a given set of root configurations
//
// space = pointer to the Hilbert where the quasihole states should be expressed
// rootConfigurations = array that contains all the root configurations
// nbrQuasiholeStates = number of root configurations
// kValue = k value of the clustered (k,r) principle
// rValue = r value of the clustered (k,r) principle
// nbrParticles = number of particles
// lzMax = number of flux quantum
// totalLz = total angular momentum along the z direction
// filePrefix = output file prefix
// return value = orthonomalized basis of quasihole states
RealMatrix FQHESphereQuasiholeMatrixElementsComputeQuasiholeStates(ParticleOnSphere* space, unsigned long** rootConfigurations, int nbrQuasiholeStates, int kValue, int rValue, int nbrParticles, int lzMax, int totalLz, char* filePrefix);


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("FQHESphereQuasiholeMatrixElements" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 9);
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "largest number of particles to consider. If zero, consider all the possible number of particles compatible with the number of flux quanta and the clustering properties ", -1);
  (*OutputGroup) += new BooleanOption  ('\n', "disable-ascii", "only store the binary files");
  (*OutputGroup) += new BooleanOption  ('\n', "disable-directory", "do not create a directory to store the data");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereQuasiholeMatrixElements -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int KValue = 1;
  int RValue = 2;
  bool Statistics = true;
  int FermionFactor = 1;
  if (Statistics == false)
    FermionFactor = 0;
  
  int LzMax = Manager.GetInteger("nbr-flux");
  
  int MaxRightStateNbrParticles = Manager.GetInteger("nbr-particles");
  int MinRightStateNbrParticles = MaxRightStateNbrParticles;
  if (MaxRightStateNbrParticles < 0)
    {
      MinRightStateNbrParticles = 0;
      if (Statistics == true)
	{
	  MaxRightStateNbrParticles = (KValue * (LzMax + 1 + RValue)) / (RValue + (KValue * FermionFactor));
	}
      else
	{
	  MaxRightStateNbrParticles = (KValue * (LzMax + RValue)) / RValue;
	}
    }

  char* TmpDirectory = 0;
  if (Manager.GetBoolean("disable-directory") == false)
    {
      TmpDirectory = new char[512];
      sprintf (TmpDirectory, "nphi_%d_sphere", LzMax); 
      CreateDirectory(TmpDirectory);
    }

  char* FilePrefix = new char[512];
  if (Statistics == true)
    {
      if (TmpDirectory != 0)
	sprintf (FilePrefix, "%s/fermions", TmpDirectory);
      else
	sprintf (FilePrefix, "fermions");
    }
  else
    {
      if (TmpDirectory != 0)
	sprintf (FilePrefix, "%s/boson", TmpDirectory);
      else
	sprintf (FilePrefix, "bosons");
    }


  for (int RightStateNbrParticles = MinRightStateNbrParticles; RightStateNbrParticles <= MaxRightStateNbrParticles; ++RightStateNbrParticles)
    {
      int LeftStateNbrParticles = RightStateNbrParticles - 1;
      int RightStateMaxTotalLz = RightStateNbrParticles * LzMax - (((RValue + (KValue * FermionFactor)) * RightStateNbrParticles * (RightStateNbrParticles - 1)));
      int LeftStateMaxTotalLz = LeftStateNbrParticles * LzMax - (((RValue + (KValue * FermionFactor)) * LeftStateNbrParticles * (LeftStateNbrParticles - 1)));
      
      for (int RightTotalLz = -RightStateMaxTotalLz; RightTotalLz <= RightStateMaxTotalLz; RightTotalLz += 2)
	{
	  ParticleOnSphere* RightSpace = 0;
	  if (Statistics == true)
	    {
	      RightSpace = new FermionOnSphere(RightStateNbrParticles, RightTotalLz, LzMax);
	    }
	  else
	    {
	      RightSpace = new BosonOnSphereShort(RightStateNbrParticles, RightTotalLz, LzMax);
	    }     
	  int RightNbrQuasiholeStates = 0;
	  for (long i = 0l; i < RightSpace->GetHilbertSpaceDimension(); ++i)
	    {
	      if (RightSpace->HasPauliExclusions(i, KValue, RValue + (KValue * FermionFactor)) == true)
		{
		  ++RightNbrQuasiholeStates;
		}
	    }
	  if (RightNbrQuasiholeStates > 0)
	    {
	      cout << "---------------------------------------" << endl;
	      cout << "---------------------------------------" << endl;
	      cout << "processing right states with N=" << RightStateNbrParticles << " LzMax=" << LzMax << " TotalLz=" << RightTotalLz << endl;
	      cout << "found " << RightNbrQuasiholeStates << " quasihole states" << endl;
	      unsigned long** RightRootConfigurations = new unsigned long*[RightNbrQuasiholeStates];
	      RightNbrQuasiholeStates = 0;
	      cout << "admissible configurations : " << endl;
	      for (long i = 0l; i < RightSpace->GetHilbertSpaceDimension(); ++i)
		{
		  if (RightSpace->HasPauliExclusions(i, KValue, RValue + (KValue * FermionFactor)) == true)
		    {
		      RightRootConfigurations[RightNbrQuasiholeStates] = new unsigned long[LzMax + 1];
		      RightSpace->GetOccupationNumber(i, RightRootConfigurations[RightNbrQuasiholeStates]);
		      for (int j = 0; j <= LzMax; ++j)
			{
			  cout << RightRootConfigurations[RightNbrQuasiholeStates][j] << " ";
			}
		      cout << endl;
		      ++RightNbrQuasiholeStates;
		    }
		}
	      RealMatrix RightVectors = FQHESphereQuasiholeMatrixElementsComputeQuasiholeStates(RightSpace, RightRootConfigurations, 
												RightNbrQuasiholeStates, KValue, RValue,
												RightStateNbrParticles, LzMax, RightTotalLz, FilePrefix);
	      for (int OperatorLzValue = -LzMax; OperatorLzValue <= LzMax; OperatorLzValue += 2)
		{
		  int ShiftedOperatorLzValue = (OperatorLzValue + LzMax) >> 1;
		  int LeftTotalLz = RightTotalLz - OperatorLzValue;
		  if ((LeftTotalLz >= -LeftStateMaxTotalLz) && (LeftTotalLz <= LeftStateMaxTotalLz))
		    {
		      ParticleOnSphere* LeftSpace = 0;
		      if (Statistics == true)
			{
			  LeftSpace = new FermionOnSphere(LeftStateNbrParticles, LeftTotalLz, LzMax);
			}
		      else
			{
			  LeftSpace = new BosonOnSphereShort(LeftStateNbrParticles, LeftTotalLz, LzMax);
			}     
		      int LeftNbrQuasiholeStates = 0;
		      for (long i = 0l; i < LeftSpace->GetHilbertSpaceDimension(); ++i)
			{
			  if (LeftSpace->HasPauliExclusions(i, KValue, RValue + (KValue * FermionFactor)) == true)
			    {
			      ++LeftNbrQuasiholeStates;
			    }
			}
		      if (LeftNbrQuasiholeStates > 0)
			{
			  cout << "  computing <Psi_{N-1}|c_{" << OperatorLzValue << "}|Psi_{N}>" << endl;
			  cout << "  found " << LeftNbrQuasiholeStates << " quasihole states with N-1=" << LeftStateNbrParticles << " LzMax=" << LzMax << " TotalLz=" << LeftTotalLz << endl;
			  unsigned long** LeftRootConfigurations = new unsigned long*[LeftNbrQuasiholeStates];
			  LeftNbrQuasiholeStates = 0;
			  cout << "  admissible configurations : " << endl;
			  for (long i = 0l; i < LeftSpace->GetHilbertSpaceDimension(); ++i)
			    {
			      if (LeftSpace->HasPauliExclusions(i, KValue, RValue + (KValue * FermionFactor)) == true)
				{
				  LeftRootConfigurations[LeftNbrQuasiholeStates] = new unsigned long[LzMax + 1];
				  LeftSpace->GetOccupationNumber(i, LeftRootConfigurations[LeftNbrQuasiholeStates]);
				  cout << "  ";
				  for (int j = 0; j <= LzMax; ++j)
				    {
				      cout << LeftRootConfigurations[LeftNbrQuasiholeStates][j] << " ";
				    }
				  cout << endl;
				  ++LeftNbrQuasiholeStates;
				}
			    }
			  RealMatrix LeftVectors = FQHESphereQuasiholeMatrixElementsComputeQuasiholeStates(LeftSpace, LeftRootConfigurations, 
													   LeftNbrQuasiholeStates, KValue, RValue,
													   LeftStateNbrParticles, LzMax, LeftTotalLz, FilePrefix);
			  RightSpace->SetTargetSpace(LeftSpace);
			  RealMatrix TmpOutputMatrix(LeftNbrQuasiholeStates, RightNbrQuasiholeStates, true);
			  for (int i = 0; i < RightNbrQuasiholeStates; ++i)
			    {
			      RealVector TmpVector (LeftSpace->GetLargeHilbertSpaceDimension(), true);
			      for (int j = 0; j < RightSpace->GetHilbertSpaceDimension(); ++j)
				{
				  double TmpCoefficient = 0.0;
				  int TmpIndex = RightSpace->A(j, ShiftedOperatorLzValue, TmpCoefficient);
				  if (TmpIndex < LeftSpace->GetHilbertSpaceDimension())
				    {
				      TmpVector[TmpIndex] += TmpCoefficient * RightVectors[i][j];
				    }
				}
			      for (int j = 0; j < LeftNbrQuasiholeStates; ++j)
				{
				  double Tmp = -(LeftVectors[j] * TmpVector);
				  cout << "    <" << j << "|c_{" << OperatorLzValue << "}|" << i << ">=" << Tmp << endl;
				  TmpOutputMatrix.SetMatrixElement(j, i, Tmp);
				}			  
			    }
			  delete[] LeftRootConfigurations;
			  char* TmpOutputFileName = new char[256 + strlen(FilePrefix)];
			  if (Statistics == true)
			    {
			      sprintf (TmpOutputFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_c_%d.mat", FilePrefix, KValue, RValue, RightStateNbrParticles, 
				       LzMax, RightTotalLz, OperatorLzValue);
			    }
			  else
			    {
			      sprintf (TmpOutputFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_c_%d.mat", FilePrefix, KValue, RValue, RightStateNbrParticles, 
				       LzMax, RightTotalLz, OperatorLzValue);
			    }
			  char* AsciiTmpOutputFileName = new char[16 + strlen(TmpOutputFileName)];
			  sprintf (AsciiTmpOutputFileName, "%s.txt", TmpOutputFileName);
			  if (Manager.GetBoolean("disable-ascii") == false)
			    TmpOutputMatrix.WriteAsciiMatrix(AsciiTmpOutputFileName, true);
			  TmpOutputMatrix.WriteMatrix(TmpOutputFileName);			  
			  delete[] TmpOutputFileName;
			}
		      delete LeftSpace;		      
		    }
		  RealSymmetricMatrix TmpOutputMatrix(RightNbrQuasiholeStates, true);
		  for (int i = 0; i < RightNbrQuasiholeStates; ++i)
		    {
		      for (int k = i; k < RightNbrQuasiholeStates; ++k)
			{
			  double TmpCoefficient = 0.0;
			  for (int j = 0; j < RightSpace->GetHilbertSpaceDimension(); ++j)
			    {
			      TmpCoefficient += RightSpace->AdA(j, ShiftedOperatorLzValue) * RightVectors[k][j] * RightVectors[i][j];
			    }
			  TmpOutputMatrix.SetMatrixElement(i, k, TmpCoefficient);
			}
		    }
		  char* TmpOutputFileName = new char[256 + strlen(FilePrefix)];
		  if (Statistics == true)
		    {
		      sprintf (TmpOutputFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d.mat", FilePrefix, KValue, RValue, RightStateNbrParticles, LzMax, RightTotalLz, OperatorLzValue);
		    }
		  else
		    {
		      sprintf (TmpOutputFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d.mat", FilePrefix, KValue, RValue, RightStateNbrParticles, LzMax, RightTotalLz, OperatorLzValue);
		    }
		  char* AsciiTmpOutputFileName = new char[16 + strlen(TmpOutputFileName)];
		  sprintf (AsciiTmpOutputFileName, "%s.txt", TmpOutputFileName);
		  if (Manager.GetBoolean("disable-ascii") == false)
		    TmpOutputMatrix.WriteAsciiMatrix(AsciiTmpOutputFileName, true);
		  TmpOutputMatrix.WriteMatrix(TmpOutputFileName);
		  delete[] TmpOutputFileName;
		}
	      delete[] RightRootConfigurations;
	    }
	  delete RightSpace;
	}
    }
  return 0;
}

// compute all the quasihole states for a given set of root configurations
//
// space = pointer to the Hilbert where the quasihole states should be expressed
// rootConfigurations = array that contains all the root configurations
// nbrQuasiholeStates = number of root configurations
// kValue = k value of the clustered (k,r) principle
// rValue = r value of the clustered (k,r) principle
// nbrParticles = number of particles
// lzMax = number of flux quantum
// totalLz = total angular momentum along the z direction
// filePrefix = output file prefix
// return value = orthonomalized basis of quasihole states

RealMatrix FQHESphereQuasiholeMatrixElementsComputeQuasiholeStates(ParticleOnSphere* space, unsigned long** rootConfigurations, int nbrQuasiholeStates, int kValue, int rValue, int nbrParticles, int lzMax, int totalLz, char* filePrefix)
{
  char* QuasiholeVectorFileName = new char[256 + strlen(filePrefix)];
  sprintf (QuasiholeVectorFileName, "%s_qh_states_k_%d_r_%d_n_%d_nphi_%d_lz_%d.mat", filePrefix, kValue, rValue, nbrParticles, lzMax, totalLz);
  if (IsFile(QuasiholeVectorFileName) == true)
    {
      RealMatrix QuasiholeVectors;
      if (QuasiholeVectors.ReadMatrix(QuasiholeVectorFileName) == false)
	{
	  cout << "error while reading " << QuasiholeVectorFileName << endl;
	}
      return QuasiholeVectors;
    }
  RealVector* QuasiholeVectors = new RealVector[nbrQuasiholeStates];
  int* ReferenceState = new int[lzMax + 1];
  double Alpha = ((double) -(kValue + 1)) / ((double) (rValue - 1));
  if (space->GetParticleStatistic() == AbstractQHEParticle::BosonicStatistic)
    {
      for (int i = 0; i < nbrQuasiholeStates; ++i)
	{	  
	  for (int j = 0; j <= lzMax; ++j)
	    {
	      ReferenceState[j] = (int) rootConfigurations[i][j];
	    }
	  double Alpha = ((double) -kValue) / ((double) (rValue + kValue));
	  if (nbrParticles >1)
	    {
	      BosonOnSphereHaldaneBasisShort SqueezedSpace (nbrParticles, totalLz, lzMax, ReferenceState);
	      RealVector TmpState(SqueezedSpace.GetLargeHilbertSpaceDimension(), true);
	      SqueezedSpace.GenerateJackPolynomial(TmpState, Alpha);
	      SqueezedSpace.ConvertFromUnnormalizedMonomial(TmpState);
	      QuasiholeVectors[i] = SqueezedSpace.ConvertToNbodyBasis(TmpState, *((BosonOnSphereShort*) space));
	    }
	  else
	    {
	      QuasiholeVectors[i] = RealVector(1);
	      QuasiholeVectors[i][0] = 1.0;
	    }
	}
    }
  else
    {
      for (int i = 0; i < nbrQuasiholeStates; ++i)
	{	  
	  for (int j = 0; j <= lzMax; ++j)
	    {
	      ReferenceState[j] = (int) rootConfigurations[i][j];
	    }
	  if (nbrParticles >1)
	    {
	      FermionOnSphereHaldaneBasis SqueezedSpace (nbrParticles, totalLz, lzMax, ReferenceState);
	      RealVector TmpState(SqueezedSpace.GetLargeHilbertSpaceDimension(), true);
	      SqueezedSpace.GenerateJackPolynomial(TmpState, Alpha);
	      SqueezedSpace.ConvertFromUnnormalizedMonomial(TmpState);
	      QuasiholeVectors[i] = SqueezedSpace.ConvertToNbodyBasis(TmpState, *((FermionOnSphere*) space));
	    }
	  else
	    {
	      QuasiholeVectors[i] = RealVector(1);
	      QuasiholeVectors[i][0] = 1.0;
	    }
	}
    }
  delete[] ReferenceState;
  char* QuasiholeDimensionFileName = new char[256 + strlen(filePrefix)];
  sprintf (QuasiholeDimensionFileName, "%s_qh_states_k_%d_r_%d_nphi_%d.dat", filePrefix, kValue, rValue, lzMax);
  if (IsFile(QuasiholeDimensionFileName) == true)
    {
      ofstream File;
      File.open(QuasiholeDimensionFileName, ios::out | ios::app);
      File << lzMax << " " << nbrParticles << " " << totalLz << " " << nbrQuasiholeStates << endl;
      File.close();
   }
  else
    {
      ofstream File;
      File.open(QuasiholeDimensionFileName, ios::out);
      File << "# N_Phi N Lz dim" << endl; 
      File << lzMax << " " << nbrParticles << " " << totalLz << " " << nbrQuasiholeStates << endl;
      File.close();
    }
  RealMatrix QuasiholeVectors2(QuasiholeVectors, nbrQuasiholeStates);
  QuasiholeVectors2.OrthoNormalizeColumns();
  if (QuasiholeVectors2.WriteMatrix(QuasiholeVectorFileName) == false)
    {
      cout << "error while writing " << QuasiholeVectorFileName << endl;
    }
  return QuasiholeVectors2;
}
