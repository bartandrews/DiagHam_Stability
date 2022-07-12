#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Options/Options.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Operator/ParticleOnSphereDensityOperator.h"
#include "Operator/ParticleOnSphereAnnihilationOperator.h"

#include "Architecture/ArchitectureOperation/MultipleVectorOperatorMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/FQHECylinderMultipleJackGeneratorOperation.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <sys/time.h>


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
// ratio = cylinder aspect ratio
// filePrefix = output file prefix
// architecture = pointer to the architecture
// return value = orthonomalized basis of quasihole states
RealMatrix FQHECylinderQuasiholeMatrixElementsComputeQuasiholeStates(ParticleOnSphere* space, unsigned long** rootConfigurations, int nbrQuasiholeStates,
								     int kValue, int rValue, int nbrParticles, int lzMax, int totalLz, double ratio, char* filePrefix, AbstractArchitecture* architecture);


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("FQHECylinderQuasiholeMatrixElements" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 9);
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "largest number of particles to consider. If negative, consider all the possible number of particles compatible with the number of flux quanta and the clustering properties ", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "fix-nbrparticles", "if positive, only consider a given number of particles", -1);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "aspect-ratio", "aspect ratio of the cylinder", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter", "if non zero, fix the cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "offdiagonal-density", "if non-zero, compute all the density terms c^+_{m+x}c_m, with x lower or equal to the provided value", 0);
  (*OutputGroup) += new BooleanOption  ('\n', "enable-ascii", "store the files both in binary format and text format");
  (*OutputGroup) += new BooleanOption  ('\n', "disable-directory", "do not create a directory to store the data");
  (*OutputGroup) += new BooleanOption  ('\n', "show-admissible", "show the admissible configurations");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderQuasiholeMatrixElements -h" << endl;
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
  double Ratio = Manager.GetDouble("aspect-ratio");
  double Perimeter = Manager.GetDouble("cylinder-perimeter");
  if (Perimeter > 0.0)
    {
      Ratio = (Perimeter * Perimeter) / (2.0 * M_PI * (LzMax + 1));
    }  
  cout << "Cylinder geometry : Length=" << (Perimeter / Ratio) <<  ", Perimeter=" << Perimeter << ", aspect ratio=" << Ratio << endl;
  
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
  if (Manager.GetInteger("fix-nbrparticles") > 0)
    {
      MinRightStateNbrParticles = Manager.GetInteger("fix-nbrparticles");
      MaxRightStateNbrParticles = Manager.GetInteger("fix-nbrparticles");
    }

  char* TmpDirectory = 0;
  if (Manager.GetBoolean("disable-directory") == false)
    {
      TmpDirectory = new char[512];
      if (Manager.GetInteger("fix-nbrparticles") > 0)
	{
	  if (Perimeter > 0.0)	
	    {
	      sprintf (TmpDirectory, "nphi_%d_n_%d_cylinder_perimeter_%.6f", LzMax, MinRightStateNbrParticles, Perimeter); 
	    }
	  else
	    {
	      sprintf (TmpDirectory, "nphi_%d_n_%d_cylinder_ratio_%.6f", LzMax, MinRightStateNbrParticles, Ratio); 
	    }
	}  
      else
	{
	  if (Perimeter > 0.0)	
	    {
	      sprintf (TmpDirectory, "nphi_%d_cylinder_perimeter_%.6f", LzMax, Perimeter); 
	    }
	  else
	    {
	      sprintf (TmpDirectory, "nphi_%d_cylinder_ratio_%.6f", LzMax, Ratio); 
	    }
	}
      CreateDirectory(TmpDirectory);
    }

  char* FilePrefix = new char[512];
  if (Perimeter > 0.0)	
    {
      if (Statistics == true)
	{
	  if (TmpDirectory != 0)
	    sprintf (FilePrefix, "%s/fermions_cylinder_perimeter_%.6f", TmpDirectory, Perimeter);
	  else
	    sprintf (FilePrefix, "fermions_cylinder_perimeter_%.6f", Perimeter);
	}
      else
	{
	  if (TmpDirectory != 0)
	    sprintf (FilePrefix, "%s/bosons_cylinder_perimeter_%.6f", TmpDirectory, Perimeter);
	  else
	    sprintf (FilePrefix, "bosons_cylinder_perimeter_%.6f", Perimeter);
	}
    }
  else
    {
      if (Statistics == true)
	{
	  if (TmpDirectory != 0)
	    sprintf (FilePrefix, "%s/fermions_cylinder_ratio_%.6f", TmpDirectory, Ratio);
	  else
	    sprintf (FilePrefix, "fermions_cylinder_ratio_%.6f", Ratio);
	}
      else
	{
	  if (TmpDirectory != 0)
	    sprintf (FilePrefix, "%s/bosons_cylinder_ratio_%.6f", TmpDirectory, Ratio);
	  else
	    sprintf (FilePrefix, "bosons_cylinder_ratio_%.6f", Ratio);
	}      
    }

  if (Manager.GetInteger("fix-nbrparticles") > 0)
    {
      for (int RightStateNbrParticles = 0; RightStateNbrParticles <= 0; ++RightStateNbrParticles)
	{
	  int RightStateMaxTotalLz = 0;
	  for (int RightTotalLz = -RightStateMaxTotalLz; RightTotalLz <= RightStateMaxTotalLz; RightTotalLz += 2)
	    {
	      timeval TotalStartingTime;
	      gettimeofday (&(TotalStartingTime), 0);
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
	      timeval TotalEndingTime;
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));               
	      cout << "finding admissible configurations done in " << Dt << "s" << endl;
	      if (RightNbrQuasiholeStates > 0)
		{
		  cout << "---------------------------------------" << endl;
		  cout << "---------------------------------------" << endl;
		  cout << "processing right states with N=" << RightStateNbrParticles << " LzMax=" << LzMax << " TotalLz=" << RightTotalLz << endl;
		  cout << "found " << RightNbrQuasiholeStates << " quasihole states" << endl;
		  unsigned long** RightRootConfigurations = new unsigned long*[RightNbrQuasiholeStates];
		  RightNbrQuasiholeStates = 0;
		  if (Manager.GetBoolean("show-admissible"))
		    {
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
		    }
		  else
		    {
		      for (long i = 0l; i < RightSpace->GetHilbertSpaceDimension(); ++i)
			{
			  if (RightSpace->HasPauliExclusions(i, KValue, RValue + (KValue * FermionFactor)) == true)
			    {
			      RightRootConfigurations[RightNbrQuasiholeStates] = new unsigned long[LzMax + 1];
			      RightSpace->GetOccupationNumber(i, RightRootConfigurations[RightNbrQuasiholeStates]);
			      ++RightNbrQuasiholeStates;
			    }
			}
		    }
		  RealMatrix RightVectors = FQHECylinderQuasiholeMatrixElementsComputeQuasiholeStates(RightSpace, RightRootConfigurations, 
												      RightNbrQuasiholeStates, KValue, RValue,
												      RightStateNbrParticles, LzMax, RightTotalLz, Ratio, FilePrefix, Architecture.GetArchitecture());
		  for (int OperatorLzValue = -LzMax; OperatorLzValue <= LzMax; OperatorLzValue += 2)
		    {
		      int ShiftedOperatorLzValue = (OperatorLzValue + LzMax) >> 1;
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
		      sprintf (TmpOutputFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d.mat", FilePrefix, KValue, RValue, RightStateNbrParticles, LzMax, RightTotalLz, OperatorLzValue);
		      char* AsciiTmpOutputFileName = new char[16 + strlen(TmpOutputFileName)];
		      sprintf (AsciiTmpOutputFileName, "%s.txt", TmpOutputFileName);
		      if (Manager.GetBoolean("enable-ascii") == true)
			TmpOutputMatrix.WriteAsciiMatrix(AsciiTmpOutputFileName, true);
		      TmpOutputMatrix.WriteMatrix(TmpOutputFileName);
		      delete[] TmpOutputFileName;
		    }
		  delete[] RightRootConfigurations;
		}
	      delete RightSpace;
	    }
	}
    }
  for (int RightStateNbrParticles = MinRightStateNbrParticles; RightStateNbrParticles <= MaxRightStateNbrParticles; ++RightStateNbrParticles)
    {
      int RightStateMaxTotalLz = RightStateNbrParticles * LzMax - (((RValue + (KValue * FermionFactor)) * RightStateNbrParticles * (RightStateNbrParticles - 1)));
      int LeftStateNbrParticles = RightStateNbrParticles - 1;
      int LeftStateMaxTotalLz = LeftStateNbrParticles * LzMax - (((RValue + (KValue * FermionFactor)) * LeftStateNbrParticles * (LeftStateNbrParticles - 1)));

      
      for (int RightTotalLz = -RightStateMaxTotalLz; RightTotalLz <= RightStateMaxTotalLz; RightTotalLz += 2)
	{
	  timeval TotalStartingTime;
	  gettimeofday (&(TotalStartingTime), 0);
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
	  timeval TotalEndingTime;
	  gettimeofday (&(TotalEndingTime), 0);
	  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));               
	  cout << "finding admissible configurations done in " << Dt << "s" << endl;
	  if (RightNbrQuasiholeStates > 0)
	    {
	      cout << "---------------------------------------" << endl;
	      cout << "---------------------------------------" << endl;
	      cout << "processing right states with N=" << RightStateNbrParticles << " LzMax=" << LzMax << " TotalLz=" << RightTotalLz << endl;
	      cout << "found " << RightNbrQuasiholeStates << " quasihole states" << endl;
	      unsigned long** RightRootConfigurations = new unsigned long*[RightNbrQuasiholeStates];
	      RightNbrQuasiholeStates = 0;
	      if (Manager.GetBoolean("show-admissible"))
		{
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
		}
	      else
		{
		  for (long i = 0l; i < RightSpace->GetHilbertSpaceDimension(); ++i)
		    {
		      if (RightSpace->HasPauliExclusions(i, KValue, RValue + (KValue * FermionFactor)) == true)
			{
			  RightRootConfigurations[RightNbrQuasiholeStates] = new unsigned long[LzMax + 1];
			  RightSpace->GetOccupationNumber(i, RightRootConfigurations[RightNbrQuasiholeStates]);
			  ++RightNbrQuasiholeStates;
			}
		    }
		}
	      RealMatrix RightVectors = FQHECylinderQuasiholeMatrixElementsComputeQuasiholeStates(RightSpace, RightRootConfigurations, 
												  RightNbrQuasiholeStates, KValue, RValue,
												  RightStateNbrParticles, LzMax, RightTotalLz, Ratio, FilePrefix, Architecture.GetArchitecture());
	      RealVector* TmpAdAInputVectors = new RealVector[RightNbrQuasiholeStates];
	      RealVector* TmpAdAOutputVectors = new RealVector[RightNbrQuasiholeStates];
	      for (int j = 0; j < RightNbrQuasiholeStates; ++j)
		{
		  TmpAdAInputVectors[j] = RightVectors[j];
		  TmpAdAOutputVectors[j] = RealVector(RightSpace->GetHilbertSpaceDimension());
		}
	      for (int OperatorLzValue = -LzMax; OperatorLzValue <= LzMax; OperatorLzValue += 2)
		{
		  int ShiftedOperatorLzValue = (OperatorLzValue + LzMax) >> 1;
		  int LeftTotalLz = RightTotalLz - OperatorLzValue;
		  if ((MinRightStateNbrParticles != MaxRightStateNbrParticles) && (LeftTotalLz >= -LeftStateMaxTotalLz) && (LeftTotalLz <= LeftStateMaxTotalLz))
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
		      timeval TotalStartingTime;
		      gettimeofday (&(TotalStartingTime), 0);
		      if (LeftNbrQuasiholeStates > 0)
			{
			  cout << "  computing <Psi_{N-1}|c_{" << OperatorLzValue << "}|Psi_{N}>" << endl;
			  cout << "  found " << LeftNbrQuasiholeStates << " quasihole states with N-1=" << LeftStateNbrParticles << " LzMax=" << LzMax << " TotalLz=" << LeftTotalLz << endl;
			  unsigned long** LeftRootConfigurations = new unsigned long*[LeftNbrQuasiholeStates];
			  LeftNbrQuasiholeStates = 0;
			  if (Manager.GetBoolean("show-admissible"))
			    {
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
			    }
			  else
			    {
			      for (long i = 0l; i < LeftSpace->GetHilbertSpaceDimension(); ++i)
				{
				  if (LeftSpace->HasPauliExclusions(i, KValue, RValue + (KValue * FermionFactor)) == true)
				    {
				      LeftRootConfigurations[LeftNbrQuasiholeStates] = new unsigned long[LzMax + 1];
				      LeftSpace->GetOccupationNumber(i, LeftRootConfigurations[LeftNbrQuasiholeStates]);
				      ++LeftNbrQuasiholeStates;
				    }
				}
			    }
			  RealMatrix LeftVectors = FQHECylinderQuasiholeMatrixElementsComputeQuasiholeStates(LeftSpace, LeftRootConfigurations, 
													     LeftNbrQuasiholeStates, KValue, RValue,
													     LeftStateNbrParticles, LzMax, LeftTotalLz, Ratio, FilePrefix, Architecture.GetArchitecture());
			  RightSpace->SetTargetSpace(LeftSpace);
			  RealVector* TmpAOutputVectors = new RealVector[RightNbrQuasiholeStates];
			  for (int j = 0; j < RightNbrQuasiholeStates; ++j)
			    {
			      TmpAOutputVectors[j] = RealVector(LeftSpace->GetHilbertSpaceDimension());
			    }
 			  RealMatrix TmpOutputMatrix(LeftNbrQuasiholeStates, RightNbrQuasiholeStates, true);
			  Architecture.GetArchitecture()->SetDimension(RightSpace->GetHilbertSpaceDimension());		  
			  ParticleOnSphereAnnihilationOperator TmpAOperator (RightSpace, ShiftedOperatorLzValue);
			  MultipleVectorOperatorMultiplyOperation TmpOperation (&TmpAOperator, TmpAdAInputVectors, TmpAOutputVectors, RightNbrQuasiholeStates);
			  TmpOperation.ApplyOperation(Architecture.GetArchitecture());
			  for (int i = 0; i < RightNbrQuasiholeStates; ++i)
			    {
			      for (int k = 0; k < LeftNbrQuasiholeStates; ++k)
				{ 
				  TmpOutputMatrix.SetMatrixElement(k, i, -(LeftVectors[k] * TmpAOutputVectors[i]));
				}
			    }
			  delete[] TmpAOutputVectors;

// 			  RealMatrix TmpVectors(LeftSpace->GetLargeHilbertSpaceDimension(), RightNbrQuasiholeStates, true);
// 			  for (int j = 0; j < RightSpace->GetHilbertSpaceDimension(); ++j)
// 			    {
// 			      double TmpCoefficient = 0.0;
// 			      int TmpIndex = RightSpace->A(j, ShiftedOperatorLzValue, TmpCoefficient);
// 			      if (TmpIndex < LeftSpace->GetHilbertSpaceDimension())
// 				{
// 				  for (int i = 0; i < RightNbrQuasiholeStates; ++i)
// 				    TmpVectors[i][TmpIndex] += TmpCoefficient * RightVectors[i][j];
// 				}
// 			    }
// 			  for (int i = 0; i < RightNbrQuasiholeStates; ++i)
// 			    {
//  			      for (int j = 0; j < LeftNbrQuasiholeStates; ++j)
//  				{
//  				  double Tmp = -(LeftVectors[j] * TmpVectors[i]);
//  				  TmpOutputMatrix.SetMatrixElement(j, i, Tmp);
//  				}			  
// 			    }
			  for (int k = 0; k < LeftNbrQuasiholeStates; ++k)
			    {
			      delete[] LeftRootConfigurations[k];
			    }
			  delete[] LeftRootConfigurations;
			  char* TmpOutputFileName = new char[256 + strlen(FilePrefix)];
			  sprintf (TmpOutputFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_c_%d.mat", FilePrefix, KValue, RValue, RightStateNbrParticles, LzMax, RightTotalLz, OperatorLzValue);
			  char* AsciiTmpOutputFileName = new char[16 + strlen(TmpOutputFileName)];
			  sprintf (AsciiTmpOutputFileName, "%s.txt", TmpOutputFileName);
			  if (Manager.GetBoolean("enable-ascii") == true)
			    TmpOutputMatrix.WriteAsciiMatrix(AsciiTmpOutputFileName, true);
			  TmpOutputMatrix.WriteMatrix(TmpOutputFileName);			  
			  delete[] TmpOutputFileName;
			}
		      delete LeftSpace;		      
		      timeval TotalEndingTime;
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));               
		      cout << "computing c matrix elements done in " << Dt << "s" << endl;
		    }
		  timeval TotalStartingTime;
		  gettimeofday (&(TotalStartingTime), 0);
		  TotalStartingTime;
		  gettimeofday (&(TotalStartingTime), 0);
		  RealSymmetricMatrix TmpOutputMatrix(RightNbrQuasiholeStates, true);
		  
		  Architecture.GetArchitecture()->SetDimension(RightSpace->GetHilbertSpaceDimension());		  
		  ParticleOnSphereDensityOperator TmpOperator (RightSpace, ShiftedOperatorLzValue);
		  MultipleVectorOperatorMultiplyOperation TmpOperation (&TmpOperator, TmpAdAInputVectors, TmpAdAOutputVectors, RightNbrQuasiholeStates);
		  TmpOperation.ApplyOperation(Architecture.GetArchitecture());
 		  for (int i = 0; i < RightNbrQuasiholeStates; ++i)
 		    {
 		      for (int k = i; k < RightNbrQuasiholeStates; ++k)
			{ 
			  TmpOutputMatrix.SetMatrixElement(i, k, RightVectors[k] * TmpAdAOutputVectors[i]);
			}
		    }
// 		  RealMatrix TmpMatrix (RightSpace->GetHilbertSpaceDimension(), RightNbrQuasiholeStates);
// 		  for (int j = 0; j < RightSpace->GetHilbertSpaceDimension(); ++j)
// 		    {
//  		      double TmpCoefficient = RightSpace->AdA(j, ShiftedOperatorLzValue);
//  		      for (int i = 0; i < RightNbrQuasiholeStates; ++i)
//  			{
//  			  TmpMatrix[i][j] = TmpCoefficient * RightVectors[i][j];
//  			}
// 		    }
//  		  for (int i = 0; i < RightNbrQuasiholeStates; ++i)
//  		    {
//  		      for (int k = i; k < RightNbrQuasiholeStates; ++k)
// 			{ 
// 			  TmpOutputMatrix.SetMatrixElement(i, k, RightVectors[k] * TmpMatrix[i]);
// 			}
// 		    }
		  char* TmpOutputFileName = new char[256 + strlen(FilePrefix)];
		  sprintf (TmpOutputFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d.mat", FilePrefix, KValue, RValue, RightStateNbrParticles, LzMax, RightTotalLz, OperatorLzValue);
		  char* AsciiTmpOutputFileName = new char[16 + strlen(TmpOutputFileName)];
		  sprintf (AsciiTmpOutputFileName, "%s.txt", TmpOutputFileName);
		  if (Manager.GetBoolean("enable-ascii") == true)
		    TmpOutputMatrix.WriteAsciiMatrix(AsciiTmpOutputFileName, true);
		  TmpOutputMatrix.WriteMatrix(TmpOutputFileName);
		  delete[] TmpOutputFileName;
		  timeval TotalEndingTime;
		  gettimeofday (&(TotalEndingTime), 0);
		  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));               
		  cout << "computing c^+c matrix elements done in " << Dt << "s" << endl;
		  if (Manager.GetInteger("offdiagonal-density") > 0)
		    {
		      for (int MomentumTransfer = 1; MomentumTransfer <= Manager.GetInteger("offdiagonal-density"); ++MomentumTransfer)
			{
			  int ShiftedOperatorLzValue = (OperatorLzValue + LzMax) >> 1;
			  int LeftStateNbrParticles = RightStateNbrParticles;
			  int LeftTotalLz = RightTotalLz - (2 * MomentumTransfer);
			  if ((LeftTotalLz >= -RightStateMaxTotalLz) && (LeftTotalLz <= RightStateMaxTotalLz) && 
			      (abs(OperatorLzValue - (2 * MomentumTransfer)) <= LzMax))
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
			      timeval TotalStartingTime;
			      gettimeofday (&(TotalStartingTime), 0);
			      if (LeftNbrQuasiholeStates > 0)
				{
				  cout << "  computing <Psi_{N}|c^+_{" << (OperatorLzValue - (2 * MomentumTransfer)) << "}c_{" << OperatorLzValue << "}|Psi_{N}>" << endl;
				  cout << "  found " << LeftNbrQuasiholeStates << " quasihole states with N=" << LeftStateNbrParticles << " LzMax=" << LzMax << " TotalLz=" << LeftTotalLz << endl;
				  unsigned long** LeftRootConfigurations = new unsigned long*[LeftNbrQuasiholeStates];
				  LeftNbrQuasiholeStates = 0;
				  if (Manager.GetBoolean("show-admissible"))
				    {
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
				    }
				  else
				    {
				      for (long i = 0l; i < LeftSpace->GetHilbertSpaceDimension(); ++i)
					{
					  if (LeftSpace->HasPauliExclusions(i, KValue, RValue + (KValue * FermionFactor)) == true)
					    {
					      LeftRootConfigurations[LeftNbrQuasiholeStates] = new unsigned long[LzMax + 1];
					      LeftSpace->GetOccupationNumber(i, LeftRootConfigurations[LeftNbrQuasiholeStates]);
					      ++LeftNbrQuasiholeStates;
					    }
					}
				    }
				  RealMatrix LeftVectors = FQHECylinderQuasiholeMatrixElementsComputeQuasiholeStates(LeftSpace, LeftRootConfigurations, 
														     LeftNbrQuasiholeStates, KValue, RValue,
														     LeftStateNbrParticles, LzMax, LeftTotalLz, Ratio, FilePrefix, Architecture.GetArchitecture());
				  RightSpace->SetTargetSpace(LeftSpace);
				  RealVector* TmpAdAOutputVectors = new RealVector[RightNbrQuasiholeStates];
				  for (int j = 0; j < RightNbrQuasiholeStates; ++j)
				    {
				      TmpAdAOutputVectors[j] = RealVector(LeftSpace->GetHilbertSpaceDimension());
				    }
				  RealMatrix TmpOutputMatrix(LeftNbrQuasiholeStates, RightNbrQuasiholeStates, true);
				  Architecture.GetArchitecture()->SetDimension(RightSpace->GetHilbertSpaceDimension());		  
				  ParticleOnSphereDensityOperator TmpOperator (RightSpace, ShiftedOperatorLzValue - MomentumTransfer, ShiftedOperatorLzValue);
				  MultipleVectorOperatorMultiplyOperation TmpOperation (&TmpOperator, TmpAdAInputVectors, TmpAdAOutputVectors, RightNbrQuasiholeStates);
				  TmpOperation.ApplyOperation(Architecture.GetArchitecture());
				  for (int i = 0; i < RightNbrQuasiholeStates; ++i)
				    {
				      for (int k = 0; k < LeftNbrQuasiholeStates; ++k)
					{ 
					  TmpOutputMatrix.SetMatrixElement(k, i, -(LeftVectors[k] * TmpAdAOutputVectors[i]));
					}
				    }
				  delete[] TmpAdAOutputVectors;
				  for (int k = 0; k < LeftNbrQuasiholeStates; ++k)
				    {
				      delete[] LeftRootConfigurations[k];
				    }
				  delete[] LeftRootConfigurations;
				  char* TmpOutputFileName = new char[256 + strlen(FilePrefix)];
				  sprintf (TmpOutputFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d_%d.mat", FilePrefix, KValue, RValue, RightStateNbrParticles, LzMax, RightTotalLz, 
					   (OperatorLzValue - (2 * MomentumTransfer)), OperatorLzValue);
				  char* AsciiTmpOutputFileName = new char[16 + strlen(TmpOutputFileName)];
				  sprintf (AsciiTmpOutputFileName, "%s.txt", TmpOutputFileName);
				  if (Manager.GetBoolean("enable-ascii") == true)
				    TmpOutputMatrix.WriteAsciiMatrix(AsciiTmpOutputFileName, true);
				  TmpOutputMatrix.WriteMatrix(TmpOutputFileName);			  
				  delete[] TmpOutputFileName;
				}
			      delete LeftSpace;		      		      
			      timeval TotalEndingTime;
			      gettimeofday (&(TotalEndingTime), 0);
			      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
						    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));               
			      cout << "computing c^+c matrix elements done in " << Dt << "s" << endl;
			    }
			  LeftTotalLz = RightTotalLz + (2 * MomentumTransfer);
			  if ((LeftTotalLz >= -RightStateMaxTotalLz) && (LeftTotalLz <= RightStateMaxTotalLz) && 
			      (abs(OperatorLzValue + (2 * MomentumTransfer)) <= LzMax))
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
			      timeval TotalStartingTime;
			      gettimeofday (&(TotalStartingTime), 0);
			      if (LeftNbrQuasiholeStates > 0)
				{
				  cout << "  computing <Psi_{N}|c^+_{" << (OperatorLzValue + (2 * MomentumTransfer)) << "}c_{" << OperatorLzValue << "}|Psi_{N}>" << endl;
				  cout << "  found " << LeftNbrQuasiholeStates << " quasihole states with N=" << LeftStateNbrParticles << " LzMax=" << LzMax << " TotalLz=" << LeftTotalLz << endl;
				  unsigned long** LeftRootConfigurations = new unsigned long*[LeftNbrQuasiholeStates];
				  LeftNbrQuasiholeStates = 0;
				  if (Manager.GetBoolean("show-admissible"))
				    {
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
				    }
				  else
				    {
				      for (long i = 0l; i < LeftSpace->GetHilbertSpaceDimension(); ++i)
					{
					  if (LeftSpace->HasPauliExclusions(i, KValue, RValue + (KValue * FermionFactor)) == true)
					    {
					      LeftRootConfigurations[LeftNbrQuasiholeStates] = new unsigned long[LzMax + 1];
					      LeftSpace->GetOccupationNumber(i, LeftRootConfigurations[LeftNbrQuasiholeStates]);
					      ++LeftNbrQuasiholeStates;
					    }
					}
				    }
				  RealMatrix LeftVectors = FQHECylinderQuasiholeMatrixElementsComputeQuasiholeStates(LeftSpace, LeftRootConfigurations, 
														     LeftNbrQuasiholeStates, KValue, RValue,
														     LeftStateNbrParticles, LzMax, LeftTotalLz, Ratio, FilePrefix, Architecture.GetArchitecture());
				  RightSpace->SetTargetSpace(LeftSpace);
				  RealVector* TmpAdAOutputVectors = new RealVector[RightNbrQuasiholeStates];
				  for (int j = 0; j < RightNbrQuasiholeStates; ++j)
				    {
				      TmpAdAOutputVectors[j] = RealVector(LeftSpace->GetHilbertSpaceDimension());
				    }
				  RealMatrix TmpOutputMatrix(LeftNbrQuasiholeStates, RightNbrQuasiholeStates, true);
				  Architecture.GetArchitecture()->SetDimension(RightSpace->GetHilbertSpaceDimension());		  
				  ParticleOnSphereDensityOperator TmpOperator (RightSpace, ShiftedOperatorLzValue + MomentumTransfer, ShiftedOperatorLzValue);
				  MultipleVectorOperatorMultiplyOperation TmpOperation (&TmpOperator, TmpAdAInputVectors, TmpAdAOutputVectors, RightNbrQuasiholeStates);
				  TmpOperation.ApplyOperation(Architecture.GetArchitecture());
				  for (int i = 0; i < RightNbrQuasiholeStates; ++i)
				    {
				      for (int k = 0; k < LeftNbrQuasiholeStates; ++k)
					{ 
					  TmpOutputMatrix.SetMatrixElement(k, i, -(LeftVectors[k] * TmpAdAOutputVectors[i]));
					}
				    }
				  delete[] TmpAdAOutputVectors;
				  for (int k = 0; k < LeftNbrQuasiholeStates; ++k)
				    {
				      delete[] LeftRootConfigurations[k];
				    }
				  delete[] LeftRootConfigurations;
				  char* TmpOutputFileName = new char[256 + strlen(FilePrefix)];
				  sprintf (TmpOutputFileName, "%s_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d_%d.mat", FilePrefix, KValue, RValue, RightStateNbrParticles, LzMax, RightTotalLz, 
					   (OperatorLzValue + (2 * MomentumTransfer)), OperatorLzValue);
				  char* AsciiTmpOutputFileName = new char[16 + strlen(TmpOutputFileName)];
				  sprintf (AsciiTmpOutputFileName, "%s.txt", TmpOutputFileName);
				  if (Manager.GetBoolean("enable-ascii") == true)
				    TmpOutputMatrix.WriteAsciiMatrix(AsciiTmpOutputFileName, true);
				  TmpOutputMatrix.WriteMatrix(TmpOutputFileName);			  
				  delete[] TmpOutputFileName;
				}
			      delete LeftSpace;		      		      
			      timeval TotalEndingTime;
			      gettimeofday (&(TotalEndingTime), 0);
			      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
						    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));               
			      cout << "computing c^+c matrix elements done in " << Dt << "s" << endl;
			    }
			}
		    }
		}
	      delete[] RightRootConfigurations;
	      delete[] TmpAdAInputVectors;
	      delete[] TmpAdAOutputVectors;
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
// ratio = cylinder aspect ratio
// filePrefix = output file prefix
// architecture = pointer to the architecture
// return value = orthonomalized basis of quasihole states

RealMatrix FQHECylinderQuasiholeMatrixElementsComputeQuasiholeStates(ParticleOnSphere* space, unsigned long** rootConfigurations, int nbrQuasiholeStates, 
								     int kValue, int rValue, int nbrParticles, int lzMax, int totalLz, double ratio, char* filePrefix, AbstractArchitecture* architecture)
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
  timeval TotalStartingTime;
  gettimeofday (&(TotalStartingTime), 0);


  FQHECylinderMultipleJackGeneratorOperation TmpOperation (space, rootConfigurations, nbrQuasiholeStates, kValue, rValue, nbrParticles, lzMax, totalLz, ratio);
  TmpOperation.ApplyOperation(architecture);
  RealMatrix QuasiholeVectors2 = TmpOperation.GetBasis();

//   RealVector* QuasiholeVectors = new RealVector[nbrQuasiholeStates];
//   int* ReferenceState = new int[lzMax + 1];
//   double Alpha = ((double) -(kValue + 1)) / ((double) (rValue - 1));
//   if (space->GetParticleStatistic() == AbstractQHEParticle::BosonicStatistic)
//     {
//       for (int i = 0; i < nbrQuasiholeStates; ++i)
// 	{	  
// 	  for (int j = 0; j <= lzMax; ++j)
// 	    {
// 	      ReferenceState[j] = (int) rootConfigurations[i][j];
// 	    }
// 	  double Alpha = ((double) -kValue) / ((double) (rValue + kValue));
// 	  if (nbrParticles >1)
// 	    {
// 	      BosonOnSphereHaldaneBasisShort SqueezedSpace (nbrParticles, totalLz, lzMax, ReferenceState);
// 	      RealVector TmpState(SqueezedSpace.GetLargeHilbertSpaceDimension(), true);
// 	      SqueezedSpace.GenerateJackPolynomial(TmpState, Alpha);
// 	      SqueezedSpace.NormalizeJackToCylinder(TmpState, ratio);
// 	      QuasiholeVectors[i] = SqueezedSpace.ConvertToNbodyBasis(TmpState, *((BosonOnSphereShort*) space));
// 	    }
// 	  else
// 	    {
// 	      QuasiholeVectors[i] = RealVector(1);
// 	      QuasiholeVectors[i][0] = 1.0;
// 	    }
// 	}
//     }
//   else
//     {
//       for (int i = 0; i < nbrQuasiholeStates; ++i)
// 	{	  
// 	  for (int j = 0; j <= lzMax; ++j)
// 	    {
// 	      ReferenceState[j] = (int) rootConfigurations[i][j];
// 	    }
// 	  if (nbrParticles >1)
// 	    {
// 	      FermionOnSphereHaldaneBasis SqueezedSpace (nbrParticles, totalLz, lzMax, ReferenceState);
// 	      RealVector TmpState(SqueezedSpace.GetLargeHilbertSpaceDimension(), true);
// 	      SqueezedSpace.GenerateJackPolynomial(TmpState, Alpha);
// 	      SqueezedSpace.NormalizeJackToCylinder(TmpState, ratio);
// 	      QuasiholeVectors[i] = SqueezedSpace.ConvertToNbodyBasis(TmpState, *((FermionOnSphere*) space));
// 	    }
// 	  else
// 	    {
// 	      QuasiholeVectors[i] = RealVector(1);
// 	      QuasiholeVectors[i][0] = 1.0;
// 	    }
// 	}
//     }
//   delete[] ReferenceState;

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
//   RealMatrix QuasiholeVectors2(QuasiholeVectors, nbrQuasiholeStates);
//   QuasiholeVectors2.OrthoNormalizeColumns();
  if (QuasiholeVectors2.WriteMatrix(QuasiholeVectorFileName) == false)
    {
      cout << "error while writing " << QuasiholeVectorFileName << endl;
    }
  timeval TotalEndingTime;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));               
  cout << "subspace generation done in " << Dt << "s" << endl;
  return QuasiholeVectors2;
}
