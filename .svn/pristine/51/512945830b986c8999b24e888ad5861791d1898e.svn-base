#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"

#include "MathTools/FactorialCoefficient.h"
#include "MathTools/LongRational.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHEMPSCreateStateOperation.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Matrix/SparseComplexMatrix.h"
#include "Matrix/SparseRealMatrix.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;



int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETorusMPSCreateState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager (false, true);

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-fluxquanta", "set the total number of flux quanta and deduce the number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption ('y', "ky-momentum", "constraint on the total momentum modulo the maximum momentum (negative it has to be guessed from the topological sector)", -1);
  (*SystemGroup) += new SingleIntegerOption ('\n', "topological-sector", "set the topological sector", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "precalculation-blocksize", " indicates the size of the block (i.e. number of B matrices) for precalculations", 1);
  (*OutputGroup) += new BooleanOption ('t', "txt-output", "output the MPS state into a text file instead of a binary file");
  (*OutputGroup) += new SingleStringOption ('\n', "alternate-output", "use an alternate output file name instead of the default one");
  (*OutputGroup) += new BooleanOption ('\n', "no-normalization", "do not normalize the final state");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusMPSCreateState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  int TotalKy = 0;
  bool TwistedTorusFlag = false;
  char* OutputFileName = 0;
  char* OutputTxtFileName = 0;

  if (Manager.GetInteger("nbr-fluxquanta") <= 0)
    {
      return -1;
    }
  NbrFluxQuanta = Manager.GetInteger("nbr-fluxquanta");
  
  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  double AspectRatio = Manager.GetDouble("aspect-ratio");
  NbrParticles = MPSMatrix->GetMatrixNaturalNbrParticles(NbrFluxQuanta, true);


  int NbrMPSSumIndices;
  int* MPSSumIndices = MPSMatrix->GetTopologicalSectorIndices(Manager.GetInteger("topological-sector"), NbrMPSSumIndices);
  SparseRealMatrix* SparseBMatrices  = MPSMatrix->GetMatrices();
  SparseComplexMatrix* SparseComplexBMatrices = MPSMatrix->GetComplexMatrices();
  if (SparseBMatrices == 0)
    {
      TwistedTorusFlag = true;
      cout << "B matrix size = " << SparseComplexBMatrices[0].GetNbrRow() << "x" << SparseComplexBMatrices[0].GetNbrColumn() << endl;
    }
  else
    {
      cout << "B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;
    }
  SparseRealMatrix StringMatrix;
  if (Manager.GetBoolean("boson") == true)
    StringMatrix = MPSMatrix->GetTorusStringMatrix(0);
  else
    StringMatrix = MPSMatrix->GetTorusStringMatrix(NbrParticles);
  int PauliKValue;
  int PauliRValue;
  MPSMatrix->GetKRExclusionPrinciple(PauliKValue, PauliRValue);
  int ReducedBrillouinZoneSize = FindGCD(NbrParticles, NbrFluxQuanta);
  TotalKy = MPSMatrix->GetTorusMinimumKyMomentum(NbrParticles, NbrFluxQuanta, !Manager.GetBoolean("boson"));

  ParticleOnTorus* Space = 0;
  RealVector State ;
  ComplexVector ComplexState ;
  if (Manager.GetBoolean("boson") == true)
    {
      int MaxOccupation = Manager.GetInteger("boson-truncation");
      bool TotalKyFlag = false;
      
      while ((TotalKyFlag == false) && (TotalKy < NbrFluxQuanta))
	{
	  if (Manager.GetInteger("ky-momentum") >= 0)
	    {
	      TotalKy = Manager.GetInteger("ky-momentum");
	      TotalKyFlag = true;
	    }
	  cout << "checking Ky=" << TotalKy << " sector" << endl;
	  if (MaxOccupation > NbrParticles)
	    Space = new BosonOnTorusShort(NbrParticles, NbrFluxQuanta, TotalKy);
	  else
	    Space = new BosonOnTorusShort(NbrParticles, NbrFluxQuanta, TotalKy, MaxOccupation);
	  if (TwistedTorusFlag == true)
	    ComplexState = ComplexVector(Space->GetHilbertSpaceDimension(), true);
	  else
	    State = RealVector(Space->GetHilbertSpaceDimension(), true);
	  for (int i = 0; (i < Space->GetHilbertSpaceDimension()) && (TotalKyFlag == false); ++i)
	    {
	      TotalKyFlag = Space->HasPauliExclusions(i, PauliKValue, PauliRValue);
	      if (TotalKyFlag == true)
		{
		  cout << "find admissible configuration at " << i << endl;
		  if (TwistedTorusFlag == true)
		    {
		      Space->CreateStateFromMPSDescription(SparseComplexBMatrices, StringMatrix, ComplexState, 
							   MPSSumIndices, NbrMPSSumIndices, 1l, (long) i, 1l);
		      if ((ComplexState[i].Re == 0.0) && (ComplexState[i].Im == 0.0))
			{
			  TotalKyFlag = false;
			  cout << ", but does not lead to a non-zero coefficient" << endl;
			}
		      else
			{
			  cout << " and leads to a non-zero coefficient (" << ComplexState[i] << ")" << endl;
			  ComplexState[i] = 0.0;
			}
		    }
		  else
		    {		
		      Space->CreateStateFromMPSDescription(SparseBMatrices, StringMatrix, State, MPSSumIndices, NbrMPSSumIndices, 1l, (long) i, 1l);
		      if (State[i] == 0.0)
			{
			  TotalKyFlag = false;
			  cout << ", but does not lead to a non-zero coefficient" << endl;
			}
		      else
			{
			  cout << " and leads to a non-zero coefficient (" << State[i] << ")" << endl;
			  State[i] = 0.0;
			}
		    }
		}
	    }
	  if (TotalKyFlag == false)
	    {
	      delete Space;
	      TotalKy += ReducedBrillouinZoneSize;
	    }
	}      
      if (TotalKy >= NbrFluxQuanta)
	{
	  cout << "error, can't find any root configuration" << endl;
	  return -1;
	}
    }
  else
    {
      bool TotalKyFlag = false;
      
      while ((TotalKyFlag == false) && (TotalKy < NbrFluxQuanta))
	{
	  if (Manager.GetInteger("ky-momentum") >= 0)
	    {
	      TotalKy = Manager.GetInteger("ky-momentum");
	      TotalKyFlag = true;
	    }
#ifdef __64_BITS__
	  if (NbrFluxQuanta <= 62)
#else
	    if (NbrFluxQuanta <= 30)
#endif
	      {
		Space = new FermionOnTorus(NbrParticles, NbrFluxQuanta, TotalKy);
	      }
	    else
	      {
#ifdef __128_BIT_LONGLONG__
		if (NbrFluxQuanta <= 126)
#else
		  if (NbrFluxQuanta <= 62)
#endif
		    {
		      Space = 0;//new FermionOnSpherePTruncatedLong(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
		    }
		  else
		    {
#ifdef __128_BIT_LONGLONG__
		      cout << "cannot generate an Hilbert space when nbr-flux > 126" << endl;
#else
		      cout << "cannot generate an Hilbert space when nbr-flux > 62" << endl;
#endif
		      return 0;
		    }
	      }
	  if (TwistedTorusFlag == true)
	    ComplexState = ComplexVector(Space->GetHilbertSpaceDimension(), true);
	  else
	    State = RealVector(Space->GetHilbertSpaceDimension(), true);
	  for (int i = 0; (i < Space->GetHilbertSpaceDimension()) && (TotalKyFlag == false); ++i)
	    {
	      TotalKyFlag = Space->HasPauliExclusions(i, PauliKValue, PauliRValue);
	      if (TotalKyFlag == true)
		{
		  cout << "find admissible configuration at " << i;
		  if (TwistedTorusFlag == true)
		    {
		      Space->CreateStateFromMPSDescription(SparseComplexBMatrices, StringMatrix, ComplexState, MPSSumIndices, NbrMPSSumIndices, 1l, (long) i, 1l);
		      Space->PrintState(cout, i) << endl;
		      if ((ComplexState[i].Re == 0.0) && (ComplexState[i].Im == 0.0))
			{
			  TotalKyFlag = false;
			  cout << ", but does not lead to a non-zero coefficient" << endl;
			}
		      else
			{
			  cout << " and leads to a non-zero coefficient (" << ComplexState[i] << ")" << endl;
			  ComplexState[i] = 0.0;
			}
		    }
		  else
		    {		
		      Space->CreateStateFromMPSDescription(SparseBMatrices, StringMatrix, State, MPSSumIndices, NbrMPSSumIndices, 1l, (long) i, 1l);
		      Space->PrintState(cout, i) << endl;
		      if (State[i] == 0.0)
			{
			  TotalKyFlag = false;
			  cout << ", but does not lead to a non-zero coefficient" << endl;
			}
		      else
			{
			  cout << " and leads to a non-zero coefficient (" << State[i] << ")" << endl;
			  State[i] = 0.0;
			}
		    }
		}
	    }
	  if (TotalKyFlag == false)
	    {
	      delete Space;
	      TotalKy += ReducedBrillouinZoneSize;
	    }
	}      
      if (TotalKy >= NbrFluxQuanta)
	{
	  cout << "error, can't find any root configuration" << endl;
	  return -1;
	}
    }
  
  cout << "Hilbert space dimension : " << Space->GetLargeHilbertSpaceDimension() << endl;

  if (Manager.GetString("alternate-output") == 0)
    {
      if (Manager.GetBoolean("boson") == true)
	{
	  if (Manager.GetBoolean("txt-output"))
	    {
	      OutputTxtFileName = new char [512 + strlen(MPSMatrix->GetName())];
	      if (Manager.GetDouble("angle") == 0.0)
		{
		  sprintf (OutputTxtFileName , "bosons_torus_kysym_mps_plevel_%ld_maxocc_%ld_%s_n_%d_2s_%d_ratio_%.6f_ky_%d.0.vec.txt", 
			   Manager.GetInteger("p-truncation"), Manager.GetInteger("boson-truncation"), MPSMatrix->GetName(), NbrParticles, NbrFluxQuanta, 
			   AspectRatio, TotalKy);
		}
	      else
		{
		  sprintf (OutputTxtFileName , "bosons_torus_kysym_mps_plevel_%ld_maxocc_%ld_%s_n_%d_2s_%d_ratio_%.6f_angle_%.6f_ky_%d.0.vec.txt", 
			   Manager.GetInteger("p-truncation"), Manager.GetInteger("boson-truncation"), MPSMatrix->GetName(), NbrParticles, NbrFluxQuanta, 
			   AspectRatio, (M_PI * Manager.GetDouble("angle")), TotalKy);
		}
	    }
	  else
	    {
	      OutputFileName = new char [512 + strlen(MPSMatrix->GetName())];
	      if (Manager.GetDouble("angle") == 0.0)
		{
		  sprintf (OutputFileName , "bosons_torus_kysym_mps_plevel_%ld_maxocc_%ld_%s_n_%d_2s_%d_ratio_%.6f_ky_%d.0.vec", 
			   Manager.GetInteger("p-truncation"), Manager.GetInteger("boson-truncation"), MPSMatrix->GetName(), NbrParticles, NbrFluxQuanta, AspectRatio, TotalKy);
		}
	      else
		{
		  sprintf (OutputFileName , "bosons_torus_kysym_mps_plevel_%ld_maxocc_%ld_%s_n_%d_2s_%d_ratio_%.6f_angle_%.6f_ky_%d.0.vec", 
			   Manager.GetInteger("p-truncation"), Manager.GetInteger("boson-truncation"), MPSMatrix->GetName(), NbrParticles, NbrFluxQuanta, 
			   AspectRatio, (M_PI * Manager.GetDouble("angle")), TotalKy);
		}
	    }
	}
      else
	{
	  if (Manager.GetBoolean("txt-output"))
	    {
	      OutputTxtFileName = new char [512 + strlen(MPSMatrix->GetName())];
	      if (Manager.GetDouble("angle") == 0.0)
		{
		  sprintf (OutputTxtFileName, "fermions_torus_kysym_mps_plevel_%ld_%s_n_%d_2s_%d_ratio_%.6f_ky_%d.0.vec.txt", Manager.GetInteger("p-truncation"), MPSMatrix->GetName(), 
			   NbrParticles, NbrFluxQuanta, AspectRatio, TotalKy);
		}
	      else
		{
		  sprintf (OutputTxtFileName, "fermions_torus_kysym_mps_plevel_%ld_%s_n_%d_2s_%d_ratio_%.6f_angle_%.6f_ky_%d.0.vec.txt", Manager.GetInteger("p-truncation"), MPSMatrix->GetName(), 
			   NbrParticles, NbrFluxQuanta, AspectRatio, (M_PI * Manager.GetDouble("angle")), TotalKy);
		}
	    }
	  else
	    {
	      OutputFileName = new char [512 + strlen(MPSMatrix->GetName())];
	      if (Manager.GetDouble("angle") == 0.0)
		{
		  sprintf (OutputFileName, "fermions_torus_kysym_mps_plevel_%ld_%s_n_%d_2s_%d_ratio_%.6f_ky_%d.0.vec", Manager.GetInteger("p-truncation"), MPSMatrix->GetName(), 
			   NbrParticles, NbrFluxQuanta, AspectRatio, TotalKy);
		}
	      else
		{
		  sprintf (OutputFileName, "fermions_torus_kysym_mps_plevel_%ld_%s_n_%d_2s_%d_ratio_%.6f_angle_%.6f_ky_%d.0.vec", Manager.GetInteger("p-truncation"), MPSMatrix->GetName(), 
			   NbrParticles, NbrFluxQuanta, AspectRatio, (M_PI * Manager.GetDouble("angle")), TotalKy);
		}
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("txt-output"))
	{
	  OutputTxtFileName = new char [strlen(Manager.GetString("alternate-output")) + 1];
	  strcpy (OutputTxtFileName, Manager.GetString("alternate-output"));
	}
      else
	{
	  OutputFileName = new char [strlen(Manager.GetString("alternate-output")) + 1];
	  strcpy (OutputFileName, Manager.GetString("alternate-output"));
	}
    }


  if (TwistedTorusFlag == true)
    {
       FQHEMPSCreateStateOperation Operation(Space, SparseComplexBMatrices, StringMatrix, &ComplexState, MPSSumIndices, NbrMPSSumIndices,
					     Manager.GetInteger("precalculation-blocksize"));
       Operation.ApplyOperation(Architecture.GetArchitecture());
    }
  else
    {
       FQHEMPSCreateStateOperation Operation(Space, SparseBMatrices, StringMatrix, &State, MPSSumIndices, NbrMPSSumIndices,
  					    Manager.GetInteger("precalculation-blocksize"));
       Operation.ApplyOperation(Architecture.GetArchitecture());
    }

  if (Architecture.GetArchitecture()->CanWriteOnDisk() == true)
    {
      if (Manager.GetBoolean("no-normalization") == false)
	{
	  if (TwistedTorusFlag == true)
	    {
	      ComplexState /= ComplexState.Norm();
	    }
	  else
	    {
	      State /= State.Norm();
	    }
	}
      if (TwistedTorusFlag == true)
	{
	  if (OutputTxtFileName != 0)
	    {
	      ofstream File;
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);	
	      for (long i = 0; i < Space->GetLargeHilbertSpaceDimension(); ++i)
		{
		  ComplexState.PrintComponent(File, i) << " ";
		  Space->PrintState(File, i) << endl;
		}
	      File.close();
	    }
	  if (OutputFileName != 0)
	    {
	      ComplexState.WriteVector(OutputFileName);
	    }
	}
      else
	{
	  if (OutputTxtFileName != 0)
	    {
	      ofstream File;
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);	
	      for (long i = 0; i < Space->GetLargeHilbertSpaceDimension(); ++i)
		{
		  State.PrintComponent(File, i) << " ";
		  Space->PrintState(File, i) << endl;
		}
	      File.close();
	    }
	  if (OutputFileName != 0)
	    {
	      State.WriteVector(OutputFileName);
	    }
	}
    }
  
  return 0;
}

