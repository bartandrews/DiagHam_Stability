#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneLargeBasis.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSzSymmetry.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>



using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::ifstream;



int main(int argc, char** argv)
{

  OptionManager Manager ("FQHESphereWithSU2SpinEntanglementEntropy" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "name of the file describing the system ground state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-la", "minimum size of the subsystem whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-la", "maximum size of the subsystem whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-alt", "use alternative Hilbert space for  bosonic states");
  (*SystemGroup) += new BooleanOption  ('\n', "spinup-spindown", "use spin up/spin down separation instead of the orbital bipartite cut");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the reduced density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "na-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed number of particles", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "lza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Lz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "sza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Sz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
  (*ToolsGroup) += new SingleDoubleOption  ('\n', "diag-precision", "convergence precision in non LAPACK mode", 1e-7);
#endif
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSU2SpinEntanglementEntropy -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  bool BipartiteFlag = !(Manager.GetBoolean("spinup-spindown"));
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  int FilterNa = Manager.GetInteger("na-eigenstate");
  int FilterLza = Manager.GetInteger("lza-eigenstate");
  int FilterSza = Manager.GetInteger("sza-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");

  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  char* FileName = Manager.GetString("input-file");
  if (FileName == 0)
    {
      cout << " an input file has to be provided" << endl;
      return -1;
    }
  bool SVDFlag = Manager.GetBoolean("use-svd");

  int NbrParticles=0;
  int LzMax=0;
  int TotalLz=0;
  int TotalSz=0;
  int LzSymmetry = 0;
  int SzSymmetry = 0;
  bool Statistics = true;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(FileName, NbrParticles, LzMax, TotalLz, TotalSz, LzSymmetry, SzSymmetry, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << FileName << endl;
      return -1;
    }
  int NbrParticlesUp = (NbrParticles + TotalSz) >> 1;
  int NbrParticlesDown = (NbrParticles - TotalSz) >> 1;
  ParticleOnSphereWithSpin* Space = 0;
  char* StatisticPrefix = new char[16];
  if (Statistics == true)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }

  if (Statistics == true)
    {  
      if (Manager.GetBoolean("haldane") == false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	    if (LzMax <= 15)
#endif
	      {
		Space = new FermionOnSphereWithSpin  (NbrParticles, TotalLz, LzMax, TotalSz);
	      }
	    else
	      {
#ifdef __128_BIT_LONGLONG__
		if (LzMax <= 63)
#else
		  if (LzMax <= 31)
#endif
		    {
		      Space = new FermionOnSphereWithSpinLong (NbrParticles, TotalLz, LzMax, TotalSz);
		    }
		  else
		    {
		      cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		      return 0;
		    }	
	      }
	}
      else
	{
	  int** ReferenceStates = 0;
	  int NbrReferenceStates;
	  if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
	    {
	      cout << "error, a reference file is needed" << endl;
	      return 0;
	    }
	  if (FQHEGetRootPartitionSU2(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceStates, NbrReferenceStates) == false)
	    {
	      cout << "error while parsing " << Manager.GetString("reference-file") << endl;	      
	      return 0;
	    }
	  Space = new FermionOnSphereWithSpinHaldaneBasis(NbrParticles, TotalLz, LzMax, TotalSz, ReferenceStates, NbrReferenceStates); 
	}
    }
  else
    {
      if (Manager.GetBoolean("haldane") == false)
	{
	  if (Manager.GetBoolean("use-alt") == false)
	    {
	      Space = new BosonOnSphereWithSpin (NbrParticles, TotalLz, LzMax, TotalSz);
	    }
	  else
	    {
	      if (LzSymmetry == 0)
		{
		  if (SzSymmetry == 0)
		    {		      
		      if (Manager.GetString("load-hilbert") == 0)
			{
			  Space = new BosonOnSphereWithSU2Spin (NbrParticles, TotalLz, LzMax, TotalSz);
			}
		      else
			{
			  Space = new BosonOnSphereWithSU2Spin (Manager.GetString("load-hilbert"));
			}
		    }
		  else
		    {		      
		      if (Manager.GetString("load-hilbert") == 0)
			{
			  Space = new BosonOnSphereWithSU2SpinSzSymmetry (NbrParticles, TotalLz, LzMax, TotalSz, (SzSymmetry == -1));
			}
		      else
			{
			  Space = new BosonOnSphereWithSU2SpinSzSymmetry (Manager.GetString("load-hilbert"));
			}
		    }		  
		}
	      else
		{
		  if (SzSymmetry == 0)
		    {		      
		      if (Manager.GetString("load-hilbert") == 0)
			{
			  Space = new BosonOnSphereWithSU2SpinLzSymmetry (NbrParticles, LzMax, TotalSz, (LzSymmetry == -1));
			}
		      else
			{
			  Space = new BosonOnSphereWithSU2SpinLzSymmetry (Manager.GetString("load-hilbert"));
			}
		    }
		  else
		    {		      
		      if (Manager.GetString("load-hilbert") == 0)
			{
			  Space = new BosonOnSphereWithSU2SpinLzSzSymmetry (NbrParticles, LzMax, TotalSz, (SzSymmetry == -1), (LzSymmetry == -1));
			}
		      else
			{
			  Space = new BosonOnSphereWithSU2SpinLzSzSymmetry (Manager.GetString("load-hilbert"));
			}
		    }
		}
	    }
	}
      else
	{
	  cout << "squeezed Hilbert space is not available for bosons" << endl;
	  return 0;
	}
    }
  
  RealVector GroundState;
  if (GroundState.ReadVector (FileName) == false)
    {
      cout << "can't open vector file " << FileName << endl;
      return -1;      
    }   
  
  if (Space->GetHilbertSpaceDimension() != GroundState.GetVectorDimension())
    {
      cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
      return -1;
    }


  if (BipartiteFlag == false)
    {
      if (Statistics == true)
	{ 
	  cout << "spin Up - spin down separation entropy is not supported for bosons " << endl;
	  return 0;
	} 
      int NbrParticlesUp= (NbrParticles+TotalSz)>>1;
      if (DensityMatrixFileName != 0)
	{
	  ofstream DensityMatrixFile;
	  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
	  DensityMatrixFile << "# Lz    lambda" << endl;
	  DensityMatrixFile.close();
	}
      
      double EntanglementEntropy = 0.0;
      double DensitySum =0.0;
      int MaxLzUp = (LzMax*NbrParticlesUp-NbrParticlesUp*(NbrParticlesUp-1)) + 1;
      for(int lzUp= (-LzMax*NbrParticlesUp+NbrParticlesUp*(NbrParticlesUp-1)); lzUp < MaxLzUp; lzUp+=2)
	{
	  RealSymmetricMatrix PartialDensityMatrix;
	  PartialDensityMatrix= Space->EvaluatePartialDensityMatrixSpinSeparation(lzUp,GroundState);
	  
	  if (PartialDensityMatrix.GetNbrRow() > 1)
	    {
	      RealDiagonalMatrix TmpDiag(PartialDensityMatrix.GetNbrRow());
	      PartialDensityMatrix.Diagonalize(TmpDiag);
	      TmpDiag.SortMatrixDownOrder();
	      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		{
		  if (TmpDiag[i] > 1e-14)
		    {
		      EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
		      DensitySum += TmpDiag[i];
		    }
		}
	      if (DensityMatrixFileName != 0)
		{
		  ofstream DensityMatrixFile;
		  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		  DensityMatrixFile.precision(14);
		  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		    DensityMatrixFile << (0.5 * ((double) lzUp)) << " " << TmpDiag[i] << endl;
		  DensityMatrixFile.close();
		}
	    }
	  else 
	    {
	      if (PartialDensityMatrix.GetNbrRow() == 1)
		{
		  if (PartialDensityMatrix(0,0) > 1e-14)
		    {
		      EntanglementEntropy += PartialDensityMatrix(0,0) * log(PartialDensityMatrix(0,0));
		      DensitySum += PartialDensityMatrix(0,0);
		      if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  DensityMatrixFile << (0.5 * ((double) lzUp)) << " " << PartialDensityMatrix(0,0) << endl;
			  DensityMatrixFile.close();
			}
		    }
		}
	    }
	}      
      cout << (-EntanglementEntropy) << " " << DensitySum;
    }
  else //standard orbital partitioning
    {
      if (DensityMatrixFileName != 0)
	{
	  ofstream DensityMatrixFile;
	  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
	  DensityMatrixFile << "# l_a    N    Lz    Sz    lambda" << endl;
	  DensityMatrixFile.close();
	}
      ofstream File;
      if (Manager.GetString("output-file") != 0)
	File.open(Manager.GetString("output-file"), ios::binary | ios::out);
      else
	{
	  char* TmpFileName;
	  TmpFileName = ReplaceExtensionToFileName(FileName, "vec", "ent");
	  if (TmpFileName == 0)
	    {
	      cout << "no vec extension was find in " << FileName << " file name" << endl;
	      return 0;
	    }
	  File.open(TmpFileName, ios::binary | ios::out);
	  delete[] TmpFileName;
	}
      File.precision(14);
      cout.precision(14);
      int MeanSubsystemSize = LzMax >> 1;
      if ((LzMax & 1) != 0)
	++MeanSubsystemSize;
      if (Manager.GetInteger("max-la") > 0)
	{
	  MeanSubsystemSize = Manager.GetInteger("max-la");
	  if (MeanSubsystemSize > LzMax)
	    MeanSubsystemSize = LzMax;
	}
      int SubsystemSize = Manager.GetInteger("min-la");
      if (SubsystemSize < 1)
	SubsystemSize = 1;
      
      for (; SubsystemSize <= MeanSubsystemSize; ++SubsystemSize)
	{
	  int ComplementarySubsystemSize = LzMax + 1 - SubsystemSize;
	  double EntanglementEntropy = 0.0;
	  double DensitySum = 0.0;
	  int MaxSubsystemNbrParticles = NbrParticles;
	  int SubsystemNbrParticles = 0;
	  if (Statistics == true)
	    { 
	      if (MaxSubsystemNbrParticles > (2 * SubsystemSize))
		MaxSubsystemNbrParticles = 2 * SubsystemSize;
	      SubsystemNbrParticles = NbrParticles - 2 * (LzMax + 1 - SubsystemSize);
	      if (SubsystemNbrParticles < 0)
		SubsystemNbrParticles = 0;
	    }
	  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	    {
	      int SubsystemTotalSz = 0;
	      int SubsystemMaxTotalSz = SubsystemNbrParticles;
	      SubsystemTotalSz = -SubsystemNbrParticles; 
	      for (; SubsystemTotalSz <= SubsystemMaxTotalSz; SubsystemTotalSz += 2)
		{
		  int SubsystemTotalLz = 0;
		  int SubsystemLzMax = SubsystemSize - 1;
		  int SubsystemNbrParticlesUp = (SubsystemNbrParticles + SubsystemTotalSz) >> 1;
		  int SubsystemNbrParticlesDown = (SubsystemNbrParticles - SubsystemTotalSz) >> 1;
		  int ComplementarySubsystemNbrParticlesUp = NbrParticlesUp - SubsystemNbrParticlesUp;
		  int ComplementarySubsystemNbrParticlesDown = NbrParticlesDown - SubsystemNbrParticlesDown;
		  if ((Statistics == false) || (((SubsystemNbrParticlesUp <= SubsystemSize) && (SubsystemNbrParticlesDown <= SubsystemSize) &&
						 (SubsystemNbrParticlesUp >= 0) && (SubsystemNbrParticlesDown >= 0) &&
						 (ComplementarySubsystemNbrParticlesUp <= ComplementarySubsystemSize) && 
						 (ComplementarySubsystemNbrParticlesDown <= ComplementarySubsystemSize) &&
						 (ComplementarySubsystemNbrParticlesUp >= 0) && (ComplementarySubsystemNbrParticlesDown >= 0))))
		    {
		      int SubsystemMaxTotalLz = 0;
		      int ComplementarySubsystemMinTotalLz = 0;
		      int ComplementarySubsystemMaxTotalLz = 0;
		      if (Statistics == false)
			{
			  SubsystemMaxTotalLz = SubsystemNbrParticles * SubsystemLzMax;
			  ComplementarySubsystemMinTotalLz = 0;
			  ComplementarySubsystemMaxTotalLz = (ComplementarySubsystemNbrParticlesUp + ComplementarySubsystemNbrParticlesDown) * (ComplementarySubsystemSize - 1);
			  ComplementarySubsystemMinTotalLz += (ComplementarySubsystemNbrParticlesUp + ComplementarySubsystemNbrParticlesDown) * SubsystemSize;
			  ComplementarySubsystemMaxTotalLz += (ComplementarySubsystemNbrParticlesUp + ComplementarySubsystemNbrParticlesDown) * SubsystemSize;
			  SubsystemTotalLz = 0;
			}
		      else
			{
			  SubsystemTotalLz = ((SubsystemNbrParticlesUp * (SubsystemNbrParticlesUp - 1))
					      + (SubsystemNbrParticlesDown * (SubsystemNbrParticlesDown - 1))) >> 1; 
			  SubsystemMaxTotalLz = ((SubsystemNbrParticlesUp + SubsystemNbrParticlesDown) * SubsystemLzMax) - SubsystemTotalLz;
			  ComplementarySubsystemMinTotalLz = ((ComplementarySubsystemNbrParticlesUp * (ComplementarySubsystemNbrParticlesUp - 1))
								  + (ComplementarySubsystemNbrParticlesDown * (ComplementarySubsystemNbrParticlesDown - 1))) >> 1;
			  ComplementarySubsystemMaxTotalLz = ((ComplementarySubsystemSize - 1) * (ComplementarySubsystemNbrParticlesUp + ComplementarySubsystemNbrParticlesDown)) - ComplementarySubsystemMinTotalLz;
			  
			  ComplementarySubsystemMinTotalLz += (ComplementarySubsystemNbrParticlesUp + ComplementarySubsystemNbrParticlesDown) * SubsystemSize;
			  ComplementarySubsystemMaxTotalLz += (ComplementarySubsystemNbrParticlesUp + ComplementarySubsystemNbrParticlesDown) * SubsystemSize;
			}
		      int ShiftedTotalLz = (TotalLz + (NbrParticles * LzMax)) >> 1;
		      for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz++)
			{
			  int SubsystemTrueTotalLz = ((SubsystemTotalLz << 1) - (SubsystemNbrParticles * SubsystemLzMax));
			  if (((ShiftedTotalLz - SubsystemTotalLz) <= ComplementarySubsystemMaxTotalLz) &&
			      ((ShiftedTotalLz - SubsystemTotalLz) >= ComplementarySubsystemMinTotalLz) && 
			      ((EigenstateFlag == false) || ((FilterNa == SubsystemNbrParticles) && (FilterLza == SubsystemTrueTotalLz) && (FilterSza == SubsystemTotalSz))))
			    {
			      cout << "processing subsystem size=" << SubsystemSize << "  subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTrueTotalLz << " subsystem total Sz=" << SubsystemTotalSz << endl;

                              RealSymmetricMatrix PartialDensityMatrix;
                              RealMatrix PartialEntanglementMatrix; 
                              if (SVDFlag == false)
                                 {
 			            RealSymmetricMatrix TmpPartialDensityMatrix = Space->EvaluatePartialDensityMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTrueTotalLz, SubsystemTotalSz, GroundState);
			            if (PartialDensityMatrix.GetNbrRow() == 0)
				       PartialDensityMatrix = TmpPartialDensityMatrix;
			            else
				       PartialDensityMatrix += TmpPartialDensityMatrix;
                                 }
                              else
                                 {
                                    RealMatrix TmpPartialEntanglementMatrix = Space->EvaluatePartialEntanglementMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTrueTotalLz, SubsystemTotalSz, GroundState);
			            if (PartialEntanglementMatrix.GetNbrRow() == 0)
				       PartialEntanglementMatrix = TmpPartialEntanglementMatrix;
			            else
				       PartialEntanglementMatrix += TmpPartialEntanglementMatrix;				
                                 }

                              if ((PartialDensityMatrix.GetNbrRow() > 1) || ((PartialEntanglementMatrix.GetNbrRow() >= 1) && (PartialEntanglementMatrix.GetNbrColumn() >= 1)))
				{
				  RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
		                  if (SVDFlag == false)
				    {
#ifdef __LAPACK__
				      if (LapackFlag == true)
					{
					  if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
					      && (FilterLza == SubsystemTrueTotalLz) && (FilterSza == SubsystemTotalSz))
					    {
					      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
								    PartialDensityMatrix.GetNbrRow(), true);
					      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
						TmpEigenstates[i][i] = 1.0;
					      PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
					      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
					      char* TmpEigenstateName = new char[512];
					      int MaxNbrEigenstates = NbrEigenstates;
					      if (NbrEigenstates == 0)
						MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
					      for (int i = 0; i < MaxNbrEigenstates; ++i)
						{
						  if (TmpDiag[i] > 1e-14)
						    {
						      sprintf (TmpEigenstateName,
							       "%s_sphere_su2_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d_sza_%d_.%d.vec",
							       StatisticPrefix, NbrParticles, LzMax, TotalLz, SubsystemSize,
							       SubsystemNbrParticles, SubsystemTrueTotalLz, SubsystemTotalSz, i);
						      TmpEigenstates[i].WriteVector(TmpEigenstateName);
						    }
						}
					      delete[] TmpEigenstateName;
					    }
					  else
					    {
					      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
					      TmpDiag.SortMatrixDownOrder();
					    }
					}
				      else
					{
					  if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
					      && (FilterLza == SubsystemTrueTotalLz ) && (FilterSza == SubsystemTotalSz))
					    {
					      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
									PartialDensityMatrix.GetNbrRow(), true);
					      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
						TmpEigenstates[i][i] = 1.0;
					      PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
					      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
					      char* TmpEigenstateName = new char[512];
					      int MaxNbrEigenstates = NbrEigenstates;
					      if (NbrEigenstates == 0)
						MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
					      for (int i = 0; i < MaxNbrEigenstates; ++i)
						{
						  if (TmpDiag[i] > 1e-14)
						    {
						      sprintf (TmpEigenstateName,
							       "%s_sphere_su2_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d_sza_%d.%d.vec",
							       StatisticPrefix, NbrParticles, LzMax, TotalLz, SubsystemSize,
							       SubsystemNbrParticles, SubsystemTrueTotalLz, SubsystemTotalSz, i);
						      TmpEigenstates[i].WriteVector(TmpEigenstateName);
						    }
						}
					      delete[] TmpEigenstateName;
					    }
					  else
					    {
					      PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
					      TmpDiag.SortMatrixDownOrder();
					    }
					}
#else
				      if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
				      && (FilterLza == SubsystemTrueTotalLz) && (FilterSza == SubsystemTotalSz))
					{
					  if (PartialDensityMatrix.GetNbrRow() == 1)
					    {
					      PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
					      TmpDiag.SortMatrixDownOrder();
					    }
					  else
					    {
					      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
									PartialDensityMatrix.GetNbrRow(), true);
					      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
						TmpEigenstates[i][i] = 1.0;
					      PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
					      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
					      char* TmpEigenstateName = new char[512];
					      int MaxNbrEigenstates = NbrEigenstates;
					      if (NbrEigenstates == 0)
						MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
					      for (int i = 0; i < MaxNbrEigenstates; ++i)
						{
						  if (TmpDiag[i] > 1e-14)
						    {
						      sprintf (TmpEigenstateName,
							       "%s_sphere_su2_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d_sza_%d.%d.vec",
							       StatisticPrefix, NbrParticles, LzMax, TotalLz, SubsystemSize,
							       SubsystemNbrParticles, SubsystemTrueTotalLz, SubsystemTotalSz, i);
						      TmpEigenstates[i].WriteVector(TmpEigenstateName);
						    }
						}
					      delete[] TmpEigenstateName;
					    }
					}
				      else
					{
					  PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
					  TmpDiag.SortMatrixDownOrder();
				    }
#endif		  
				      
				    } //SVD
				  else
				    {
				      cout<<"Using SVD. "<<endl;  
				      if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
					{	
					  cout << "PartialEntanglementMatrix = " << PartialEntanglementMatrix.GetNbrRow() << " x " << PartialEntanglementMatrix.GetNbrColumn() << endl;
					  double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
					  int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
					  if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
					    {
					      TmpDimension = PartialEntanglementMatrix.GetNbrRow();
					    }
					  for (int i = 0; i < TmpDimension; ++i)
					    {
					      TmpValues[i] *= TmpValues[i];
					    }
					  TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);
					  TmpDiag.SortMatrixDownOrder();
					}
				      else
					{
					  double TmpValue = 0.0;
					  if (PartialEntanglementMatrix.GetNbrRow() == 1)
					    {
					      for (int i = 0; i < PartialEntanglementMatrix.GetNbrColumn(); ++i)
						TmpValue += PartialEntanglementMatrix[i][0] * PartialEntanglementMatrix[i][0];
					    }
					  else
					    {
					      for (int i = 0; i < PartialEntanglementMatrix.GetNbrRow(); ++i)
						TmpValue += PartialEntanglementMatrix[0][i] * PartialEntanglementMatrix[0][i];				  
					    }
					  TmpDiag = RealDiagonalMatrix(1, 1);
					  TmpDiag[0] = TmpValue;
					}
				    }   
				  
				  
				  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
					  DensitySum += TmpDiag[i];
					}
				    }
				  if (DensityMatrixFileName != 0)
				    {
				      ofstream DensityMatrixFile;
				      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
				      DensityMatrixFile.precision(14);
				      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
					DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTrueTotalLz << " " << SubsystemTotalSz << " " << TmpDiag[i] << endl;
				      DensityMatrixFile.close();
				    }
				}
			      else
				{
				  if (PartialDensityMatrix.GetNbrRow() == 1)
				    {
				      double TmpValue = PartialDensityMatrix(0,0);
				      if (TmpValue > 1e-14)
					{
					  EntanglementEntropy += TmpValue * log(TmpValue);
					  DensitySum += TmpValue;
					}
				      if (DensityMatrixFileName != 0)
					{
					  ofstream DensityMatrixFile;
					  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
					  DensityMatrixFile.precision(14);
					  DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTrueTotalLz << " " << SubsystemTotalSz << " " << TmpValue << endl;
					  DensityMatrixFile.close();
					}		  
				    }
				}
			    }
			}
		    }
		}
	    }
	  File << SubsystemSize << " " << (-EntanglementEntropy) << " " << DensitySum << endl;
	}
    }
  return 0;
}
