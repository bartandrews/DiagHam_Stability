#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHESpectrum/LatticePhases.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"

#include "HilbertSpace/FermionOnLattice.h"
#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/BosonOnLatticeGeneric.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLatticeGeneric.h"


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

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHELatticeEntanglementEntropyParticlePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  LatticePhases::AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum size of the particles whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  int NbrSpaces = 1;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  ParticleOnLattice** Spaces = 0;
  double Interaction=-1.0;
  int TmpI=-1;
  bool HardCore=false;
  bool GenericLattice=false;
  int NbrSites=0;
  int NbrSubLattices=1;
  bool HaveContFlux = false;
  double ContFlux =0.0;
  LatticePhases *Lattice = NULL;  
  ComplexVector* GroundStates = 0;
  char** GroundStateFiles = 0;
  int Lx = 0;
  int Ly = 0;
  int NbrFluxQuanta;
  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  bool Statistics = true;
  double* Coefficients = 0;
  bool ShowTimeFlag = Manager.GetBoolean("show-time");

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHETopInsulatorEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) && 
      (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }
  if ((Manager.GetString("degenerated-groundstate") != 0) && 
      (IsFile(Manager.GetString("degenerated-groundstate")) == false))
    {
      cout << "can't open file " << Manager.GetString("degenerated-groundstate") << endl;
      return -1;
    }


  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      Coefficients = new double[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));
      Coefficients[0] = 1.0;
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
       NbrSpaces = DegeneratedFile.GetNbrLines();
       GroundStateFiles = new char* [NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
       if (DegeneratedFile.GetNbrColumns() == 1)
	 {
	   Coefficients = new double[NbrSpaces];
	   for (int i = 0; i < NbrSpaces; ++i)
	     Coefficients[i] = 1.0 / ((double) NbrSpaces);
	 }
       else
	 {
	   double TmpSum = 0.0;
	   Coefficients = DegeneratedFile.GetAsDoubleArray(1);
	   for (int i = 0; i < NbrSpaces; ++i)
	     TmpSum += Coefficients[i];
	   TmpSum = 1.0 / TmpSum;
	   for (int i = 0; i < NbrSpaces; ++i)
	     Coefficients[i] *= TmpSum;
	 }
    }

  // read in files, and check consistency
  // first check if generic lattice is used
  char* LatticeName = NULL;
  if ((Manager.GetString("lattice-definition")!=NULL)||(FQHEOnLatticeHaveGeneralLattice(GroundStateFiles[0])))
    {
      GenericLattice=true;
      if (Manager.GetString("lattice-definition")==NULL)
	{
	  cout << "Please indicate the file with the lattice-definition for this vector"<<endl;
	  exit(1);
	}
      // get the lattice geometry
      Lattice = new LatticePhases();
      NbrSites = Lattice->GetNbrSites();
      Lx = Lattice->GetLatticeLength(0);
      Ly = Lattice->GetLatticeLength(1);
      NbrSubLattices = Lattice->GetNbrSubLattices();
      LatticeName = Lattice->GeometryString();
      if (strstr(GroundStateFiles[0], LatticeName)==0)
	{
	  cout << "The given lattice parameters do not coincide with the filename, verify lattice definition, and repetition of unit cells"<<endl;
	}
      if (FQHEOnLatticeFindSystemInfoFromGeneralVectorFileName(GroundStateFiles[0], NbrParticles, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore, HaveContFlux, ContFlux) == false)
	{
	  cout<<"Please use standard file-names, or indicate all necessary system parameters!"<<endl;
	  exit(1);
	}
    }
  else
    {
      if (FQHEOnLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrParticles, Lx, Ly, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore) == false)
	{
	  cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
	  exit(1);
	}
      NbrSites = Lx*Ly;
    }
  int TmpNbrParticles;
  bool TmpStatistics, TmpHardCore;
  for (int i = 1; i < NbrSpaces; ++i)
    {
      TmpStatistics=true;
			TmpNbrParticles =0;
      if (GenericLattice)
	{
	  if (strstr(GroundStateFiles[i], LatticeName)==0)
	    {
	      cout << "File number "<<i<<" not consistent with given lattice definition"<<endl;
	      return -1;
	    }
	  if (FQHEOnLatticeFindSystemInfoFromGeneralVectorFileName(GroundStateFiles[i], TmpNbrParticles, Interaction, NbrFluxQuanta, TmpI, TmpStatistics, TmpHardCore, HaveContFlux, ContFlux) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
      else
	{
	  if (FQHEOnLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i], TmpNbrParticles, Lx, Ly, Interaction, NbrFluxQuanta, TmpI, TmpStatistics, TmpHardCore) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
      if ((TmpNbrParticles!=NbrParticles)||(TmpHardCore!=HardCore)||(TmpStatistics!=Statistics))
	{
		cout <<"i = " <<i<<" " << TmpNbrParticles << " " <<NbrParticles<< " " <<TmpHardCore <<" "<< HardCore<<endl;
	  cout << "Error : Particle number and hardcore condition need to be the same for all vectors in multiplet."<<endl;
	  return -1;
	}
    }
  delete [] LatticeName;

  GroundStates = new ComplexVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
      {
	cout << "can't open vector file " << GroundStateFiles[i] << endl;
	return -1;      
      }

  Spaces = new ParticleOnLattice*[NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (Statistics)
	{
	  if (GenericLattice)
	    {
	      cout << "Generic lattice for fermions not yet defined!"<<endl;
	      return -1;
	    }
	  else
	    Spaces[i] = new FermionOnLattice(NbrParticles, Lx, Ly, NbrFluxQuanta, MemorySpace);
	}
      else
	{
	  if (GenericLattice)
	    {
	      if (HardCore)
		Spaces[i] = new HardCoreBosonOnLatticeGeneric(NbrParticles, Lattice, NbrFluxQuanta, MemorySpace);
	      else Spaces[i] = new BosonOnLatticeGeneric(NbrParticles, Lattice, NbrFluxQuanta, MemorySpace);
	    }
	  else
	    {
	      if (HardCore)
		Spaces[i] =new HardCoreBosonOnLattice(NbrParticles, Lx, Ly, NbrFluxQuanta, MemorySpace);
	      else Spaces[i] = new BosonOnLattice(NbrParticles, Lx, Ly, NbrFluxQuanta, MemorySpace);
	    }
	}
      if (Spaces[i]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return -1;
	}
    }


  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "#  N    lambda";
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }

  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "partent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);

  int MaxSubsystemNbrParticles = (NbrParticles >> 1) + (NbrParticles & 1);
  if (Manager.GetInteger("max-na") > 0)
    MaxSubsystemNbrParticles = Manager.GetInteger("max-na");
  int SubsystemNbrParticles = Manager.GetInteger("min-na");

  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
//       for (int SubsystemTotalKx = 0; SubsystemTotalKx < NbrSiteX; ++SubsystemTotalKx)
// 	{
// 	  for (int SubsystemTotalKy = 0; SubsystemTotalKy < NbrSiteY; ++SubsystemTotalKy)
// 	    {
	  
      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << endl; // << " subsystem total Kx=" << SubsystemTotalKx << " Ky=" << SubsystemTotalKy << endl;
	      timeval TotalStartingTime;
	      timeval TotalEndingTime;
	      if (ShowTimeFlag == true)
		{
		  gettimeofday (&(TotalStartingTime), 0);
		}
	      HermitianMatrix PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStates[0], Architecture.GetArchitecture());
	      PartialDensityMatrix *= Coefficients[0];
	      for (int i = 1; i < NbrSpaces; ++i)
		{
		  HermitianMatrix TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, GroundStates[i], Architecture.GetArchitecture());
		  TmpMatrix *= Coefficients[i];
		  PartialDensityMatrix += TmpMatrix;
		}
	      if (ShowTimeFlag == true)
		{
		  gettimeofday (&(TotalEndingTime), 0);
		  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		  cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
		}
	      if (PartialDensityMatrix.GetNbrRow() > 1)
		{
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalStartingTime), 0);
		    }
		  RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
		  if (LapackFlag == true)
		    {
		      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		    }
		  else
		    {
		      PartialDensityMatrix.Diagonalize(TmpDiag);
		    }
#else
		  PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
		  TmpDiag.SortMatrixDownOrder();
		  if (DensityMatrixFileName != 0)
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			DensityMatrixFile << SubsystemNbrParticles << " " << TmpDiag[i] << endl;
		      //DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpDiag[i] << endl;
		      DensityMatrixFile.close();
		    }
		  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		    {
		      if (TmpDiag[i] > 1e-14)
			{
			  EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
			  DensitySum +=TmpDiag[i];
			}
		    }
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << "diagonalization done in " << Dt << "s" << endl;
		    }
		}
	      else
		if (PartialDensityMatrix.GetNbrRow() == 1)
		  {
		    double TmpValue = PartialDensityMatrix(0,0);
		    if (DensityMatrixFileName != 0)
		      {
			ofstream DensityMatrixFile;
			DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			DensityMatrixFile.precision(14);
			//DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalKx << " " << SubsystemTotalKy << " " << TmpValue << endl;
			DensityMatrixFile << SubsystemNbrParticles << " " << TmpValue << endl;
			DensityMatrixFile.close();
		      }		  
		    if (TmpValue > 1e-14)
		      {
			EntanglementEntropy += TmpValue * log(TmpValue);
			DensitySum += TmpValue;
		      }
		  }
// 	    }
// 	}
      File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << endl;
    }
  File.close();
  return 0;
}
