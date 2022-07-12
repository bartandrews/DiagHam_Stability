#include "Options/Options.h"

#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU4SpinMomentumSpace.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong.h"



#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeTwoBandHofstadterHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFourBandHofstadterHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterTriangular.h"
#include "Tools/FTITightBinding/TightBindingInteractionCoulombForHofstadterLattice.h"
#include "Tools/FTITightBinding/TightBindingInteractionCoulomb2DEwald.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"

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
  OptionManager Manager ("FCIHofstadterModel" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  
  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-cellx", "number of unit cells along the x direction", 5);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-celly", "number of unit cells along the y direction", 1);
  
  (*SystemGroup) += new SingleIntegerOption  ('X', "unit-cellx", "number of sites in unit cell along the x direction", 1);
  (*SystemGroup) += new SingleIntegerOption  ('Y', "unit-celly", "number of sites in unit cell along the y direction", 7);
  
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux-per-cell", "number of flux quanta per unit cell", 1);
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-sz", "only evalute a given spin sector (negative if all sz sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-parity", "only evalute a given sz parity sector (0 if all sz parity sectors have to be computed)", 0);
  (*SystemGroup) += new  BooleanOption  ('\n', "redundant-kx", "Calculate all kx subspaces", false);
  (*SystemGroup) += new  BooleanOption  ('\n', "redundant-ky", "Calculate all ky subspaces", false);
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "spin", "use fermions with spin");
  (*SystemGroup) += new BooleanOption  ('\n', "triangular", "use the Hofstadter model for a triangular lattice");

  (*SystemGroup) += new BooleanOption  ('\n', "coulomb", "assume Coulomb interactions");
  (*SystemGroup) += new BooleanOption  ('\n', "generic-coulomb", "use the generic interaction class to generate Coulomb interactions");
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive onsite(boson) or NN (fermion) potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive NN(boson) or NNN (fermion) potential strength", 0.0);
  
  (*SystemGroup) += new SingleDoubleOption  ('e', "periodic-potential", "strength of an additional periodic potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "on-site chemical potential (enlarges unit cell by factor of 2 in X direction)", 0.0);
  (*SystemGroup) += new  SingleIntegerOption  ('\n', "enlarge-cellx", "enlarge unit cell by factor * in X direction without changing the flux density", 1);
  (*SystemGroup) += new  SingleIntegerOption  ('\n', "enlarge-celly", "enlarge unit cell by factor * in Y direction without changing the flux density", 1);
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "tunneling for lattice enclosing one magnetic flux per cell (period is UnitCellX / 2 * UnitCellY)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "tunnelling for lattice enclosing half a magnetic flux per cell (period is UnitCellX / 2 * UnitCellY)", 0.0);
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "delta", "on-site kinetic energy term (does not break translation symmetry) from Motruk paper", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "M", "on-site kinetic energy term (breaks translation symmetry) from Motruk paper", 0.0);
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-min", "lowest band to be populated (-1=highest band)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-max", "highest band to be populated", 0);
  
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the chern number (only in singleparticle-spectrum mode)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "singleparticle-precision", "length of floating point representations used for single-particle diagonalization, in bits", 64);
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  //  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new BooleanOption  ('\n', "landau-x", "Use Landau gauge along the x-axis within unit cell");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n', "embedding", "compute the band structure witht the embedding");
  
  (*SystemGroup) += new BooleanOption  ('\n', "real-space", "use the real space representation when considering the system with all bands");
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "use the real space representation when considering the system with all bandswithout the translations");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-szsymmetry", "disable the Sz<->-Sz symmetry");
  (*SystemGroup) += new BooleanOption  ('\n', "hardcore", "consider hardcore bosons (oly valid in real space mode)");
  
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleIntegerOption ('\n', "min-nbrsinglets", "mininum number of on-site singlet that defines the projected space", -1);  


  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleDoubleOption  ('\n',"testhermitian-error", "precision of the hermeticity test",0);


  
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrCellX = Manager.GetInteger("nbr-cellx"); 
  int NbrCellY = Manager.GetInteger("nbr-celly");
  
  int UnitCellX = Manager.GetInteger("unit-cellx"); 
  int UnitCellY = Manager.GetInteger("unit-celly");
  
  int FluxPerCell = Manager.GetInteger("flux-per-cell");
  
  int FullUnitCellX = Manager.GetInteger("enlarge-cellx") * UnitCellX;
  int FullUnitCellY = Manager.GetInteger("enlarge-celly") * UnitCellY;
  
  double PeriodicPotentialStrength = Manager.GetDouble("periodic-potential");
  double MuPotential = Manager.GetDouble("mu-s");
  double T1 = Manager.GetDouble("t1");
  double T2 = Manager.GetDouble("t2");  
  double delta = Manager.GetDouble("delta");
  double M = Manager.GetDouble("M");
  
  if ((T2 != 0) && ((UnitCellX % 2) != 0))
  {
    cout << "Error: UnitCellX must be even to add a potential with half a flux per cell" << endl;
    return -1;
  }
  
  bool EnlargeCellXFlag = false;
  bool EnlargeCellYFlag = false;
  if (Manager.GetInteger("enlarge-cellx") != 1)
    EnlargeCellXFlag = true;
  if (Manager.GetInteger("enlarge-celly") != 1)
    EnlargeCellYFlag = true;
  
  int MinBand = Manager.GetInteger("band-min");
  int MaxBand = Manager.GetInteger("band-max");
  bool EmbeddingFlag = Manager.GetBoolean("embedding");
  
  int  MinNbrSinglets = Manager.GetInteger("min-nbrsinglets");

  if ((MinNbrSinglets >= 0 )&&(Manager.GetBoolean("spin") == false) )
    {
      cout <<" One cannot constain the number of singlets if there is no spin boulet!"<<endl;
      exit(1);
    }
  
  char Axis ='y';
  
  if (Manager.GetBoolean("landau-x"))
    Axis ='x';
  
  if ((MaxBand<0)||(MaxBand >= FullUnitCellX * FullUnitCellY))
    {
      cout << "max-band out of range"<<endl;
      exit(1);
    }
  if (MinBand > MaxBand)
    {
      cout << "min-band out of range"<<endl;
      exit(1);
    }
  if (MinBand<0)
    MinBand=MaxBand;
  
  int NbrBands = MaxBand-MinBand+1;
  
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  bool SzSymmetryFlag = false;
  if ((Manager.GetBoolean("disable-szsymmetry") == false) && (Manager.GetBoolean("spin") == true))
    {
      if (Manager.GetBoolean("boson") == false)
	{
	  SzSymmetryFlag = true;
	}
      else
	{
	  cout << "warning, Sz<->-Sz is not implemented for bosons" << endl;	  
	}	
    }

  char* StatisticPrefix = new char [50];
  sprintf (StatisticPrefix, "fermions");
  
  if (Manager.GetBoolean("boson") == false)
    {
      if (Manager.GetBoolean("spin") == false)
	{
	  sprintf (StatisticPrefix, "fermions");
	}
      else
	{
	  if (MinNbrSinglets <= 0 )
	    {
	      sprintf (StatisticPrefix, "fermions_su2"); 
	    }
	  else
	    {
	      sprintf (StatisticPrefix, "fermions_su2_minnbrsinglet_%d",MinNbrSinglets); 
	    }
	}
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }
  
  
  char* FilePrefix = new char [512];
  int lenFilePrefix=0;
  
  
  if (Manager.GetBoolean("triangular")==false)
    {
      if (Manager.GetBoolean("real-space") == false)
	{
	  lenFilePrefix += sprintf (FilePrefix, "%s_hofstadter_X_%d_Y_%d_q_%d", StatisticPrefix, UnitCellX, UnitCellY, FluxPerCell);
	  
	  if ((NbrBands>1)||(MaxBand>0))
	    {
	      if (NbrBands==1)
		lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_b_%d", MaxBand);
	      else
		lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_b_%d-%d", MinBand, MaxBand);
	    }
	  if (EnlargeCellXFlag == true)
	    lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_Xeff_%d", (FullUnitCellX));
	  if (EnlargeCellYFlag == true)
	    lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_Yeff_%d", (FullUnitCellY));
	  if (Manager.GetBoolean("landau-x"))
	    lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_landau-x");
	}
      else
	{
	  if ( Manager.GetBoolean("no-translation") == false)
	    {
	      if ( Manager.GetBoolean("hardcore") == false)
		lenFilePrefix += sprintf (FilePrefix, "%s_realspace_hofstadter_X_%d_Y_%d_q_%d", StatisticPrefix, UnitCellX, UnitCellY, FluxPerCell);
	      else
		lenFilePrefix += sprintf (FilePrefix, "%s_realspace_gutzwiller_hofstadter_X_%d_Y_%d_q_%d", StatisticPrefix, UnitCellX, UnitCellY, FluxPerCell);
	      
	      if (EnlargeCellXFlag == true)
		lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_Xeff_%d", (FullUnitCellX));
	      if (EnlargeCellYFlag == true)
		lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_Yeff_%d", (FullUnitCellY));
	    }
	  else
	    {
	      if ( Manager.GetBoolean("hardcore") == false)
		lenFilePrefix += sprintf (FilePrefix, "%s_realspace_notranslation_hofstadter_X_%d_Y_%d_q_%d", StatisticPrefix, UnitCellX, UnitCellY, FluxPerCell);
	      else
		lenFilePrefix += sprintf (FilePrefix, "%s_realspace_notranslation_gutzwiller_hofstadter_X_%d_Y_%d_q_%d", StatisticPrefix, UnitCellX, UnitCellY, FluxPerCell);
	    }
	}
      
    }
  else
    {
      // only quarter flux density implemented for the moment:
      lenFilePrefix += sprintf (FilePrefix, "%s_hofstadter-tri_X_1_Y_2_q_1", StatisticPrefix);
      
      if ((NbrBands>1)||(MaxBand>0))
	{
	  if (NbrBands==1)
	    lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_b_%d", MaxBand);
	  else
	    lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_b_%d-%d", MinBand, MaxBand);
	}
      
    }
  // precision
  int Precision = Manager.GetInteger("singleparticle-precision");
  if (Precision!=64)
    {
      if (Precision <= 32)
	lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_sp-float");
      else if (Precision > 64)
	lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_sp-gmp_%d", Precision);
    }
  // common naming options:
  lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_n_%d_x_%d_y_%d", NbrParticles, NbrCellX, NbrCellY);
  
  char* CommentLine = new char [256];
  if(Manager.GetBoolean("spin") == false) 
    sprintf (CommentLine, "eigenvalues\n# kx ky ");
  else
    {
      if (SzSymmetryFlag == false)
	sprintf (CommentLine, "eigenvalues\n# sz kx ky ");
      else
	sprintf (CommentLine, "eigenvalues\n# sz szsym kx ky ");
    }
  
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    {
      strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
      delete [] FilePrefix;
      FilePrefix = RemoveExtensionFromFileName(EigenvalueOutputFile,".dat");
      if (FilePrefix==0)
	strcpy(FilePrefix, EigenvalueOutputFile);
    }
  else
    {
      if (((Manager.GetBoolean("flat-band") == false)&&(Manager.GetBoolean("hardcore") == false ))||(Manager.GetDouble("v-potential")!=0.0)||(PeriodicPotentialStrength != 0.0))
	lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_u_%g",Manager.GetDouble("u-potential"));
      if (Manager.GetDouble("v-potential")!=0.0)
	lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_v_%g",Manager.GetDouble("v-potential"));
      if (PeriodicPotentialStrength!=0.0)
	lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_epsilon_%g",PeriodicPotentialStrength);
      if (MuPotential != 0.0)
	lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_mus_%g", MuPotential);
      if (T1 != 0.0)
	lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_t1_%g", T1);
      if (T2 != 0.0)
	lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_t2_%g", T2);
      
      if (delta != 0.0)
	lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_delta_%g", M);
      if (M != 0.0)
	lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_M_%g", M);
      
      lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_gx_%g_gy_%g", Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
      if (EmbeddingFlag)
	lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_emb");
      
      sprintf (EigenvalueOutputFile,"%s.dat",FilePrefix);
    }
  
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
	ExportOneBody = true;
      
      Abstract2DTightBindingModel *TightBindingModel;
      if (Manager.GetBoolean("triangular")==false)
      {
	if ((MuPotential == 0.0) && (EnlargeCellXFlag == false) && (EnlargeCellYFlag == false) && (T1 == 0.0) && (T2 == 0.0) && (delta == 0.0) && (M == 0.0))
	  TightBindingModel= new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody, EmbeddingFlag, Precision);
	else
	  if ((EnlargeCellYFlag == false) && (delta == 0.0) && (M == 0.0))
	    TightBindingModel= new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, MuPotential, FullUnitCellX, T1, T2, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody, EmbeddingFlag, Precision);
	  else
	  {
	    if ((MuPotential != 0.0) || (T1 != 0.0) || (T2 != 0.0))
	    {
	      cout << "Set of parameters not implemented" << endl;
	      return -1;
	    }
	    TightBindingModel= new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, FullUnitCellX, FullUnitCellY, delta, M, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody, EmbeddingFlag, Precision);
	  }
      }
      else
      {
	TightBindingModel= new TightBindingModelHofstadterTriangular(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody, EmbeddingFlag);
      }
      
      TightBindingModel->WriteAsciiSpectrum(EigenvalueOutputFile);
      
      for (int n=0; n<TightBindingModel->GetNbrBands()-1; ++n)
	{
	  double BandSpread = TightBindingModel->ComputeBandSpread(n);
	  double DirectBandGap = TightBindingModel->ComputeDirectBandGap(n);
	  cout << "Spread("<<n<<") = " << BandSpread << "  Direct Gap = " << DirectBandGap  << "  Flattening = " << (BandSpread / DirectBandGap) << endl;
	  if (Manager.GetBoolean("singleparticle-chernnumber") == true)
	    {
	      cout << "Chern number("<<n<<") = " << TightBindingModel->ComputeChernNumber(n) << endl;
	    }
	}
      double BandSpread = TightBindingModel->ComputeBandSpread(TightBindingModel->GetNbrBands()-1);
      cout << "Spread("<<TightBindingModel->GetNbrBands()-1<<") = " << BandSpread << endl;
      
      if (Manager.GetBoolean("singleparticle-chernnumber") == true)
	{
	  cout << "Chern number("<<TightBindingModel->GetNbrBands()-1<<") = " << TightBindingModel->ComputeChernNumber(TightBindingModel->GetNbrBands()-1) << endl;
	}
      if (ExportOneBody == true)
	{
	  char* BandStructureOutputFile = new char [512];
	  if (Manager.GetString("export-onebodyname") != 0)
	    strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	  else
	    sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
	  if (Manager.GetBoolean("export-onebody") == true)
	    {
	      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
	    }
	  else
	    {
	      TightBindingModel->WriteBandStructureASCII(BandStructureOutputFile);
	    }
	  delete[] BandStructureOutputFile;
	}
      delete TightBindingModel;
      return 0;
    }
  
  
  int MinKx = 0;
  int MaxKx = NbrCellX - 1;
  if ((Manager.GetBoolean("redundant-kx")==false) && (fabs(Manager.GetDouble("gamma-x"))<1e-12)) // want to reduce zone, and no offset?
    MaxKx = NbrCellX/2;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  int MinKy = 0;
  int MaxKy = NbrCellY - 1;
  if ((Manager.GetBoolean("redundant-ky")==false) && (fabs(Manager.GetDouble("gamma-y"))<1e-12)) // want to reduce zone, and no offset?
    MaxKy = NbrCellY/2;
  if (Manager.GetInteger("only-ky") >= 0)
    {						
      MinKy = Manager.GetInteger("only-ky");
      MaxKy = MinKy;
    }
  
  if(Manager.GetBoolean("no-translation") == true)
    {  
      MaxKx = 0;
      MaxKy = 0;
    }

  int MinSz = 0;
  int MaxSz = 2*NbrParticles;

  if(Manager.GetBoolean("spin") == false)
    {  
      MaxSz = 0; 
    }
  
  if (Manager.GetInteger("only-sz") >= 0)
    {						
      MinSz = Manager.GetInteger("only-sz");
      MaxSz = MinSz;
    }
  

  Abstract2DTightBindingModel *TightBindingModel;
  if (Manager.GetBoolean("triangular")==false)
      {
	if ((MuPotential == 0.0) && (EnlargeCellXFlag == false) && (EnlargeCellYFlag == false) && (T1 == 0.0) && (T2 == 0.0))
	  TightBindingModel= new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true, EmbeddingFlag, Precision);
	else
	{
	  if ((EnlargeCellYFlag == false) && (delta == 0.0) && (M == 0.0))
	    TightBindingModel= new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, MuPotential, FullUnitCellX, T1, T2, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true, EmbeddingFlag, Precision);
	  else
	  {
	    if ((MuPotential != 0.0) || (T1 != 0.0) || (T2 != 0.0))
	    {
	      cout << "Set of parameters not implemented" << endl;
	      return -1;
	    }
	    TightBindingModel= new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, FullUnitCellX, FullUnitCellY, delta, M, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true, EmbeddingFlag, Precision);
	  }
	}	  
      }
  else
    TightBindingModel= new TightBindingModelHofstadterTriangular(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis,
								 Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true, EmbeddingFlag);

  if (Precision > 64) // always save tightbinding model eigenstates for expensive arbitrary precision calculations
    {
      char* BandStructureOutputFile = new char [512];
      if (Manager.GetString("export-onebodyname") != 0)
	strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
      else
	sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
      delete[] BandStructureOutputFile;
    }
  
//  HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
//   cout <<   TightBindingMatrix <<endl;
/*   RealDiagonalMatrix TmpDiag (TightBindingMatrix.GetNbrRow(),true);
   TightBindingMatrix.LapackDiagonalize(TmpDiag) ;
      for (int i = 0 ; i < TightBindingMatrix.GetNbrRow(); i++)
     cout <<TmpDiag[i]<<" ";
   cout <<endl;*/

  if (Manager.GetBoolean("boson") == false)
    {
      int FilledNbrBands=-1;
      
      double E=TightBindingModel->ComputeGroundstateEnergy(NbrParticles,FilledNbrBands, true);

      cout << "Total energy of groundstate: "<<E<<" ("<<FilledNbrBands<<" filled bands)"<<endl;
    }

  bool FirstRunFlag = true;
  for (int Sz = MinSz; Sz <= MaxSz; Sz+=2)
    {
      int SzSymmetrySector = 0;
      int MaxSzSymmetrySector = 0;
      if ((SzSymmetryFlag == true) && (Sz == 0))
	{
	  SzSymmetrySector = -1;
	  MaxSzSymmetrySector = 1;
	  if (Manager.GetInteger("sz-parity") != 0)
	  {
	    SzSymmetrySector = Manager.GetInteger("sz-parity");
	    MaxSzSymmetrySector = Manager.GetInteger("sz-parity");
	  }
	}
      for (; SzSymmetrySector <= MaxSzSymmetrySector; SzSymmetrySector += 2)
	{
	  for (int i = MinKx; i <= MaxKx; ++i)
	    {
	      for (int j = MinKy; j <= MaxKy; ++j)
		{
		  if (Manager.GetBoolean("spin") == true) 
		    {
		      if ((SzSymmetryFlag == true) && (Sz == 0))
			{
			  cout << "(kx=" << i << ",ky=" << j << "), Sz= " << Sz << ", Sz<->-Sz=" << SzSymmetrySector << " : " << endl;			  
			}
		      else
			{
			  cout << "(kx=" << i << ",ky=" << j << "), Sz= " << Sz << " : " << endl;
			}
		    }
		  else
		    {
		      cout << "(kx=" << i << ",ky=" << j << ") : " << endl;
		    }
		  ParticleOnSphere* Space = 0;
		  AbstractQHEHamiltonian* Hamiltonian = 0;
		  
		  if (Manager.GetBoolean("real-space") == false)
		    {
		      if (NbrBands==1)
			{
			  if (Manager.GetBoolean("boson") == false)
			    {
			      if ((NbrCellX * NbrCellY) <= 63)
				{
				  Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrCellX, NbrCellY, i, j);
				}
			      else
				{
				  Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrCellX, NbrCellY, i, j);
				}
			    }
			  else
			    {
			      Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrCellX, NbrCellY, i, j);
			    }
			  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
			  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
			    Memory = Architecture.GetArchitecture()->GetLocalMemory();
			  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
			  // assign Hamiltonian:
			  if (Manager.GetBoolean("coulomb") == true)
			    {
			      std::cout << "Using Coulomb interactions"<<std::endl;
			      AbstractTightBindingInteraction *CoulombInteraction; 
			      if (Manager.GetBoolean("generic-coulomb") == true)
				CoulombInteraction = new TightBindingInteractionCoulombForHofstadterLattice(TightBindingModel);
			      else
				CoulombInteraction = new TightBindingInteractionCoulomb2DEwald(TightBindingModel);
			      Hamiltonian = new ParticleOnLatticeHofstadterSingleBandGenericHamiltonian(Space, NbrParticles, NbrCellX, NbrCellY, MaxBand, CoulombInteraction, TightBindingModel, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
			    }
			  else
			  {
			    if (PeriodicPotentialStrength == 0)
			      Hamiltonian = new ParticleOnLatticeHofstadterSingleBandHamiltonian(Space, NbrParticles, NbrCellX, NbrCellY, MaxBand, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), TightBindingModel, Manager.GetBoolean("flat-band"),Architecture.GetArchitecture(), Memory);
			    else
			    {
			      double** PeriodicOneBodyPotential = new double* [NbrCellX];
			      double KxFactor = 2.0 * M_PI / ((double) (NbrCellX));
			      double KyFactor = 2.0 * M_PI / ((double) (NbrCellY));
			      for (int kx = 0; kx < NbrCellX; ++kx)
			      {
				PeriodicOneBodyPotential[kx] = new double[NbrCellY];
				for (int ky = 0; ky < NbrCellY; ++ky)
				  PeriodicOneBodyPotential[kx][ky] = - PeriodicPotentialStrength * (cos(((double) kx  + Manager.GetDouble("gamma-x"))* KxFactor) + cos(((double) ky + Manager.GetDouble("gamma-y"))* KyFactor));
			      }
			      Hamiltonian = new ParticleOnLatticeHofstadterSingleBandHamiltonian(Space, NbrParticles, NbrCellX, NbrCellY, MaxBand, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), TightBindingModel, PeriodicOneBodyPotential, Manager.GetBoolean("flat-band"),Architecture.GetArchitecture(), Memory);
			      for (int kx = 0; kx < NbrCellX; ++kx)
				delete[] PeriodicOneBodyPotential[kx];
			      delete[] PeriodicOneBodyPotential;
			    }      
			  }
			}
		      else
			{
			  if (NbrBands==2)
			    {
			      
			      if (Manager.GetBoolean("boson") == false)
				{
				  if ((NbrCellX * NbrCellY) <= 63)
				    {
				      Space = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrCellX, NbrCellY, i, j);
				    }
				  else
				    {
				      Space = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrCellX, NbrCellY, i, j);
				    }
				}
			      else
				{
				  Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrCellX, NbrCellY, i, j);
				}
			      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
			      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
				Memory = Architecture.GetArchitecture()->GetLocalMemory();
			      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
			      // assign Hamiltonian:
			      Hamiltonian = new ParticleOnLatticeTwoBandHofstadterHamiltonian((ParticleOnSphereWithSpin*)Space, NbrParticles, NbrCellX, NbrCellY, MinBand, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), TightBindingModel,Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
			    }
			  else
			    {
			      if (NbrBands==3)
				{
				  cout << "Three-band case not implemented, yet"<<endl;
				  exit(1);
				}
			      if (NbrBands==4)
				{
				  
				  if (Manager.GetBoolean("boson") == false)
				    {
				      if ((NbrCellX * NbrCellY) <= 63)
					{
					  Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrCellX, NbrCellY, i, j);
					}
				      else
					{
					  Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrCellX, NbrCellY, i, j);
					}
				    }
				  else
				    {
				      Space = new BosonOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrCellX, NbrCellY, i, j);
				    }
				  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
				  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
				    Memory = Architecture.GetArchitecture()->GetLocalMemory();
				  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
				  // assign Hamiltonian:
				  Hamiltonian = new ParticleOnLatticeFourBandHofstadterHamiltonian((ParticleOnSphereWithSU4Spin*)Space, NbrParticles, NbrCellX, NbrCellY, MinBand, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), TightBindingModel,Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
				}
			      else
				{
				  cout << "Multi-band n>4 situations not implemented, yet"<<endl;
				  exit(1);
				}
			    }
			}
		    }
		  else
		    {
		      RealSymmetricMatrix DensityDensityInteraction(TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), true);
		      if (Manager.GetBoolean("boson") == true)
			{
			  if(Manager.GetBoolean("hardcore") == false)
			    {
			      if(Manager.GetBoolean("no-translation") == true)
				Space = new BosonOnLatticeRealSpace(NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand());
			      else
				Space = new BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation(NbrParticles, UnitCellX*NbrCellX, UnitCellY*NbrCellY, i, NbrCellX, j,  NbrCellY); 
			      
			      double UPotential = Manager.GetDouble("u-potential");
			      for (int x = 0; x <  NbrCellX; ++x)
				{
				  for (int y = 0; y <  NbrCellY; ++y)
				    {
				      for (int k = 0; k < TightBindingModel->GetNbrBands(); ++k)
					{
					  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y, k),TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y, k), UPotential);
					}
				    }
				}
			    }
			  else
			    {
			      if(Manager.GetBoolean("no-translation") == true)
				{
				  Space = new BosonOnLatticeGutzwillerProjectionRealSpace(NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand());
				}		      
			      else
				{
				  if (UnitCellX*NbrCellX*  UnitCellY*NbrCellY <= 63 ) 
				    {
				      Space = new BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation(NbrParticles, UnitCellX*NbrCellX, UnitCellY*NbrCellY, i, NbrCellX, j,  NbrCellY);
				    }
				  else
				    {
				      Space = new BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong(NbrParticles, UnitCellX*NbrCellX, UnitCellY*NbrCellY, i, NbrCellX, j,  NbrCellY);
				      
				    }

				}
			    }
			  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
			  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
			  HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
			  
			  if(Manager.GetBoolean("no-translation") == false)
			    {
			      double FluxDensity =  (((double) FluxPerCell)/( (double) (UnitCellX*UnitCellY)));
			      
			      double PhaseTranslationX = 2.0* M_PI * FluxDensity * UnitCellX;
			      double PhaseTranslationY = 0.0;
			      Hamiltonian = new ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian (Space, NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), i,  NbrCellX, j,  NbrCellY, PhaseTranslationX, PhaseTranslationY,TightBindingMatrix, DensityDensityInteraction, Architecture.GetArchitecture(), Memory);
			    }
			  else
			    {
			      Hamiltonian = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), TightBindingMatrix, DensityDensityInteraction, Architecture.GetArchitecture(), Memory);
			    }
			}
		      else 
			{
			  if (Manager.GetBoolean("no-translation") == false)
			    {
			      if (Manager.GetBoolean("spin") == true) 
				{
				  if ((SzSymmetryFlag == true) && (Sz == 0))
				    {				      
#ifdef __64_BITS__
				      if (((int) TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()) < 31)
#else
					if (((int) TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()) < 15)
#endif
					  {
					    if (MinNbrSinglets <= 0 )
					      {
						Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, Sz,
															 ((int) TightBindingModel->GetNbrBands() 
															  * TightBindingModel->GetNbrStatePerBand()),
															 (SzSymmetrySector == -1),
															 i, NbrCellX, j,  NbrCellY );	  
					      }
					    else
					      {
						Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, MinNbrSinglets, Sz, 													   ((int) TightBindingModel->GetNbrBands()  * TightBindingModel->GetNbrStatePerBand()), (SzSymmetrySector == -1),   i, NbrCellX, j,  NbrCellY, 10000000);
					      } 
					  }
					else
					  {
					    if (MinNbrSinglets <= 0 )
					      {
						Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong (NbrParticles, Sz,
															 ((int) TightBindingModel->GetNbrBands() 
															  * TightBindingModel->GetNbrStatePerBand()),
															 (SzSymmetrySector == -1),
															 i, NbrCellX, j,  NbrCellY );	  
					      }
					    else
					      {
						Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, MinNbrSinglets, Sz, 													   ((int) TightBindingModel->GetNbrBands()  * TightBindingModel->GetNbrStatePerBand()), (SzSymmetrySector == -1),   i, NbrCellX, j,  NbrCellY);
						
					      }
					    
					  }
				    }
				  else
				    {
#ifdef __64_BITS__
				      if (((int) TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()) < 31)
#else
					if (((int) TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()) < 15)
#endif
					  {
					    if (MinNbrSinglets <= 0 )
					      {
						Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, Sz, 
													       ((int) TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()),
													       i, NbrCellX, j,  NbrCellY);	  
					      }
					    else
					      {
						Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, MinNbrSinglets, Sz, 													   ((int) TightBindingModel->GetNbrBands()  * TightBindingModel->GetNbrStatePerBand()), i, NbrCellX, j,  NbrCellY, 10000000ul);
					      }
					  }
					else
					  {
					    if (MinNbrSinglets <= 0 )
					      {
						Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong (NbrParticles, Sz, 
														   ((int) TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()),
														   i, NbrCellX, j,  NbrCellY);	  
						
					      }
					    else
					      {
						Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, MinNbrSinglets, Sz, 													   ((int) TightBindingModel->GetNbrBands()  * TightBindingModel->GetNbrStatePerBand()), i, NbrCellX, j,  NbrCellY, 10000000ul);
					      }
					  }
				    }
				      
				      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
					Memory = Architecture.GetArchitecture()->GetLocalMemory();
				  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
				  
				  HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
				  int NbrSites =   TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand();
				  RealSymmetricMatrix DensityDensityInteractionupup(NbrSites, true);
				  RealSymmetricMatrix DensityDensityInteractiondowndown(NbrSites, true);
				  RealSymmetricMatrix DensityDensityInteractionupdown(NbrSites, true);
				  if (Manager.GetDouble("u-potential") != 0.0)
				    {
				      double UPotential = Manager.GetDouble("u-potential");
				      for (int i = 0; i < NbrSites; ++i)
					{
					  DensityDensityInteractionupdown.SetMatrixElement(i, i, UPotential);
					}
				    }
				  
				  Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian((ParticleOnSphereWithSpin  *)  Space, NbrParticles, NbrSites, i, NbrCellX, j,  NbrCellY,
														  TightBindingMatrix, TightBindingMatrix,
														  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
														  DensityDensityInteractionupdown, 
														  Architecture.GetArchitecture(), Memory);
				  
				}
			      else
				{
				  cout <<"case not implented yet"<<endl;
				}
			    }
			  else
			    {
			      if(Manager.GetBoolean("spin") == false) 
				{
				  Space =  new FermionOnLatticeRealSpace (NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand());
				  double UPotential = Manager.GetDouble("u-potential");
				  int  NbrTranslationX;
				  int  NbrTranslationY;
				  int  NbrTranslation2X;
				  int  NbrTranslation2Y;
				  Complex TranslationPhase = 0.0;
				  for (int x = 0; x <  NbrCellX; ++x)
				    {
				      for (int y = 0; y <  NbrCellY; ++y)
					{
					  for (int X = 0; X <  UnitCellX; X++)
					    {
					      for (int Y = 0; Y <  UnitCellY; Y++)
						{
						  int OrbitalIndex =  TightBindingModel->EncodeSublatticeIndex(X,Y,NbrTranslationX,NbrTranslationY,TranslationPhase);
						  int OrbitalIndex2 =  TightBindingModel->EncodeSublatticeIndex(X+1,Y,NbrTranslation2X,NbrTranslation2Y,TranslationPhase);
						  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-NbrTranslationX, y-NbrTranslationY,OrbitalIndex),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-NbrTranslation2X, y-NbrTranslation2Y, OrbitalIndex2), UPotential);
						  OrbitalIndex2 =  TightBindingModel->EncodeSublatticeIndex(X-1,Y,NbrTranslation2X,NbrTranslation2Y,TranslationPhase);
						  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-NbrTranslationX, y-NbrTranslationY,OrbitalIndex),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-NbrTranslation2X, y-NbrTranslation2Y, OrbitalIndex2), UPotential);
						  OrbitalIndex2 =  TightBindingModel->EncodeSublatticeIndex(X,Y+1,NbrTranslation2X,NbrTranslation2Y,TranslationPhase);
						  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-NbrTranslationX, y-NbrTranslationY,OrbitalIndex),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-NbrTranslation2X, y-NbrTranslation2Y, OrbitalIndex2), UPotential);
						  OrbitalIndex2 =  TightBindingModel->EncodeSublatticeIndex(X,Y-1,NbrTranslation2X,NbrTranslation2Y,TranslationPhase);
						  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-NbrTranslationX, y-NbrTranslationY,OrbitalIndex),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-NbrTranslation2X, y-NbrTranslation2Y, OrbitalIndex2), UPotential);
						}
					    }
					}
				    }
				  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
				  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
				  HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
				  //			      cout <<"TightBindingMatrix = "<<TightBindingMatrix<<endl;
				  Hamiltonian = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), TightBindingMatrix, DensityDensityInteraction,Architecture.GetArchitecture(), Memory);
				}
			      else
				{
				  if ((SzSymmetryFlag == true) && (Sz == 0))
				    {
				      Space = new FermionOnLatticeWithSpinSzSymmetryRealSpace (NbrParticles, Sz,
											       ((int) TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()),
											       (SzSymmetrySector == -1));	  
				    }
				  else
				    {
				      Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, Sz, ((int) TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()),
										     10000000);	  
				    }
				  
				  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
				    Memory = Architecture.GetArchitecture()->GetLocalMemory();
				  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
				  
				  HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
				  int NbrSites =   TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand();
				  RealSymmetricMatrix DensityDensityInteractionupup(NbrSites, true);
				  RealSymmetricMatrix DensityDensityInteractiondowndown(NbrSites, true);
				  RealSymmetricMatrix DensityDensityInteractionupdown(NbrSites, true);
				  if (Manager.GetDouble("u-potential") != 0.0)
				    {
				      double UPotential = Manager.GetDouble("u-potential");
				      for (int i = 0; i < NbrSites; ++i)
					{
					  DensityDensityInteractionupdown.SetMatrixElement(i, i, UPotential);
					}
				    }
				  
				  Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceHamiltonian((ParticleOnSphereWithSpin  *)  Space, NbrParticles, NbrSites,
												  TightBindingMatrix, TightBindingMatrix,
												  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
												  DensityDensityInteractionupdown, 
												  Architecture.GetArchitecture(), Memory);
				}
			    }
			}
		    }
		  
		  double Shift = 0.0;
		  Hamiltonian->ShiftHamiltonian(Shift);
		  
		  if (Manager.GetString("energy-expectation") != 0 )
		    {
		      char* StateFileName = Manager.GetString("energy-expectation");
		      if (IsFile(StateFileName) == false)
			{
			  cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
			  return -1;           
			}
		      ComplexVector State;
		      if (State.ReadVector(StateFileName) == false)
			{
			  cout << "error while reading " << StateFileName << endl;
			  return -1;
			}
		      if (State.GetVectorDimension() != Space->GetHilbertSpaceDimension())
			{
			  cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
			  return -1;
			}
		      ComplexVector TmpState(Space->GetHilbertSpaceDimension());
		      VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
		      Operation.ApplyOperation(Architecture.GetArchitecture());
		      Complex EnergyValue = State * TmpState;
		      cout << "< Energy > = " << (EnergyValue.Re - Shift) << " " << EnergyValue.Im << endl;
		      return 0; 
		    }
		  char* ContentPrefix = new char[256];
		  if(Manager.GetBoolean("spin") == false) 
		    sprintf (ContentPrefix, "%d %d", i, j);
		  else
		    {
		      if (SzSymmetryFlag == true)
			{
			  sprintf (ContentPrefix, "%d %d %d %d", Sz, SzSymmetrySector, i, j);
			}
		      else
			{
			  sprintf (ContentPrefix, "%d %d %d", Sz, i, j);
			}
		    }
		  char* EigenstateOutputFile = new char [512];
		  
		  if(Manager.GetBoolean("spin") == false) 
		    {
		      sprintf (EigenstateOutputFile,"%s_kx_%d_ky_%d",FilePrefix, i, j);
		    }
		  else
		    {
		      if ((SzSymmetryFlag == true) && (Sz == 0))
			{
			  sprintf (EigenstateOutputFile,"%s_sz_%d_szsym_%d_kx_%d_ky_%d",FilePrefix, Sz, SzSymmetrySector, i, j);
			}
		      else
			{
			  sprintf (EigenstateOutputFile,"%s_sz_%d_kx_%d_ky_%d",FilePrefix, Sz, i, j);
			}
		    }
		  
		  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
		  FirstRunFlag = false;
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  cout << "------------------------------------" << endl;
		  delete Hamiltonian;
		  delete Space;
		  delete[] EigenstateOutputFile;
		  delete[] ContentPrefix;
		}
	    }
	}
    }
  delete TightBindingModel;
  return 0;
}

