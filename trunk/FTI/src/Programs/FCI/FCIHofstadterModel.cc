#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMultiBandMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU4SpinMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandHamiltonian.h"

#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeHofstadterMultiBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeTwoBandHofstadterHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFourBandHofstadterHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterTriangularQuarter.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

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
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "triangular", "use the Hofstadter model for a triangular lattice");
    
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive onsite(boson) or NN (fermion) potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "bosonic NN potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "bosonic NNN potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "three-body-potential", "bosonic three body onsite interaction", 0.0);
  //(*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "tunnelling along a selected direction (t=1 along others)", 1.0);
  (*SystemGroup) += new  SingleIntegerOption  ('\n', "enlarge-cellx", "enlarge unit cell by factor * in X direction without changing the flux density", 1);
  (*SystemGroup) += new  SingleIntegerOption  ('\n', "enlarge-celly", "enlarge unit cell by factor * in Y direction without changing the flux density", 1);

  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup)+=new SingleDoubleOption('\n',"t2","straight-line NNN hopping t2",0.0);
  (*SystemGroup)+=new SingleDoubleOption('\n',"t3","diagonal NNN hopping t3",0.0);
  (*SystemGroup)+=new SingleDoubleOption('\n',"alpha","Anisotropy factor",1.0);

  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-min", "lowest band to be populated (-1=highest band)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-max", "highest band to be populated", 0);
  
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the chern number (only in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  //  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new BooleanOption  ('\n', "landau-x", "Use Landau gauge along the x-axis within unit cell");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  (*SystemGroup) += new BooleanOption ('\n',"symmetric-gauge","use the symmetric gauge A=(By/2,-Bx/2)");
  (*SystemGroup) += new BooleanOption ('\n',"multi-band","project onto the lowest two bands");

  Manager.StandardProceedings(argv, argc, cout);

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrCellX = Manager.GetInteger("nbr-cellx"); 
  int NbrCellY = Manager.GetInteger("nbr-celly");

  int UnitCellX = Manager.GetInteger("unit-cellx"); 
  int UnitCellY = Manager.GetInteger("unit-celly");

  int FluxPerCell = Manager.GetInteger("flux-per-cell");

  int FullUnitCellX = Manager.GetInteger("enlarge-cellx") * UnitCellX;
  int FullUnitCellY = Manager.GetInteger("enlarge-celly") * UnitCellY;

  bool SymmetricGauge = Manager.GetBoolean("symmetric-gauge");

  int MinBand = Manager.GetInteger("band-min");
  int MaxBand = Manager.GetInteger("band-max");
  bool MultiBand = Manager.GetBoolean("multi-band");
  double ThreeBody = Manager.GetDouble("three-body-potential");
  double TTwo = Manager.GetDouble("t2");
  double TThree = Manager.GetDouble("t3");
  double Alpha = Manager.GetDouble("alpha");
  if (Alpha <= 0)
  {
    std::cout << "Anisotropy factor should be positive and nonzero.";
    exit(1);
  }
  char Axis ='y';

  if (Manager.GetBoolean("landau-x"))
    Axis ='x';
  if (Manager.GetBoolean("symmetric-gauge"))
    Axis = 's';
  
  if ((MaxBand<0)||(MaxBand >= FullUnitCellX*FullUnitCellY))
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

  char* StatisticPrefix = new char [16];
  sprintf (StatisticPrefix, "fermions");
  
  if (Manager.GetBoolean("boson") == false)
  {
		sprintf (StatisticPrefix, "fermions");
  }
  else
  {
		sprintf (StatisticPrefix, "bosons");
  }

  char* FilePrefix = new char [512];
  int lenFilePrefix=0;

  if (Manager.GetBoolean("triangular")==false)
  {
		lenFilePrefix += sprintf (FilePrefix, "%s_hofstadter_X_%d_Y_%d_q_%d", StatisticPrefix, UnitCellX, UnitCellY, FluxPerCell);

	  if ((NbrBands>1)||(MaxBand>0))
	  {
			if (NbrBands==1)
				lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_b_%d", MaxBand);
			else
				lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_b_%d-%d", MinBand, MaxBand);
	  }
	  if(SymmetricGauge)
	  {
	    lenFilePrefix += sprintf(FilePrefix+lenFilePrefix,"_sym");
	  }
	  if (Manager.GetBoolean("landau-x"))
	    lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_landau-x");
	  if(MultiBand == true)
			lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_multiband");
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
  // common naming options:
  lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_n_%d_x_%d_y_%d", NbrParticles, NbrCellX, NbrCellY);
  lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_t2_%g_t3_%g", TTwo,TThree);
  lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_alpha_%g", Alpha);

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
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
		if (Manager.GetBoolean("flat-band") == false)
			lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_u_%g",Manager.GetDouble("u-potential"));

		//lenFilePrefix += sprintf(FilePrefix+lenFilePrefix, "_3b_%g",Manager.GetDouble("three-body-potential"));
		lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_gx_%g_gy_%g", Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		//lenFilePrefix += sprintf(FilePrefix+lenFilePrefix, "_v_%g_w_%g",Manager.GetDouble("v-potential"),Manager.GetDouble("w-potential"));
		sprintf (EigenvalueOutputFile,"%s.dat",FilePrefix);
	}
  
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
	{
		bool ExportOneBody = false;
		if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
			ExportOneBody = true;

	// Abstract2DTightBindingModel *TightBindingModel;
		TightBindingModelHofstadterSquare *TightBindingModel;
		if (Manager.GetBoolean("triangular")==false)
			TightBindingModel= new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody,false,TTwo,TThree,Alpha);
//      else
//	TightBindingModel= new TightBindingModelHofstadterTriangularQuarter(NbrCellX, NbrCellY, Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);

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
  if (Manager.GetInteger("only-kx") >= 0)
	{
		MinKx = Manager.GetInteger("only-kx");
		MaxKx = MinKx;
	}
  int MinKy = 0;
  int MaxKy = NbrCellY - 1;
  if (Manager.GetInteger("only-ky") >= 0)
	{
		MinKy = Manager.GetInteger("only-ky");
		MaxKy = MinKy;
	}

 // Abstract2DTightBindingModel *TightBindingModel;
  TightBindingModelHofstadterSquare *TightBindingModel;
  if (Manager.GetBoolean("triangular")==false)
  {
    //std::cout << "Now creating TB model." << endl;
    TightBindingModel= new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis,
							     Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(),true,false,TTwo,TThree,Alpha);
    //std::cout << "Done with TB model." << endl;
  }
  //else
  //  TightBindingModel= new TightBindingModelHofstadterTriangularQuarter(NbrCellX, NbrCellY, Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());


  if (Manager.GetBoolean("boson") == false)
	{
		int FilledNbrBands=-1;

		double E=TightBindingModel->ComputeGroundstateEnergy(NbrParticles,FilledNbrBands, true);

		cout << "Total energy of groundstate: "<<E<<" ("<<FilledNbrBands<<" filled bands)"<<endl;
	}
  
  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
  {
		for (int j = MinKy; j <= MaxKy; ++j)
		{
			cout << "(kx=" << i << ",ky=" << j << ") : " << endl;
			ParticleOnSphere* Space = 0;
			AbstractQHEHamiltonian* Hamiltonian = 0;
			if (NbrBands==1)
			{
				if(MultiBand==true)
				{
					std::cout << "Don't use multiband.\n";
					exit(1);
					if(Manager.GetBoolean("boson") == false)
					{
						std::cout << "Only bosons implemented in multiband case currently.\n";
						exit(1);
					}
					else
					{
						//Space = new BosonOnSquareLatticeMultiBandMomentumSpace (NbrParticles, NbrCellX, NbrCellY, i, j, 2);
					}
					if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
						Memory = Architecture.GetArchitecture()->GetLocalMemory();
					Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
				// assign Hamiltonian:
					Hamiltonian = new ParticleOnLatticeHofstadterMultiBandHamiltonian(Space, NbrParticles, NbrCellX, NbrCellY, MaxBand,
																																					 Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"),
																																					 TightBindingModel, Manager.GetBoolean("flat-band"),MultiBand, Architecture.GetArchitecture(), Memory);
				}
				else
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
					if(Manager.GetDouble("three-body-potential") != 0)
					{
						Hamiltonian = new ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian(Space, NbrParticles, NbrCellX, NbrCellY, MaxBand,
																																					 Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"),
																																					 Manager.GetDouble("three-body-potential"),
																																					 TightBindingModel, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
					}
					else
					{
						std::cout << "Now creating MB hamiltonian." << endl;
						Hamiltonian = new ParticleOnLatticeHofstadterSingleBandHamiltonian(Space, NbrParticles, NbrCellX, NbrCellY, MaxBand,
																																					 Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"),
																																					 TightBindingModel, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
						std::cout << "Done with MB hamiltonian." << endl;
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
					Hamiltonian = new ParticleOnLatticeTwoBandHofstadterHamiltonian((ParticleOnSphereWithSpin*)Space, NbrParticles,
					NbrCellX, NbrCellY, MinBand, Manager.GetDouble("u-potential"), 0.0 /*Manager.GetDouble("v-potential")*/, TightBindingModel,
													Manager.GetBoolean("flat-band"),
													Architecture.GetArchitecture(), Memory);
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
						Hamiltonian = new ParticleOnLatticeFourBandHofstadterHamiltonian((ParticleOnSphereWithSU4Spin*)Space, NbrParticles, NbrCellX, NbrCellY, MinBand, Manager.GetDouble("u-potential"), 0.0 /*Manager.GetDouble("v-potential")*/, TightBindingModel,
														 Manager.GetBoolean("flat-band"),
														 Architecture.GetArchitecture(), Memory);
					}
					else
					{
						cout << "Multi-band n>4 situations not implemented, yet"<<endl;
						exit(1);
					}
				}
			}

			char* ContentPrefix = new char[256];
			sprintf (ContentPrefix, "%d %d", i, j);
			char* EigenstateOutputFile = new char [512];

			sprintf (EigenstateOutputFile,"%s_kx_%d_ky_%d",FilePrefix, i, j);

			GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
			FirstRunFlag = false;
			MainTaskOperation TaskOperation (&Task);
			std::cout << "Now applying main operation." <<endl;
			TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			cout << "------------------------------------" << endl;
			delete Hamiltonian;
			delete Space;
			delete[] EigenstateOutputFile;
			delete[] ContentPrefix;
		}
  }
  delete TightBindingModel;
  return 0;
}
