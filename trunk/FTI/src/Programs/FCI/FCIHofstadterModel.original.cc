#include "Options/Options.h"

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
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"

#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeTwoBandHofstadterHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFourBandHofstadterHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"


#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterTriangularQuarter.h"

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
    (*SystemGroup) += new  BooleanOption  ('\n', "redundant-kx", "Calculate all kx subspaces", false);
    (*SystemGroup) += new  BooleanOption  ('\n', "redundant-ky", "Calculate all ky subspaces", false);
    (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
    (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
    (*SystemGroup) += new BooleanOption  ('\n', "triangular", "use the Hofstadter model for a triangular lattice");

    (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive onsite(boson) or NN (fermion) potential strength", 1.0);
    (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive NN(boson) or NNN (fermion) potential strength", 0.0);

    (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "tunnelling along a selected direction (t=1 along others)", 1.0);

    (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
    (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);

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
    (*SystemGroup) += new BooleanOption  ('\n', "embedding", "compute the band structure witht the embedding");

    (*SystemGroup) += new BooleanOption  ('\n', "real-space", "use the real space representation when considering the system with all bands");
    (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "use the real space representation when considering the system with all bandswithout the translations");
    (*SystemGroup) += new BooleanOption  ('\n', "hardcore", "consider hardcore bosons (oly valid in real space mode)");

    (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
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

    int MinBand = Manager.GetInteger("band-min");
    int MaxBand = Manager.GetInteger("band-max");
    bool EmbeddingFlag = Manager.GetBoolean("embedding");

    char Axis ='y';

    if (Manager.GetBoolean("landau-x"))
        Axis ='x';

    if ((MaxBand<0)||(MaxBand >= UnitCellX*UnitCellY))
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
                    lenFilePrefix += sprintf (FilePrefix, "%s_realspace_hardcore_hofstadter_X_%d_Y_%d_q_%d", StatisticPrefix, UnitCellX, UnitCellY, FluxPerCell);
            }
            else
            {
                if ( Manager.GetBoolean("hardcore") == false)
                    lenFilePrefix += sprintf (FilePrefix, "%s_realspace_notranslation_hofstadter_X_%d_Y_%d_q_%d", StatisticPrefix, UnitCellX, UnitCellY, FluxPerCell);
                else
                    lenFilePrefix += sprintf (FilePrefix, "%s_realspace_notranslation_hardcore_hofstadter_X_%d_Y_%d_q_%d", StatisticPrefix, UnitCellX, UnitCellY, FluxPerCell);
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
    // common naming options:
    lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_n_%d_x_%d_y_%d", NbrParticles, NbrCellX, NbrCellY);



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
        if (((Manager.GetBoolean("flat-band") == false)&&(Manager.GetBoolean("hardcore") == false ))||(Manager.GetDouble("v-potential")!=0.0))
            lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_u_%g",Manager.GetDouble("u-potential"));
        if (Manager.GetDouble("v-potential")!=0.0)
            lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_v_%g",Manager.GetDouble("v-potential"));

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
            TightBindingModel= new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody, EmbeddingFlag);
        else
            TightBindingModel= new TightBindingModelHofstadterTriangularQuarter(NbrCellX, NbrCellY, Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);

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

    Abstract2DTightBindingModel *TightBindingModel;
    if (Manager.GetBoolean("triangular") == false)
        TightBindingModel= new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true, EmbeddingFlag);
    else
        TightBindingModel= new TightBindingModelHofstadterTriangularQuarter(NbrCellX, NbrCellY, Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());

    if (Manager.GetBoolean("boson") == false)
    {
        int FilledNbrBands = -1;
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
                    Hamiltonian = new ParticleOnLatticeHofstadterSingleBandHamiltonian(Space, NbrParticles, NbrCellX, NbrCellY, MaxBand, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), TightBindingModel, Manager.GetBoolean("flat-band"),Architecture.GetArchitecture(), Memory);
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
                            Space = new BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation(NbrParticles, UnitCellX*NbrCellX, UnitCellY*NbrCellY, i, NbrCellX, j,  NbrCellY);
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
                        //		      cout <<  TightBindingMatrix <<endl;
                    }
                }
                else
                {
                    if(Manager.GetBoolean("no-translation") == false)
                    {
                        cout << "Fermion in real space with translation not yet supported"<<endl;
                    }
                    else
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

                        Hamiltonian = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), TightBindingMatrix, DensityDensityInteraction,Architecture.GetArchitecture(), Memory);
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
            sprintf (ContentPrefix, "%d %d", i, j);
            char* EigenstateOutputFile = new char [512];

            sprintf (EigenstateOutputFile,"%s_kx_%d_ky_%d",FilePrefix, i, j);

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
    delete TightBindingModel;
    return 0;
}

