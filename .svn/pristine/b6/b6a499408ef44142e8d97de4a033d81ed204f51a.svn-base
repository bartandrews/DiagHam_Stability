#include "Options/Options.h"

#include "HilbertSpace/BosonOnCubicLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeMomentumSpace.h"
#include "Hamiltonian/ParticleOnCubicLatticeSingleBandHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelHopfCubicLattice.h"
#include "Tools/FTITightBinding/TightBindingModel3DAtomicLimitLattice.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"

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
    OptionManager Manager ("FTI3DHopf" , "0.01");
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

    (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
    (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
    (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 2);
    (*SystemGroup) += new SingleIntegerOption  ('z', "nbr-sitez", "number of sites along the y direction", 2);
    (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
    (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
    (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kz", "only evalute a given z momentum sector (negative if all ky sectors have to be computed)", -1);
    (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
    (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
    (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site potential between like orbitals", 1.0);
    (*SystemGroup) += new SingleDoubleOption  ('\n', "uab-potential", "repulsive on-site potential between different orbitals", 1.0);
    (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive nearest neighbor potential strength", 0.0);
    (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three body interaction instead of a two body interaction");
    (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four body interaction instead of a two body interaction");
    (*SystemGroup) += new SingleDoubleOption  ('\n', "lambda", "hopping strength", 1.0);
    (*SystemGroup) += new BooleanOption  ('\n', "atomic", "take atomic limit, with chemical potential given by (-mus, 0, mus) for the three orbitals");
    (*SystemGroup) += new BooleanOption  ('\n', "self-energy", "include the self-energy term from non-normal ordering");
    (*SystemGroup) += new BooleanOption  ('\n', "intra", "include on-site intra-orbital repulsion for fermions");
    (*SystemGroup) += new BooleanOption('\n', "shift", "shift energy by +1.0 to help convergence");
    (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
    (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
    (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-z", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
    (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
    (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the Chern number of the fully filled band (only available in singleparticle-spectrum mode)");
    (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
    (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
    (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
    (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
    (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest energy band");
    (*SystemGroup) += new SingleIntegerOption  ('\n', "band-index", "index of the band that has to be partially filled, should be 0 (lower band), 1 (middle band) or 2 (upper band)", 1);
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

    if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
        cout << "see man page for option syntax or type FTI3DHopf -h" << endl;
        return -1;
    }
    if (Manager.GetBoolean("help") == true)
    {
        Manager.DisplayHelp (cout);
        return 0;
    }

    int NbrParticles = Manager.GetInteger("nbr-particles");
    int NbrSiteX = Manager.GetInteger("nbr-sitex");
    int NbrSiteY = Manager.GetInteger("nbr-sitey");
    int NbrSiteZ = Manager.GetInteger("nbr-sitez");
    long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

    char* StatisticPrefix = new char [16];
    if (Manager.GetBoolean("boson") == false)
    {
        sprintf (StatisticPrefix, "fermions");
    }
    else
    {
        sprintf (StatisticPrefix, "bosons");
    }

    char* FilePrefix = new char[256];
    sprintf(FilePrefix, "%s_singleband_hopf_p_%d_x_%d_y_%d_z_%d_l_%g", StatisticPrefix, NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ,
            Manager.GetDouble("lambda"));

    char* CommentLine = new char [256];
    sprintf (CommentLine, "eigenvalues\n# kx ky kz ");
    char* EigenvalueOutputFile = new char [512];
    if (Manager.GetString("eigenvalue-file")!=0)
        strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
    else
        sprintf(EigenvalueOutputFile, "%s.dat", FilePrefix);

    Abstract3DTightBindingModel *TightBindingModel = 0;
    if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
        bool ExportOneBody = false;
        if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true)
                || (Manager.GetBoolean("singleparticle-chernnumber") == true))
            ExportOneBody = true;
        if (Manager.GetBoolean("atomic") == true)
        {
            double chemicalPotentials[2] = {0.0, 1.0};
            TightBindingModel = new TightBindingModel3DAtomicLimitLattice(NbrSiteX, NbrSiteY, NbrSiteZ, 2, chemicalPotentials,
                    Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), Architecture.GetArchitecture(), ExportOneBody);
        }
        else
            TightBindingModel = new TightBindingModelHopfCubicLattice(NbrSiteX, NbrSiteY, NbrSiteZ, Manager.GetDouble("lambda"), Manager.GetInteger("band-index"),
                    Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), Architecture.GetArchitecture(), ExportOneBody);
        if (Manager.GetBoolean("singleparticle-chernnumber") == true)
        {
            cout << "Chern number = " << TightBindingModel->ComputeChernNumber(Manager.GetInteger("nbr-layers") - 1) << endl;
        }
        TightBindingModel->WriteAsciiSpectrum(EigenvalueOutputFile);
        double BandSpread = TightBindingModel->ComputeBandSpread(Manager.GetInteger("band-index"));
        double DirectBandGap = TightBindingModel->ComputeDirectBandGap(Manager.GetInteger("band-index"));
        cout << "Spread = " << BandSpread << "  Direct Gap = " << DirectBandGap  << "  Flattening = " << (BandSpread / DirectBandGap) << endl;
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
    int MaxKx = NbrSiteX - 1;
    if (Manager.GetInteger("only-kx") >= 0)
    {
        MinKx = Manager.GetInteger("only-kx");
        MaxKx = MinKx;
    }
    int MinKy = 0;
    int MaxKy = NbrSiteY - 1;
    if (Manager.GetInteger("only-ky") >= 0)
    {
        MinKy = Manager.GetInteger("only-ky");
        MaxKy = MinKy;
    }
    int MinKz = 0;
    int MaxKz = NbrSiteZ - 1;
    if (Manager.GetInteger("only-kz") >= 0)
    {
        MinKz = Manager.GetInteger("only-kz");
        MaxKz = MinKz;
    }

    TightBindingModel = new TightBindingModelHopfCubicLattice(NbrSiteX, NbrSiteY, NbrSiteZ, Manager.GetDouble("lambda"), Manager.GetInteger("band-index"),
            Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), Architecture.GetArchitecture());

    bool FirstRunFlag = true;
    for (int i = MinKx; i <= MaxKx; ++i)
    {
        for (int j = MinKy; j <= MaxKy; ++j)
        {
            for (int k = MinKz; k <= MaxKz; ++k)
            {
                cout << "(kx=" << i << ",ky=" << j << ",kz=" << k << ") : " << endl;

                ParticleOnSphere* Space = 0;
                if (Manager.GetBoolean("boson") == false)
                    Space = new FermionOnCubicLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, i, j, k);
                else
                    Space = new BosonOnCubicLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, i, j, k);
                cout << "dim = " << Space->GetHilbertSpaceDimension() << endl;
                if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
                    Memory = Architecture.GetArchitecture()->GetLocalMemory();
                Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
                AbstractQHEHamiltonian* Hamiltonian = 0;

                Hamiltonian = new ParticleOnCubicLatticeSingleBandHamiltonian(Space, NbrParticles,
                        NbrSiteX, NbrSiteY, NbrSiteZ, Manager.GetDouble("u-potential"), Manager.GetDouble("uab-potential"), Manager.GetDouble("v-potential"),
                        TightBindingModel, Manager.GetInteger("band-index"), Manager.GetBoolean("flat-band"), Manager.GetBoolean("self-energy"), Manager.GetBoolean("intra"), Architecture.GetArchitecture(), Memory);

                if (Manager.GetBoolean("shift"))
                    Hamiltonian->ShiftHamiltonian(1.0);

                char* ContentPrefix = new char[256];
                sprintf (ContentPrefix, "%d %d %d", i, j, k);
                char* EigenstateOutputFile = new char [512];
                if (Manager.GetString("eigenstate-file")!=0)
                    sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d_kz_%d", Manager.GetString("eigenstate-file"), i, j, k);
                else
                {
                    char* TmpExtention = new char [512];
                    sprintf (TmpExtention, "_kx_%d_ky_%d_kz_%d", i, j, k);
                    EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
                }

                GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian,
                        ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
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
    delete TightBindingModel;
    delete[] StatisticPrefix;
    delete[] FilePrefix;
    delete[] CommentLine;
    delete[] EigenvalueOutputFile;
    return 0;
}
