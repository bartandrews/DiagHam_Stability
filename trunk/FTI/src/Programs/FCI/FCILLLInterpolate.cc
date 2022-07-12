#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "Hamiltonian/ParticleOnLatticeAlternativeKagomeLatticeLLLInterpolateSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeChern3TwoOrbitalTriangularLatticeLLLInterpolateSingleBandHamiltonian.h"

#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModelHaldaneHoneycombLattice.h"
#include "Tools/FTITightBinding/TightBindingModelAlternativeKagomeLattice.h"
#include "Tools/FTITightBinding/TightBindingModelRubyLattice.h"
#include "Tools/FTITightBinding/TightBindingModelTwoOrbitalSquareLattice.h"
#include "Tools/FTITightBinding/TightBindingModelPyrochloreSlabLattice.h"
#include "Tools/FTITightBinding/TightBindingModelChern3TwoOrbitalTriangularLattice.h"

#include "LanczosAlgorithm/LanczosManager.h"
#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
    cout.precision(14);

    OptionManager Manager("FCIColorfulLLLBlochBasis" , "0.01");
    OptionGroup* ToolsGroup = new OptionGroup("tools options");
    OptionGroup* MiscGroup = new OptionGroup("misc options");
    OptionGroup* SystemGroup = new OptionGroup("system options");
    OptionGroup* PrecalculationGroup = new OptionGroup("precalculation options");

    ArchitectureManager Architecture;
    LanczosManager Lanczos(true); // Complex
    Manager += SystemGroup;
    Architecture.AddOptionGroup(&Manager);
    Lanczos.AddOptionGroup(&Manager);
    Manager += PrecalculationGroup;
    Manager += ToolsGroup;
    Manager += MiscGroup;

    (*SystemGroup) += new SingleIntegerOption('p', "nbr-particles", "number of particles", 8);
    (*SystemGroup) += new SingleIntegerOption('x', "nbr-sitex", "number of sites along the x direction", 6);
    (*SystemGroup) += new SingleIntegerOption('y', "nbr-sitey", "number of sites along the y direction", 4);
    (*SystemGroup) += new SingleDoubleOption('\n', "interpolate", "interpolation parameter", 0.5);
    (*SystemGroup) += new SingleStringOption('\n', "model", "model name", "kagome");
    (*SystemGroup) += new SingleDoubleOption('\n', "angle", "angle between the two fundamental cycles of the torus divided by Pi. set to lattice model value if negative", -1);
    (*SystemGroup) += new SingleDoubleOption('\n', "ratio", "aspect ratio Lx / Ly of the torus. set to Nx / Ny if negative", -1);
    (*SystemGroup) += new SingleIntegerOption('\n', "only-kx", "constraint on the total momentum in the x direction (negative if none)", -1);
    (*SystemGroup) += new SingleIntegerOption('\n', "only-ky", "constraint on the total momentum in the y direction (negative if none)", -1);
    (*SystemGroup) += new BooleanOption('\n', "inversion", "use inversion symmetry to reduce calculation");
    (*SystemGroup) += new BooleanOption('\n', "boson", "use bosonic statistics instead of fermionic statistics");
    (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive two-body nearest neighbor potential strength", 1.0);
    (*SystemGroup) += new SingleStringOption('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
    (*SystemGroup) += new SingleStringOption('\n', "eigenvalue-file", "filename for eigenvalues output");
    (*SystemGroup) += new SingleStringOption('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
    (*SystemGroup) += new BooleanOption('\n', "shift", "shift energy by +1.0 to help convergence");

    (*PrecalculationGroup) += new SingleIntegerOption('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
    (*PrecalculationGroup) += new SingleStringOption('\n', "load-precalculation", "load precalculation from a file",0);
    (*PrecalculationGroup) += new SingleStringOption('\n', "save-precalculation", "save precalculation in a file",0);

#ifdef __LAPACK__
    (*ToolsGroup) += new BooleanOption('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
    (*MiscGroup) += new BooleanOption('h', "help", "display this help");

    if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
        cout << "see man page for option syntax or type FCIColorfulLLLBlochBasis -h" << endl;
        return -1;
    }
    if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
        Manager.DisplayHelp (cout);
        return 0;
    }


    int NbrParticles = Manager.GetInteger("nbr-particles");
    int NbrSiteX = Manager.GetInteger("nbr-sitex");
    int NbrSiteY = Manager.GetInteger("nbr-sitey");

    Abstract2DTightBindingModel* TB;
    if (!strcmp(Manager.GetString("model"), "kagome"))
        TB = new TightBindingModelAlternativeKagomeLattice(NbrSiteX, NbrSiteY, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, Architecture.GetArchitecture());
    else if (!strcmp(Manager.GetString("model"), "haldane"))
        TB = new TightBindingModelHaldaneHoneycombLattice(NbrSiteX, NbrSiteY, 1.0, 1.0, 0.13, 0.0, 0.0, 0.0, Architecture.GetArchitecture());
    else if (!strcmp(Manager.GetString("model"), "ruby"))
        TB = new TightBindingModelRubyLattice(NbrSiteX, NbrSiteY, 1.0, 1.2, -1.2, 2.6, -1.2, 0.0, 0.0, 0.0, Architecture.GetArchitecture());
    else if (!strcmp(Manager.GetString("model"), "square"))
        TB = new TightBindingModelTwoOrbitalSquareLattice(NbrSiteX, NbrSiteY, 1.0, 1.0, 1.0, 1, 2.0, 0.0, 0.0, Architecture.GetArchitecture());
    else if (!strcmp(Manager.GetString("model"), "pyrochlore2"))
        TB = new TightBindingModelPyrochloreSlabLattice(NbrSiteX, NbrSiteY, 2, 1.0, 0.0, 0.83, 0.0, -1.03, 0.0, 0.0, 0.0, Architecture.GetArchitecture());
    else if (!strcmp(Manager.GetString("model"), "triangle2"))
        TB = new TightBindingModelChern3TwoOrbitalTriangularLattice(NbrSiteX, NbrSiteY, 1.0, 0.28, -0.22, 0.0, 0.0, 0.0, Architecture.GetArchitecture());
    else
    {
        cout << "please provide a VALID model name" << endl;
        return -1;
    }

    double TwistAngle = M_PI * Manager.GetDouble("angle");
    if (TwistAngle < 0)
        TwistAngle = TB->GetTwistAngle();

    double AspectRatio = Manager.GetDouble("ratio");
    if (AspectRatio < 0)
        AspectRatio = ((double) NbrSiteX) / NbrSiteY;

    char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
    char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
    long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

    int NbrPseudoPotentials = 0;
    double* PseudoPotentials = 0;

    char* InteractionName = 0;
    if (Manager.GetString("interaction-file") != 0)
    {
        ConfigurationParser InteractionDefinition;
        if (InteractionDefinition.Parse(Manager.GetString("interaction-file")) == false)
        {
            InteractionDefinition.DumpErrors(cout) << endl;
            exit(-1);
        }
        InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials, NbrPseudoPotentials);
        cout << "loaded pseudopotentials:" << endl;
        for (int i = 0; i < NbrPseudoPotentials; ++i)
            cout << "  V_" << i << " = " << PseudoPotentials[i] << endl;
        cout << endl;
    }
    else
    {
        if (Manager.GetBoolean("boson"))
        {
            cout << "use default pseudopotentials (Halperin 221)" << endl;
            NbrPseudoPotentials = 1;
            PseudoPotentials = new double[1];
            PseudoPotentials[0] = 1.0;
        }
        else
        {
            cout << "use default pseudopotentials (Halperin 332)" << endl;
            NbrPseudoPotentials = 2;
            PseudoPotentials = new double[2];
            PseudoPotentials[0] = 1.0; // needed for inter-layer repulsion, due to the absence of identical particle minus sign
            PseudoPotentials[1] = 1.0;
        }
    }

    char* StatisticPrefix = new char[16];
    if (Manager.GetBoolean("boson"))
        sprintf(StatisticPrefix, "bosons");
    else
        sprintf(StatisticPrefix, "fermions");
    char* FilePrefix = new char[256];
    sprintf(FilePrefix, "%s_%s-LLL_n_%d_x_%d_y_%d_u_%f_angle_%f_interpolate_%f", StatisticPrefix, Manager.GetString("model"),
            NbrParticles, NbrSiteX, NbrSiteY, Manager.GetDouble("u-potential"), TwistAngle, Manager.GetDouble("interpolate"));

    char* CommentLine = new char[256];
    sprintf (CommentLine, "eigenvalues\n# kx ky ");
    char* EigenvalueOutputFile = new char[512];
    if (Manager.GetString("eigenvalue-file") != 0)
        strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
    else
        sprintf(EigenvalueOutputFile, "%s.dat", FilePrefix);


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
    bool FirstRunFlag = true;
    for (int i = MinKx; i <= MaxKx; ++i)
    {
        for (int j = MinKy; j <= MaxKy; ++j)
        {
            int k = TB->GetLinearizedMomentumIndexSafe(i, j);
            int m = TB->GetLinearizedMomentumIndexSafe(-i, -j);
            if (Manager.GetBoolean("inversion") && (k > m))
                continue;

            cout << "(kx=" << i << ",ky=" << j << ") : " << endl;
            ParticleOnSphere* Space;
            if (Manager.GetBoolean("boson"))
            {
                Space = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, i, j);
            }
            else
            {
                if ((NbrSiteX * NbrSiteY) <= 63)
                    Space = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, i, j);
                else
                    Space = new FermionOnSquareLatticeMomentumSpaceLong(NbrParticles, NbrSiteX, NbrSiteY, i, j);
            }

            if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
                Memory = Architecture.GetArchitecture()->GetLocalMemory();
            Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

            AbstractQHEHamiltonian* Hamiltonian;

            if (!strcmp(Manager.GetString("model"), "kagome"))
                Hamiltonian = new ParticleOnLatticeAlternativeKagomeLatticeLLLInterpolateSingleBandHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, TB, 
                        Manager.GetDouble("u-potential"), 0.0, 0, true,
                        Manager.GetDouble("interpolate"), TwistAngle, AspectRatio, NbrPseudoPotentials, PseudoPotentials,
                        Architecture.GetArchitecture(), Memory);
            else if (!strcmp(Manager.GetString("model"), "triangle2"))
                Hamiltonian = new ParticleOnLatticeChern3TwoOrbitalTriangularLatticeLLLInterpolateSingleBandHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, TB, 
                        Manager.GetDouble("u-potential"), true,
                        Manager.GetDouble("interpolate"), TwistAngle, AspectRatio, NbrPseudoPotentials, PseudoPotentials,
                        Architecture.GetArchitecture(), Memory);
            else
            {
                cout << "please provide a VALID model name for Hamiltonian" << endl;
                return -1;
            }

            if (Manager.GetBoolean("shift"))
                Hamiltonian->ShiftHamiltonian(1.0);

            char* ContentPrefix = new char[256];
            sprintf (ContentPrefix, "%d %d", i, j);
            char* EigenstateOutputFile = new char[512];
            if (Manager.GetString("eigenstate-file")!=0)
                sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), i, j);
            else
            {
                char* TmpExtention = new char [512];
                sprintf(TmpExtention, "_kx_%d_ky_%d", i, j);
                EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
            }

            GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian,
                    ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
            FirstRunFlag = false;
            MainTaskOperation TaskOperation(&Task);
            TaskOperation.ApplyOperation(Architecture.GetArchitecture());
            cout << "------------------------------------" << endl;
            delete Hamiltonian;
            delete Space;
            delete[] EigenstateOutputFile;
            delete[] ContentPrefix;
        }
    }

    return 0;
}


