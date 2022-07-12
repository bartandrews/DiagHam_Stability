#include "Architecture/ArchitectureManager.h"

#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ArrayTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"

#include "Tools/FTITightBinding/TightBindingModelHaldaneHoneycombLattice.h"
#include "Tools/FTITightBinding/TightBindingModelAlternativeKagomeLattice.h"
#include "Tools/FTITightBinding/TightBindingModelRubyLattice.h"
#include "Tools/FTITightBinding/TightBindingModelTwoOrbitalSquareLattice.h"
#include "Tools/FTITightBinding/TightBindingModelPyrochloreSlabLattice.h"
#include "Tools/FTITightBinding/TightBindingModelChern3TwoOrbitalTriangularLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::setw;


// sort in ascending order
// returns the perumtation sign
double SortAscending(int Ne, int* list)
{
    int count = 0;

    for (int i = 0; i < Ne - 1; ++i)
    {
        int m = i + 1;
        for (int j = i + 2; j < Ne; ++j)
        {
            if (list[j] < list[m])
                m = j;
            else if (list[j] == list[m])
                return 0.;
        }
        if (list[i] == list[m])
            return 0.;
        if (list[m] < list[i])
        {
            ++count;
            int temp = list[m];
            list[m] = list[i];
            list[i] = temp;
        }
    }

    for (int i = 0; i < Ne - 1; ++i)
        if (list[i] >= list[i + 1])
        {
            cout << "sorting failed" << endl;
            exit(1);
        }
    return 1.0 - 2.0 * (count & 1);
}

int main(int argc, char** argv)
{
    ArchitectureManager Architecture;

    OptionManager Manager ("FCIBlochConstruction" , "0.01");
    OptionGroup* MiscGroup = new OptionGroup ("misc options");
    OptionGroup* SystemGroup = new OptionGroup ("system options");
    OptionGroup* OutputGroup = new OptionGroup ("output options");
    Manager += SystemGroup;
    Manager += OutputGroup;
    Manager += MiscGroup;
    (*SystemGroup) += new SingleIntegerOption('p', "nbr-particles", "number of particles", 8);
    (*SystemGroup) += new SingleIntegerOption('x', "nbr-sitex", "number of sites along the x direction", 6);
    (*SystemGroup) += new SingleIntegerOption('y', "nbr-sitey", "number of sites along the y direction", 4);
    (*SystemGroup) += new SingleIntegerOption('\n', "kx", "total momentum along the x direction of the FQH state", 0);
    (*SystemGroup) += new SingleIntegerOption('\n', "ky", "total momentum along the y direction of the FQH state", 0);
    (*SystemGroup) += new SingleStringOption('\n', "model", "model name", "triangle2");
    (*SystemGroup) += new SingleStringOption('\n', "dump", "dump the band structure used in this calculation to a file", "bandstructure");
    (*SystemGroup) += new BooleanOption('\n', "boson", "use bosonic statistics");
    (*SystemGroup) += new BooleanOption('\n', "pt", "use Gamma-shaped parallel-transport gauge");
    (*SystemGroup) += new BooleanOption('\n', "csign", "assume the FQH state is calculated using C with correct sign");
    (*SystemGroup) += new SingleIntegerOption('\n', "band", "index of the chern band", 0);
    (*SystemGroup) += new SingleStringOption('\n', "fqh", "filename of the FQH state vector (produced by FQHColorfulLLLBlochBasis)");
    (*SystemGroup) += new SingleStringOption('\n', "output", "filename template of the output vector. _kx_#_ky_#.#.vec will be appended", "output");
    (*SystemGroup) += new SingleStringOption('\n', "exact", "(optional) filename template of the exact diagonalization vector");
    (*MiscGroup) += new BooleanOption('h', "help", "display this help");

    if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
        cout << "see man page for option syntax or type FCIBlochConstruction -h" << endl;
        return -1;
    }
    if (Manager.GetBoolean("help") == true)
    {
        Manager.DisplayHelp (cout);
        return 0;
    }

    int Ne = Manager.GetInteger("nbr-particles");
    int Nx = Manager.GetInteger("nbr-sitex");
    int Ny = Manager.GetInteger("nbr-sitey");
    int Ns = Nx * Ny;
    int Band = Manager.GetInteger("band");
    int Kx = Manager.GetInteger("kx");
    int Ky = Manager.GetInteger("ky");
    bool NegativeCHack = !Manager.GetBoolean("csign");

    if (Manager.GetString("model") == 0)
    {
        cout << "please provide a model name" << endl;
        return -1;
    }

    Abstract2DTightBindingModel* tb;
    if (!strcmp(Manager.GetString("model"), "kagome"))
        tb = new TightBindingModelAlternativeKagomeLattice(Nx, Ny, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, Architecture.GetArchitecture());
    else if (!strcmp(Manager.GetString("model"), "haldane"))
        tb = new TightBindingModelHaldaneHoneycombLattice(Nx, Ny, 1.0, 1.0, 0.13, 0.0, 0.0, 0.0, Architecture.GetArchitecture());
    else if (!strcmp(Manager.GetString("model"), "ruby"))
        tb = new TightBindingModelRubyLattice(Nx, Ny, 1.0, 1.2, -1.2, 2.6, -1.2, 0.0, 0.0, 0.0, Architecture.GetArchitecture());
    else if (!strcmp(Manager.GetString("model"), "square"))
        tb = new TightBindingModelTwoOrbitalSquareLattice(Nx, Ny, 1.0, 1.0, 1.0, 1, 2.0, 0.0, 0.0, Architecture.GetArchitecture());
    else if (!strcmp(Manager.GetString("model"), "pyrochlore2"))
        tb = new TightBindingModelPyrochloreSlabLattice(Nx, Ny, 2, 1.0, 0.0, 0.83, 0.0, -1.03, 0.0, 0.0, 0.0, Architecture.GetArchitecture());
    else if (!strcmp(Manager.GetString("model"), "triangle2"))
        tb = new TightBindingModelChern3TwoOrbitalTriangularLattice(Nx, Ny, 1.0, 0.28, -0.22, 0.0, 0.0, 0.0, Architecture.GetArchitecture());
    else
    {
        cout << "please provide a VALID model name" << endl;
        return -1;
    }
    tb->WriteBandStructure(Manager.GetString("dump"));
    delete tb;

    Generic2DTightBindingModel TB(Manager.GetString("dump"));

    // Gamma = LLL boundary condition twisting angle
    double GammaX = TB.GetLLLGammaX(Band);
    double GammaY = TB.GetLLLGammaY(Band);
    int Chern = TB.GetChernNumber(Band);

    cout << "Curvature fluctuations / average" << endl;
    for (int ky = Ny - 1; ky >= 0; --ky)
    {
        for (int kx = 0; kx < Nx; ++kx)
            cout << "  " << (TB.GetCurvature(kx, ky, Band) / Chern - 1.0 / (Nx * Ny)) << "\t";
        cout << endl;
    }
    cout << endl;

    ComplexMatrix Gauge(Nx, Ny, true);
    if (Manager.GetBoolean("pt"))
        TB.BuildParallelTransportGauge(Band, Gauge);
    else
    {
        if (TB.BuildGeneralizedLandauGauge(Band, Gauge))
        {
            cout << "failed to build the generalized Landau gauge!";
            return -1;
        }
    }

    // lattice total momentum
    int KxLattice = Kx;
    int KyLattice = Ky;
    if (NegativeCHack && Chern < 0)
        KxLattice = (Nx - KxLattice) % Nx;

    double Nf = ((double)(Nx * Ny)) / Chern;

    // twist angle. just for show here
    double TwistAngle = TB.GetTwistAngle();
    if (NegativeCHack && Chern < 0)
        TwistAngle = M_PI - TwistAngle;

    if (NegativeCHack && Chern < 0)
        cout << "Use C>0 FQH state for both signs" << endl;

    cout << "Chern = " << Chern << endl;
    cout << "Gamma = (" << fabs(GammaX) << "," << GammaY << ")" << endl; // fabs needed for Chern < 0 with csign off
    cout << "K = (" << Kx << "," << Ky << ")" << endl;
    cout << "Klat = (" << KxLattice << "," << KyLattice << ")" << endl;
    cout << "TwistAngle = " << TwistAngle / M_PI << " Pi" << endl;

    cout << "Phase(W_x) = " << Arg(TB.GetAbelianWilsonLoopX(0, Band)) / (2 * M_PI) << " * 2 Pi" << endl;
    cout << "Phase(W_y) = " << Arg(TB.GetAbelianWilsonLoopY(0, Band)) / (2 * M_PI) << " * 2 Pi" << endl;

    ComplexVector input;
    if (input.ReadVector(Manager.GetString("fqh")) == false)
        return -1;


    ParticleOnSphere* CI = 0;
    ParticleOnSphere* LL = 0;
    if (Manager.GetBoolean("boson"))
    {
        CI = new BosonOnSquareLatticeMomentumSpace(Ne, Nx, Ny, KxLattice, KyLattice);
        LL = new BosonOnSquareLatticeMomentumSpace(Ne, Nx, Ny, Kx, Ky);
    }
    else
    {
        CI = new FermionOnSquareLatticeMomentumSpace(Ne, Nx, Ny, KxLattice, KyLattice);
        LL = new FermionOnSquareLatticeMomentumSpace(Ne, Nx, Ny, Kx, Ky);
    }
    if (LL->GetHilbertSpaceDimension() != CI->GetHilbertSpaceDimension())
    {
        cout << "dimension mismatch between the LLL Hilbert space (" << LL->GetHilbertSpaceDimension() << ") and the Chern insulator Hilbert space (" << CI->GetHilbertSpaceDimension() << ")" << endl;
        return -1;
    }
    if (LL->GetHilbertSpaceDimension() != input.GetVectorDimension())
    {
        cout << "dimension mismatch between the FQH state (" << input.GetVectorDimension() << ") and the LLL Hilbert space (" << LL->GetHilbertSpaceDimension() << ")" << endl;
        return -1;
    }

    int* Orbitals = new int[Ne];
    ComplexVector output(CI->GetHilbertSpaceDimension(), true);

    int percent = LL->GetHilbertSpaceDimension() / 100;
    if (percent == 0)
        percent = 1;

    for (int i = 0; i < LL->GetHilbertSpaceDimension(); ++i)
    {
        if (percent > 2000 && i % percent == 0)
            cout << i / percent << "%" << endl;

        Complex Amplitude = input[i];
        LL->GetOccupied(i, Orbitals);

        for (int l = 0; l < Ne; ++l)
        {
            int kx, ky;
            TB.GetLinearizedMomentumIndexSafe(Orbitals[l], kx, ky);
            if (NegativeCHack && Chern < 0)
                kx = (Nx - kx) % Nx;

            Amplitude *= Phase(- 2 * M_PI * ky * GammaX / Nf);

            Amplitude *= Gauge[ky][kx];
            Orbitals[l] = TB.GetLinearizedMomentumIndexSafe(kx, ky);
        }
        if (!Manager.GetBoolean("boson"))
            Amplitude *= SortAscending(Ne, Orbitals);
        int j = CI->FindStateIndex(Orbitals);
        if (j == -1)
        {
            cout << "cannot find the shifted state" << endl;
            exit(1);
        }
        output[j] = Amplitude;
    }

    // find the .?.vec tail
    char* tail = Manager.GetString("fqh");
    tail += strlen(tail) - 1;
    while (((*tail) != '.') && tail >= Manager.GetString("fqh"))
        --tail;
    if (tail >= Manager.GetString("fqh"))
    {
        --tail;
        while (((*tail) != '.') && tail >= Manager.GetString("fqh"))
            --tail;
    }

    if (tail < Manager.GetString("fqh"))
    {
        tail = Manager.GetString("fqh");
        tail += strlen(tail);
    }

    char* OutputFilename = new char[200];
    sprintf(OutputFilename, "%s_kx_%d_ky_%d%s", Manager.GetString("output"), KxLattice, KyLattice, tail);
    output.WriteVector(OutputFilename);
    delete[] OutputFilename;

    if (Manager.GetString("exact") != 0)
    {
        ComplexVector exact;
        if (exact.ReadVector(Manager.GetString("exact")) == true)
            cout << SqrNorm(exact * output) << endl;
    }

    delete[] Orbitals;
    delete LL;
    delete CI;
    return 0;
}
