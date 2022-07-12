#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iomanip>

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

#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "Hamiltonian/ParticleOnLatticeHaldaneModelSingleBandHamiltonian.h"


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::setw;

static bool YBasis = true;
static bool Unitary = true;
static bool FlatCurvature = false;
static bool XLQ = false;
static bool Extreme = false;
static bool Bold = false;

static int Ne;
static int Nx;
static int Ny;
static int Nb;
static int Band;
static int Kpx;
static int Kpy;
static int Nf;
static int N;
static int P;
static int Q;
static int N0x;
static int N0y;
static int Px;
static int Py;
static int Qx;
static int Qy;

static int Delta;
static int C;

static ComplexVector ShiftX;
static ComplexVector ShiftY;
static ComplexMatrix* Bloch = NULL;

// Perform Wannier construction
void Construct(OptionManager& Manager, FermionOnTorusWithMagneticTranslations& torus, ComplexVector& input, ComplexMatrix& xGauge, ComplexMatrix& yGauge, Complex* expN, Complex* expNx);
Complex SumKxPermutations(int* xlist, int* kxlist, int* kylist, Complex* expNx); // sum over all allowed permutations of kx

// utils for the manybody construction
inline unsigned long Shift(unsigned long state, int shift, int length);
inline int Count(unsigned long state);
void Print(unsigned long state, int length);
double SortDescending(int* list);
inline void Unflatten(unsigned long state, int* kxlist, int* kylist);
inline void UnflattenJ(unsigned long state, int* xlist, int* kylist); // obtain the list of (x,ky) from the binary representation of j, in descending order of j
inline double SortYMajor(int* kxlist, int* kylist); // sort {kx,ky} according to ky * Nx + kx. return the permutation sign
inline unsigned long HistSum(int* kylist); // calculate sum(Ne ** ky for ky in kylist)
inline void Print2DList(int* xlist, int* ylist); // print 2D list
inline void Swap(int* list, int i, int j); // swap two elements in a list

Complex PowInt(Complex z, int k);
Complex Dot(ComplexVector bra, ComplexVector ket, ComplexVector shift); // < bra | ket > with sublattice shift
Complex Ax(int kx, int ky); // < kx, ky | x | kx+1, ky >
Complex Ay(int kx, int ky); // < kx, ky | y | kx, ky+1 >
Complex Wx(int ky); // Wilson loop along fixed ky
Complex Wy(int kx); // Wilson loop along fixed kx
Complex Wsq(int kx, int ky); // Wilson loop around (0,ky)(kx,ky)(kx,ky+1)(0,ky+1)(0,ky)
Complex Wpq(int kx, int ky); // Wilson loop around a plaquette
Complex Lx(int ky); // eigenvalue of projected x operator in X=0 unit cell
Complex TildeLy(void); // WyU() averaged to Ny bonds, with minimal abs(arg angle)
int Principal(int ky); // pull ky back to the principal BZ, defined by C * ky + Delta in [0, Ny)
int J(int X, int ky); // j = X * Ny + C * ky + delta
void FindXky(int j, int& X, int& ky); // inverse map of j = X * Ny + C * ky + delta, C * ky + delta \in [0, Ny). pull ky back to [0, Ny) after the mapping
Complex Wsqbar(int kx, int ky); // "average" Wilson loop around (0,ky)(kx,ky)(kx,ky+1)(0,ky+1)(0,ky)
Complex Ux(int ky); // unitary measure of curvature fluctuation at ky
Complex WyU(void); // modified Wilson loop along kx = 0
Complex Phi(int X, int ky); // phase in front of |X,ky>

// build wannier state
// yGauge[ky][x] = e^iPhi(x,ky)
// xGauge[ky][kx] = lambda^kx / prod_kappa^kx Ax(kappa,ky)
void BuildWannier(char* onebody, ComplexMatrix& xGauge, ComplexMatrix& yGauge);
void CheckParallelTransport(ComplexMatrix xGauge); // check parallel transport along each ky
void CheckOrthonormality(ComplexMatrix xGauge); // check orthonormality of the Wannier states
void CheckWannierConnection(ComplexMatrix xGauge, ComplexMatrix yGauge); // check <X,ky|Y|X',ky'>
void BuildSewing(char *inversion, RealVector sublattice, ComplexMatrix& sewing); // build sewing matrix e^{i xi_{kx,ky}}
void CheckWannierInversion(ComplexMatrix sewing, ComplexMatrix xGauge, ComplexMatrix yGauge); // check inversion symmetry of |X,ky>

int main(int argc, char** argv)
{
    ArchitectureManager Architecture;

    OptionManager Manager ("FQHETopInsulatorWannierConstruction" , "0.01");
    OptionGroup* MiscGroup = new OptionGroup ("misc options");
    OptionGroup* SystemGroup = new OptionGroup ("system options");
    OptionGroup* OutputGroup = new OptionGroup ("output options");
    Manager += SystemGroup;
    Manager += OutputGroup;
    Manager += MiscGroup;
    (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
    (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
    (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 4);
    (*SystemGroup) += new SingleIntegerOption  ('\n', "torus-kx", "total momentum along the x direction for the state on torus", 2);
    (*SystemGroup) += new SingleIntegerOption  ('\n', "torus-ky", "total momentum along the y direction for the state on torus", 6);
    (*SystemGroup) += new SingleIntegerOption ('\n', "band", "index of the chern band", 0);
    (*SystemGroup) += new BooleanOption  ('\n', "x-basis", "use Wannier states localized along a line parallel to the x direction");
    (*SystemGroup) += new BooleanOption  ('\n', "non-unitary", "use non-unitary connection");
    (*SystemGroup) += new BooleanOption  ('\n', "xlq", "use XL Qi's recipe (parallel transport in y without considering Wy)");
    (*SystemGroup) += new BooleanOption  ('\n', "flat", "assume flat curvature");
    (*SystemGroup) += new BooleanOption  ('\n', "extreme", "make the Wannier connection as close to the positive real axis as possible");
    (*SystemGroup) += new BooleanOption  ('\n', "bold", "do not choke on inversion breaking, etc.");
    (*SystemGroup) += new BooleanOption  ('\n', "test", "dry run, without constructing the many-body wave functions");
    (*SystemGroup) += new SingleStringOption ('\n', "onebody", "filename template of the onebody wave functions. _kx_#_ky_#.onebody will be appended");
    (*SystemGroup) += new SingleStringOption ('\n', "sublattice", "filename of the sublattice shifts");
    (*SystemGroup) += new SingleStringOption ('\n', "inversion", "filename of the inversion matrix");
    (*SystemGroup) += new SingleStringOption ('\n', "input", "filename of the input state vector (FQH on torus, with manybody translational symmetry fully reduced)", "input");
    (*SystemGroup) += new SingleStringOption ('\n', "output", "filename template of the output vector. _kx_#_ky_#.#.vec will be appended.", "output");
    (*SystemGroup) += new SingleStringOption ('\n', "exact", "filename template of the exact diagonalization vector. _kx_#_ky_#.#.vec will be appended.");
    (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

    if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
        cout << "see man page for option syntax or type FQHETopInsulatorWannierConstruction -h" << endl;
        return -1;
    }
    if (Manager.GetBoolean("help") == true)
    {
        Manager.DisplayHelp (cout);
        return 0;
    }
    if (Manager.GetString("sublattice") == 0)
    {
        cout << "please provide --sublattice parameter." << endl;
        return 1;
    }
    if (Manager.GetString("onebody") == 0)
    {
        cout << "please provide --onebody parameter." << endl;
        return 1;
    }

    // FIXME: This code has only been checked for the Laughlin state, which always have Kpx = Kpy, and they are always inversion symmetric, i.e. Kpx = -Kpx mod N
    // It's quite possible that there are bugs that will only show up from non-Abelian states.

    YBasis = !Manager.GetBoolean("x-basis");
    Unitary = !Manager.GetBoolean("non-unitary");
    FlatCurvature = Manager.GetBoolean("flat");
    XLQ = Manager.GetBoolean("xlq");
    Extreme = Manager.GetBoolean("extreme");
    Bold = Manager.GetBoolean("bold");

    Ne = Manager.GetInteger("nbr-particles");
    Band = Manager.GetInteger("band");

    if (YBasis)
    {
        Nx = Manager.GetInteger("nbr-sitex");
        Ny = Manager.GetInteger("nbr-sitey");
        Kpx = Manager.GetInteger("torus-kx");
        Kpy = Manager.GetInteger("torus-ky");
    }
    else
    {
        Nx = Manager.GetInteger("nbr-sitey");
        Ny = Manager.GetInteger("nbr-sitex");
        Kpx = Manager.GetInteger("torus-ky");
        Kpy = Manager.GetInteger("torus-kx");
    }

    if (pow((double) Ne, Ny) > 1.9 * (1ul << 63))
    {
        cout << "HistSum will fail!" << endl;
        exit(1);
    }

    Bloch = new ComplexMatrix[Nx * Ny];
    Nf = Nx * Ny;
    N = FindGCD(Ne, Nf);
    P = Ne / N;
    Q = Nf / N;
    N0x = FindGCD(Ne, Nx);
    N0y = FindGCD(Ne, Ny);
    Px = Ne / N0x;
    Py = Ne / N0y;
    Qx = Nx / N0x;
    Qy = Ny / N0y;

    cout << "Ne, Nx, Ny = " << Ne << ", " << Nx << ", " << Ny << endl;
    cout << "Kappa_x, Kappa_y = " << Kpx << ", " << Kpy << endl;
    cout << "N, Nf = " << N << ", " << Nf << endl;
    cout << "N0x, N0y = " << N0x << ", " << N0y << endl;
    cout << "p, q = " << P << ", " << Q << endl;
    cout << "px, qx = " << Px << ", " << Qx << endl;
    cout << "py, qy = " << Py << ", " << Qy << endl;

    ComplexMatrix xGauge(Nx, Ny);
    ComplexMatrix yGauge(Nx, Ny);

    RealVector sublattice;
    sublattice.ReadVector(Manager.GetString("sublattice"));
    Nb = sublattice.GetVectorDimension() / 2;
    ShiftX.ResizeAndClean(Nb);
    ShiftY.ResizeAndClean(Nb);
    for (int i = 0; i < Nb; ++i)
    {
        ShiftX[i] = Phase(- 2.0 * M_PI * ((double)sublattice[2 * i + !YBasis]) / Nx);
        ShiftY[i] = Phase(- 2.0 * M_PI * ((double)sublattice[2 * i + YBasis]) / Ny);
    }

    BuildWannier(Manager.GetString("onebody"), xGauge, yGauge);
    if (!Unitary)
        CheckOrthonormality(xGauge);

    cout << "C = " << C << endl;
    cout << "delta = " << Delta << endl;

    cout << "Wx(0) = " << Wx(0) << endl;
    cout << "Wy(0) = " << Wy(0) << endl;

    cout << "xGauge" <<endl;
    for (int kx = 0; kx < Nx; ++kx)
    {
        for (int ky = 0; ky < Ny; ++ky)
        {
            cout << Arg(xGauge[ky][kx]) / (2.0 * M_PI) << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "yGauge" <<endl;
    for (int x = 0; x < Nx; ++x)
    {
        for (int ky = 0; ky < Ny; ++ky)
        {
            cout << Arg(yGauge[ky][x]) / (2.0 * M_PI) << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "TildeLy: " << TildeLy() << endl << endl;

    cout << "Ay(0,ky)" << endl;
    for (int ky = 0; ky < Ny; ++ky)
        cout << Arg(Ay(0, ky)) / (2.0 * M_PI) << endl;
    cout << endl;

    cout << "Ux(ky)" << endl;
    for (int ky = 0; ky < Ny; ++ky)
        cout << Arg(Ux(ky)) / (2.0 * M_PI) << endl;
    cout << endl;

    cout << "Wx(ky)" << endl;
    for (int ky = 0; ky < Ny; ++ky)
        cout << Arg(Wx(ky)) / (2.0 * M_PI) << endl;
    cout << endl;

    cout << "Wilson around plaquette. flat = " << 1. / Nf << endl;
    for (int kx = 0; kx < Nx; ++kx)
    {
        for (int ky = 0; ky < Ny; ++ky)
        {
            cout << Arg(Wpq(kx, ky)) / (2.0 * M_PI) << " ";
        }
        cout << endl;
    }
    cout << endl;

    CheckParallelTransport(xGauge);
    CheckWannierConnection(xGauge, yGauge);

    if (Manager.GetString("inversion") != 0)
    {
        ComplexMatrix sewing(Nx, Ny);
        BuildSewing(Manager.GetString("inversion"), sublattice, sewing);
        CheckWannierInversion(sewing, xGauge, yGauge);

        if (Nx % 2 == 0)
            if (Norm(sewing[0][0] * sewing[0][Nx / 2] - Wx(0)) > 1e-14)
            {
                cout << "Sewing matrix incompatible with Wx" << endl;
                if (!XLQ && Unitary && !Bold)
                    exit(1);
            }
        if (Ny % 2 == 0)
            if (Norm(sewing[0][0] * sewing[Ny / 2][0] - Wy(0)) > 1e-14)
            {
                cout << "Sewing matrix incompatible with Wy" << endl;
                if (!XLQ && Unitary && !Bold)
                    exit(1);
            }
    }

    // phase tables
    Complex* expN = new Complex[N];
    Complex* expNx = new Complex[Nx];
    for (int i = 0; i < N; ++i)
        expN[i] = Phase(- 2.0 * M_PI * ((double) i) / N);
    for (int i = 0; i < Nx; ++i)
        expNx[i] = Phase(- 2.0 * M_PI * ((double) i) / Nx);


    FermionOnTorusWithMagneticTranslations torus(Ne, Nf, Kpx, YBasis ? Kpy : ((Nf - Kpy) % Nf));
    ComplexVector input;
    if (input.ReadVector(Manager.GetString("input")) == false)
        return -1;
    if (torus.GetHilbertSpaceDimension() != input.GetVectorDimension())
    {
        cout << "dimension mismatch between the state (" << input.GetVectorDimension() << ") and the Hilbert space (" << torus.GetHilbertSpaceDimension() << ")" << endl;
        return -1;
    }

    if (Manager.GetBoolean("test") == false)
        Construct(Manager, torus, input, xGauge, yGauge, expN, expNx);

    delete[] expN;
    delete[] expNx;
    return 0;
}

// Perform Wannier construction
void Construct(OptionManager& Manager, FermionOnTorusWithMagneticTranslations& torus, ComplexVector& input, ComplexMatrix& xGauge, ComplexMatrix& yGauge, Complex* expN, Complex* expNx)
{
    int* kxlist = new int[Ne];
    int* kylist = new int[Ne];
    int* jlist = new int[Ne];
    int* xlist = new int[Ne];
    int* lxlist = new int[Ne];
    int* lylist = new int[Ne];
    char* outputname = new char[512];

    for (int s = 0; s < Qx; ++s)
    {
        int Kx = Kpx + s * Ne;
        while (Kx < 0)
            Kx += Nx;
        Kx %= Nx;
        for (int r = 0; r < Q / Qx ; ++r)
        {
            int ysign = YBasis ? 1 : (-1);
            int Ky = C * ysign * (Kpy + (r - Delta * ysign) * Ne);
            while (Ky < 0)
                Ky += Ny;
            Ky %= Ny;

            FermionOnSquareLatticeMomentumSpace space(Ne, YBasis?Nx:Ny, YBasis?Ny:Nx, YBasis?Kx:Ky, YBasis?Ky:Kx);
            ComplexVector output(space.GetHilbertSpaceDimension(), true);

            cout << "============ (Kx, r, Ky) = (" << Kx << ", " << r << ", ";
            cout << Ky << "), Dimension = " << space.GetHilbertSpaceDimension() << endl;

            for (int index = 0; index < space.GetHilbertSpaceDimension(); ++index)
            {
                if (index % 100 == 0)
                    cout << index << "/" << space.GetHilbertSpaceDimension() << endl;

                Complex amp_k = 0.;

                unsigned long state_k = space.GetStateDescription(index);
                Unflatten(state_k, kxlist, kylist);
                unsigned long hist = HistSum(kylist);
                double sign_k = SortYMajor(kxlist, kylist);

                for (int orbit = 0; orbit < torus.GetHilbertSpaceDimension(); ++orbit)
                {
                    unsigned long reorder = torus.GetReorderingSign(orbit);
                    int z = torus.GetNbrStateInOrbit(orbit);
                    unsigned long canonical = torus.GetStateDescription(orbit);

                    Complex amp_o = 0.;
                    for (int m = 0; m < z; ++m)
                    {
                        unsigned long state_m = Shift(canonical, (m * Q) % Nf, Nf);

                        int shift = YBasis ? ((Nf - r) % Nf) : r;
                        int count = Count(state_m & ~((1 << shift) - 1)); // # of particles that are shifted across the BZ boundary
                        double sign_mr = 1.0 - 2.0 * ((count * (Ne - count)) & 1);
                        state_m = Shift(state_m, shift, Nf);

                        UnflattenJ(state_m, xlist, lylist);
                        if (hist != HistSum(lylist))
                            continue;

                        double sign_l = SortYMajor(xlist, lylist); // then lylist is the same as kylist

                        for (int i = 0; i < Ne; ++i)
                            lxlist[i] = kxlist[i];
                        Complex sum_lx = SumKxPermutations(xlist, lxlist, lylist, expNx); // only lxlist is changed.
                        Complex phi = 1.;
                        for (int i = 0; i < Ne; ++i)
                            phi *= yGauge[lylist[i]][xlist[i]];

                        Complex amp_m = sign_mr * sign_l * sum_lx * phi;
                        amp_m *= 1.0 - 2.0 * ((reorder >> m) & 1ul);
                        amp_m *= expN[(Kpx * m) % N];

                        amp_o += amp_m;
                    }
                    amp_o *= input[orbit] / sqrt(z);

                    amp_k += amp_o;
                }
                amp_k *= sign_k * sqrt(Qx) / pow(sqrt(Nx), Ne);
                output[index] = amp_k;
            }

            for (int index = 0; index < space.GetHilbertSpaceDimension(); ++index)
            {
                unsigned long state_k = space.GetStateDescription(index);
                Unflatten(state_k, kxlist, kylist);
                for (int i = 0; i < Ne; ++i)
                    output[index] *= xGauge[kylist[i]][kxlist[i]];
            }
            cout << "Norm = "<< output.Norm() << endl;

            if (!Unitary)
            {
                double norm = output.Norm();
                for (int index = 0; index < space.GetHilbertSpaceDimension(); ++index)
                    output[index] = output[index] / norm;
                cout << "Manually normalized!" << endl;
            }

            if (YBasis)
                sprintf(outputname, "%s_kx_%d_ky_%d.%d.vec", Manager.GetString("output"), Kx, Ky, r / Qy);
            else
                sprintf(outputname, "%s_kx_%d_ky_%d.%d.vec", Manager.GetString("output"), Ky, Kx, r / Qy);
            output.WriteVector(outputname);
            
            if (Manager.GetString("exact") != NULL)
            {
                ComplexVector exact;
                char exactname[512];
                if (YBasis)
                    sprintf(exactname, "%s_kx_%d_ky_%d.%d.vec", Manager.GetString("exact"), Kx, Ky, r / Qy);
                else
                    sprintf(exactname, "%s_kx_%d_ky_%d.%d.vec", Manager.GetString("exact"), Ky, Kx, r / Qy);

                exact.ReadVector(exactname);
                double overlap = SqrNorm(exact * output);
                int precision = cout.precision();
                cout.precision(14);
                cout << "Overlap with ED = " << overlap << endl;
                cout.precision(precision);
            }
        }
    }

    cout << endl << "Overlap summary: ";
    if (!Unitary)
        cout << " Nonunitary";
    if (FlatCurvature)
        cout << " Flat";
    cout << endl;

    int precision = cout.precision();
    cout.precision(14);
    for (int s = 0; s < Qx; ++s)
    {
        int Kx = Kpx + s * Ne;
        while (Kx < 0)
            Kx += Nx;
        Kx %= Nx;
        for (int r = 0; r < Qy ; ++r)
        {
            int ysign = YBasis ? 1 : (-1);
            int Ky = C * ysign * (Kpy + (r - Delta * ysign) * Ne);
            while (Ky < 0)
                Ky += Ny;
            Ky %= Ny;

            cout << "(Kx, Ky) = (" << Kx << ", " << Ky << ")  Norm = ";

            ComplexVector v1, v2;
            double overlap = 0.;
            for (int d1 = 0; d1 < Q / (Qx * Qy); ++d1)
            {
                if (YBasis)
                    sprintf(outputname, "%s_kx_%d_ky_%d.%d.vec", Manager.GetString("output"), Kx, Ky, d1);
                else
                    sprintf(outputname, "%s_kx_%d_ky_%d.%d.vec", Manager.GetString("output"), Ky, Kx, d1);
                v1.ReadVector(outputname);
                cout << v1.Norm() << " ";
                for (int d2 = 0; d2 < Q / (Qx * Qy); ++d2)
                {
                    if (YBasis)
                        sprintf(outputname, "%s_kx_%d_ky_%d.%d.vec", Manager.GetString("exact"), Kx, Ky, d2);
                    else
                        sprintf(outputname, "%s_kx_%d_ky_%d.%d.vec", Manager.GetString("exact"), Ky, Kx, d2);
                    v2.ReadVector(outputname);

                    overlap += SqrNorm(v1 * v2);
                }
            }
            cout << "  Overlaps = " << overlap / (Q / (Qx * Qy)) << endl;
        }
    }
    cout << endl;
    cout.precision(precision);

    delete[] kxlist;
    delete[] kylist;
    delete[] xlist;
    delete[] jlist;
    delete[] lxlist;
    delete[] lylist;
    delete[] outputname;
}

// sum over all allowed permutations of kx
Complex SumKxPermutations(int* xlist, int* kxlist, int* kylist, Complex* expNx)
{
    Complex sum = 0.;
    int sign = 0;
//    int count = 0;
//    for (int ii = 0; ii < Ne; ++ii)
//        cout << kylist[ii];
//    cout << endl;
//    for (int ii = 0; ii < Ne; ++ii)
//        cout << kxlist[ii];
    while (1)
    {
//        ++count;
        int kxX = 0;
        for (int ii = 0; ii < Ne; ++ii)
            kxX += kxlist[ii] * xlist[ii];
        sum += expNx[kxX % Nx] * (1.0 - 2.0 * (sign & 1));

        int i;
        for (i = Ne - 2; i >= 0; --i)
        {
            if (kylist[i] == kylist[i + 1])
            {
                if (kxlist[i] > kxlist[i + 1])
                    break;
            }
        }
        if (i < 0)
            break;

        int maxind = i + 1;
        int max = kxlist[maxind];
        int j = i + 2;
        for (; (j < Ne) && (kylist[j] == kylist[i]); ++j)
            if ((kxlist[i] > kxlist[j]) && (kxlist[j] > max))
            {
                max = kxlist[j];
                maxind = j;
            }

        Swap(kxlist, i, maxind);
        ++sign;
        for (j = i + 1; j < Ne;)
        {
            int k;
            for (k = j + 1; (k < Ne) && (kylist[k] == kylist[j]); ++k);
            // need to reverse [i+1 .. k)
            for (int l = j; l < k - 1 - (l - j); ++l)
            {
                Swap(kxlist, l, k - 1 - (l - j));
                ++sign;
            }
            j = k;
        }
    }
//    cout << endl << count << endl << endl;
    return sum;
}

// build wannier state
// yGauge[ky][x] = e^iPhi(x,ky)
// xGauge[ky][kx] = lambda^kx / prod_kappa^kx Ax(kappa,ky)
void BuildWannier(char* onebody, ComplexMatrix& xGauge, ComplexMatrix& yGauge)
{
    // ComplexMatrix, stored as column vectors
    // m[column][row]
    // Bloch[(kx * Ny) + ky][band][orbital] = Conj(<orbital|band>)
    char* filename = new char[512];
    for (int kx = 0; kx < Nx; ++kx)
    {
        for (int ky = 0; ky < Ny; ++ky)
        {
            if (YBasis)
                sprintf(filename, "%s_kx_%d_ky_%d.onebody", onebody, kx, ky);
            else
                sprintf(filename, "%s_kx_%d_ky_%d.onebody", onebody, ky, kx);

            Bloch[(kx * Ny) + ky].ReadMatrix(filename);
            if (Nb != Bloch[(kx * Ny) + ky].GetNbrRow())
            {
                cout << "sublattice and onebody matrix has inconsistent dimension!" << endl;
                exit(1);
            }
        }
    }

    double* center = new double[Ny];
    for (int ky = 0; ky < Ny; ++ky)
    {
        center[ky] = - Arg(Lx(ky));
        if (center[ky] < 0)
            center[ky] += 2.0 * M_PI;
        center[ky] *= Nx / (2.0 * M_PI);
    }

    cout << "Wannier centers: ";
    for (int ky = 0; ky < Ny; ++ky)
        cout << center[ky] << " ";
    cout << endl;

    // find Delta
    Delta = 0;
    for (int ky = 1; ky < Ny; ++ky)
        if (center[ky] < center[0])
            ++Delta;

    // find sign of Chern #. If !YBasis, this actually gives -C, as desired.
    int count = 0;
    for (int ky = 0; ky < Ny; ++ky)
        if (center[ky] < center[(ky + 1) % Ny])
            ++count;
    if (count == Ny - 1)
        C = 1;
    else if (count == 1)
        C = -1;
    else
    {
        cout << "Wannier center locations do not resemble the lowest Landau level. " << count << endl;
        cout << "Fall back to C = 1, Delta = 0." << endl;
        C = 1;
        Delta = 0;
    }

    for (int ky = 0; ky < Ny; ++ky)
        for (int x = 0; x < Nx; ++x)
            yGauge[ky][x] = Phi(x, ky);

    for (int ky = 0; ky < Ny; ++ky)
    {
        for (int kx = 0; kx < Nx; ++kx)
        {
            Complex z = PowInt(Lx(ky), kx);
            for (int i = 0; i < kx; ++i)
                z /= Ax(i, ky);
            xGauge[ky][kx] = z;
        }
    }

    if (!Unitary)
    {
        for (int ky = 0; ky < Ny; ++ky)
        {
            double normalization = 0.;
            for (int kx = 0; kx < Nx; ++kx)
                normalization += SqrNorm(xGauge[ky][kx]);
            normalization = sqrt(Nx) / sqrt(normalization);
            for (int kx = 0; kx < Nx; ++kx)
                xGauge[ky][kx] *= normalization;
        }
    }

    delete[] filename;
    delete[] center;
}

// check parallel transport along each ky
void CheckParallelTransport(ComplexMatrix xGauge)
{
    for (int ky = 0; ky < Ny; ++ky)
    {
        Complex z = Ax(0, ky) * xGauge[ky][1] / xGauge[ky][0];
        for (int kx = 0; kx < Nx; ++kx)
        {
            int kx1 = (kx + 1) % Nx;
            if (Norm(Ax(kx, ky) * xGauge[ky][kx1] / xGauge[ky][kx] - z) > 1e-14)
                cout << "Parallel transport failed!" << endl;
        }
    }
}

// check orthonormality of the Wannier states
void CheckOrthonormality(ComplexMatrix xGauge)
{
    cout << "orthogonality breaking:" << endl;
    for (int ky = 0; ky < Ny; ++ky)
    {
        ComplexMatrix braket(Nx, Nx, true);
        double max = 0.;
        int m1, m2;
        for (int x1 = 0; x1 < Nx; ++x1)
        {
            for (int x2 = 0; x2 < Nx; ++x2)
            {
                Complex sum = 0.;
                for (int kx = 0; kx < Nx; ++kx)
                    sum += Phase(2.0 * M_PI * kx * ((double)(x1 - x2)) / Nx) * SqrNorm(xGauge[ky][kx]);
                braket[x2][x1] = sum;
                if ((x1 == x2) && fabs(Real(sum) - Nx) > 1e-14)
                {
                    cout << "Wannier states not normalized!" << endl;
                    exit(1);
                }
                if ((x1 != x2) && Norm(sum) > max)
                {
                    m1 = x1;
                    m2 = x2;
                    max = Norm(sum);
                }
            }
        }
        cout << "ky = " << ky << ", " << max << endl;
    }
}

// check <X,ky|Y|X',ky'>
void CheckWannierConnection(ComplexMatrix xGauge, ComplexMatrix yGauge)
{
    for (int j1 = 0; j1 < Nx * Ny; ++j1)
    {
        int j2 = j1 + C;
        while (j2 < 0)
            j2 += Nf;
        j2 %= Nf;

        int x1, ky1, x2, ky2;
        FindXky(j1, x1, ky1);
        FindXky(j2, x2, ky2);

        if ((ky2 - ky1 + Ny) % Ny != 1)
            cout << "ky Check FAILED" <<endl;

        Complex sum = 0.;
        for (int kx = 0; kx < Nx; ++kx)
        {
            Complex term = Phase(- 2.0 * M_PI * kx * ((double)(x2 - x1)) / Nx);
            term *= Ay(kx, ky1);
            term *= xGauge[ky2][kx] / xGauge[ky1][kx];
            sum += term;
        }
        sum *= yGauge[ky2][x2] / yGauge[ky1][x1];
        sum /= Nx;
        if (Norm(sum / Norm(sum) - TildeLy()) > 1e-14)
            cout << "reality condition broken: <" << x1 << "," << ky1 << "|Y|" << x2 << "," << ky2 << "> = " << sum << endl;
    }
}

// build sewing matrix
// sewring = e^{i xi_{kx,ky}}
void BuildSewing(char *inversion, RealVector sublattice, ComplexMatrix& sewing)
{
    ComplexMatrix p;
    p.ReadMatrix(inversion);
    if (Nb != p.GetNbrRow())
    {
        cout << "sublattice and inversion matrix has inconsistent dimension!" << endl;
        exit(1);
    }

    // inversion maps |x,y,a> to |-x-dx[a],-y-dy[a],b>
    RealVector dx(Nb, true);
    RealVector dy(Nb, true);

    for (int a = 0; a < Nb; ++a)
    {
        for (int b = 0; b < Nb; ++b)
        {
            if (Norm(p[b][a]) > 1e-14)
            {
                dx[a] = sublattice[2 * b + !YBasis] + sublattice[2 * a + !YBasis];
                dy[a] = sublattice[2 * b + YBasis] + sublattice[2 * a + YBasis];
            }
        }
    }

    for (int ky = 0; ky < Ny; ++ky)
    {
        for (int kx = 0; kx < Nx; ++kx)
        {
            int k = kx * Ny + ky;
            int mx = (Nx - kx) % Nx;
            int my = (Ny - ky) % Ny;
            int m = mx * Ny + my;

            Complex sum = 0.;
            for (int a = 0; a < Nb; ++a)
                for (int b = 0; b < Nb; ++b)
                    sum += Phase(2.0 * M_PI * (kx * ((double)dx[a]) / Nx + ky * ((double)dy[a]) / Ny)) * Bloch[k][Band][b] * p[b][a] * Conj(Bloch[m][Band][a]);
            sewing[ky][kx] = sum;
        }
    }

    // check sewing
    for (int ky = 0; ky < Ny; ++ky)
    {
        for (int kx = 0; kx < Nx; ++kx)
        {
            int mx = (Nx - kx) % Nx;
            int my = (Ny - ky) % Ny;
            if (Norm(sewing[ky][kx] * sewing[my][mx] - 1) > 1e-14)
                cout << "bad sewing matrix" << endl;
            if (fabs(Norm(sewing[ky][kx]) - 1) > 1e-14)
                cout << "bad sewing matrix" << endl;
        }
    }
}

// check inversion symmetry of |X,ky>
void CheckWannierInversion(ComplexMatrix sewing, ComplexMatrix xGauge, ComplexMatrix yGauge)
{
    Complex lambda = TildeLy();
    for (int x = 0; x < Nx; ++x)
    {
        for (int ky = 0; ky < Ny; ++ky)
        {
            int my = (Ny - ky) % Ny;
            int m;
            if ((Delta == 0 || 2 * Delta == Ny) && ky == Delta) // if Wx(ky) == 1
                m = (Nx - x) % Nx;
            else
                m = (2 * Nx - x - 1) % Nx;

            for (int kx = 0; kx < Nx; ++kx)
            {
                int mx = (Nx - kx) % Nx;
                Complex amp1 = yGauge[ky][x] * Phase(- 2.0 * M_PI * kx * ((double)x) / Nx) * xGauge[ky][kx] / sewing[ky][kx];
                Complex amp2 = yGauge[my][m] * Phase(- 2.0 * M_PI * mx * ((double)m) / Nx) * xGauge[my][mx];
                // sewing[ky][kx]
                if (Norm(sewing[0][0] / PowInt(lambda, 2 * ky) * amp1 - amp2) > 1e-14)
                {
                    cout << "Inversion broken!" << endl;
                    cout << sewing[0][0] / PowInt(lambda, 2 * ky) * amp1 << "  " << amp2 << endl;
                    if (!XLQ && Unitary && !Bold)
                        exit(1);
                }
            }
        }
    }
}

Complex PowInt(Complex z, int k)
{
    Complex prod = 1.;
    for (; k > 0; --k)
        prod *= z;
    return prod;
}

// < bra | ket > with sublattice shift
Complex Dot(ComplexVector bra, ComplexVector ket, ComplexVector shift)
{
    Complex sum = 0.;
    for (int i = 0; i < shift.GetVectorDimension(); ++i)
        sum += shift[i] * Conj(bra[i]) * ket[i];
    if (Norm(sum) < 1e-13)
    {
        cout << "Cannot make connection unitary due to orthogonality!" << endl;
        exit(1);
    }
    if (Unitary)
        sum /= Norm(sum);
    return sum;
}

// < kx, ky | x | kx+1, ky >
Complex Ax(int kx, int ky)
{
    while (kx < 0)
        kx += Nx;
    while (ky < 0)
        ky += Ny;
    int k = (kx % Nx) * Ny + (ky % Ny);
    int k1 = ((kx + 1) % Nx) * Ny + (ky % Ny);
    return Dot(Bloch[k1][Band], Bloch[k][Band], ShiftX); // Bloch[(kx * Ny) + ky][Band][orbital] = Conj(<orbital|Band>)
}

// < kx, ky | y | kx, ky+1 >
Complex Ay(int kx, int ky)
{
    while (kx < 0)
        kx += Nx;
    while (ky < 0)
        ky += Ny;
    int k = (kx % Nx) * Ny + (ky % Ny);
    int k1 = (kx % Nx) * Ny + ((ky + 1) % Ny);
    return Dot(Bloch[k1][Band], Bloch[k][Band], ShiftY); // Bloch[(kx * Ny) + ky][Band][orbital] = Conj(<orbital|Band>)
}

// Wilson loop along fixed ky
Complex Wx(int ky)
{
    Complex prod = 1.;
    for (int i = 0; i < Nx; ++i)
        prod *= Ax(i, ky);
    return prod;
}

// Wilson loop along fixed kx
Complex Wy(int kx)
{
    Complex prod = 1.;
    for (int i = 0; i < Ny; ++i)
        prod *= Ay(kx, i);
    return prod;
}

// Wilson loop around (0,ky)(kx,ky)(kx,ky+1)(0,ky+1)(0,ky)
Complex Wsq(int kx, int ky)
{
    Complex prod = Ay(kx, ky) / Ay(0, ky);
    for (int i = 0; i < kx; ++i)
        prod *= Ax(i, ky) / Ax(i, ky + 1);
    return prod;
}

// Wilson loop around a plaquette
Complex Wpq(int kx, int ky)
{
    return Ax(kx, ky) * Ay(kx + 1, ky) / (Ax(kx, ky + 1) * Ay(kx, ky));
}

// eigenvalue of projected x operator in X=0 unit cell
Complex Lx(int ky)
{
    // pow(z, 1.0 / Nx) has argument angle in (- pi / Nx, pi / Nx].
    // reason: in Complex.h, pow(Complex,double) is implemented as 
    //         pow(z,y) = exp(y*ln(|z|)) * (cos(y*Imag(ln(z))) + i*sin(y*Imag(ln(z))))
    //         and Imag(ln(z)) is implemented by atan2(Imag(z), Real(z)), which takes range (-pi,pi]
    Complex lambda = pow(Wx(ky), 1.0 / Nx);
    if (fabs(Imag(lambda)) < 1e-14) // avoid inversion symmetry breaking from round-off
        lambda = Real(lambda);
    if (Imag(lambda) > 0)
        lambda *= Phase(- 2.0 * M_PI / Nx);
    return lambda;
}

// WyU() averaged to Ny bonds, with minimal abs(arg angle)
Complex TildeLy(void)
{
    Complex lambda = Extreme ? pow(PowInt(WyU(), Nx), 1.0 / Nf) : pow(WyU(), 1.0 / Ny);
    if (fabs(Imag(lambda)) < 1e-14)
        lambda = Real(lambda);
    return lambda;
}

// pull ky back to the principal BZ, defined by C * ky + Delta in [0, Ny)
int Principal(int ky)
{
    while (C * ky + Delta < 0)
        ky += C * Ny;
    ky = ((C * ky + Delta) % Ny - Delta) * C;
    return ky;
}

// j = X * Ny + C * ky + delta
int J(int X, int ky)
{
    while (X < 0)
        X += Nx;
    X %= Nx;
    ky = Principal(ky);
    return X * Ny + C * ky + Delta;
}

// inverse map of j = X * Ny + C * ky + delta, C * ky + delta \in [0, Ny)
// pull ky back to [0, Ny) after the mapping
void FindXky(int j, int& X, int& ky)
{
    while (j < 0)
        j += Nx * Ny;
    j %= Nx * Ny;
    X = j / Ny;
    ky = (j % Ny - Delta) * C;

    while (ky < 0)
        ky += Ny;
    ky %= Ny;
}

// "average" Wilson loop around (0,ky)(kx,ky)(kx,ky+1)(0,ky+1)(0,ky)
Complex Wsqbar(int kx, int ky)
{
    ky = Principal(ky);
    Complex prod = PowInt(Lx(ky) / Lx(ky + 1), kx);
    if ((C * ky + Delta + C < 0) || (C * ky + Delta + C >= Ny))
        prod *= Phase(2.0 * M_PI * C * ((double) kx) / Nx);
    return prod;
}

// unitary measure of curvature fluctuation at ky
Complex Ux(int ky)
{
    if (FlatCurvature)
        return 1.;
    Complex sum = 0.;
    for (int i = 0; i < Nx; ++i)
        sum += Wsq(i, ky) / Wsqbar(i, ky);
    sum /= Nx;
    sum /= Norm(sum);
    return sum;
}

// modified Wilson loop along kx = 0
Complex WyU(void)
{
    Complex prod = Wy(0);
    for (int i = 0; i < Ny; ++i)
        prod *= Ux(i);
    return prod;
}

// phase in front of |X,ky>
Complex Phi(int X, int ky)
{
    // The recursive construction starts from ky = -C delta
    // This only differs from the convention in paper (starts from ky = 0) by a trivial overall phase factor.
    Complex prod = 1.;
    Complex lambda = XLQ ? 1. : TildeLy();
    ky = Principal(ky);
    if (C > 0)
    {
        for (int i = - Delta; i < ky; ++i)
            prod *= lambda / (Ay(0, i) * Ux(i));
        prod *= Extreme ? PowInt(PowInt(lambda, Ny) / WyU(), X) : 1.;
        prod /= Norm(prod); // in case Unitary = False
        return prod;
    }
    else
    {
        for (int i = ky; i < Delta; ++i)
            prod *= lambda / (Ay(0, i) * Ux(i));
        prod *= Extreme ? PowInt(PowInt(lambda, Ny) / WyU(), Nx - X - 1) : 1.;
        prod /= Norm(prod); // in case Unitary = False
        return 1. / prod;
    }
}

// circular shift. corresponds to |{j - shift}> 
unsigned long Shift(unsigned long state, int shift, int length)
{
    // shift must belong to [0, length) !!!
    unsigned long mask = (1ul << shift) - 1;
    return (state >> shift) | ((state & mask) << (length - shift));
}

// print binary rep of state
void Print(unsigned long state, int length)
{
    for (int i = length - 1; i >= 0; --i) // programmer-friendly
        cout << ((state >> i) & 1ul);
//    for (int i = 0; i < length; ++i) // human-friendly
//        cout << ((state >> i) & 1ul);
}

// count the number of 1's
int Count(unsigned long state)
{
    int count;
    for (count = 0; state; ++count)
        state &= state - 1ul;
    return count;
}

// sort in descending order
// returns the perumtation sign
double SortDescending(int* list)
{
    // FIXME: plain stupid algorithm... But Ne is small anyways..
    int count = 0;

    for (int i = 0; i < Ne - 1; ++i)
    {
        int m = i + 1;
        for (int j = i + 2; j < Ne; ++j)
        {
            if (list[j] > list[m])
                m = j;
            else if (list[j] == list[m])
                return 0.;
        }
        if (list[i] == list[m])
            return 0.;
        if (list[m] > list[i])
        {
            ++count;
            int temp = list[m];
            list[m] = list[i];
            list[i] = temp;
        }
    }

    for (int i = 0; i < Ne - 1; ++i)
        if (list[i] <= list[i + 1])
        {
            cout << "sorting failed" << endl;
            exit(1);
        }
    return 1.0 - 2.0 * (count & 1);
}

// obtain the list of (kx,ky) from the binary representation, in descending order of kxNy+ky
// don't check anything. assume Count(state) is the length of klists
void Unflatten(unsigned long state, int* kxlist, int* kylist)
{
    if (YBasis)
    {
        for (int i = Nf - 1; i >= 0; --i)
        {
            if ((state >> i) & 1ul)
            {
                *(kxlist++) = i / Ny;
                *(kylist++) = i % Ny;
            }
        }
    }
    else
    {
        // kx, ky are not swapped in the binary representation
        for (int i = Nf - 1; i >= 0; --i)
        {
            if ((state >> i) & 1ul)
            {
                *(kylist++) = i / Nx;
                *(kxlist++) = i % Nx;
            }
        }
    }
}

// obtain the list of (x,ky) from the binary representation of j, in descending order of j
// don't check anything. assume Count(state) is the length of klists
void UnflattenJ(unsigned long state, int* xlist, int* kylist)
{
    for (int j = Nf - 1; j >= 0; --j)
    {
        if ((state >> j) & 1ul)
        {
            *(xlist++) = j / Ny;
            *(kylist++) = ((j % Ny - Delta) * C + 10 * Ny) % Ny;
        }
    }
}

// sort {kx,ky} according to ky * Nx + kx. return the permutation sign
double SortYMajor(int* kxlist, int* kylist)
{
    int *klist = new int[Ne];
    for (int i = 0; i < Ne; ++i)
        klist[i] = kylist[i] * Nx + kxlist[i];
    double sign = SortDescending(klist);
    for (int i = 0; i < Ne; ++i)
    {
        kylist[i] = klist[i] / Nx;
        kxlist[i] = klist[i] % Nx;
    }
    delete[] klist;
    return sign;
}

// calculate sum(Ne ** ky for ky in kylist)
unsigned long HistSum(int* kylist)
{
    unsigned long sum = 0ul;
    for (int i = 0; i < Ne; ++i)
    {
        unsigned long term = 1ul;
        for (int j = 0; j < kylist[i]; ++j)
            term *= (unsigned long)Ne;
        sum += term;
    }
    return sum;
}

// print 2D list
void Print2DList(int* xlist, int* ylist)
{
    for (int i = 0; i < Ne; ++i)
        cout << "(" << xlist[i] << "," << ylist[i] << ")";
}

// swap two elements in a list
void Swap(int* list, int i, int j)
{
    int tmp = list[i];
    list[i] = list[j];
    list[j] = tmp;
}

