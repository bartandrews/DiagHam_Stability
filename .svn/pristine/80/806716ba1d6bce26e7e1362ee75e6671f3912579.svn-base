#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"

#include "HilbertSpace/FermionOnSphereMPSAlternativeWrapper.h"

#include "FunctionBasis/ParticleOnDiskFunctionBasis.h"
#include "FunctionBasis/ParticleOnCylinderFunctionBasis.h"

#include "Matrix/SparseComplexMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Options/Options.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

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
    cout.precision(14);
    OptionManager Manager ("FQHESphereMPSDensity" , "0.01");
    OptionGroup* MiscGroup = new OptionGroup ("misc options");

    ArchitectureManager Architecture;
    FQHEMPSMatrixManager MPSMatrixManager;

    MPSMatrixManager.AddOptionGroup(&Manager);
    OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
    OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
    OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
    Architecture.AddOptionGroup(&Manager);
    Manager += MiscGroup;

    (*SystemGroup) += new SingleIntegerOption('p', "nbr-particles", "number of particles", 2);
    (*SystemGroup) += new SingleStringOption('\n', "location-file", "evaluate density at locations given in this text file");
    (*SystemGroup) += new SingleStringOption('\n', "output-file", "write density values to this text file");
    (*MiscGroup) += new BooleanOption('h', "help", "display this help");

    if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
        cout << "see man page for option syntax or type FQHESphereMPSDensity -h" << endl;
        return -1;
    }
    if (Manager.GetBoolean("help") == true)
    {
        Manager.DisplayHelp (cout);
        return 0;
    }

    int NbrParticles = Manager.GetInteger("nbr-particles"); 
    bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");

    int NbrQuasiholes = 0;
    Complex* QuasiholePositions = 0;
    if (Manager.GetString("with-quasiholes") != 0)
    {
        MultiColumnASCIIFile InputQuasiholePosition;
        if (InputQuasiholePosition.Parse(Manager.GetString("with-quasiholes")) == false)
        {
            InputQuasiholePosition.DumpErrors(cout) << endl;
            return -1;
        }
        QuasiholePositions = InputQuasiholePosition.GetAsComplexArray(0);
        NbrQuasiholes = InputQuasiholePosition.GetNbrLines();

        cout << "coordinates of " << NbrQuasiholes << " quasihole(s): (before cylinder exp, if at all)" << endl;
        for (int i = 0; i < NbrQuasiholes; ++i)
            cout << QuasiholePositions[i] << endl;
    }

    int NbrFluxQuanta = Manager.GetInteger("laughlin-index") * (NbrParticles - 1) + NbrQuasiholes;

    double AspectRatio = Manager.GetDouble("aspect-ratio");
    double Kappa = 0.0;
    double Perimeter = 0.0;
    double H;
    if (Manager.GetDouble("cylinder-perimeter") > 0.0)
    {
        Kappa = (2.0 * M_PI) / Manager.GetDouble("cylinder-perimeter");
        Perimeter = Manager.GetDouble("cylinder-perimeter");
        H =  2.0 * M_PI * (NbrFluxQuanta + 1.0)/Perimeter;
        AspectRatio = Perimeter/H;
    }
    else
    {
        Kappa = (2.0 * M_PI)/sqrt(2.0 * M_PI * (NbrFluxQuanta + 1) * AspectRatio);
        Perimeter = sqrt(2.0 * M_PI * (NbrFluxQuanta + 1) * AspectRatio); 
        H = sqrt(2.0 * M_PI * (NbrFluxQuanta + 1.0))/sqrt(AspectRatio);
    }
    cout << "Cylinder geometry: perimeter = " << Perimeter << " H = " << H << ", kappa = " << Kappa << endl;


    int LandauLevel = 0;
    AbstractFunctionBasis* Basis;
    if (CylinderFlag)
        Basis = new ParticleOnCylinderFunctionBasis(NbrFluxQuanta, LandauLevel, AspectRatio, false); // do not add shift to index
    else
        Basis = new ParticleOnDiskFunctionBasis(NbrFluxQuanta);

    AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta); 
    SparseRealMatrix* SparseBMatrices = MPSMatrix->GetMatrices();
    int NbrQuasiholeBMatrices = 0;
    SparseComplexMatrix* SparseQuasiholeBMatrices = 0;
    if (NbrQuasiholes > 0)
    {
        SparseQuasiholeBMatrices = MPSMatrix->GetQuasiholeMatrices(NbrQuasiholes, QuasiholePositions);
        NbrQuasiholeBMatrices = 1; // single edge matrix
    }

    int MPSRowIndex = 0;
    int MPSColumnIndex = 0;
    MPSMatrix->GetMatrixBoundaryIndices(MPSRowIndex, MPSColumnIndex, 0); // never use padding
    cout << "MPS boundary indices (" << MPSRowIndex << "," << MPSColumnIndex << ")" << endl;

    FermionOnSphereMPSAlternativeWrapper* SpaceWrapper = 0;

    SpaceWrapper = new FermionOnSphereMPSAlternativeWrapper(NbrParticles, NbrFluxQuanta,
            MPSRowIndex, MPSColumnIndex, !CylinderFlag,
            SparseBMatrices, SparseQuasiholeBMatrices, NbrQuasiholeBMatrices, Architecture.GetArchitecture());

    Complex** OrbitalDensityMatrix = new Complex*[NbrFluxQuanta + 1];	  
    for (int i = 0; i <= NbrFluxQuanta; ++i)
        OrbitalDensityMatrix[i] = new Complex[NbrFluxQuanta + 1];

    // square norm of the unnormalized one-body wave-function, 2 * pi * 2^j * j!
    double* SqrNormAnnulus = new double[NbrFluxQuanta + 1];
    SqrNormAnnulus[0] = 2 * 3.1415926535897932385;
    for (int j = 1; j <= NbrFluxQuanta; ++j)
        SqrNormAnnulus[j] = SqrNormAnnulus[j - 1] * 2 * j;

    cout << endl;
    cout << "density matrix in orbital basis: < MPS | c^_m c_n | MPS >" << endl;
    for (int i = 0; i <= NbrFluxQuanta; ++i)
    {
        for (int j = 0; j <= NbrFluxQuanta; ++j)
        {
            Complex TmpCoef;
            SpaceWrapper->AdA(0, i, j, TmpCoef);
            OrbitalDensityMatrix[i][j] = TmpCoef;

            if (!CylinderFlag) // unnormalized
                OrbitalDensityMatrix[i][j] /= sqrt(SqrNormAnnulus[i] * SqrNormAnnulus[j]);

            cout << "(" << i << "," << j << ") = " << OrbitalDensityMatrix[i][j] << endl;
        }
    }
    cout << "trace should be equal to " << NbrParticles << endl;
    cout << endl;

    int NbrLocations = 0;
    Complex* Locations;
    if (Manager.GetString("location-file") == 0)
    {
        NbrLocations = 1;
        Locations = new Complex;
        Locations[0] = 0;
    }
    else
    {
        MultiColumnASCIIFile LocationFile;
        if (LocationFile.Parse(Manager.GetString("location-file")) == false)
        {
            LocationFile.DumpErrors(cout) << endl;
            return -1;
        }
        Locations = LocationFile.GetAsComplexArray(0);
        NbrLocations = LocationFile.GetNbrLines();
    }

    ofstream File;
    if (Manager.GetString("output-file") != 0)
    {
        File.open(Manager.GetString("output-file"), ios::out);
        File.precision(14);
    }

    cout << "evaluate density at " << NbrLocations << " location(s): (before cylinder exp, if at all)" << endl;
    for (int i = 0; i < NbrLocations; ++i)
    {
        Complex Sum = 0;
        double x = Locations[i].Re;
        double y = Locations[i].Im;
        for (int m = 0; m <= NbrFluxQuanta; ++m)
        {
            Complex bra = Conj(Basis->GetFunctionValue(x, y, m));
            for (int n = 0; n <= NbrFluxQuanta; ++n)
            {
                Complex ket = Basis->GetFunctionValue(x, y, n);
                Sum += bra * OrbitalDensityMatrix[m][n] * ket;
            }
        }
        if (!CylinderFlag) // ParticleOnDiskFunctionBasis does not provide the Gaussian factor
            Sum *= exp(- 0.5 * (x * x + y * y));
        cout << "(" << x << "," << y << ") " << Sum << endl;

        if (Manager.GetString("output-file") != 0)
            File << Sum.Re << endl;
    }
    File.close();

    delete Basis;
    delete SpaceWrapper;
    for (int i = 0; i <= NbrFluxQuanta; ++i)
        delete[] OrbitalDensityMatrix[i];
    delete[] OrbitalDensityMatrix;
    return 0;
}
