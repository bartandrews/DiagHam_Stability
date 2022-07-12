#include "Options/Options.h"
#include "Architecture/ArchitectureManager.h"
#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"

#include "MathTools/LongRational.h"
#include "Matrix/RealMatrix.h"

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
    OptionManager Manager("FQHEMPSCFTMatrix" , "0.01");

    OptionGroup* SystemGroup = new OptionGroup("system options");
    Manager += SystemGroup;
    (*SystemGroup) += new SingleIntegerOption('k', "k-index", "k as in (k, r)", 2);
    (*SystemGroup) += new SingleIntegerOption('r', "r-index", "r as in (k, r)", 2);
    (*SystemGroup) += new SingleIntegerOption('\n', "p-truncation", "truncation level", 1);
    (*SystemGroup) += new BooleanOption('\n', "use-nonrational", "use double numbers instead of rational numbers for CFT calcaultions");
    (*SystemGroup) += new SingleStringOption('\n', "matrices-cft", "optional directory where the geomerty independent CFT matrices are stored");

    ArchitectureManager Architecture;
    Architecture.AddOptionGroup(&Manager);

    OptionGroup* MiscGroup = new OptionGroup("misc options");
    Manager += MiscGroup;
    (*MiscGroup) += new BooleanOption('h', "help", "display this help");

    if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
        cout << "see man page for option syntax or type FQHEMPSCFTMatrix -h" << endl;
        return -1;
    }
    if (Manager.GetBoolean("help") == true)
    {
        Manager.DisplayHelp(cout);
        return 0;
    }

    int k = Manager.GetInteger("k-index");
    int r = Manager.GetInteger("r-index");
    int pLevel = Manager.GetInteger("p-truncation");
    bool useRational = !Manager.GetBoolean("use-nonrational");
    char* cftDirectory = Manager.GetString("matrices-cft");

    if (k == 2)
    {
        LongRational centralCharge((r + 2l) - (2l * (r - 1l) * (r - 1l)), r + 2l);
        if (r == 2)
        {
            FQHEMPSClustered2RMatrix M(pLevel, centralCharge, "pfaffian", useRational);

            int nbrSectors = 3;
            char** sectorNames = new char*[3];
	    sectorNames[0] = new char[16];
	    sprintf (sectorNames[0], "identity");
	    sectorNames[1] = new char[16];
	    sprintf (sectorNames[1], "psi");
	    sectorNames[2] = new char[16];
	    sprintf (sectorNames[2], "sigma");
            LongRational weights[3] = {LongRational(0l, 1l), LongRational(1l, 2l), LongRational(1l, 16l)};

            M.ComputeMatrixElements(cftDirectory, Architecture.GetArchitecture(), "identity", LongRational(0l, 1l), nbrSectors, sectorNames, weights, RealMatrix());

            RealMatrix psi(nbrSectors, nbrSectors, true);
            psi.SetMatrixElement(1, 0, 1.0); // <ψ|ψ|1>
            psi.SetMatrixElement(0, 1, 1.0); // <1|ψ|ψ>
            psi.SetMatrixElement(2, 2, 1 / sqrt(2.0)); // <σ|ψ|σ>
            M.ComputeMatrixElements(cftDirectory, Architecture.GetArchitecture(), "psi", LongRational(1l, 2l), nbrSectors, sectorNames, weights, psi);

            RealMatrix sigma(nbrSectors, nbrSectors, true);
            sigma.SetMatrixElement(2, 0, 1.0); // <σ|σ|1>
            sigma.SetMatrixElement(2, 1, 1 / sqrt(2.0)); // <σ|σ|ψ>
            sigma.SetMatrixElement(0, 2, 1.0); // <1|σ|σ>
            sigma.SetMatrixElement(1, 2, 1 / sqrt(2.0)); // <ψ|σ|σ>
            M.ComputeMatrixElements(cftDirectory, Architecture.GetArchitecture(), "sigma", LongRational(1l, 16l), nbrSectors, sectorNames, weights, sigma);

            return 0;
        }
        else if (r == 3)
        {
            FQHEMPSClustered2RMatrix M(pLevel, centralCharge, "gaffnian", useRational);
            int nbrSectors = 4;
            char** sectorNames = new char*[4];
	    sectorNames[0] = new char[16];
	    sprintf (sectorNames[0], "identity");
	    sectorNames[1] = new char[16];
	    sprintf (sectorNames[1], "psi");
	    sectorNames[2] = new char[16];
	    sprintf (sectorNames[2], "sigma");
 	    sectorNames[3] = new char[16];
	    sprintf (sectorNames[3], "phi");
            LongRational weights[4] = {LongRational(0l, 1l), LongRational(3l, 4l), LongRational(-1l, 20l), LongRational(1l, 5l)};

            M.ComputeMatrixElements(cftDirectory, Architecture.GetArchitecture(), "identity", LongRational(0l, 1l), nbrSectors, sectorNames, weights, RealMatrix());

            RealMatrix psi(nbrSectors, nbrSectors, true);
            psi.SetMatrixElement(1, 0, 1.0); // <ψ|ψ|1>
            psi.SetMatrixElement(0, 1, 1.0); // <1|ψ|ψ>
            psi.SetMatrixElement(3, 2, 1 / sqrt(2.0)); // <φ|ψ|σ>
            psi.SetMatrixElement(2, 3, 1 / sqrt(2.0)); // <σ|ψ|φ>
            M.ComputeMatrixElements(cftDirectory, Architecture.GetArchitecture(), "psi", LongRational(3l, 4l), nbrSectors, sectorNames, weights, psi);

            RealMatrix sigma(nbrSectors, nbrSectors, true);
            sigma.SetMatrixElement(2, 0, 1.0); // <σ|σ|1>
            sigma.SetMatrixElement(3, 1, 1 / sqrt(2.0)); // <φ|σ|ψ>
            sigma.SetMatrixElement(0, 2, 1.0); // <1|σ|σ>
            sigma.SetMatrixElement(3, 2, 1.0); // <φ|σ|σ> // FIXME
            sigma.SetMatrixElement(1, 3, 1 / sqrt(2.0)); // <ψ|σ|φ>
            sigma.SetMatrixElement(2, 3, 1.0); // <σ|σ|φ> // FIXME
            M.ComputeMatrixElements(cftDirectory, Architecture.GetArchitecture(), "sigma", LongRational(-1l, 20l), nbrSectors, sectorNames, weights, sigma);

            return 0;
        }
    }
    else if (r == 2)
    {
        FQHEMPSClustered2RMatrix M(pLevel, LongRational(4l, 5l), "readrezayi", useRational);

        int nbrSectors = 6;
        //                          0         1     2      3          4       5
	char** sectorNames = new char*[6];
	sectorNames[0] = new char[16];
	sprintf (sectorNames[0], "identity");
	sectorNames[1] = new char[16];
	sprintf (sectorNames[1], "psi");
	sectorNames[2] = new char[16];
	sprintf (sectorNames[2], "w");
	sectorNames[3] = new char[16];
	sprintf (sectorNames[3], "epsilon");
	sectorNames[4] = new char[16];
	sprintf (sectorNames[4], "sigma");
	sectorNames[5] = new char[16];
	sprintf (sectorNames[5], "phi");
        LongRational weights[6] = {LongRational(0l, 1l), LongRational(2l, 3l), LongRational(3l, 1l),
            LongRational(2l, 5l), LongRational(1l, 15l), LongRational(7l, 5l)};

        double sqrtc = sqrt(0.5) * exp(0.25 * (lgamma(0.2) + 3 * lgamma(0.6) - lgamma(0.8) - 3 * lgamma(0.4)));

        M.ComputeMatrixElements(cftDirectory, Architecture.GetArchitecture(), "identity", LongRational(0l, 1l), nbrSectors, sectorNames, weights, RealMatrix());

        RealMatrix psi(nbrSectors, nbrSectors, true);
        psi.SetMatrixElement(1, 0, 1.0); // <ψ₁|ψ₁|1>
        psi.SetMatrixElement(1, 1, 2 / sqrt(3.0)); // <ψ₂|ψ₁|ψ₁>
        psi.SetMatrixElement(0, 1, 1.0); // <1|ψ₁|ψ₂>
        psi.SetMatrixElement(2, 1, - sqrt(26.0) / 9); // <W|ψ₁|ψ₂>
        psi.SetMatrixElement(1, 2, sqrt(26.0) / 9); // <ψ₁|ψ₁|W>
        psi.SetMatrixElement(4, 3, sqrt(2.0 / 3)); // <σ₂|ψ₁|ε>
        psi.SetMatrixElement(3, 4, sqrt(2.0 / 3)); // <ε|ψ₁|σ₁>
        psi.SetMatrixElement(5, 4, - sqrt(7.0 / 2) / 3); // <φ|ψ₁|σ₁>
        psi.SetMatrixElement(4, 4, 1.0 / sqrt(3.0)); // <σ₁|ψ₁|σ₂>
        psi.SetMatrixElement(4, 5, sqrt(7.0 / 2) / 3); // <σ₂|ψ₁|φ>
        M.ComputeMatrixElements(cftDirectory, Architecture.GetArchitecture(), "psi", LongRational(2l, 3l), nbrSectors, sectorNames, weights, psi);

        RealMatrix sigma(nbrSectors, nbrSectors, true);
        sigma.SetMatrixElement(4, 0, 1.0); // <σ₁|σ₁|1>
        sigma.SetMatrixElement(3, 1, sqrt(2.0 / 3)); // <ε|σ₁|ψ₁>
        sigma.SetMatrixElement(5, 1, - sqrt(7.0 / 2) / 3); // <φ|σ₁|ψ₁>
        sigma.SetMatrixElement(4, 1, 1 / sqrt(3.0)); // <σ₂|σ₁|ψ₂>
        sigma.SetMatrixElement(4, 2, 1 / (9 * sqrt(26.0))); // <σ₁|σ₁|W>
        sigma.SetMatrixElement(1, 3, sqrt(2.0 / 3)); // <ψ₂|σ₁|ε>
        sigma.SetMatrixElement(4, 3, sqrtc); // <σ₁|σ₁|ε>
        sigma.SetMatrixElement(1, 4, 1 / sqrt(3.0)); // <ψ₁|σ₁|σ₁>
        sigma.SetMatrixElement(4, 4, sqrt(2.0) * sqrtc); // <σ₂|σ₁|σ₁>
        sigma.SetMatrixElement(0, 4, 1.0); // <1|σ₁|σ₂>
        sigma.SetMatrixElement(2, 4, - 1 / (9 * sqrt(26.0))); // <W|σ₁|σ₂>
        sigma.SetMatrixElement(3, 4, sqrtc); // <ε|σ₁|σ₂>
        sigma.SetMatrixElement(5, 4, sqrtc / sqrt(21.0)); // <φ|σ₁|σ₂>
        sigma.SetMatrixElement(1, 5, sqrt(7.0 / 2) / 3); // <ψ₂|σ₁|φ>
        sigma.SetMatrixElement(4, 5, - sqrtc / sqrt(21.0)); // <σ₁|σ₁|φ>
        M.ComputeMatrixElements(cftDirectory, Architecture.GetArchitecture(), "sigma", LongRational(1l, 15l), nbrSectors, sectorNames, weights, sigma);

        return 0;
    }

    FQHEMPSMatrixManager MPSMatrixManager(false); // need this dummy call to link properly...
    return -1;
}

