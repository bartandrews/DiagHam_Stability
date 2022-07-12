#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <cassert>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;
//new namespaces (added by ba340)
using std::setw;

int main ( int argc, char** argv )
{
    cout.precision ( 14 );

    // some running options and help
    OptionManager Manager ( "FCIHofstadterCorrelation" , "0.01" );
    OptionGroup* SystemGroup = new OptionGroup ( "system options" );
    OptionGroup* PlotOptionGroup = new OptionGroup ( "plot options" );
    OptionGroup* PrecalculationGroup = new OptionGroup ( "precalculation options" );
    OptionGroup* MiscGroup = new OptionGroup ( "misc options" );

    ArchitectureManager Architecture;

    Manager += SystemGroup;
    Manager += PlotOptionGroup;
    Architecture.AddOptionGroup ( &Manager );
    Manager += PrecalculationGroup;
    Manager += MiscGroup;

    ( *SystemGroup ) += new SingleStringOption ( '\0', "state", "name of the vector file describing the state whose density has to be plotted" );
    ( *SystemGroup ) += new BooleanOption ( '\n', "density", "plot density instead of density-density correlation", false );
    ( *SystemGroup ) += new BooleanOption ( '\n', "k-space", "compute the density/correlation in momentum space", false );
    ( *SystemGroup ) += new BooleanOption ( '\n', "subtract-density", "subtract the average density from the density-density correlation", false );
    ( *SystemGroup ) += new SingleIntegerOption ( 'x', "reference-x", "x-coordinate or reference site", 0 );
    ( *SystemGroup ) += new SingleIntegerOption ( 'y', "reference-y", "y-coordinate or reference site", 0 );

    ( *PlotOptionGroup ) += new SingleStringOption ( '\n', "output", "output file name (default output name replace the .vec extension of the input file with .rho or .rhorho)", 0 );
    ( *PlotOptionGroup ) += new SingleIntegerOption ( '\n', "nbr-samplesx", "number of samples along the x direction", 30, true, 3 );
    ( *PlotOptionGroup ) += new SingleIntegerOption ( '\n', "nbr-samplesy", "number of samples along the y direction", 30, true, 3 );

    ( *MiscGroup ) += new BooleanOption ( 'h', "help", "display this help" );

    if ( Manager.ProceedOptions ( argv, argc, cout ) == false )
    {
        cout << "see man page for option syntax or type FCICorrelation -h" << endl;
        return -1;
    }
    if ( Manager.GetBoolean ( "help" ) == true )
    {
        Manager.DisplayHelp ( cout );
        return 0;
    }
    if ( Manager.GetString ( "state" ) == 0 )
    {
        cout << "FCICorrelation requires an input state" << endl;
        return -1;
    }
    if ( IsFile ( Manager.GetString ( "state" ) ) == false )
    {
        cout << "can't find vector file " << Manager.GetString ( "state" ) << endl;
        return -1;
    }

    int NbrParticles = 0;
    int NbrFluxQuanta = 0;
    int Momentum = 0;
    double Ratio = 0;
    bool Statistics = false;
    
    int NbrSamplesX = Manager.GetInteger ( "nbr-samplesx" );
    int NbrSamplesY = Manager.GetInteger ( "nbr-samplesy" );
    bool DensityFlag = Manager.GetBoolean ( "density" );
    bool SubtractDensityFlag = Manager.GetBoolean ( "subtract-density" );
    
    if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, NbrFluxQuanta, Momentum, Ratio, Statistics)==false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
      return -1;
    }
   
    cout << setw ( 20 ) << std::left << "NbrParticles" << setw ( 20 ) << std::left << NbrParticles << endl;
    cout << setw ( 20 ) << std::left << "NbrFluxQuanta" << setw ( 20 ) << std::left << NbrFluxQuanta << endl;
    cout << setw ( 20 ) << std::left << "Momentum" << setw ( 20 ) << std::left << Momentum << endl;
    cout << setw ( 20 ) << std::left << "Ratio" << setw ( 20 ) << std::left << Ratio << endl;
    cout << setw ( 20 ) << std::left << "Statistics" << setw ( 20 ) << std::left << Statistics << endl;
    
    ParticleOnTorus* Space = 0;
    if ( Statistics == true )
        Space = new FermionOnTorus ( NbrParticles, NbrFluxQuanta, Momentum );
    else
    {

#ifdef  __64_BITS__
        if ( ( NbrFluxQuanta + NbrParticles - 1 ) < 63 )
#else
        if ( ( NbrFluxQuanta + NbrParticles - 1 ) < 31 )
#endif
        {
            Space = new BosonOnTorusShort ( NbrParticles, NbrFluxQuanta, Momentum );
        }
        else
        {
            Space = new BosonOnTorus ( NbrParticles, NbrFluxQuanta, Momentum );
        }
        cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;
    }
    
    RealVector* RealState = new RealVector();
    
    if ( RealState->ReadVector ( Manager.GetString ( "state" ) ) == false )
    {
        cout << "can't open vector file " << Manager.GetString ( "state" ) << endl;
        return -1;
    }
    Complex* PrecalculatedValues_rho = 0;
    Complex* PrecalculatedValues_rhorho = 0;
    int* PrecalculatedIndices = 0;
    int NbrPrecalculatedValues = 0;

    const double L_x = sqrt(2*M_PI*NbrFluxQuanta / Ratio );
    const double invL_x = 1.0 / L_x;
    const double L_y = sqrt(2*M_PI*NbrFluxQuanta*Ratio);
    const double invL_y = 1.0 / L_y;
    const double L_x_scaled = L_x * NbrSamplesX;
    const double L_y_scaled = L_y * NbrSamplesY;
    const double invNbrFluxQuanta=1.0/(double)NbrFluxQuanta;
    
    if ( DensityFlag == false )
    {
        for ( int j1 =0; j1 < NbrFluxQuanta; ++j1 )
            for ( int j2 =0; j2 < NbrFluxQuanta; ++j2 )
                for ( int j3 =0; j3 < NbrFluxQuanta; ++j3 )
                    for ( int j4 =0; j4 < NbrFluxQuanta; ++j4 )
                    {
                        if ( ( ( j1 + j2 - j3 - j4 ) % NbrFluxQuanta ) == 0 )
                        {
                            ++NbrPrecalculatedValues;
                        }
                    }
                    
        PrecalculatedValues_rhorho = new Complex [NbrPrecalculatedValues];
        PrecalculatedIndices = new int [4 * NbrPrecalculatedValues];
        NbrPrecalculatedValues = 0;
        for ( int j1 =0; j1 < NbrFluxQuanta; ++j1 )
            for ( int j2 =0; j2 < NbrFluxQuanta; ++j2 )
                for ( int j3 =0; j3 < NbrFluxQuanta; ++j3 )
                    for ( int j4 =0; j4 < NbrFluxQuanta; ++j4 )
                    {
                        if ( ( ( j1 + j2 - j3 - j4 ) % NbrFluxQuanta ) == 0 )
                        {
                            ParticleOnSphereDensityDensityOperator Operator ( Space, j1, j2, j3, j4 );
                            PrecalculatedValues_rhorho[NbrPrecalculatedValues] = Operator.MatrixElement ( *RealState, *RealState );
                            PrecalculatedIndices[4*NbrPrecalculatedValues] = j1;
                            PrecalculatedIndices[4*NbrPrecalculatedValues + 1] = j2;
                            PrecalculatedIndices[4*NbrPrecalculatedValues + 2] = j3;
                            PrecalculatedIndices[4*NbrPrecalculatedValues + 3] = j4;
                            ++NbrPrecalculatedValues;
                        }
                    }
    }
    
    delete RealState;
    
    delete Space;
    ofstream File;
    File.precision ( 14 );
    RealVector Position ( 2, true );
    RealVector Position_start ( 2, true );
    RealVector q(2,true);
    if ( Manager.GetString ( "output" ) != 0 )
        File.open ( Manager.GetString ( "output" ), ios::binary | ios::out );
    else
    {
        char* TmpFileName = 0;
        if ( DensityFlag == false )
        {
            TmpFileName = ReplaceExtensionToFileName ( Manager.GetString ( "state" ), "vec", "rhorho" );
        }
        else
        {
            TmpFileName = ReplaceExtensionToFileName ( Manager.GetString ( "state" ), "vec", "rho" );
        }
        if ( TmpFileName == 0 )
        {
            cout << "no vec extension was find in " << Manager.GetString ( "state" ) << " file name" << endl;
            return 0;
        }
        File.open ( TmpFileName, ios::binary | ios::out );
        delete[] TmpFileName;
    }

//     if ( Manager.GetBoolean ( "k-space" ) == true )
//     {
//         if ( DensityFlag == true )
//         {
//             File << "# j n(j)" << endl;
//             for ( int j =0; j < NbrCellX; ++j )
//             {
//                 File << j << " " << PrecalculatedValues_rhorho[j].Re << endl;
//             }
//             File.close();
//         }
//         return 0;
//     }

    //
    Complex* Coefficients = new Complex[NbrFluxQuanta];
    Complex* Coefficients2 = new Complex[NbrFluxQuanta];
    Position_start[0] = Manager.GetInteger ( "reference-x" );
    Position_start[1] = Manager.GetInteger ( "reference-y" );

    //routine to calculate the two-particle correlation function of the form <psi|n_i n_0|psi>
    //
    const double Normalisation = ( 1.0/ ( double ) ( NbrFluxQuanta ) ); //normalisation factor for the <c^+ c> term i.e. 1/N_c
    const double Normalisation2 = ( 1.0/ ( double ) ( NbrFluxQuanta*NbrFluxQuanta ) ); //normalisation factor for the <c^+ c^+ c c> term i.e. 1/N_c^2
    const double DensityPrefactor = ( ( double ) ( NbrFluxQuanta*NbrFluxQuanta ) / ( double ) NbrParticles );
    const double CorrelationTorusPrefactor = (1.0 / (( double ) (NbrParticles)*(double)(NbrParticles-1)));
    const double AverageDensity = ( ( double ) NbrParticles / ( double ) ( NbrFluxQuanta ) );
    
    double Correlations[NbrSamplesX][NbrSamplesY];
    Complex summand = 0;
    int j1=0,j2=0,j3=0,j4=0;
    
    const int s_cutoff = ( int ) ( ( L_x / 2*M_PI ) *sqrt ( 2*34.5 ) );
    const int t_cutoff = ( int ) ( ( L_y / 2*M_PI ) *sqrt ( 2*34.5 ) );
    
    for ( int x=0; x<NbrSamplesX; ++x )
    {
        for ( int y=0; y<NbrSamplesY; ++y )
        {

            Position[0]=fmod ( Position_start[0] + ( double ) x*L_x/double ( NbrSamplesX),L_x );
            Position[1]=fmod ( Position_start[1] + ( double ) y*L_y/double ( NbrSamplesY),L_y );
	  
            Complex TmpValue = 0.0; //<psi|n_i n_0|psi>

            if ( DensityFlag == false )
            {
                for ( int s=-s_cutoff; s<s_cutoff; ++s )
                {
                    for ( int t=-t_cutoff; t<t_cutoff; ++t )
                    {

                        q[0]=2*M_PI*s*invL_x;
                        q[1]=2*M_PI*t*invL_y;

                        for ( int i = 0; i < NbrPrecalculatedValues; ++i )
                        {
                            j1 = PrecalculatedIndices[4*i];
                            j2 = PrecalculatedIndices[4*i + 1];
                            j3 = PrecalculatedIndices[4*i + 2];
                            j4 = PrecalculatedIndices[4*i + 3];

                            if ( ( ( ( j1-j4 ) - t )  % NbrFluxQuanta )  == 0 )
                            {
                                if ( q.SqrNorm() <=2*34.5 )
                                {
                                    //cout << "q = " << q << endl;
                                    //cout <<"q.SqrNorm() = " << q.SqrNorm() << endl;
                                    summand = Polar ( q*Position - ( j1-j3 ) * 2* M_PI * s *invNbrFluxQuanta ) *exp ( - 0.5 * q.SqrNorm() );
                                    //cout << "summand = " << summand << endl;
				    //cout << "CorrelationTorusPrefactor = " << CorrelationTorusPrefactor << endl;
				    //cout << "PrecalculatedValues_rhorho["<<i<<"] = " << PrecalculatedValues_rhorho[i] << endl;
                                    TmpValue += CorrelationTorusPrefactor * PrecalculatedValues_rhorho[i] * summand ;
                                }
                            }
                        }
                    }
                }
                if ( SubtractDensityFlag == true )
		{
		  TmpValue-=AverageDensity;
		}
                
                //TmpValue-=0.5*NbrParticles*NbrFluxQuanta*CorrelationTorusPrefactor;
            }
            else
            {
                cout << "Density not implemented yet for the continuous torus." << endl;
                for ( int i = 0; i < NbrFluxQuanta; ++i )
                {
                    TmpValue += PrecalculatedValues_rho[i] * SqrNorm ( Coefficients[i] );
                }
                TmpValue *= Normalisation;
            }
            cout << Position[0] << " " << Position[1] << " " << TmpValue.Re << endl;
	    Correlations[x][y]=TmpValue.Re;
        }
    }

    // output in format suitable for Gnuplot:
    for ( int x=0; x< NbrSamplesX; ++x )
    {
        Position[0]=fmod(Position_start[0] + (double)x*L_x/double(NbrSamplesX),L_x);
        for ( int y=0; y< NbrSamplesY; ++y )
        {
            Position[1]=fmod(Position_start[1] + (double)y*L_y/double(NbrSamplesY),L_y);
            File << Position[0] << " " << Position[1] << " " << Correlations[x][y] << endl;
        }
        File << endl;
    }

    File << endl;
    File.close();
    if ( PrecalculatedValues_rhorho!=0 ) delete [] PrecalculatedValues_rhorho;
    delete[] PrecalculatedValues_rho;
    delete[] Coefficients;
    delete[] Coefficients2;
    return 0;
}
