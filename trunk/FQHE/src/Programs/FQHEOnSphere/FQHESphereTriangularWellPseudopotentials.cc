#include "Options/Options.h"

#include "Tools/FQHESpectrum/PseudoPotentials.h"
#include "Tools/FQHESpectrum/AbstractZDensityProfile.h"

#include "Vector/RealVector.h"
 
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

using std::ofstream;
using std::ios;
using std::cout;
using std::endl;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereTriangularWellPseudopotentials" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (i.e. twice the maximum momentum for a single particle)", 8);
  (*SystemGroup) += new SingleIntegerOption  ('l', "landau-level", "index of the Landau level (0 for the lowest Landau level)", 0, true, 0);
  
  (*SystemGroup) += new  SingleDoubleOption ('t', "layer-thickness", "assume finite layer thickness",1.0);
  (*SystemGroup) += new  SingleDoubleOption ('b', "bias", "energy bias of the triangular quantum well, compared to the square well",1.0);
  (*SystemGroup) += new  SingleDoubleOption ('e', "epsilon", "energy correction for subband energy gap",0.0);
  (*SystemGroup) += new  SingleDoubleOption ('d', "delta", "manual entry of the energy difference between the two pseudospin states",-1.0);
  (*SystemGroup) += new  SingleDoubleOption ('m', "magnetic-field", "magnetic field for energy unit conversion between bias of the quantum well and coulomb energy",7.0);
  (*SystemGroup) += new  SingleDoubleOption ('c', "converter", "specify conversion parameter from coulomb energy to units of hbar/(2m l_B^2) manually",-1.0);

  (*SystemGroup) += new MultipleIntegerOption ('\n', "rescale", "rescale separation and thickness for a particular number of particles, filling factor and shift (4 arguments N,p,q,Sigma)",',');
  
  (*SystemGroup) += new  SingleIntegerOption ('\n', "profile-points", "number of points where profile is evaluated", 50);

  (*SystemGroup) += new BooleanOption  ('\n', "tab", "output a table with the results in columns (extension .tab)");
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output file name (default is pseudopotential_coulomb_l_x_2s_y[_v_xxx_w_yyy].dat)");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereTriangularWellPseudopotentials -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
// ***************************************** RESCALE ********************************************
    bool Rescale=false;
    double NewScale=1.0;
    int NbrParticles=1, FillingP=1, FillingQ=1, ShiftSigma=0;
    int TmpLength;
    int *ScaleParams=Manager.GetIntegers("rescale",TmpLength);
    if (TmpLength>0)
      {
	if (TmpLength!=4)
	  {
	    cout << "error: --rescale requires 4 parameters: NbrParticles,FillingP,FillingQ,ShiftSigma"<<endl;
	    exit(1);
	  }
	Rescale=true;
	NbrParticles=ScaleParams[0];
	FillingP=ScaleParams[1];
	FillingQ=ScaleParams[2];
	ShiftSigma=ScaleParams[3];
	NewScale=sqrt((double)NbrParticles/((double)NbrParticles-(double)FillingP/(double)FillingQ*ShiftSigma));
	cout << "Rescaling layer separation and thickness with scale factor: "<<NewScale<<endl;
      }
   
// ********************************* EXTRACT PARAMETERS ******************************************   
    int NbrFlux = Manager.GetInteger("nbr-flux");
    int landauLevel = Manager.GetInteger("landau-level");
    int MaxMomentum = NbrFlux + (landauLevel << 1);
    double width = Manager.GetDouble("layer-thickness")*NewScale;
    double bias = Manager.GetDouble("bias")/NewScale;
    double MagField = Manager.GetDouble("magnetic-field");
    double EnergyUnitsConv = Manager.GetDouble("converter");
 
// ************************** ENERGY UNITS CONVERSION *************************************

    cout << " ------  Triangular well pseudopotentials  ----- " << endl;
    cout << "w=" << width << endl << "V=" << bias << endl;
    
    if( EnergyUnitsConv == -1 )
    {
	const double BohrRadius = 52.917721092e-12;
	const double hbar = 1.054571726e-34;
	const double ElementaryCharge = 1.602176565e-19;
	double Permittivity = 12.9;
	double EffectiveMass = 0.067;
	double EffectiveBohrRadius = Permittivity*BohrRadius/EffectiveMass;
	cout << "Effective Bohr Radius : " << EffectiveBohrRadius << endl;
	double MagLength = sqrt(hbar/ElementaryCharge/MagField);
	cout << "Magnetic Length : " << MagLength << endl;
	EnergyUnitsConv=(2*MagLength)/EffectiveBohrRadius;	
	cout << "EnergyUnitsConv = " << EnergyUnitsConv << endl;
    }

    bias *= EnergyUnitsConv; // transforms the given bias ( in e^2/(\epsilon l_B) units) in the right units for TriangularWell subroutines (\hbar/(2m l_B)=1)
    int nbrPointsInteg = Manager.GetInteger("profile-points");
    double EnergyDifference = Manager.GetDouble("delta");
    double EnergyCorrection = Manager.GetDouble("epsilon");
    double E0,E1;

    if ( bias < 0 ) 
    {
          cout << "error: bias has to be a positive number" << endl;
          exit(1);
    }

    bool ManualEnergyDifference = ( EnergyDifference != -1.0 );
    if ( ManualEnergyDifference )
    {  
	 E0=-EnergyDifference/2.;
	 E1=EnergyDifference/2.;
    }
    else
    {
      if ( bias == 0 )
      {
	E0=M_PI*M_PI/width/width;
	E1=4*E0;
      }
      else
      {
	TriangularWellEigenfunction Well(width,bias);
	E0=Well.GetEnergy(0);
	E1=Well.GetEnergy(1);
      }
      E0/=EnergyUnitsConv;// convert energies given by TriangularWell subroutine in e^2/(\epsilon l_B) units
      E1/=EnergyUnitsConv;
      E1+=EnergyCorrection;
      //EnergyDifference=E1-E0;
    }

// ******************************* COMPUTE PSEUDOPOTENTIALS *************************************
     int pseudospins[4];
     pseudospins[0] = 0; pseudospins[1] = 0; pseudospins[2] = 0; pseudospins[3] = 0;
     double* PseudopotentialsUpUpUpUp = EvaluateTriangularWellPseudopotentials(NbrFlux, landauLevel, width, bias, pseudospins, nbrPointsInteg, true);
     pseudospins[0] = 0; pseudospins[1] = 1; pseudospins[2] = 0; pseudospins[3] = 0;
     double* PseudopotentialsUpDownUpUp = EvaluateTriangularWellPseudopotentials(NbrFlux, landauLevel, width, bias, pseudospins, nbrPointsInteg, true);
     pseudospins[0] = 1; pseudospins[1] = 1; pseudospins[2] = 0; pseudospins[3] = 0;
     double* PseudopotentialsDownDownUpUp = EvaluateTriangularWellPseudopotentials(NbrFlux, landauLevel, width, bias, pseudospins, nbrPointsInteg, true);
     pseudospins[0] = 0; pseudospins[1] = 1; pseudospins[2] = 0; pseudospins[3] = 1;
     double* PseudopotentialsUpDownUpDown = EvaluateTriangularWellPseudopotentials(NbrFlux, landauLevel, width, bias, pseudospins, nbrPointsInteg, true);
     pseudospins[0] = 0; pseudospins[1] = 1; pseudospins[2] = 1; pseudospins[3] = 1;
     double* PseudopotentialsUpDownDownDown = EvaluateTriangularWellPseudopotentials(NbrFlux, landauLevel, width, bias, pseudospins, nbrPointsInteg, true);
     pseudospins[0] = 1; pseudospins[1] = 1; pseudospins[2] = 1; pseudospins[3] = 1;
     double* PseudopotentialsDownDownDownDown = EvaluateTriangularWellPseudopotentials(NbrFlux, landauLevel, width, bias, pseudospins, nbrPointsInteg, true);

// ********************************** WRITE OUTPUT FILE ****************************************
     char* OutputFile;
     ((SingleDoubleOption*)Manager["layer-thickness"])->SetStringFormat("%g");
     ((SingleDoubleOption*)Manager["bias"])->SetStringFormat("%g");
     ((SingleDoubleOption*)Manager["delta"])->SetStringFormat("%g");
     if (Manager.GetString("output") == 0l) 
        {
         if ( ManualEnergyDifference )
                OutputFile = Manager.GetFormattedString("pseudopotential_triangular_width_%layer-thickness%_bias_%bias%_delta_%delta%_l_%landau-level%_2s_%nbr-flux%.dat");
         else
                OutputFile = Manager.GetFormattedString("pseudopotential_triangular_width_%layer-thickness%_bias_%bias%_l_%landau-level%_2s_%nbr-flux%.dat");
        }
     else
       {
         OutputFile = new char [strlen(Manager.GetString("output")) + 1];
         strcpy (OutputFile, Manager.GetString("output"));
       }

      ofstream File;
      File.open(OutputFile, ios::binary | ios::out);
      File.precision(14);
      File << "PseudopotentialsUpUpUpUp=";        
      for (int i = 0; i <= MaxMomentum; ++i)
        File << " " << PseudopotentialsUpUpUpUp[i];
      File << endl;    

      File << "PseudopotentialsUpDownUpUp=";        
      for (int i = 0; i <= MaxMomentum; ++i)
        File << " " << PseudopotentialsUpDownUpUp[i];
      File << endl;    
 
      File << "PseudopotentialsDownDownUpUp=";        
      for (int i = 0; i <= MaxMomentum; ++i)
        File << " " << PseudopotentialsDownDownUpUp[i];
      File << endl;    
 
      File << "PseudopotentialsUpUpUpDown=";        
      for (int i = 0; i <= MaxMomentum; ++i)
        File << " " << PseudopotentialsUpDownUpUp[i];
      File << endl;    
 
      File << "PseudopotentialsUpDownUpDown=";        
      for (int i = 0; i <= MaxMomentum; ++i)
        File << " " << PseudopotentialsUpDownUpDown[i];
      File << endl;    
 
      File << "PseudopotentialsDownUpUpDown=";        
      for (int i = 0; i <= MaxMomentum; ++i)
        File << " " << PseudopotentialsDownDownUpUp[i];
      File << endl;    
 
      File << "PseudopotentialsDownDownUpDown=";        
      for (int i = 0; i <= MaxMomentum; ++i)
        File << " " << PseudopotentialsUpDownDownDown[i];
      File << endl;    
 
      File << "PseudopotentialsUpUpDownDown=";        
      for (int i = 0; i <= MaxMomentum; ++i)
        File << " " << PseudopotentialsDownDownUpUp[i];
      File << endl;    
 
      File << "PseudopotentialsUpDownDownDown=";        
      for (int i = 0; i <= MaxMomentum; ++i)
        File << " " << PseudopotentialsUpDownDownDown[i];
      File << endl;    
 
      File << "PseudopotentialsDownDownDownDown=";        
      for (int i = 0; i <= MaxMomentum; ++i)
        File << " " << PseudopotentialsDownDownDownDown[i];
      File << endl;    
 
      if ( EnergyDifference != 0 )
      {
        File << "OneBodyPotentialUpUp=";        
        for (int i = 0; i <= MaxMomentum; ++i)
          File << " " << E0;
        File << endl;    
      
        File << "OneBodyPotentialDownDown=";        
        for (int i = 0; i <= MaxMomentum; ++i)
          File << " " << E1;
        File << endl;    
        File.close();
      }

    delete[] OutputFile;

    delete[] PseudopotentialsUpUpUpUp;
    delete[] PseudopotentialsUpDownUpUp;
    delete[] PseudopotentialsDownDownUpUp;
    delete[] PseudopotentialsUpDownUpDown;
    delete[] PseudopotentialsUpDownDownDown;
    delete[] PseudopotentialsDownDownDownDown;

   return 0;
}
