#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereLong.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "GeneralTools/ConfigurationParser.h"

#include "Operator/ParticleOnSphereDensityOperator.h"
#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "FunctionBasis/ParticleOnCylinderFunctionBasis.h"
#include "Hamiltonian/ParticleOnCylinderStructureFactor.h"
#include "Hamiltonian/ParticleOnCylinderPairAmplitude.h"
#include "Hamiltonian/ParticleOnCylinderDensityDensity.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term

//Complex EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int MaxMomentum, double kappa, double x, double y);

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHECylinderCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('e', "eigenstate", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "ky-max", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "total-y", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "landau-level", "index of the Landau level (0 being the LLL)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "hopping-cutoff", "cutoff on the hopping processes (=|m1-m3|)", -1);
  (*SystemGroup) += new SingleDoubleOption  ('r', "ratio", "aspect ratio of the cylinder", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of points to evaluate", 50);
  (*SystemGroup) += new BooleanOption  ('\n', "rho-rho","evaluate rho-rho correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "structure-factor","evaluate the guiding center structure factor", false);

  (*SystemGroup) += new BooleanOption  ('\n', "pair-amplitude","evaluate the pair amplitude", false);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pair-center", "twice the momentum of center of the pair", -1);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "pair-anisotropy", "anisotropy parameter for the pair", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pair-minmomentum", "cutoff for the minimal momentum of a pair", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pair-maxmomentum", "cutoff for the maximal momentum of a pair", -1);

  (*SystemGroup) += new SingleIntegerOption  ('\n', "x-points", "number of points along the cylinder", 100);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsCorrelation -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int KyMax = Manager.GetInteger("ky-max");
  int TotalKy = Manager.GetInteger("total-y");
  int LandauLevel = Manager.GetInteger("landau-level");
  int NbrPoints = Manager.GetInteger("nbr-points");
  int XPoints = Manager.GetInteger("x-points");
  double Ratio = Manager.GetDouble("ratio");
  bool EvaluateS0Q = Manager.GetBoolean("structure-factor");
  bool EvaluateRhoRho = Manager.GetBoolean("rho-rho");
  bool EvaluatePairAmplitude = Manager.GetBoolean("pair-amplitude");
  int PairCenter = Manager.GetInteger("pair-center");
  double PairAnisotropy = Manager.GetDouble("pair-anisotropy");
  int MinPairMomentumCutoff = Manager.GetInteger("pair-minmomentum");
  int MaxPairMomentumCutoff = Manager.GetInteger("pair-maxmomentum");

  int HoppingCutoff = Manager.GetInteger("hopping-cutoff");

  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  bool Statistics = true;
  if (Manager.GetString("eigenstate") == 0)
    {
      cout << "FQHECylinderFermionsCorrelation requires a state" << endl;
      return -1;
    }

  //if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("eigenstate"),
//						  NbrParticles, KyMax, TotalKy, Statistics) == false)
//    {
//      cout << "error while retrieving system parameters from file name " << Manager.GetString("eigenstate") << endl;
//      return -1;
//    }

  if (IsFile(Manager.GetString("eigenstate")) == false)
    {
      cout << "can't find vector file " << Manager.GetString("eigenstate") << endl;
      return -1;      
    }

  ParticleOnSphere* Space = 0;
#ifdef __64_BITS__
	  if (TotalKy <= 62)
#else
	  if (TotalKy <= 30)
#endif
  	    Space = new FermionOnSphere(NbrParticles, TotalKy, KyMax);

	  else
#ifdef __128_BIT_LONGLONG__
	    if (TotalKy <= 126)
#else
	      if (TotalKy <= 62)
#endif
 	        Space = new FermionOnSphereLong(NbrParticles, TotalKy, KyMax);
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, TotalKy, KyMax);


  cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;

  ofstream File;
  File.precision(14);
  if (Manager.GetString("output-file") != 0)
     File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
   {
    cout << "Enter output file! " << endl;
    exit(1);
   }
  
  ComplexVector State;

  if (State.ReadVector (Manager.GetString("eigenstate")) == false)
   {
     cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
     return -1;      
   }

  if (EvaluatePairAmplitude == true)
   {
      long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;

      if (PairCenter < 0)
        PairCenter = KyMax + 1;

      cout << "Evaluating pair amplitudes. Pair centered at " << PairCenter << " , pair anisotropy = " << PairAnisotropy << endl;

      cout << "Pair amplitude is defined as P^+P (M, a, p) = sum_r,rp exp(-Xr*Xr-Xrp*Xrp) * Hermite(Xr * math.sqrt(2.0)) * Hermit(Xrp * math.sqrt(2.0)) /(2^M M!) c_{p+r}^+ c_{p-r}^+ c_{p+rp} c_{p-rp} " << endl;

      double Length = sqrt(2.0 * M_PI * Ratio * (KyMax + 1));
      double kappaL = 2.0 * M_PI/Length;
      double Height = Length/Ratio;
      double kappaH = 2.0 * M_PI/Height;
      double MinDimension = Length;
      if (Height < Length)
         MinDimension = Height;

      if (MaxPairMomentumCutoff < 0)
         MaxPairMomentumCutoff = pow(MinDimension, 2.0)/8;

      for (int TmpMomentum = MinPairMomentumCutoff; TmpMomentum <= MaxPairMomentumCutoff; TmpMomentum++)
        {	
           AbstractHamiltonian* Hamiltonian = new ParticleOnCylinderPairAmplitude (Space, NbrParticles, KyMax, Ratio, PairCenter, PairAnisotropy, TmpMomentum, Architecture.GetArchitecture(), Memory);

           ComplexVector TmpState(Space->GetHilbertSpaceDimension(), true);
           VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
           Operation.ApplyOperation(Architecture.GetArchitecture());
           Complex TmpPairAmp = State * TmpState;
           cout << "M= " << TmpMomentum << " "<<TmpPairAmp.Re << " "<<TmpPairAmp.Im << endl;
           File << TmpMomentum << " "<<TmpPairAmp.Re << " " << TmpPairAmp.Im << endl;
        }

      return 0;
   }

  if (EvaluateS0Q == true)
    {
      cout << "Evaluating structure factor...."<<endl;
      double Length = sqrt(2.0 * M_PI * Ratio * (KyMax + 1));
      double kappaL = 2.0 * M_PI/Length;
      double Height = Length/Ratio;
      double kappaH = 2.0 * M_PI/Height;
      long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;

      int counter = 1;
      for (int s = 0; s <= KyMax; s++)
        for (int t = 0; t <= KyMax; t++)
          {
            cout << "Doing sector " << counter << " out of " << (KyMax+1)*(KyMax+1) << endl;
            counter++;

            int mint = t;
            double QyValue = kappaL * t;
            double QyValueRescaled = kappaL * (t - (KyMax + 1));
            if (QyValue*QyValue > QyValueRescaled*QyValueRescaled)
              {
               mint = t - (KyMax+1);
               QyValue = QyValueRescaled;
              }
   

            int mins = s;
            double QxValue = kappaH * s;
            double QxValueRescaled = kappaH * (s - (KyMax + 1));
            if (QxValue*QxValue > QxValueRescaled*QxValueRescaled)
              {
               mins = s - (KyMax+1);
               QxValue = QxValueRescaled;
              }
   
            double Q = sqrt(QxValue * QxValue + QyValue * QyValue);
       
            AbstractHamiltonian* Hamiltonian = new ParticleOnCylinderStructureFactor (Space, NbrParticles, KyMax, Ratio, QxValue, mint, Architecture.GetArchitecture(), Memory);

           ComplexVector TmpState(Space->GetHilbertSpaceDimension(), true);
           VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
           Operation.ApplyOperation(Architecture.GetArchitecture());
           Complex S0Q = State*TmpState;

           //Add the term \sum_i <exp(iq.R_i)exp(-iq.R_i)>=N_e 
           S0Q += (double)NbrParticles/(double)(KyMax+1);

           if (Q == 0.0)
           {  
            //Subtract the term 1/N_phi \sum_i <exp(iq.R_i)> \sum_j <exp(-iq.R_j)>
            S0Q -= (double)(NbrParticles*NbrParticles)/(double)(KyMax+1);
           }

           delete Hamiltonian;

           cout << "Qx= "<<QxValue<<" , Qy= "<<QyValue<<" "<<S0Q<<endl; 
           File << Q << " " << QxValue << " " << QyValue <<" "<<S0Q.Re<<" "<<S0Q.Im<<endl; 
         }  
	
      return 0;
    }

  cout << "NbrParticles "<<NbrParticles<<" "<<KyMax<<" "<<LandauLevel<<" "<<Ratio<<endl;
  cout << "NbrPoints "<<NbrPoints<<endl;
  
  ParticleOnCylinderFunctionBasis Basis (KyMax, LandauLevel, Ratio);

  double H = sqrt(2.0 * M_PI * (KyMax + 1.0))/sqrt(Ratio);
  double Length = sqrt(2.0 * M_PI * Ratio * (KyMax + 1));
  cout << "Cylinder L= " << Length<<" H= "<<H<<endl;

  Complex* Occupations= new Complex [KyMax+1];
      
  cout<<"Orbital occupations: "<<endl;
  Complex CheckOccupations(0,0);
  for (int i = 0; i <= KyMax; ++i)
   {
//     AbstractOperator* Operator = new ParticleOnSphereDensityOperator(Space, i);
//     ComplexVector TmpState(Space->GetHilbertSpaceDimension(), true);
//     VectorOperatorMultiplyOperation Operation (Operator, &State, &TmpState);
//     Operation.ApplyOperation(Architecture.GetArchitecture());
//     Occupations[i] = State*TmpState;
//     delete Operator;
//     CheckOccupations += Occupations[i];
//     cout << i << " " << Occupations[i].Re << " " << Occupations[i].Im << endl;

     ParticleOnSphereDensityOperator Operator (Space, i);
     Occupations[i] = Operator.MatrixElement(State, State);
     CheckOccupations += Occupations[i];
     cout << i << " " << Occupations[i].Re << " " << Occupations[i].Im << endl;
   }
  cout<<endl;
  cout<<"Sum of all occupations: "<<CheckOccupations<<endl;

  double Offset = 10.0;
  Complex DensityIntegral (0.0,0.0);
  double XPosition = -0.5*(H + Offset);
  double Step = (H + Offset)/1000.0;
  while (XPosition <= 0.5*(H + Offset))
    {
      Complex Density (0.0, 0.0);        
      for (int i = 0; i <= KyMax; ++i)
        {
          Density += Conj(Basis.GetFunctionValue(XPosition, 0.0, (double)i-0.5*KyMax)) * Basis.GetFunctionValue(XPosition, 0.0, (double)i-0.5*KyMax) * Occupations[i];
	  //cout<<"i= "<<i<<" "<<DensityMatEl[i]<<" ";
        }

      DensityIntegral += Density * Step * Length;
      if (EvaluateRhoRho == false)
         File << XPosition <<" "<< Density.Re << " " << Density.Im << endl;
      XPosition += Step;
    }

  cout << "Integrated density along the cylinder: "<<DensityIntegral <<endl;




  if (EvaluateRhoRho)
   {

/*
        long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
        double XInc = 0.8*H/(double)NbrPoints;
        double X0 = 0.0;
        double X; 
        Complex RhoRho; 

        int m4;
        int NbrIndices = 0;
        for (int m1 = 0; m1 <= KyMax; ++m1)
  	  for (int m2 = 0; m2 < m1; ++m2)
	    for (int m3 = 0; m3 <= KyMax; ++m3)
 	      {
	        m4 = m1 + m2 - m3;
	        if ((m4 >= 0) && (m4 <= KyMax))
                  if (m3 > m4)
  		    {
                       NbrIndices++;
                    }
               }
         cout << "nbr j1,j2,j3,j4 = " << NbrIndices << endl; 

         int* M1Values = new int [NbrIndices];
         int* M2Values = new int [NbrIndices];
         int* M3Values = new int [NbrIndices];
         int* M4Values = new int [NbrIndices];
         Complex* CoeffValues = new Complex [NbrIndices];

        NbrIndices = 0;
        for (int m1 = 0; m1 <= KyMax; ++m1)
  	  for (int m2 = 0; m2 < m1; ++m2)
	    for (int m3 = 0; m3 <= KyMax; ++m3)
	      {
	        m4 = m1 + m2 - m3;
	        if ((m4 >= 0) && (m4 <= KyMax))
                  if (m3 > m4)
  		    {
                       M1Values[NbrIndices] = m1;
                       M2Values[NbrIndices] = m2;
                       M3Values[NbrIndices] = m3;
                       M4Values[NbrIndices] = m4;
                       ParticleOnSphereDensityDensityOperator Operator (Space, m1, m2, m3, m4);
                       CoeffValues[NbrIndices] = Operator.MatrixElement(State, State);
                       NbrIndices++;
                    }
               }
         cout << "Finished precalculation." << endl;       

        Complex RhoRhoIntegrated(0,0);
        double kappa = 2.0 * M_PI/Length;
        for (int k = 0; k < NbrPoints; ++k)
	  {
            X = X0 + k * XInc;

            RhoRho.Re = 0.0;
            RhoRho.Im = 0.0;
            for (int i = 0; i < NbrIndices; i++)
              {
                RhoRho += CoeffValues[i] * (EvaluateInteractionCoefficient(M1Values[i], M2Values[i], M3Values[i], M4Values[i], KyMax, kappa, X, 0.0)
                                            -EvaluateInteractionCoefficient(M2Values[i], M1Values[i], M3Values[i], M4Values[i], KyMax, kappa, X, 0.0)
                                            -EvaluateInteractionCoefficient(M1Values[i], M2Values[i], M4Values[i], M3Values[i], KyMax, kappa, X, 0.0)
                                            +EvaluateInteractionCoefficient(M2Values[i], M1Values[i], M4Values[i], M3Values[i], KyMax, kappa, X, 0.0)  
                                            );
              }

           //RhoRho /= (NbrParticles * (NbrParticles - 1)); 
           RhoRhoIntegrated += Length * RhoRho * XInc;

           File << X << " " << RhoRho.Re << " " << RhoRho.Im << endl;
           cout << X << " " << RhoRho.Re << " " << RhoRho.Im << endl;
	  }
           
      cout<<"Integral of rho rho = " << RhoRhoIntegrated << endl;

      delete[] M1Values;
      delete[] M2Values;
      delete[] M3Values;
      delete[] M4Values;
      delete[] CoeffValues;
*/
        long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
        double XInc = 0.99*H/(double)NbrPoints;
        double X0 = -0.49 * H;
        double X; 
        Complex RhoRho; 

        Complex Density0 (0.0, 0.0);        
        for (int i = 0; i <= KyMax; ++i)
          {
            Density0 += Conj(Basis.GetFunctionValue(X0, 0.0, (double)i-0.5*KyMax)) * Basis.GetFunctionValue(X0, 0.0, (double)i-0.5*KyMax) * Occupations[i];
	  //cout<<"i= "<<i<<" "<<DensityMatEl[i]<<" ";
          }

        double kappa = 2.0 * M_PI/Length;
        for (int k = 0; k < NbrPoints; ++k)
	  {
            X = X0 + k * XInc;

            AbstractHamiltonian* Hamiltonian = new ParticleOnCylinderDensityDensity (Space, NbrParticles, KyMax, Ratio, LandauLevel, X0, 0.0, X, 0.0, HoppingCutoff, Architecture.GetArchitecture(), Memory);

           ComplexVector TmpState(Space->GetHilbertSpaceDimension(), true);
           VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
           Operation.ApplyOperation(Architecture.GetArchitecture());
           RhoRho = State*TmpState;

           Complex Density (0.0, 0.0);        
           for (int i = 0; i <= KyMax; ++i)
             {
                Density += Conj(Basis.GetFunctionValue(X, 0.0, (double)i-0.5*KyMax)) * Basis.GetFunctionValue(X, 0.0, (double)i-0.5*KyMax) * Occupations[i];
	        //cout<<"i= "<<i<<" "<<DensityMatEl[i]<<" ";
             }
           
           //RhoRho *= 1.0 /(NbrParticles * Density0);
           RhoRho /= (Density0 * Density);

           File << X << " " << RhoRho.Re << " " << RhoRho.Im << endl;
           cout << X << " " << RhoRho.Re << " " << RhoRho.Im << endl;
         }

  }

  if (HoppingCutoff >= 0)
    cout << "Hopping was truncated at |m1-m3|<= " << HoppingCutoff << endl;

  delete[] Occupations;

  File.close();
 
  return 0;
}

/*
// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term

Complex EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int MaxMomentum, double kappa, double x, double y)
{
  double Xm1 = kappa * m1;
  double Xm2 = kappa * m2;
  double Xm3 = kappa * m3;
  double Xm4 = kappa * m4;  

  Complex Phase;
  Phase.Re = cos(y * (Xm1 - Xm4));
  Phase.Im = sin(y * (Xm1 - Xm4));

  Complex Res;
  Res.Re = exp(-0.5*pow(Xm1-Xm4,2.0));
  Res.Im = 0.0;
  Res *= exp(-0.5*pow(x - (Xm1-Xm3),2.0));

  return (Res * Phase);
}
*/
