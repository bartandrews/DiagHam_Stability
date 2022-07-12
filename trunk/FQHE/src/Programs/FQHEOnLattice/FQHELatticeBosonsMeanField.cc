#include "Tools/FQHESpectrum/LatticePhases.h"
#include "Tools/FQHEWaveFunction/GrossPitaevskiiOnLatticeState.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeBosonsMeanField" , "0.01");  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OptimizationGroup = new OptionGroup ("optimization options");
  OptionGroup* ExportGroup = new OptionGroup ("export options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  LatticePhases::AddOptionGroup(&Manager);
  Manager += OptimizationGroup;
  Manager += ExportGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleDoubleOption  ('m', "mu", "chemical potential", 1.0);
  (*SystemGroup) += new SingleStringOption  ('e', "interaction-file", "use definition of two-body interactions from a file");
  (*SystemGroup) += new SingleStringOption  ('E', "interaction-name", "descriptor of external interaction (if in use)","ext");
  (*SystemGroup) += new SingleStringOption  ('f', "potential-file", "use definition of one-body interactions from a file");
  (*SystemGroup) += new SingleStringOption  ('F', "potential-name", "descriptor of external single particle potential (if in use)","");

  (*OptimizationGroup) += new BooleanOption('\n', "gradient", "Use optimization based on gradients");
  (*OptimizationGroup) += new BooleanOption('\n', "simplex", "Use optimization based on a simplex method");
  (*OptimizationGroup) += new SingleDoubleOption('\n', "tolerance", "tolerance for variational parameters in condensate",1e-6);
  (*OptimizationGroup) += new SingleIntegerOption('i', "nbr-iter", "number of iterations for optimization",10000);
  (*OptimizationGroup) += new SingleStringOption('\n', "parameters", "file with initial parameters");
  (*OptimizationGroup) += new BooleanOption('u', "uniform", "start from a uniform vector field configuration (first attempt only)");
  (*OptimizationGroup) += new BooleanOption('r', "random-amplitude", "randomize phases and amplitudes");

  (*OptimizationGroup) += new SingleIntegerOption('a', "nbr-attempts", "number of separate attempts to optimize a configuration",1);
  (*OptimizationGroup) += new SingleIntegerOption('s', "nbr-save", "maximum number of (distinct) configurations to be saved",10);
  (*OptimizationGroup) += new SingleDoubleOption('\n', "overlap-same", "overlap above which two conf's are considered identical", 0.99);

  (*ExportGroup) += new BooleanOption  ('\n', "export", "export a given parameter file (input: parameters)");
  (*ExportGroup) += new BooleanOption  ('\n', "no-rectify", "take bare one-particle wavefunctions for a given trapping potential");
  (*ExportGroup) += new SingleStringOption  ('\n', "basis-path", "path to basis file","./");
  (*ExportGroup) += new SingleStringOption  ('\n', "basis-pattern", "pattern of basis file names (variable site number given as XX)");
  
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file", NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  double ChemicalPotential = Manager.GetDouble("mu");

  int NbrAttempts = Manager.GetInteger("nbr-attempts");
  int NbrToSave = Manager.GetInteger("nbr-save");
  double IdentityThreshold = Manager.GetDouble("overlap-same");
  bool Uniform = Manager.GetBoolean("uniform");
  
  if ((Manager.GetBoolean("export")==false)&&(Manager.GetString("interaction-file")==NULL)&&(Manager.GetString("potential-file")==NULL))
    {
      cout << "An external definition of the hopping or interaction matrix elements is required. Use option -e and/or -f"<<endl;
      exit(1);
    }

  // get the lattice geometry
  LatticePhases *Lattice = new LatticePhases();
  char* LatticeName = Lattice->GeometryString();
  int NbrSites = Lattice->GetNbrSites();


  if (Manager.GetBoolean("export"))
    {
      if (Manager.GetString("parameters")==NULL)
	{
	  cout << "Please indicate a wavefunction to be expanded using option --parameters"<<endl;
	  exit(1);
	}
      ComplexVector InputParameters;
      if (InputParameters.ReadVector(Manager.GetString("parameters"))==false)
	{
	  cout << "Could not read vector of initial parameters" <<Manager.GetString("parameters")<<endl;
	  exit(1);
	}
      int NbrCells = Lattice->GetNbrCells();
      int NbrParameters = InputParameters.GetVectorDimension();
      if (NbrParameters%NbrCells!=0)
	{
	  cout << "Lattice geometry and number of parameters do not match up"<<endl;
	  exit(1);
	}
      int NbrStatesPerCell = NbrParameters/NbrCells;
      cout << "System has "<<NbrStatesPerCell<<" basis states per unit cell of the lattice"<<endl;
      char *BasisPath = Manager.GetString("basis-path");
      if (BasisPath==0)
	{
	  BasisPath = new char[3];
	  sprintf(BasisPath,"./");
	}
      else if (BasisPath[strlen(BasisPath)-1]!='/')
	{
	  BasisPath = new char[strlen(BasisPath)+2];
	  sprintf(BasisPath,"%s/",Manager.GetString("basis-path"));
	}
      else
	{
	  BasisPath = new char[strlen(BasisPath)+1];
	  sprintf(BasisPath,"%s",Manager.GetString("basis-path"));
	}
      char *BasisPattern = Manager.GetString("basis-pattern");
      if (BasisPattern==0)
	{
	  cout << "A pattern for the names of the basis file names is required"<<endl;
	  exit(1);
	}
      char *SiteLabel = strstr(BasisPattern, "XX");
      if (SiteLabel ==0)
	{
	  cout << "Error: Pattern does not contain place-holder for site number 'XX'"<<endl;
	  exit(1);
	}
      char *BasisInitial = new char[SiteLabel-BasisPattern+1];
      strncpy(BasisInitial, BasisPattern, SiteLabel-BasisPattern);
      BasisInitial[SiteLabel-BasisPattern]='\0'; // not required!
      char *BasisFinal = new char[strlen(BasisPattern)-strlen(BasisInitial)];
      strcpy(BasisFinal, SiteLabel+2);

      cout << "Reading basis from files "<<BasisPath<<BasisInitial<<" XX "<<BasisFinal<<endl;
      
      int NbrBasisStates=NbrParameters;
      ComplexVector *BasisStates = new ComplexVector[NbrBasisStates];
      char *TmpVectorName = new char[strlen(BasisPath)+strlen(BasisInitial)+strlen(BasisFinal)+6];
      for (int i=0; i<NbrBasisStates; ++i)
	{
	  sprintf(TmpVectorName,"%s%s%d%s", BasisPath, BasisInitial, i, BasisFinal);
	  if (BasisStates[i].ReadVector(TmpVectorName)==false)
	    {
	      cout << "Error: could not read vector "<<TmpVectorName<<endl;
	      exit(1);
	    }
	  if (BasisStates[i].GetVectorDimension() != NbrSites)
	    {
	      cout << "Number of sites on lattice does not correspond to dimension of wavefunction "<<TmpVectorName<<endl;
	      exit(1);
	    }
	  // choose gauge such that largest entry is real (as requested from MatrixElement in DiceLatticeModel.pl
	  int maxI=0;
	  double maxN=0.0;
	  for (int n=0; n<NbrSites; ++n)
	    if (Norm(BasisStates[i][n])>maxN)
	      {
		maxI=n;
		maxN=Norm(BasisStates[i][n]);
	      }
	  Complex Phase=Polar(1.0,-Arg(BasisStates[i][maxI]));
	  BasisStates[i]*=Phase;
	  // rectify absolute values of entries to the case of vanishing trapping potential
	  if (!Manager.GetBoolean("no-rectify"))
	    {
	      for (int n=0; n<NbrSites; ++n)
		{
		  if (fabs(Norm(BasisStates[i][n])-1.0/sqrt(2.0))<0.05)
		    {
		      BasisStates[i][n]*=1.0/sqrt(2.0)/Norm(BasisStates[i][n]);
		    }
		  else if (fabs(Norm(BasisStates[i][n])-1.0/sqrt(12.0))<0.05)
		    {
		      BasisStates[i][n]*=1.0/sqrt(12.0)/Norm(BasisStates[i][n]);
		    }
		  else if (Norm(BasisStates[i][n])>1e-6)
		    {
		      cout << "Unknown amplitude for single-particle wavefunction!"<<endl;
		      exit(1);
		    }
		}
	    }
	}

      ComplexVector ResultingState;
      ResultingState.Copy(BasisStates[0], InputParameters[0]);
      for (int i=1; i<NbrBasisStates; ++i)
	ResultingState.AddLinearCombination(InputParameters[i], BasisStates[i]);

      char* OutputName = ReplaceExtensionToFileName(Manager.GetString("parameters"),"par","full.vec");
      
      ResultingState.WriteVector(OutputName);

      char* FieldName = ReplaceExtensionToFileName(Manager.GetString("parameters"),"par","full.dat");
      char* CurrentFieldName = ReplaceExtensionToFileName(Manager.GetString("parameters"),"par","j.dat");
      ofstream VectorField;
      VectorField.open(FieldName,ios::out);
      RealVector SitePosition, SitePosition2, RelativePosition;
      int CellPosition[2];
      int Sub;
      for (int s=0; s<NbrSites; ++s)
	{
	  Lattice->GetSiteCoordinates(s, CellPosition, Sub);
	  SitePosition = Lattice->GetSitePosition(CellPosition,Sub);
	  VectorField << SitePosition[0] << "\t" << SitePosition[1]
		      << "\t" << ResultingState[NbrSites-s-1].Re
		      << "\t" << ResultingState[NbrSites-s-1].Im << endl;
	}
      VectorField.close();
      VectorField.open(CurrentFieldName,ios::out);
      cout << "Successfully wrote effective model wavefunction to "<<OutputName<<" / dat"<<endl;

      double FluxDensity=0.5;
      cout << "Assuming flux density of n_phi="<<FluxDensity<<endl;
      // calculate current operator
      int NbrNeighbors;
      int *Neighbors;
      double *Phases;
      double *Amplitudes;
      int **PeriodicTranslations;
      for (int s=0; s<NbrSites; ++s)
	{
	  Lattice->GetSiteCoordinates(s, CellPosition, Sub);
	  SitePosition = Lattice->GetSitePosition(CellPosition,Sub);
	  Lattice->GetNeighbors(s, NbrNeighbors, Neighbors, Phases, PeriodicTranslations, Amplitudes);
	  for (int n=0; n<NbrNeighbors; ++n)
	    if (Neighbors[n]>s)
	      {
		double Current = -2.0*Imag(Conj(ResultingState[NbrSites-s-1])*ResultingState[NbrSites-Neighbors[n]-1]*Polar(1.0, -2.0*M_PI*FluxDensity*Phases[n]));
		if (Amplitudes!=NULL)
		  Current *= Amplitudes[n];
		//cout << "J_"<<s<<","<<Neighbors[n]<<" = "<<Current<<endl;
		Lattice->GetSiteCoordinates(Neighbors[n], CellPosition, Sub);
		SitePosition2 = Lattice->GetSitePosition(CellPosition,Sub);
		
		for (int d=0; d<2; ++d)
		  SitePosition2.AddLinearCombination((double)(-PeriodicTranslations[n][d]*Lattice->GetLatticeLength(d)),Lattice->GetLatticeVector(d));
		//cout << "Connection "<<endl<<SitePosition<<"to"<<endl<<SitePosition2<< "(shift "<<PeriodicTranslations[n][0]<<","<< PeriodicTranslations[n][1]<<")"<<endl;
		RelativePosition.Copy(SitePosition2);
		RelativePosition.AddLinearCombination(-1.0,SitePosition);
		RelativePosition/=RelativePosition.Norm();
		SitePosition2.AddLinearCombination(1.0,SitePosition);
		SitePosition2*=0.5;

		RelativePosition*=Current;
		
		VectorField << SitePosition2[0] << "\t" << SitePosition2[1]
			    << "\t" << RelativePosition[0] << "\t" << RelativePosition[1]<<endl;

		cout << s << " " << Neighbors[n] << " " << Arg(Conj(ResultingState[NbrSites-s-1])*ResultingState[NbrSites-Neighbors[n]-1]*Polar(1.0, -2.0*M_PI*FluxDensity*Phases[n])) << " "<< -PeriodicTranslations[n][0] << " " << -PeriodicTranslations[n][1] << endl;

		
	      }
	}
      VectorField.close();
      
      delete [] BasisInitial;
      delete [] BasisFinal;
      delete [] BasisPath;
      delete [] TmpVectorName;
      delete [] OutputName;
      delete [] FieldName;
      delete [] CurrentFieldName;

      exit(0);
    }
  
  char* OutputName;
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      char* TmpChar = new char[1024];
      OutputName = new char [1024];
      sprintf(OutputName,"bosons_lattice_mf_%s",LatticeName);
      if (Manager.GetString("interaction-file")!=NULL)
	{
	  if (Manager.GetString("interaction-name")!=NULL)
	    sprintf(TmpChar,"%s_%s",OutputName,Manager.GetString("interaction-name"));
	  else
	    sprintf(TmpChar,"%s_%s",OutputName,Manager.GetString("interaction-file"));
	  strcpy(OutputName,TmpChar);
	}
      if (Manager.GetString("potential-file")!=NULL)
	{
	  if (Manager.GetString("potential-name")!=NULL)
	    sprintf(TmpChar,"%s_%s",OutputName,Manager.GetString("potential-name"));
	  else
	    sprintf(TmpChar,"%s_%s",OutputName,Manager.GetString("potential-file"));
	  strcpy(OutputName,TmpChar);
	}
      sprintf(TmpChar,"%s_mu_%g.dat",OutputName,ChemicalPotential);
      strcpy(OutputName,TmpChar);
      delete [] TmpChar;
    }

  char* BaseName = RemoveExtensionFromFileName(OutputName,".dat");

  ComplexVector *InitialParameters = NULL;

  if ((Manager.GetString("parameters")!=NULL)&&(NbrAttempts==1))
    {
      InitialParameters = new ComplexVector;
      if (InitialParameters->ReadVector(Manager.GetString("parameters"))==false)
	{
	  cout << "Could not read vector of initial parameters" <<Manager.GetString("parameters")<<endl;
	  exit(1);
	}
    }

  ofstream File;
  ifstream TestFile;
  TestFile.open(OutputName, ios::in);
  if (TestFile.is_open())
    {
      TestFile.close();
      File.open(OutputName, ios::app);
    }
  else
    {
      File.open(OutputName, ios::out );
      File << "#E_site\tDensity\tE_int\tParameters"<<endl;
    }
    
  ComplexVector *OptimalWaveFunctions = new ComplexVector[NbrToSave];
  int NbrFound=0;
  double *LowestEnergies = new double[NbrToSave];
  timeval RandomTime;
  gettimeofday (&(RandomTime), 0);
  NumRecRandomGenerator *RandomNumberGenerator = new NumRecRandomGenerator(RandomTime.tv_sec);

  for (int i=0; i<NbrAttempts; ++i)
    {
      GrossPitaevskiiOnLatticeState MeanFieldState(NbrSites, Manager.GetString("potential-file"), Manager.GetString("interaction-file"), Lattice, InitialParameters, RandomNumberGenerator);
      MeanFieldState.SetChemicalPotential(ChemicalPotential);
      if (InitialParameters==NULL)
	{
	  if ((Uniform)&&(i==0))
	    MeanFieldState.SetToUniformState();
	  else	  
	    MeanFieldState.SetToRandomPhase(1.0,Manager.GetBoolean("random-amplitude"));
	}
      int MaxEval = 2*NbrSites*Manager.GetInteger("nbr-iter");
      double Energy;
      if (Manager.GetBoolean("gradient"))
	Energy=MeanFieldState.GradientOptimize(Manager.GetDouble("tolerance"), MaxEval, /*initialStep*/ 0.01, /*lineMinPar*/ Manager.GetDouble("tolerance")/10.0);
      else
	
	if (Manager.GetBoolean("simplex"))
	  Energy=MeanFieldState.GradientOptimize(Manager.GetDouble("tolerance"), MaxEval, 1.0);
	else
	  Energy=MeanFieldState.Optimize(Manager.GetDouble("tolerance"), MaxEval);
      RealVector Parameters=MeanFieldState.GetVariationalParameters();
      bool Recorded=false;
      ComplexVector TmpWaveFunction=MeanFieldState.GetWaveFunction();
      for (int k=0; (k<NbrFound) && (Recorded==false); ++k)
	if (Energy <= LowestEnergies[k])
	  {
	    if (Norm(TmpWaveFunction*OptimalWaveFunctions[k])>IdentityThreshold*TmpWaveFunction.Norm()*OptimalWaveFunctions[k].Norm())
	      { // same configuration: simply replace with the one of lower energy
		OptimalWaveFunctions[k]=TmpWaveFunction;
		LowestEnergies[k] = Energy;
		Recorded=true;
	      }
	    else
	      {
		int UpperLimit;
		if (NbrFound < NbrToSave)
		  {
		    LowestEnergies[NbrFound]=LowestEnergies[NbrFound-1];
		    OptimalWaveFunctions[NbrFound]=OptimalWaveFunctions[NbrFound-1];
		    UpperLimit=NbrFound-1;
		    NbrFound++;
		  }
		else UpperLimit = NbrToSave-1;
		for (int s=UpperLimit; s>k; --s)
		  {
		    LowestEnergies[s]=LowestEnergies[s-1];
		    OptimalWaveFunctions[s]=OptimalWaveFunctions[s-1];
		  }
		OptimalWaveFunctions[k]=TmpWaveFunction;
		LowestEnergies[k] = Energy;
		Recorded=true;
	      }
	    }
	else
	  {
	    if (Norm(TmpWaveFunction*OptimalWaveFunctions[k])>IdentityThreshold*TmpWaveFunction.Norm()*OptimalWaveFunctions[k].Norm())
	      Recorded=true;
	  }
      if ((Recorded==false)&&(NbrFound<NbrToSave))
	{
	  OptimalWaveFunctions[NbrFound]=MeanFieldState.GetWaveFunction();
	  LowestEnergies[NbrFound] = Energy;
	  ++NbrFound;
	}
      
      cout << "Found mean field state with energy: "<<Energy<<" and density "<< MeanFieldState.GetNbrParticles()/NbrSites<<endl;

      File.precision(10);
      File << Energy/NbrSites << "\t" << MeanFieldState.GetNbrParticles()/NbrSites << "\t" << Energy/MeanFieldState.GetNbrParticles() + ChemicalPotential;
      File.precision(nearbyint(-log(Manager.GetDouble("tolerance"))/log((double)10.0)));
      for (int i=0; i<Parameters.GetVectorDimension(); ++i)
	File << "\t" << Parameters[i];
      File << endl;
    }
  File.close();

  if (NbrAttempts>1)
    {
      char* SelectedName = AddExtensionToFileName(BaseName,"select.dat");
      ofstream SelectFile ( SelectedName, ios::app );
      int CellPosition[2];
      int Sub;
      RealVector SitePosition;

      for (int i=0; i<NbrFound; ++i)
	{
	  int Counter=0;
	  char* ParameterName = GetUniqueFileName(BaseName,Counter,".par");
	  char* FieldName = ReplaceExtensionToFileName(ParameterName,"par","wf");
	  ofstream VectorField;
	  VectorField.open(FieldName,ios::out);
	  if (Counter==0)
	    SelectFile << "#Count\tE_site\tFilename"<<endl;
	  SelectFile.precision(12);
	  SelectFile << Counter << "\t" << LowestEnergies[i]/NbrSites << "\t" << ParameterName << endl;
	  OptimalWaveFunctions[i].WriteVector(ParameterName);
	  for (int s=0; s<NbrSites; ++s)
	    {
	      Lattice->GetSiteCoordinates(s, CellPosition, Sub);
	      SitePosition = Lattice->GetSitePosition(CellPosition,Sub);
	      VectorField << SitePosition[0] << "\t" << SitePosition[1]
			  << "\t" << (OptimalWaveFunctions[i])[s].Re
			  << "\t" << (OptimalWaveFunctions[i])[s].Im << endl;
	    }
	  VectorField.close();
	  delete [] ParameterName;
	}
      SelectFile.close();
      delete [] SelectedName;
    }

  delete [] BaseName;
  delete [] OutputName;
  delete Lattice;
  delete [] LatticeName;
  if (InitialParameters != NULL)
    delete InitialParameters;
  delete RandomNumberGenerator;
  return 0;
}
