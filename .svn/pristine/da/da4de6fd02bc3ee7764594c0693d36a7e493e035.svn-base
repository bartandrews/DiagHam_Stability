#include "SphereWithSpinGeneralEnergy.h"
#include "GeneralTools/ConfigurationParser.h"

#include <cmath>
using std::sqrt;
using std::cout;
using std::endl;


double dsqrarg1;   // shift declaration to nrutil.cc to suppress warnings if unused in module that includes nrutil.h
#define DSQR(a) ((dsqrarg1=(a)) == 0.0 ? 0.0 : dsqrarg1*dsqrarg1)


// default constructor
SphereWithSpinGeneralEnergy::SphereWithSpinGeneralEnergy()
{
  this->InteractionType=SphereWithSpinGeneralEnergy::Unknown;
}

// constructor
// nbrUp = Number of particles with up spin
// nbrFlux = Number of Flux piercing sphere
// parametersInter = file describing parameters of the interaction (inter-spin)
// parametersIntra = file describing parameters of the interaction (intra-spin)
// parametersIntra2 = file describing parameters of the interaction (intra-spin, other spin-channel)
SphereWithSpinGeneralEnergy::SphereWithSpinGeneralEnergy(int nbrUp, int nbrFlux, const char* parametersInter, const char* parametersIntra, const char* parametersIntra2)
{
  this->Type=AbstractObservable::RealObservable;
  this->NbrUp = nbrUp;
  this->NbrFlux = nbrFlux;
  this->NbrParticles = 0;
  this->Values = new WeightedRealObservable();
  this->NbrObservations=0;

  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(parametersInter) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      exit(-1);
    }
  char **TmpType;
  int TmpLength=0;
  this->InteractionType=0;
  if (InteractionDefinition.GetAsStringArray("Type", ' ', TmpType, TmpLength) == false)
    {
      cout << "Assuming polynomial effective interaction!"<<endl;
      this->InteractionType=SphereWithSpinGeneralEnergy::Polynomial;
    }
  else
    {
      if (TmpLength == 1)
	{
	  if (strcmp(TmpType[0], "Polynomial") == 0)
	    {
	      this->InteractionType=SphereWithSpinGeneralEnergy::Polynomial;	      
	    }
	  else
	    {
	      cout << "For other interaction types, please use a single configuration file"<<endl;
	      exit(1);
	    }
	}
      else
	{
	  cout<<"Assuming interaction to be of polynomial type"<<endl;
	  this->InteractionType=SphereWithSpinGeneralEnergy::Polynomial;	      
	}
    }
  if (this->InteractionType!=SphereWithSpinGeneralEnergy::Polynomial)
    {
      cout << "Attention, this interaction is not implemented when using multiple configuration files"<<endl;
      exit(1);
    }
  for (int i=0; i<TmpLength; ++i)
    delete [] TmpType[0];
  if (TmpLength>0) delete [] TmpType;

  int TmpNbrParameters;
  if (InteractionDefinition.GetAsDoubleArray("Parameters", ' ', this->CoefficientsInter, TmpNbrParameters) == false)
    {
      cout << "Parameters are not defined or has a wrong value in " << parametersInter << endl;
      exit(-1);
    }
  if (InteractionDefinition.GetAsSingleInteger("NumCoefficients", this->NbrParametersInter) == false)
    {
      cout << "NumCoefficients are not defined or has a wrong value in " << parametersInter << endl;
      exit(-1);
    }  
  if (NbrParametersInter!=TmpNbrParameters)
    {
      cout << "Values not consistent in " << parametersInter << endl;
      exit(-1);
    }
  int TmpNphi;
  if (InteractionDefinition.GetAsSingleInteger("Nphi", TmpNphi) == false)
    {
      cout << "Nphi is not defined or has a wrong value in " << parametersInter << endl;
      exit(-1);
    }
  if (this->NbrFlux==0)
    this->NbrFlux=TmpNphi;
  if (NbrFlux!=TmpNphi)
    {
      cout << "Interspin interactions have inconsistent flux value"<<endl;
      exit(-1);
    }
  this->Radius = sqrt(0.5*(double)TmpNphi); // the radius is also the inverse magnetic length


  // Intra-spin interactions:
  ConfigurationParser InteractionDefinition2;
  if (InteractionDefinition2.Parse(parametersIntra) == false)
    {
      InteractionDefinition2.DumpErrors(cout) << endl;
      exit(-1);
    }  
  if (InteractionDefinition2.GetAsDoubleArray("Parameters", ' ', this->CoefficientsIntra, TmpNbrParameters) == false)
	{
	  cout << "Parameters are not defined or has a wrong value in " << parametersIntra << endl;
	  exit(-1);
    }
  if (InteractionDefinition2.GetAsSingleInteger("NumCoefficients", this->NbrParametersIntra) == false)
    {
      cout << "NumCoefficients are not defined or has a wrong value in " << parametersIntra << endl;
      exit(-1);
    }  
  if (NbrParametersIntra!=TmpNbrParameters)
    {
      cout << "Values not consistent in " << parametersIntra << endl;
      exit(-1);
    }
  if (InteractionDefinition2.GetAsSingleInteger("Nphi", TmpNphi) == false)
    {
      cout << "Nphi is not defined or has a wrong value in " << parametersIntra << endl;
      exit(-1);
    }
  if (NbrFlux!=TmpNphi)
    {
      cout << "Intraspin interactions have inconsistent flux value"<<endl;
      exit(-1);
    }

  if (parametersIntra2!=NULL) // have asymmetric interactions?
    {
      ConfigurationParser InteractionDefinition3;
      this->HaveIntra2=true;
      // Intra-spin interactions:
      if (InteractionDefinition3.Parse(parametersIntra2) == false)
	{
	  InteractionDefinition3.DumpErrors(cout) << endl;
	  exit(-1);
	}
      int TmpNbrParameters;
      if (InteractionDefinition3.GetAsDoubleArray("Parameters", ' ', this->CoefficientsIntra2, TmpNbrParameters) == false)
	{
	  cout << "Parameters are not defined or has a wrong value in " << parametersIntra2 << endl;
	  exit(-1);
	}
      if (InteractionDefinition3.GetAsSingleInteger("NumCoefficients", this->NbrParametersIntra2) == false)
	{
	  cout << "NumCoefficients are not defined or has a wrong value in " << parametersIntra2 << endl;
	  exit(-1);
	}  
      if (NbrParametersIntra2!=TmpNbrParameters)
	{
	  cout << "Values not consistent in " << parametersIntra2 << endl;
	  exit(-1);
	}
      if (InteractionDefinition3.GetAsSingleInteger("Nphi", TmpNphi) == false)
	{
	  cout << "Nphi is not defined or has a wrong value in " << parametersIntra2 << endl;
	  exit(-1);
	}
      if (NbrFlux!=TmpNphi)
	{
	  cout << "Second  intraspin interactions have inconsistent flux value"<<endl;
	  exit(-1);
	}
    }
  else
    {
      this->HaveIntra2=false;
      CoefficientsIntra2=CoefficientsIntra;
      NbrParametersIntra2=NbrParametersIntra;
    }
    
}


// constructor (new style, single parameter file)
// nbrUp = Number of particles with up spin
// nbrFlux = Number of Flux piercing sphere
// parameters = file describing parameters of the interaction (all spin channels UpUp, UpDown, DownDown)
SphereWithSpinGeneralEnergy::SphereWithSpinGeneralEnergy(int nbrUp, int nbrFlux, const char* parameters)
{
  this->Type=AbstractObservable::RealObservable;
  this->NbrUp = nbrUp;
  this->NbrParticles = 0;
  this->NbrFlux = nbrFlux;
  this->Values = new WeightedRealObservable();
  this->NbrObservations=0;

  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(parameters) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      exit(-1);
    }
  char **TmpType;
  int TmpLength;
  this->InteractionType=0;
  if (InteractionDefinition.GetAsStringArray("Type", ' ', TmpType, TmpLength) == false)
    {
      cout<<"Please indicate an interaction type in parameter file "<<parameters<<endl;
      exit(1);
    }
  else
    {
      if (TmpLength == 1)
	{
	  if (strcmp(TmpType[0], "Polynomial") == 0)
	    this->InteractionType=SphereWithSpinGeneralEnergy::Polynomial;	      
	  if (strcmp(TmpType[0], "AsymptoticExp") == 0)
	    this->InteractionType=SphereWithSpinGeneralEnergy::AsymptoticExp;
	}
      else
	{
	  cout<<"Interaction name must be a single word."<<endl;
	  exit(1);
	}
    }
  delete [] TmpType[0];
  delete [] TmpType;
  
  if (this->InteractionType==0)
    {
      cout << "Attention, the interaction requested in "<<parameters<<" is not defined, yet."<<endl;
      exit(1);
    }
  if (this->InteractionType==SphereWithSpinGeneralEnergy::Polynomial)
    {
      int TmpNbrParameters;
      if ((InteractionDefinition.GetAsDoubleArray("ParametersUpDown", ' ', this->CoefficientsInter, TmpNbrParameters) == false)&&(InteractionDefinition.GetAsDoubleArray("ParametersInter", ' ', this->CoefficientsInter, TmpNbrParameters) == false))
	{
	  cout << "ParametersUpDown are not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      if ((InteractionDefinition.GetAsSingleInteger("NumCoefficientsUpDown", this->NbrParametersInter) == false)&&
	  (InteractionDefinition.GetAsSingleInteger("NumCoefficientsInter", this->NbrParametersInter) == false))
	{
	  cout << "NumCoefficients are not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}  
      if (NbrParametersInter!=TmpNbrParameters)
	{
	  cout << "Number of parameters (updown) not consistent in " << parameters << endl;
	  exit(-1);
	}
      int TmpNphi;
      if (InteractionDefinition.GetAsSingleInteger("Nphi", TmpNphi) == false)
	{
	  cout << "Nphi is not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      if (this->NbrFlux==0)
	this->NbrFlux=TmpNphi;
      if (NbrFlux!=TmpNphi)
	{
	  cout << "Interactions have inconsistent flux value"<<endl;
	  exit(-1);
	}
      this->Radius = sqrt(0.5*(double)TmpNphi); // the radius is also the inverse magnetic length

      // UpUp/Intra Interactions
      if ((InteractionDefinition.GetAsDoubleArray("ParametersUpUp", ' ', this->CoefficientsIntra, TmpNbrParameters) == false)&&
	  (InteractionDefinition.GetAsDoubleArray("ParametersIntra", ' ', this->CoefficientsIntra, TmpNbrParameters) == false))
	{
	  cout << "Parameters (upup/intra) are not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      if ((InteractionDefinition.GetAsSingleInteger("NumCoefficientsUpUp", this->NbrParametersIntra) == false)&&
	  (InteractionDefinition.GetAsSingleInteger("NumCoefficientsIntra", this->NbrParametersIntra) == false))
	{
	  cout << "NumCoefficients (upup/intra) are not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	} 
      if (NbrParametersIntra!=TmpNbrParameters)
	{
	  cout << "Values not consistent in " << parameters << endl;
	  exit(-1);
	}

      // DownDown interactions:
      if (InteractionDefinition.GetAsDoubleArray("ParametersDownDown", ' ', this->CoefficientsIntra2, TmpNbrParameters) == false)
	{
	  this->HaveIntra2=false;   // UpUp and DownDown identical
	  CoefficientsIntra2=CoefficientsIntra;
	  NbrParametersIntra2=NbrParametersIntra;
	}
      else
	{
	  
	  if (InteractionDefinition.GetAsSingleInteger("NumCoefficientsDownDown", this->NbrParametersIntra2) == false)
	    {
	      cout << "NumCoefficients are not defined or has a wrong value in " << parameters << endl;
	      exit(-1);
	    }  
	  if (NbrParametersIntra2!=TmpNbrParameters)
	    {
	      cout << "Values not consistent in " << parameters << endl;
	      exit(-1);
	    }
	}
    }
  else if (this->InteractionType==SphereWithSpinGeneralEnergy::AsymptoticExp)
    {
      if (this->NbrFlux<=0)
	{
	  cout << "Please indicate the appropriate flux with AsymptoticExp interactions!"<<endl;
	  exit(-1);
	}
      this->Radius = sqrt(0.5*(double)NbrFlux); // the radius is also the inverse magnetic length
      
      /* Coding format for Asymptotic expansion:
	 ParametersUpDown = C0 C1 [...] C_NumCoefficientsUpDown | B1 B3 ... B_(2*NumAsymptoticsUpDown+1) | A1 A3 ... A_(2*NumAsymptoticsUpDown+1)
	 with A0 === d^2 for double layer form of interaction
      */
      int TmpNbrParameters;
      int MaxNbrParameters=0;
      int MaxNbrAsymptotics=0;
      if ((InteractionDefinition.GetAsDoubleArray("ParametersUpDown", ' ', this->CoefficientsInter, TmpNbrParameters) == false)&&
	  (InteractionDefinition.GetAsDoubleArray("ParametersInter", ' ', this->CoefficientsInter, TmpNbrParameters) == false))
	{
	  cout << "Parameters (UpDown) are not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      if ((InteractionDefinition.GetAsSingleInteger("NumCoefficientsUpDown", this->NbrParametersInter) == false)&&
	  (InteractionDefinition.GetAsSingleInteger("NumCoefficientsInter", this->NbrParametersInter) == false))
	{
	  cout << "NumCoefficients (UpDown) is not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      MaxNbrParameters = NbrParametersInter;
      if ((InteractionDefinition.GetAsSingleInteger("NumAsymptoticsUpDown", this->NbrAsymptoticsInter) == false)&&
	  (InteractionDefinition.GetAsSingleInteger("NumAsymptoticsInter", this->NbrAsymptoticsInter) == false))
	{
	  cout << "NumAsymptotics (UpDown) is not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      MaxNbrAsymptotics = NbrAsymptoticsInter;
      if (NbrParametersInter+2*NbrAsymptoticsInter!=TmpNbrParameters)
	{
	  cout << "Total number of parameters (updown) not consistent in " << parameters << endl;
	  exit(-1);
	}
      this->AsymptoticsInter=this->CoefficientsInter+this->NbrParametersInter;
      this->AsymptoticsInterReg = this->CoefficientsInter+this->NbrParametersInter+this->NbrAsymptoticsInter;

      // UpUp/Intra Interactions
      if ((InteractionDefinition.GetAsDoubleArray("ParametersUpUp", ' ', this->CoefficientsIntra, TmpNbrParameters) == false)&&
	  (InteractionDefinition.GetAsDoubleArray("ParametersIntra", ' ', this->CoefficientsIntra, TmpNbrParameters) == false))
	{
	  cout << "Parameters (upup/intra) are not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      if ((InteractionDefinition.GetAsSingleInteger("NumCoefficientsUpUp", this->NbrParametersIntra) == false)&&
	  (InteractionDefinition.GetAsSingleInteger("NumCoefficientsIntra", this->NbrParametersIntra) == false))
	{
	  cout << "NumCoefficients (upup/intra) are not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      if (MaxNbrParameters < NbrParametersIntra) MaxNbrParameters = NbrParametersIntra;
      if ((InteractionDefinition.GetAsSingleInteger("NumAsymptoticsUpUp", this->NbrAsymptoticsIntra) == false)&&
	  (InteractionDefinition.GetAsSingleInteger("NumAsymptoticsIntra", this->NbrAsymptoticsIntra) == false))
	{
	  cout << "NumAsymptotics (upup/intra) is not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      if (MaxNbrAsymptotics < NbrAsymptoticsIntra) MaxNbrAsymptotics = NbrAsymptoticsIntra;
      if (NbrParametersIntra+2*NbrAsymptoticsIntra!=TmpNbrParameters)
	{
	  cout << "Total number of parameters (upup/intra) not consistent in " << parameters << endl;
	  exit(-1);
	}
      this->AsymptoticsIntra=this->CoefficientsIntra+this->NbrParametersIntra;
      this->AsymptoticsIntraReg = this->CoefficientsIntra+this->NbrParametersIntra+this->NbrAsymptoticsIntra;
      // DownDown interactions:
      if (InteractionDefinition.GetAsDoubleArray("ParametersDownDown", ' ', this->CoefficientsIntra2, TmpNbrParameters) == false)
	{
	  this->HaveIntra2=false;   // UpUp and DownDown identical
	  CoefficientsIntra2=CoefficientsIntra;
	  NbrParametersIntra2=NbrParametersIntra;
	  NbrAsymptoticsIntra2=NbrAsymptoticsIntra;
	  AsymptoticsIntra2=AsymptoticsIntra;
	  AsymptoticsIntra2Reg=AsymptoticsIntraReg;
	}
      else
	{	  
	  if (InteractionDefinition.GetAsSingleInteger("NumCoefficientsDownDown", this->NbrParametersIntra2) == false)
	    {
	      cout << "NumCoefficients (downdown/intra2) are not defined or has a wrong value in " << parameters << endl;
	      exit(-1);
	    }
	  if (MaxNbrParameters < NbrParametersIntra2) MaxNbrParameters = NbrParametersIntra2;
	  if (InteractionDefinition.GetAsSingleInteger("NumAsymptoticsDownDown", this->NbrAsymptoticsIntra2) == false)
	    {
	      cout << "NumAsymptotics (downdown/intra2) is not defined or has a wrong value in " << parameters << endl;
	      exit(-1);
	    }
	  if (NbrParametersIntra2+2*NbrAsymptoticsIntra2!=TmpNbrParameters)
	    {
	      cout << "Total number of parameters (downdown/intra2) not consistent in " << parameters << endl;
	      exit(-1);
	    }
	  if (MaxNbrAsymptotics < NbrAsymptoticsIntra2) MaxNbrAsymptotics = NbrAsymptoticsIntra2;
	  this->AsymptoticsIntra2=this->CoefficientsIntra2+this->NbrParametersIntra2;
	  this->AsymptoticsIntra2Reg = this->CoefficientsIntra2+this->NbrParametersIntra2+this->NbrAsymptoticsIntra2;
	}
      this->RijSq=NULL;
      MaxNbrAsymptotics = 2*MaxNbrAsymptotics+1;
      MaxNbrParameters = MaxNbrParameters - 1;
      this->NumSqPowers = (MaxNbrAsymptotics>MaxNbrParameters ? MaxNbrAsymptotics : MaxNbrParameters);
    }

}


// destructor
SphereWithSpinGeneralEnergy::~SphereWithSpinGeneralEnergy()
{
  if (this->InteractionType!=SphereWithSpinGeneralEnergy::Unknown)
    {
      delete Values;
      delete [] CoefficientsInter;
      delete [] CoefficientsIntra;
      if (this->HaveIntra2)
	delete [] CoefficientsIntra2;
      if (this->InteractionType==SphereWithSpinGeneralEnergy::AsymptoticExp)
	{
	  for (int j=0; j<this->NbrParticles; ++j)
	    {
	      delete [] this->RijSq[j];
	      delete [] this->GaussianIJ[j];
	    }
	  delete [] this->RijSq;
	  delete [] this->GaussianIJ;
	  
	  for (int i=1; i<this->NumSqPowers; ++i)
	    {	      
	      for (int j=0; j<this->NbrParticles; ++j)	
		delete [] this->RijSqPowers[i][j];
	      delete [] this->RijSqPowers[i];
	    }
	  delete [] this->RijSqPowers;
	}
    }
}

// call to make an observation
void SphereWithSpinGeneralEnergy::RecordValue(double weight)
{
  if (this->InteractionType==SphereWithSpinGeneralEnergy::Polynomial)
    {
      int N = this->NbrParticles;
      ++NbrObservations;
      double rst, dij, sum=0.0;
      for (int i=1;i<this->NbrUp;i++)
	{
	  for(int j=0;j<i;j++)
	    {
	      dij = 2.0*Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	      rst = this->CoefficientsIntra[0] / dij;
	      double p = this->CoefficientsIntra[this->NbrParametersIntra-1];
	      for (int k=this->NbrParametersIntra-2; k>0; --k)
		{
		  p=p*dij + this->CoefficientsIntra[k];
		}
	      rst+=p;
	      sum+=rst;
	    }
	}
      for (int i=this->NbrUp+1;i<N;i++)
	{
	  for(int j=this->NbrUp;j<i;j++)
	    {
	      dij = 2.0*Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	      rst = this->CoefficientsIntra2[0] / dij;
	      double p = this->CoefficientsIntra2[this->NbrParametersIntra2-1];
	      for (int k=this->NbrParametersIntra2-2; k>0; --k)
		{
		  p=p*dij + this->CoefficientsIntra2[k];
		}
	      rst+=p;
	      sum+=rst;
	    }
	}
      for (int i=0;i<this->NbrUp;i++)
	{
	  for(int j=this->NbrUp;j<N;j++)
	    {
	      dij = 2.0*Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	      rst = this->CoefficientsInter[0]/ dij;
	      double p = this->CoefficientsInter[this->NbrParametersInter-1];
	      for (int k=this->NbrParametersInter-2; k>0; --k)
		{
		  p=p*dij + this->CoefficientsInter[k];
		}
	      rst+=p;
	      sum+=rst;
	    }
	}
      this->Values->Observe(sum/Radius, weight);
    }
  else if (this->InteractionType==SphereWithSpinGeneralEnergy::AsymptoticExp)
    {
      this->EvaluateGaussianTables();
      int N = this->NbrParticles;
      ++NbrObservations;
      double rst, sum=0.0;
      for (int i=0;i<this->NbrUp;++i)
	for(int j=i+1;j<this->NbrUp;++j)
	  {
	    rst = AsymptoticsIntra[0]/sqrt(RijSq[i][j]+AsymptoticsIntraReg[0]);  // B1/sqrt(r^2+d^2)
	    rst += CoefficientsIntra[0] * GaussianIJ[i][j]; // gaussian term C0 exp(-r^2)
	    // gaussian terms Ck r^2k exp(-r^2)
	    for (int k=1; k<NbrParametersIntra; ++k)
	      {
		rst += CoefficientsIntra[k] * RijSqPowers[k-1][i][j] * GaussianIJ[i][j];
	      }
	    for (int k=1; k<NbrAsymptoticsIntra; ++k)
	      {
		  rst += AsymptoticsIntra[k] / sqrt(RijSqPowers[k<<1][i][j]+AsymptoticsIntraReg[k]);
	      }
	    sum+=rst;
	  }
      for (int i=this->NbrUp;i<N;++i)
	for(int j=i+1;j<N;++j)
	  {
	    rst = AsymptoticsIntra2[0]/sqrt(RijSq[i][j]+AsymptoticsIntra2Reg[0]);  // B1/sqrt(r^2+d^2)
	    rst += CoefficientsIntra2[0] * GaussianIJ[i][j]; // gaussian term C0 exp(-r^2)
	    // gaussian terms Ck r^2k exp(-r^2)
	    for (int k=1; k<NbrParametersIntra2; ++k)
	      {
		rst += CoefficientsIntra2[k] * RijSqPowers[k-1][i][j] * GaussianIJ[i][j];
	      }
	    for (int k=1; k<NbrAsymptoticsIntra2; ++k)
	      {
		  rst += AsymptoticsIntra2[k] / sqrt(RijSqPowers[k<<1][i][j]+AsymptoticsIntra2Reg[k]);
	      }
	    sum+=rst;
	  }
      for (int i=0;i<this->NbrUp;i++)
	for(int j=this->NbrUp;j<N;j++)
	  {
	    rst = AsymptoticsInter[0]/sqrt(RijSq[i][j]+AsymptoticsInterReg[0]);  // B1/sqrt(r^2+d^2)
	    rst += CoefficientsInter[0] * GaussianIJ[i][j]; // gaussian term C0 exp(-r^2)
	    // gaussian terms Ck r^2k exp(-r^2)
	    for (int k=1; k<NbrParametersInter; ++k)
	      {
		rst += CoefficientsInter[k] * RijSqPowers[k-1][i][j] * GaussianIJ[i][j];
	      }
	    for (int k=1; k<NbrAsymptoticsInter; ++k)
	      {
		  rst += AsymptoticsInter[k] / sqrt(RijSqPowers[k<<1][i][j]+AsymptoticsInterReg[k]);
	      }
	    sum+=rst;
	  }
      this->Values->Observe(sum, weight);
    }
}

// print legend to the given stream
void SphereWithSpinGeneralEnergy::PrintLegend(std::ostream &output, bool all)
{
  if (all)
    {
      output << "E\t+/-";
    }
  else
    {
      output << "E\t+/-";
    }
}

// print status to the given stream
void SphereWithSpinGeneralEnergy::PrintStatus(std::ostream &output, bool all)
{
  if (NbrObservations>0)
    {
      if (all)
	{
	  output << this->Values->Average()<<"\t"<<this->Values->ErrorEstimate();
	}
      else
	{
	  int tmp=output.precision();
	  output.precision(6);
	  output << this->Values->Average()<<"\t"<<this->Values->ErrorEstimate();
	  output.precision(tmp);
	}
    }
}

// print formatted data suitable for plotting
// ouput = the target stream
void SphereWithSpinGeneralEnergy::WriteDataFile(std::ostream &output)
{
  output << "#  E  \t err  "<<endl;
  output << this->Values->Average()
	 <<"\t"<<this->Values->ErrorEstimate()<<"\t"<<this->GetTotalBackgroundEnergy()<<endl;  
}

// set particle collection that the observable operates on
// system = particle collection
void SphereWithSpinGeneralEnergy::SetParticleCollection(AbstractParticleCollection *system)
{
  if (system->GetCollectionType()!=AbstractParticleCollection::OnSphereCollection)
    {
      cout << "Need a particle collection on the sphere for SphereCoulombEnergy"<<endl;
      exit(1);
    }
  this->System = (ParticleOnSphereCollection*) system;
  // take care of internal tables which depend on particle number
  if ((this->InteractionType == SphereWithSpinGeneralEnergy::AsymptoticExp)
      &&(this->RijSq!=NULL))
    {
      for (int j=0; j<this->NbrParticles; ++j)
	{
	  delete [] this->RijSq[j];
	  delete [] this->GaussianIJ[j];
	}
      delete [] this->RijSq;
      delete [] this->GaussianIJ;
      
      for (int i=1; i<this->NumSqPowers; ++i)
	{	      
	  for (int j=0; j<this->NbrParticles; ++j)	
	    delete [] this->RijSqPowers[i][j];
	  delete [] this->RijSqPowers[i];
	}
      delete [] this->RijSqPowers;
    }
  this->NbrParticles = System->GetNbrParticles();
  if (this->NbrParticles<this->NbrUp)
    {
      cout << "The number of particles needs to be superiour to the assigned value of up-spins in SphereWithSpinGeneralEnergy" << endl;
      exit(1);
    }
  this->System->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);

  if (this->InteractionType == SphereWithSpinGeneralEnergy::AsymptoticExp)
    {
      this->RijSq=new double*[NbrParticles];
      this->GaussianIJ=new double*[NbrParticles];
      this->RijSqPowers=new double**[this->NumSqPowers];
      for (int i=1; i<this->NumSqPowers; ++i)
	{
	  this->RijSqPowers[i]=new double*[NbrParticles];
	  for (int j=0; j<this->NbrParticles; ++j)	
	    this->RijSqPowers[i][j]=new double[NbrParticles];
	}
      for (int j=0; j<this->NbrParticles; ++j)
	{
	  this->RijSq[j]=new double[NbrParticles];
	  this->GaussianIJ[j]=new double[NbrParticles];
	}
      this->RijSqPowers[0]=this->RijSq;
    }
}

// additional routines for energy observables:
// returns the total background energy
double SphereWithSpinGeneralEnergy::GetTotalBackgroundEnergy()
{
  double PowerTwo=1.0;
  double ResultUU=0.0;
  double ResultUD=0.0;
  double ResultDD=0.0;
  for (int i=0; i<this->NbrParametersIntra; ++i)
    {
      ResultUU+= this->CoefficientsIntra[i]*PowerTwo/(1.0+(double)i);
      PowerTwo*=2.0;
    }
  ResultUU*=this->NbrUp*this->NbrUp/(2.0*Radius);
  
  PowerTwo=1.0;
  for (int i=0; i<this->NbrParametersIntra2; ++i)
    {
      ResultDD+= this->CoefficientsIntra2[i]*PowerTwo/(1.0+(double)i);
      PowerTwo*=2.0;
    }
  ResultDD*=(this->NbrParticles-this->NbrUp)*(this->NbrParticles-this->NbrUp)/(2.0*Radius);

  PowerTwo=1.0;
  for (int i=0; i<this->NbrParametersInter; ++i)
    {
      ResultUD+= this->CoefficientsInter[i]*PowerTwo/(1.0+(double)i);
      PowerTwo*=2.0;
    }
  ResultUD*=this->NbrUp*(this->NbrParticles-this->NbrUp)/Radius;

  return ResultUU+ResultDD+ResultUD;
}

// evaluate exponentials and powers of r^2
void SphereWithSpinGeneralEnergy::EvaluateGaussianTables()
{
  double factor=DSQR(2.0*Radius);
  for (int i=0; i<NbrParticles; ++i)
    {
      for (int j=i+1; j<NbrParticles; ++j)
	{
	  RijSq[i][j] = factor*SqrNorm(SpinorUCoordinates[i]*SpinorVCoordinates[j]
				    -SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	  for (int k=1; k<NumSqPowers; ++k)
	    RijSqPowers[k][i][j]=RijSqPowers[k-1][i][j]*RijSq[i][j];
	  GaussianIJ[i][j]=exp(-RijSq[i][j]);
	}
    }
}
