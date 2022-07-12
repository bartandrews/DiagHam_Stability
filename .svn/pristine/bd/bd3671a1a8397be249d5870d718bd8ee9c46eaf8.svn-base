#include "SphereGeneralEnergy.h"
#include "GeneralTools/ConfigurationParser.h"

#include <cmath>
using std::sqrt;
using std::cout;
using std::endl;


double dsqrarg1;   // shift declaration to nrutil.cc to suppress warnings if unused in module that includes nrutil.h
#define DSQR(a) ((dsqrarg1=(a)) == 0.0 ? 0.0 : dsqrarg1*dsqrarg1)


// default constructor
SphereGeneralEnergy::SphereGeneralEnergy()
{
  this->InteractionType=SphereGeneralEnergy::Unknown;
}

// constructor
// nbrFlux = Number of Flux piercing sphere
// parameters = file describing parameters of the interaction
SphereGeneralEnergy::SphereGeneralEnergy(int nbrFlux, const char* parameters)
{
  this->Type=AbstractObservable::RealObservable;
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
  int TmpLength=0;
  this->InteractionType=0;
  if (InteractionDefinition.GetAsStringArray("Type", ' ', TmpType, TmpLength) == false)
    {
      cout << "Assuming polynomial effective interaction!"<<endl;
      this->InteractionType=SphereGeneralEnergy::Polynomial;
    }
  else
    {
      if (TmpLength == 1)
	{
	  if (strcmp(TmpType[0], "Polynomial") == 0)
	    this->InteractionType=SphereGeneralEnergy::Polynomial;	      
	  if (strcmp(TmpType[0], "AsymptoticExp") == 0)
	    this->InteractionType=SphereGeneralEnergy::AsymptoticExp;
	}
      else
	{
	  cout<<"Interaction name must be a single word."<<endl;
	  exit(1);
	}
    }
  for (int i=0; i<TmpLength; ++i)
    delete [] TmpType[0];
  if (TmpLength>0) delete [] TmpType;
  
  if (this->InteractionType==0)
    {
      cout << "Attention, the interaction requested in "<<parameters<<" is not defined, yet."<<endl;
      exit(1);
    }
  if (this->InteractionType==SphereGeneralEnergy::Polynomial)
    {
      int TmpNbrParameters;
      if (InteractionDefinition.GetAsDoubleArray("Parameters", ' ', this->Coefficients, TmpNbrParameters) == false)
	{
	  cout << "Parameters are not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      if (InteractionDefinition.GetAsSingleInteger("NumCoefficients", this->NbrParameters) == false)
	{
	  cout << "NumCoefficients are not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}  
      if (NbrParameters!=TmpNbrParameters)
	{
	  cout << "Values not consistent in " << parameters << endl;
	  exit(-1);
	}
      int TmpNphi;
      if (InteractionDefinition.GetAsSingleInteger("Nphi", TmpNphi) == false)
	{
	  cout << "Nphi is not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      cout << "Using TmpNphi="<<TmpNphi<<endl;
      this->Radius = sqrt(0.5*(double)TmpNphi); // the radius is also the inverse magnetic length

      if (this->NbrParameters==1)
	{
	  double *tmpC = new double[2];
	  tmpC[0]= this->Coefficients[0];
	  tmpC[1]= 0.0;
	  this->NbrParameters=2;
	  delete [] this->Coefficients;
	  this->Coefficients=tmpC;
	}

      if (this->NbrFlux==0)
	this->NbrFlux=TmpNphi;
    }
  else if (this->InteractionType==SphereGeneralEnergy::AsymptoticExp)
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
      if (InteractionDefinition.GetAsDoubleArray("Parameters", ' ', this->Coefficients, TmpNbrParameters) == false)
	{
	  cout << "Parameters are not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      if (InteractionDefinition.GetAsSingleInteger("NumCoefficients", this->NbrParameters) == false)
	{
	  cout << "NumCoefficients is not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      if (InteractionDefinition.GetAsSingleInteger("NumAsymptotics", this->NbrAsymptotics) == false)
	{
	  cout << "NumAsymptotics is not defined or has a wrong value in " << parameters << endl;
	  exit(-1);
	}
      if (NbrParameters+2*NbrAsymptotics!=TmpNbrParameters)
	{
	  cout << "Total number of parameters not consistent in " << parameters << endl;
	  exit(-1);
	}
      this->Asymptotics=this->Coefficients+this->NbrParameters;
      this->AsymptoticsReg = this->Coefficients+this->NbrParameters+this->NbrAsymptotics;
      this->RijSq=NULL;
      int MaxNbrAsymptotics = 2*NbrAsymptotics+1;
      int MaxNbrParameters = NbrParameters - 1;
//       cout << "Parameters used in interaction:"<<endl;
//       for (int i=0; i<NbrParameters; ++i)
// 	cout << "C["<<i<<"]="<<this->Coefficients[i]<<endl;
//       for (int i=0; i<NbrAsymptotics; ++i)
// 	cout << "B["<<i<<"]="<<this->Asymptotics[i]<<", A["<<i<<"]="<<AsymptoticsReg[i]<<endl;
      this->NumSqPowers = (MaxNbrAsymptotics>MaxNbrParameters ? MaxNbrAsymptotics : MaxNbrParameters);
    }
}


// destructor
SphereGeneralEnergy::~SphereGeneralEnergy()
{
  if (this->InteractionType!=SphereGeneralEnergy::Unknown)
    {
      delete Values;
      delete [] Coefficients;
      if (this->InteractionType==SphereGeneralEnergy::AsymptoticExp)
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
void SphereGeneralEnergy::RecordValue(double weight)
{
  if (this->InteractionType==SphereGeneralEnergy::Polynomial)
    {
      int N = this->NbrParticles;
      ++NbrObservations;
      double rst, dij, sum=0.0;
      for (int i=1;i<N;i++)
	{
	  for(int j=0;j<i;j++)
	    {
	      dij = 2.0*Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	      rst = this->Coefficients[0]/ dij;
	      double p = this->Coefficients[this->NbrParameters-1];
	      for (int k=this->NbrParameters-2; k>0; --k)
		{
		  p=p*dij + this->Coefficients[k];
		}
	      rst+=p;
	      sum+=rst;
	    }
	}
      this->Values->Observe(sum/Radius, weight);
    }
  else if (this->InteractionType==SphereGeneralEnergy::AsymptoticExp)
    {
      this->EvaluateGaussianTables();
      int N = this->NbrParticles;
      ++NbrObservations;
      double rst, sum=0.0;
      for (int i=0;i<N;++i)
	for(int j=i+1;j<N;++j)
	  {
	    rst = Asymptotics[0]/sqrt(RijSq[i][j]+AsymptoticsReg[0]);  // B1/sqrt(r^2+d^2)
	    rst += Coefficients[0] * GaussianIJ[i][j]; // gaussian term C0 exp(-r^2)
	    // gaussian terms Ck r^2k exp(-r^2)
	    for (int k=1; k<NbrParameters; ++k)
	      {
		rst += Coefficients[k] * RijSqPowers[k-1][i][j] * GaussianIJ[i][j];
	      }
	    for (int k=1; k<NbrAsymptotics; ++k)
	      {
		  rst += Asymptotics[k] / sqrt(RijSqPowers[k<<1][i][j]+AsymptoticsReg[k]);
	      }
	    sum+=rst;
	  }
      this->Values->Observe(sum, weight);
    }
}

// old version:
// // call to make an observation
// void SphereGeneralEnergy::RecordValue(double weight)
// {
//   int N = this->NbrParticles;
//   ++NbrObservations;
//   double rst, dij, sum=0.0;
//   for (int i=1;i<N;i++)
//     {
//       for(int j=0;j<i;j++)
// 	{
// 	  dij = 2.0*Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
// 	  rst = 1.0/ dij;
// 	  double p = this->Coefficients[this->NbrParameters-1];
// 	  for (int k=this->NbrParameters-2; k>-1; --k)
// 	    {
// 	      p=p*dij + this->Coefficients[k];
// 	    }
// 	  rst+=p;
// 	  sum+=rst;
// 	}
//     }
//   this->Values->Observe(sum/Radius, weight);
// }

// print legend to the given stream
void SphereGeneralEnergy::PrintLegend(std::ostream &output, bool all)
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
void SphereGeneralEnergy::PrintStatus(std::ostream &output, bool all)
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
void SphereGeneralEnergy::WriteDataFile(std::ostream &output)
{
  output << "#  E  \t err  "<<endl;
  output << this->Values->Average()
	 <<"\t"<<this->Values->ErrorEstimate()<<"\t"<<this->GetTotalBackgroundEnergy()<<endl;  
}

// set particle collection that the observable operates on
// system = particle collection
void SphereGeneralEnergy::SetParticleCollection(AbstractParticleCollection *system)
{
  if (system->GetCollectionType()!=AbstractParticleCollection::OnSphereCollection)
    {
      cout << "Need a particle collection on the sphere for SphereCoulombEnergy"<<endl;
      exit(1);
    }
  this->System = (ParticleOnSphereCollection*) system;
  // take care of internal tables which depend on particle number
  if ((this->InteractionType == SphereGeneralEnergy::AsymptoticExp)
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
  this->System->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);

  if (this->InteractionType == SphereGeneralEnergy::AsymptoticExp)
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


#ifdef HAVE_GSL  

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace SphereGeneralEnergyBG
{
  // inner integration of density profile
  double Integrand (double theta, void * params)
  {
    SphereGeneralEnergy *Energy=(SphereGeneralEnergy *)params;
    return sin(theta)*Energy->GetPotentialValue(2.0*Energy->GetRadius()*sin(theta/2.0));
  }
}

#endif


// additional routines for energy observables:
// returns the total background energy
double SphereGeneralEnergy::GetTotalBackgroundEnergy()
{
  if (this->InteractionType == SphereGeneralEnergy::Polynomial)
    {
      double PowerTwo=1.0;
      double Result=0.0;
      for (int i=0; i<this->NbrParameters; ++i)
	{
	  Result+= this->Coefficients[i]*PowerTwo/(1.0+(double)i);
	  PowerTwo*=2.0;
	}
      return Result*this->NbrParticles*this->NbrParticles/(2.0*Radius);
    }
  else if (this->InteractionType == SphereGeneralEnergy::AsymptoticExp)
    {
      double rst=0.0;
#ifdef HAVE_GSL  
      // integration workspaces
      int numInt=1000;
      gsl_integration_workspace * IntW= gsl_integration_workspace_alloc (numInt);
      gsl_function TheIntegrand;
      TheIntegrand.function = &SphereGeneralEnergyBG::Integrand;
      TheIntegrand.params = this;
      double error;
      gsl_integration_qag (&TheIntegrand, 1e-8, M_PI, 0.0, 1e-8, numInt, /* KEY */ GSL_INTEG_GAUSS41,
			   IntW, &rst, &error);
#endif
      return rst*this->NbrParticles*this->NbrParticles/4.0;
    }
  else
    return 0.0;
}
  
// evaluate exponentials and powers of r^2
void SphereGeneralEnergy::EvaluateGaussianTables()
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

// plot effective interaction
// str = stream to write to
ostream & SphereGeneralEnergy::PlotPotential(ostream &str, int numpoints)
{
  str << "# R\tValue\n";
  double Rmax = 2.0*this->Radius;
  double delta = Rmax/(double)numpoints;
  for (int i=1; i<=numpoints; ++i)
    {
      double R=i*delta;
      str << R << "\t" << this->GetPotentialValue(R)<<endl;
    }
  return str;
}

// obtain value of the interaction for a given separation theta between particles
// theta = separation
double SphereGeneralEnergy::GetPotentialValue(double R)
{
  if (this->InteractionType==SphereGeneralEnergy::Polynomial)
    {
      double rst, dij=R/this->Radius;
      rst = this->Coefficients[0]/ dij;
      double p = this->Coefficients[this->NbrParameters-1];
      for (int k=this->NbrParameters-2; k>0; --k)
	{
	  p=p*dij + this->Coefficients[k];
	}
      rst+=p;
      return rst/Radius;
    }
  else if (this->InteractionType==SphereGeneralEnergy::AsymptoticExp)
    {      
      double *RijSqPowers = new double[NumSqPowers];
      double *RijSq = RijSqPowers;
      *RijSq = R*R;
      for (int k=1; k<NumSqPowers; ++k)
	RijSqPowers[k]=RijSqPowers[k-1]* *RijSq;
      double GaussianIJ=std::exp(-*RijSq);
      double rst = Asymptotics[0]/sqrt(*RijSq+AsymptoticsReg[0]);  // B1/sqrt(r^2+d^2)
      rst += Coefficients[0] * GaussianIJ; // gaussian term C0 exp(-r^2)
      // gaussian terms Ck r^2k exp(-r^2)
      for (int k=1; k<NbrParameters; ++k)
	{
	  rst += Coefficients[k] * RijSqPowers[k-1] * GaussianIJ;
	}
      for (int k=1; k<NbrAsymptotics; ++k)
	{
	  rst += Asymptotics[k] / sqrt(RijSqPowers[k<<1]+AsymptoticsReg[k]);
	}
      return rst;
    }
  else return 0.0;
}
