#include "SphereBilayerCoulombEnergy.h"

#include <cmath>
using std::sqrt;
using std::cout;
using std::endl;


double dsqrarg;   // shift declaration to nrutil.cc to suppress warnings if unused in module that includes nrutil.h
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)


// default constructor
SphereBilayerCoulombEnergy::SphereBilayerCoulombEnergy()
{
  this->NbrSeparations=0;
}

// constructor
// nbrFlux = Number of Flux piercing sphere
// nbrSeparations = number of layer separations where observations should be made
// lowestSeparation = value of the lowest layer separations
// spacing = spacing of further layer separations
SphereBilayerCoulombEnergy::SphereBilayerCoulombEnergy(int nbrFlux, int nbrSeparations, double lowestSeparation, double spacing)
{
  this->Type=AbstractObservable::VectorValued|AbstractObservable::RealObservableT;
  if (nbrSeparations <= 0)
    {
      cout << "Number of layer separation must be > 0 in SphereBilayerCoulombEnergy" << endl;
      exit(1);
    }
  this->NbrFlux = nbrFlux;
  this->Radius = sqrt(0.5*(double)NbrFlux); // the radius is also the inverse magnetic length
  this->NbrSeparations = nbrSeparations;
  this->SqrSeparations = new double [nbrSeparations];
  this->Separations = new double [nbrSeparations];
  this->Energies = new double [nbrSeparations];
  this->Separations[0]=lowestSeparation;
  for (int i=1; i<NbrSeparations; ++i)
    this->Separations[i] = this->Separations[i-1] + spacing;
  for (int i=0; i<NbrSeparations; ++i) // assign squared distances in magnetic lengths
    {
      this->SqrSeparations[i] = DSQR(this->Separations[i]/Radius);
    }
  this->Values = new WeightedRealVectorObservable(NbrSeparations);
  this->NbrObservations=0;
}


// destructor
SphereBilayerCoulombEnergy::~SphereBilayerCoulombEnergy()
{
  if (NbrSeparations>0)
    {
      delete Values;
      delete [] SqrSeparations;
      delete [] Separations;
      delete [] Energies;
    }
}

// call to make an observation
void SphereBilayerCoulombEnergy::RecordValue(double weight)
{
  int N1 = this->NbrParticles>>1;
  ++NbrObservations;
  double E1=0.0, E2=0.0, dij, dsqr;
  for (int i=1;i<N1;i++)
    for(int j=0;j<i;j++)
      {
	E1+=0.5/Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	E2+=0.5/Norm(SpinorUCoordinates[i+N1]*SpinorVCoordinates[j+N1]-SpinorUCoordinates[j+N1]*SpinorVCoordinates[i+N1]);
      }
  for (int s=0; s<NbrSeparations; ++s) Energies[s]=E1+E2;
  for (int i=0; i<N1; ++i)
    for (int j=N1; j<this->NbrParticles; ++j)
      {
	dij=SqrNorm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	dsqr = 4.0*dij;
	for (int s=0; s<NbrSeparations; ++s)
	  {
	    Energies[s]+=1.0/sqrt( dsqr  + SqrSeparations[s] );
	  }
      }
  for (int s=0; s<NbrSeparations; ++s) Energies[s]/=Radius;
  this->Values->Observe(Energies, weight);
}

// print legend to the given stream
void SphereBilayerCoulombEnergy::PrintLegend(std::ostream &output, bool all)
{
  if ((all)||(NbrSeparations<3))
    {
      if (this->NbrSeparations>0)
	output << "d="<<this->Separations[0]<<"\t+/-";
      for (int i=1; i<this->NbrSeparations; ++i)
	output << "\t" << this->Separations[i]<<"\t+/-";  
    }
  else
    {      
      output << "d=" << this->Separations[0]<<"\t+/-";
      output << "\t" << this->Separations[NbrSeparations/2]<<"\t+/-";
      output << "\t" << this->Separations[NbrSeparations-1]<<"\t+/-";
    }
}

// print status to the given stream
void SphereBilayerCoulombEnergy::PrintStatus(std::ostream &output, bool all)
{
  if (NbrObservations>0)
    {
      if ((all)||(NbrSeparations<3))
	{
	  if (this->NbrSeparations>0)
	    output << this->Values->Average(0)<<"\t"<<this->Values->ErrorEstimate(0);
	  for (int i=1; i<this->NbrSeparations; ++i)
	    output << "\t" << this->Values->Average(i)<<"\t"<<this->Values->ErrorEstimate(i);
	}
      else
	{
	  int tmp=output.precision();
	  output.precision(6);
	  output << this->Values->Average(0)<<"\t"<<this->Values->ErrorEstimate(0);
	  output << "\t" << this->Values->Average(NbrSeparations/2)<<"\t"
		 <<this->Values->ErrorEstimate(NbrSeparations/2);
	  output << "\t" << this->Values->Average(NbrSeparations-1)<<"\t"
		 <<this->Values->ErrorEstimate(NbrSeparations-1);
	  output.precision(tmp);
	}
    }
}

// print formatted data suitable for plotting
// ouput = the target stream
void SphereBilayerCoulombEnergy::WriteDataFile(std::ostream &output)
{
  output << "#  d\t  E  \t err  "<<endl;
  int OldPrecision=output.precision();
  output.precision(8);
  for (int i=0; i<this->NbrSeparations; ++i)
    output << this->Separations[i] << "\t" << this->Values->Average(i)
	   <<"\t"<<this->Values->ErrorEstimate(i)<<endl;
  output.precision(OldPrecision);
}

// set particle collection that the observable operates on
// system = particle collection
void SphereBilayerCoulombEnergy::SetParticleCollection(AbstractParticleCollection *system)
{
  this->System = (ParticleOnSphereCollection*) system;
  this->NbrParticles = System->GetNbrParticles();
  if (NbrParticles & 1)
    {
      cout << "SphereBilayerCoulombEnergy requires an even number of particles"<<endl;
      exit(1);
    }  
  this->System->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);
}

// accessor function to return the legend and numerical values for legend
void SphereBilayerCoulombEnergy::GetVectorLegend(std::string &legendParameters, std::string &legendValue, RealVector &parameterValues)
{
  legendParameters = std::string("Separation 'd'");
  legendValue = std::string("Energy 'E'");
  parameterValues.Resize(this->NbrSeparations);
  for (int i=0; i<this->NbrSeparations; ++i)
    parameterValues[i] = this->Separations[i];
}

// accessor function for average and error for variables with real measurements
void SphereBilayerCoulombEnergy::GetRealVectorMeasurement(RealVector &values, RealVector &errors)
{
  for (int i=0; i<this->NbrSeparations; ++i)
    {
      values[i]= this->Values->Average(i);
      errors[i]= this->Values->ErrorEstimate(i);
    }
}

// additional routines for energy observables:
// sep = layer separation
// returns the total background energy for a layer separation
//
// attention: in previous work, convention of E_bg ~ N1(N1-1) for intra-layer background energy
double SphereBilayerCoulombEnergy::GetTotalBackgroundEnergy(double sep)
{
  int N1 = this->NbrParticles>>1;
  return 0.5*(N1*N1+N1*N1)/this->Radius + N1*N1/this->NbrFlux*(sqrt(2.0*this->NbrFlux+sep*sep)-sep); // not tested...
}
