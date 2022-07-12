#include "SphereCoulombEnergy.h"

#include <cmath>
using std::sqrt;
using std::cout;
using std::endl;


double dsqrarg;   // shift declaration to nrutil.cc to suppress warnings if unused in module that includes nrutil.h
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)


// default constructor
SphereCoulombEnergy::SphereCoulombEnergy()
{
  this->NbrFlux=0;
}

// constructor
// nbrFlux = Number of Flux piercing sphere
// width = simple model of finite layer width for 1/sqrt(w^2+r^2)
SphereCoulombEnergy::SphereCoulombEnergy(int nbrFlux, double width)
{
  this->Type=AbstractObservable::RealObservable;
  this->NbrFlux = nbrFlux;
  this->WidthSqr=width*width;
  cout << "nbrFlux="<<nbrFlux<<endl;
  this->Radius = sqrt(0.5*(double)NbrFlux); // the radius is also the inverse magnetic length
  this->Values = new WeightedRealObservable();
  this->NbrObservations=0;
}


// destructor
SphereCoulombEnergy::~SphereCoulombEnergy()
{
  if (NbrFlux>0)
    {
      delete Values;
    }
}

// call to make an observation
void SphereCoulombEnergy::RecordValue(double weight)
{
  int N = this->NbrParticles;
  ++NbrObservations;
  double E=0.0;
  if (WidthSqr>0.0)
    {
      for (int i=1;i<N;i++)
	for(int j=0;j<i;j++)
	  {
	    E+=1.0/sqrt(WidthSqr+SqrNorm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]));
	  }
    }
  else
    {
      for (int i=1;i<N;i++)
	for(int j=0;j<i;j++)
	  {
	    E+=1.0/Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	  }
    }
  
  this->Values->Observe(0.5*E/Radius, weight);
}

// print legend to the given stream
void SphereCoulombEnergy::PrintLegend(std::ostream &output, bool all)
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
void SphereCoulombEnergy::PrintStatus(std::ostream &output, bool all)
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
void SphereCoulombEnergy::WriteDataFile(std::ostream &output)
{
  output << "#  E  \t err  "<<endl;
  output << this->Values->Average()
	 <<"\t"<<this->Values->ErrorEstimate()<<endl;  
}

// set particle collection that the observable operates on
// system = particle collection
void SphereCoulombEnergy::SetParticleCollection(AbstractParticleCollection *system)
{
  if (system->GetCollectionType()!=AbstractParticleCollection::OnSphereCollection)
    {
      cout << "Need a particle collection on the sphere for SphereCoulombEnergy"<<endl;
      exit(1);
    }
  this->System = (ParticleOnSphereCollection*) system;
  this->NbrParticles = System->GetNbrParticles();
  this->System->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);
}
