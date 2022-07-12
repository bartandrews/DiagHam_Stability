#include "SphereBlockCoulombEnergy.h"

#include <cmath>
using std::sqrt;
using std::cout;
using std::endl;


double dsqrarg;   // shift declaration to nrutil.cc to suppress warnings if unused in module that includes nrutil.h
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)


// default constructor
SphereBlockCoulombEnergy::SphereBlockCoulombEnergy()
{
  this->NbrFlux=0;
}

// constructor
// nbrFlux = Number of Flux piercing sphere
// nbrBlockEntries = total number of entries for block vector observable  
// width = simple model of finite layer width for 1/sqrt(w^2+r^2)
SphereBlockCoulombEnergy::SphereBlockCoulombEnergy(int nbrBlockEntries, int nbrFlux, double width, bool useComplex)
{
  this->Type=AbstractObservable::RealObservableT;
  this->NbrFlux = nbrFlux;
  this->WidthSqr=width*width;
  cout << "nbrFlux="<<nbrFlux<<endl;
  this->Radius = sqrt(0.5*(double)NbrFlux); // the radius is also the inverse magnetic length
  this->NbrObservations=0;
  this->NbrBlocks = nbrBlockEntries;
  this->UseComplex = useComplex;
  if(this->UseComplex)
    {
      this->ComplexValues = new ComplexVectorObservable(nbrBlockEntries);
      this->TempComplex = new Complex[nbrBlockEntries];
    }
  else
    {
      this->RealValues = new RealObservable[nbrBlockEntries*2];
    }
}


// destructor
SphereBlockCoulombEnergy::~SphereBlockCoulombEnergy()
{
  if (NbrFlux>0)
    {
      if(this->UseComplex)
	{
	  delete ComplexValues;
	  delete [] this->TempComplex;
	}
      else
	{
	  delete [] this->RealValues;
	}
    }
}

// call to make an observation
// weight = relative weight of this sample
void SphereBlockCoulombEnergy::RecordRealValue(Complex *weights)
{
  double E=this->GetEnergy();
  if (this->UseComplex)
    {
      for (int i=0; i<NbrBlocks; ++i)
	this->TempComplex[i]=E*weights[i].Re;
      this->ComplexValues->Observe(TempComplex);
    }
  else
    {
      for (int i=0; i<NbrBlocks; ++i)
	{
	  this->RealValues[i]<<(E*weights[i].Re);
	}
    }
}

// call to make an observation
// weight = relative weight of this sample
void SphereBlockCoulombEnergy::RecordValue(Complex *weights)
{
  double E=this->GetEnergy();
  if (this->UseComplex)
    {
      for (int i=0; i<NbrBlocks; ++i)
	this->TempComplex[i]=E*weights[i];
      this->ComplexValues->Observe(TempComplex);
    }
  else
    {
      for (int i=0, j=NbrBlocks; i<NbrBlocks; ++i, ++j)
	{
	  this->RealValues[i]<<(E*weights[i].Re);
	  this->RealValues[j]<<(E*weights[i].Im);
	}
    }
}
  

// call to make an observation
void SphereBlockCoulombEnergy::RecordValue(double weight)
{
  cout << "Function SphereBlockCoulombEnergy::RecordValue(double weight) is not implemented"<<endl;
}


// print legend to the given stream
void SphereBlockCoulombEnergy::PrintLegend(std::ostream &output, bool all)
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
void SphereBlockCoulombEnergy::PrintStatus(std::ostream &output, bool all)
{
  if (NbrObservations>0)
    {
      if (all)
	{
	  if (UseComplex)
	    output << this->ComplexValues->Average(NbrBlocks-1)<<"\t"
		   <<this->ComplexValues->ErrorEstimateRe(NbrBlocks-1)<<"+I*"<<this->ComplexValues->ErrorEstimateIm(NbrBlocks-1);
	  else
	    output << this->RealValues[NbrBlocks-1].Average()<<"+I*"<<this->RealValues[2*NbrBlocks-1].Average()
	     <<"\t"<<this->RealValues[NbrBlocks-1].ErrorEstimate()<<"+I*"<<this->RealValues[2*NbrBlocks-1].ErrorEstimate();
	  cout << endl;
	}
      else
	{
	  int tmp=output.precision();
	  output.precision(6);
	  if (UseComplex)
	    output << this->ComplexValues->Average(NbrBlocks-1)<<"\t"
		   <<this->ComplexValues->ErrorEstimateRe(NbrBlocks-1)<<"+I*"<<this->ComplexValues->ErrorEstimateIm(NbrBlocks-1);
	  else
	    output << this->RealValues[NbrBlocks-1].Average()<<"+I*"<<this->RealValues[2*NbrBlocks-1].Average()
		   <<"\t"<<this->RealValues[NbrBlocks-1].ErrorEstimate()<<"+I*"<<this->RealValues[2*NbrBlocks-1].ErrorEstimate();
	  cout << endl;
	  output.precision(tmp);
	}
    }
}

// print formatted data suitable for plotting
// ouput = the target stream
void SphereBlockCoulombEnergy::WriteDataFile(std::ostream &output)
{
  output << "#  E  \t err  \t  Blocks"<<endl;
  if (UseComplex)
    {
      output << this->ComplexValues->Average(NbrBlocks-1)
	     <<"\t"<<this->ComplexValues->ErrorEstimateRe(NbrBlocks-1)<<"+I*"<<this->ComplexValues->ErrorEstimateIm(NbrBlocks-1);
      for (int i=0; i<NbrBlocks-1; ++i)
	output << this->ComplexValues->Average(i)
	       <<"\t"<<this->ComplexValues->ErrorEstimateRe(i)<<"+I*"<<this->ComplexValues->ErrorEstimateIm(i);
      output <<endl;
    }
  else
    {
      output << this->RealValues[NbrBlocks-1].Average()<<"+I*"<<this->RealValues[2*NbrBlocks-1].Average()
	     <<"\t"<<this->RealValues[NbrBlocks-1].ErrorEstimate()<<"+I*"<<this->RealValues[2*NbrBlocks-1].ErrorEstimate();
      for (int i=0; i<NbrBlocks-1; ++i)
	output << this->RealValues[i].Average()<<"+I*"<<this->RealValues[NbrBlocks+i].Average()
	     <<"\t"<<this->RealValues[i].ErrorEstimate()<<"+I*"<<this->RealValues[NbrBlocks+i].ErrorEstimate();
      output <<endl;
    }
}

// set particle collection that the observable operates on
// system = particle collection
void SphereBlockCoulombEnergy::SetParticleCollection(AbstractParticleCollection *system)
{
  if (system->GetCollectionType()!=AbstractParticleCollection::OnSphereCollection)
    {
      cout << "Need a particle collection on the sphere for SphereBlockCoulombEnergy"<<endl;
      exit(1);
    }
  this->System = (ParticleOnSphereCollection*) system;
  this->NbrParticles = System->GetNbrParticles();
  this->System->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);
  cout << "Set particle collection for Energy observable with N="<<this->NbrParticles<<endl;
}

// additional routines for energy observables:
// returns the total background energy
double SphereBlockCoulombEnergy::GetTotalBackgroundEnergy()
{
  return 0.5*this->NbrParticles*this->NbrParticles/this->Radius; // not tested...
}
