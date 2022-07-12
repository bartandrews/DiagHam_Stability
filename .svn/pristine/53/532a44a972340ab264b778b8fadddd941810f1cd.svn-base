#include "ParticleOnSphereCollection.h"

// switch for debugging output:
//#define DEBUG_OUTPUT


ParticleOnSphereCollection::ParticleOnSphereCollection()
{
  this->NbrParticles = 0;
}

ParticleOnSphereCollection::ParticleOnSphereCollection(int N, long seed)
{
  this->NbrParticles = N;
  this->SpinorUCoordinates = new Complex[N];
  this->SpinorVCoordinates = new Complex[N];
  this->LastMoved=-1;
  this->Flag.Initialize();
  this->Generator = new NumRecRandomGenerator(seed);
  ExternalGenerator = false;
  this->Theta0 = (0.2*M_PI/sqrt((double)N/10.0));
  this->ThetaPhi.Resize(2*N);  
  double phi, theta,  c, s;
  for (int i=0; i<N; ++i)
    {
      phi=Generator->GetRealRandomNumber()*2.0*M_PI;
      theta=acos(1.0-2.0*Generator->GetRealRandomNumber());
      s=sin(theta/2.); c = cos(theta/2.);
      this->SpinorUCoordinates[i].Re =c*cos(phi/2.0);
      this->SpinorUCoordinates[i].Im =-c*sin(phi/2.0);      
      this->SpinorVCoordinates[i].Re = s*cos(phi/2.0);
      this->SpinorVCoordinates[i].Im = s*sin(phi/2.0);
      this->ThetaPhi[i<<1] = (2.0*acos(Norm(SpinorUCoordinates[i])));
      this->ThetaPhi[(i<<1)+1] = (Arg(SpinorVCoordinates[i])-Arg(SpinorUCoordinates[i]));
    }
}

ParticleOnSphereCollection::ParticleOnSphereCollection(int N, AbstractRandomNumberGenerator *generator)
{
  this->NbrParticles = N;
  this->SpinorUCoordinates = new Complex[N];
  this->SpinorVCoordinates = new Complex[N];
  this->LastMoved=-1;
  this->Flag.Initialize();
  this->Generator = generator;
  this->ExternalGenerator = true;
  this->Theta0 = (0.2*M_PI/sqrt((double)N/10.0));
  this->ThetaPhi.Resize(2*N);  
  double phi, theta,  c, s;
  for (int i=0; i<N; ++i)
    {
      phi=Generator->GetRealRandomNumber()*2.0*M_PI;
      theta=acos(1.0-2.0*Generator->GetRealRandomNumber());
      s=sin(theta/2.); c = cos(theta/2.);
      this->SpinorUCoordinates[i].Re =c*cos(phi/2.0);
      this->SpinorUCoordinates[i].Im =-c*sin(phi/2.0);      
      this->SpinorVCoordinates[i].Re = s*cos(phi/2.0);
      this->SpinorVCoordinates[i].Im = s*sin(phi/2.0);
      this->ThetaPhi[i<<1] = (2.0*acos(Norm(SpinorUCoordinates[i])));
      this->ThetaPhi[(i<<1)+1] = (Arg(SpinorVCoordinates[i])-Arg(SpinorUCoordinates[i]));
    }
}


ParticleOnSphereCollection::ParticleOnSphereCollection(const ParticleOnSphereCollection &tocopy)
{
  this->NbrParticles = tocopy.NbrParticles;
  this->LastMoved = tocopy.LastMoved;
  this->LastU = tocopy.LastU;
  this->LastV = tocopy.LastV;
  this->LastTheta = tocopy.LastTheta;
  this->LastPhi = tocopy.LastPhi;
  this->SpinorUCoordinates = tocopy.SpinorUCoordinates;
  this->SpinorVCoordinates = tocopy.SpinorVCoordinates;
  this->Flag  = tocopy.Flag;
  this->Generator = tocopy.Generator;
  this->ExternalGenerator = tocopy.ExternalGenerator;
  this->ThetaPhi = tocopy.ThetaPhi;
  this->Theta0 = tocopy.Theta0;
}

ParticleOnSphereCollection::~ParticleOnSphereCollection()
{
  if ((this->Flag.Used() == true)&&(this->Flag.Shared() == false))
    {
      delete [] SpinorUCoordinates;
      delete [] SpinorVCoordinates;
      if (!this->ExternalGenerator)
	delete Generator;
    }
}

// randomly moves particle number nbrParticle
void ParticleOnSphereCollection::Move(int nbrParticle)
{
  // store old positions
  this->LastMoved = nbrParticle;
  this->LastU = SpinorUCoordinates[nbrParticle];
  this->LastV = SpinorVCoordinates[nbrParticle];
  this->LastTheta = this->ThetaPhi[nbrParticle<<1];
  this->LastPhi = this->ThetaPhi[(nbrParticle<<1)+1];
  // make a random move with gaussian distribution around the initial position
  double theta = this->Theta0*fabs(Generator->GetGaussianRandomNumber());
  double phi = 2.0*M_PI*Generator->GetRealRandomNumber();
  Complex move_v = sin(theta/2.0)*Polar(1.0, phi/2.0);
  Complex move_u = cos(theta/2.0)*Polar(1.0,-phi/2.0);  
  SpinorUCoordinates[nbrParticle] =LastU*move_u-Conj(LastV)*move_v;
  SpinorVCoordinates[nbrParticle] =LastV*move_u+Conj(LastU)*move_v;
  this->ThetaPhi[nbrParticle<<1] = (2.0*acos(Norm(SpinorUCoordinates[nbrParticle])));
  this->ThetaPhi[(nbrParticle<<1)+1] = (Arg(SpinorVCoordinates[nbrParticle])-Arg(SpinorUCoordinates[nbrParticle]));
  // cout << "new coordinates: ("<<this->ThetaPhi[nbrParticle<<1]<<", " << this->ThetaPhi[(nbrParticle<<1)+1] << ")" <<endl;
}

// randomly select a particle and move it
int ParticleOnSphereCollection::Move()
{
  this->LastMoved = (int) (((double) NbrParticles) * Generator->GetRealRandomNumber());
  if (LastMoved == NbrParticles) --LastMoved;
  this->Move(LastMoved);
  return LastMoved;
}

// restore last move
void ParticleOnSphereCollection::RestoreMove()
{
  if (LastMoved>-1)
    {      
      SpinorUCoordinates[LastMoved] = this->LastU;
      SpinorVCoordinates[LastMoved] = this->LastV;
      this->ThetaPhi[LastMoved<<1] = this->LastTheta;
      this->ThetaPhi[(LastMoved<<1)+1] = this->LastPhi;
      this->LastMoved = -1;
    }
}


void ParticleOnSphereCollection::MultiplyStepLength(double multiplier)
{
  this->Theta0*=multiplier;
}

// get single position theta
double ParticleOnSphereCollection::Theta(int nbrParticle)
{
  return this->ThetaPhi[nbrParticle<<1];
}


// get single position phi
double ParticleOnSphereCollection::Phi(int nbrParticle)
{
  return this->ThetaPhi[(nbrParticle<<1)+1];
}

// set single position
void ParticleOnSphereCollection::SetPosition(int nbrParticle, double theta, double phi)
{
  double s2, s=sin(theta/2.);
  double c2, c = cos(theta/2.);
  this->SpinorUCoordinates[nbrParticle].Re =c*(c2=cos(phi/2.0));
  this->SpinorUCoordinates[nbrParticle].Im =-c*(s2=sin(phi/2.0));      
  this->SpinorVCoordinates[nbrParticle].Re = s*c2;
  this->SpinorVCoordinates[nbrParticle].Im = s*s2;  
}

RealVector& ParticleOnSphereCollection::GetPositions()
{
  return this->ThetaPhi;
}

double ParticleOnSphereCollection::GetRandomNumber()
{
  return Generator->GetRealRandomNumber();
}

// randomize particle positions
void ParticleOnSphereCollection::Randomize()
{
  double phi, theta,  c, s;
  for (int i=0; i<NbrParticles; ++i)
    {
      phi=Generator->GetRealRandomNumber()*2.0*M_PI;
      theta=acos(1.0-2.0*Generator->GetRealRandomNumber());
      s=sin(theta/2.); c = cos(theta/2.);
      this->SpinorUCoordinates[i].Re =c*cos(phi/2.0);
      this->SpinorUCoordinates[i].Im =-c*sin(phi/2.0);      
      this->SpinorVCoordinates[i].Re = s*cos(phi/2.0);
      this->SpinorVCoordinates[i].Im = s*sin(phi/2.0);
      this->ThetaPhi[i<<1] = (2.0*acos(Norm(SpinorUCoordinates[i])));
      this->ThetaPhi[(i<<1)+1] = (Arg(SpinorVCoordinates[i])-Arg(SpinorUCoordinates[i]));
    }
}
