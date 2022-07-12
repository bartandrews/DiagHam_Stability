#include "ParticleOnSphereCollectionGauged.h"
#include <iostream>
using std::cout;
using std::endl;
// switch for debugging output:
//#define DEBUG_OUTPUT


ParticleOnSphereCollectionGauged::ParticleOnSphereCollectionGauged()
{
  this->NbrParticles = 0;
}

ParticleOnSphereCollectionGauged::ParticleOnSphereCollectionGauged(int N, long seed)
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
  this->Randomize();

  this->Distances = new RealSymmetricMatrix(NbrParticles, true);
  this->DistancesUpToDate = false;

}

ParticleOnSphereCollectionGauged::ParticleOnSphereCollectionGauged(int N, AbstractRandomNumberGenerator *generator)
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
  this->Randomize();

  this->Distances = new RealSymmetricMatrix(NbrParticles, true);
  this->DistancesUpToDate = false;

}


ParticleOnSphereCollectionGauged::ParticleOnSphereCollectionGauged(const ParticleOnSphereCollectionGauged &tocopy)
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

ParticleOnSphereCollectionGauged::~ParticleOnSphereCollectionGauged()
{
  if ((this->Flag.Used() == true)&&(this->Flag.Shared() == false))
    {
      delete [] SpinorUCoordinates;
      delete [] SpinorVCoordinates;
      if (!this->ExternalGenerator)
	delete Generator;
      delete [] Distances;
    }
}

// randomly moves particle number nbrParticle
void ParticleOnSphereCollectionGauged::Move(int nbrParticle)
{
  // store old positions
  this->LastMoved = nbrParticle;
  this->LastU = SpinorUCoordinates[nbrParticle].Re;
  this->LastV = SpinorVCoordinates[nbrParticle];
  this->LastTheta = this->ThetaPhi[nbrParticle<<1];
  this->LastPhi = this->ThetaPhi[(nbrParticle<<1)+1];
  // make a random move with gaussian distribution around the initial position
  double theta = this->Theta0*fabs(Generator->GetGaussianRandomNumber());
  double phi = 2.0*M_PI*Generator->GetRealRandomNumber();
  Complex move_v = sin(theta/2.0)*Polar(phi);
  double move_u = cos(theta/2.0);
  SpinorUCoordinates[nbrParticle] =LastU*move_u-Conj(LastV)*move_v;
  SpinorVCoordinates[nbrParticle] =LastV*move_u+LastU*move_v;
  double n=Norm(SpinorUCoordinates[nbrParticle]);
  if (n>1e-13)
    {
      Complex Phase = Conj(SpinorUCoordinates[nbrParticle]);
      SpinorUCoordinates[nbrParticle] = n;
      n=1.0/n;
      Phase *= n;
      SpinorVCoordinates[nbrParticle] *= Phase;
    }
  else
    {
      double phase = -Arg(SpinorUCoordinates[nbrParticle]);
      Complex Phase = Polar(phase);
      SpinorUCoordinates[nbrParticle] *= Phase;
      SpinorUCoordinates[nbrParticle].Im = 0.0;
      SpinorVCoordinates[nbrParticle] *= Phase;
    }
  this->ThetaPhi[nbrParticle<<1] = (2.0*acos(SpinorUCoordinates[nbrParticle].Re));
  this->ThetaPhi[(nbrParticle<<1)+1] = Arg(SpinorVCoordinates[nbrParticle]);
  // cout << "new coordinates: ("<<this->ThetaPhi[nbrParticle<<1]<<", " << this->ThetaPhi[(nbrParticle<<1)+1] << ")" <<endl;
}

// randomly select a particle and move it
int ParticleOnSphereCollectionGauged::Move()
{
  this->LastMoved = (int) (((double) NbrParticles) * Generator->GetRealRandomNumber());
  if (LastMoved == NbrParticles) --LastMoved;
  this->Move(LastMoved);
  return LastMoved;
}


// randomly moves all particles 
void ParticleOnSphereCollectionGauged::RotateAll()
{
  double theta = this->Theta0*fabs(Generator->GetGaussianRandomNumber());
  double phi = 2.0*M_PI*Generator->GetRealRandomNumber();
  this->RotateAll(theta, phi);
}

// moves all particles by angles theta and phi
void ParticleOnSphereCollectionGauged::RotateAll(double theta, double phi)
{
  Complex move_v = sin(theta/2.0)*Polar(phi);
  double move_u = cos(theta/2.0);
  Complex lastV;
  double lastU;
  for (int nbrParticle=0; nbrParticle<NbrParticles; ++nbrParticle)
    {
      lastU=SpinorUCoordinates[nbrParticle].Re;
      lastV=SpinorVCoordinates[nbrParticle];
      SpinorUCoordinates[nbrParticle] =lastU*move_u+lastV*Conj(move_v);
      SpinorVCoordinates[nbrParticle] =-lastU*move_v+lastV*move_u;

      double phase = -Arg(SpinorUCoordinates[nbrParticle]);
      Complex Phase = Polar(phase);
      SpinorUCoordinates[nbrParticle] *= Phase;
      SpinorUCoordinates[nbrParticle].Im = 0.0;
      SpinorVCoordinates[nbrParticle] *= Phase;

      this->ThetaPhi[nbrParticle<<1] = (2.0*acos(SpinorUCoordinates[nbrParticle].Re));
      this->ThetaPhi[(nbrParticle<<1)+1] = Arg(SpinorVCoordinates[nbrParticle]);
      cout << "new coordinates: ("<<this->ThetaPhi[nbrParticle<<1]<<", " << this->ThetaPhi[(nbrParticle<<1)+1] << ") Norm" << SqrNorm(SpinorUCoordinates[nbrParticle])+SqrNorm(SpinorVCoordinates[nbrParticle])<<endl;
    }

  // rotate old positions, as well
  lastU=this->LastU;
  lastV=this->LastV;
  Complex TmpC=lastU*move_u+lastV*Conj(move_v);
  this->LastV=-lastU*move_v+lastV*move_u;
  double phase = -Arg(TmpC);
  Complex Phase = Polar(phase);
  TmpC *= Phase;
  this->LastU = TmpC.Re;
  this->LastV *= Phase; 
  this->LastTheta = (2.0*acos(LastU));
  this->LastPhi = Arg(LastV);
  
}


// restore last move
void ParticleOnSphereCollectionGauged::RestoreMove()
{
  if (LastMoved>-1)
    {      
      SpinorUCoordinates[LastMoved].Re = this->LastU;
      SpinorVCoordinates[LastMoved] = this->LastV;
      this->ThetaPhi[LastMoved<<1] = this->LastTheta;
      this->ThetaPhi[(LastMoved<<1)+1] = this->LastPhi;
      this->LastMoved = -1;
    }
}


void ParticleOnSphereCollectionGauged::MultiplyStepLength(double multiplier)
{
  this->Theta0*=multiplier;
}

// set single position
void ParticleOnSphereCollectionGauged::SetPosition(int nbrParticle, double theta, double phi)
{
  double s=sin(theta/2.);
  double c = cos(theta/2.);
  this->SpinorUCoordinates[nbrParticle].Re =c;
  this->SpinorUCoordinates[nbrParticle].Im =0.0;
  this->SpinorVCoordinates[nbrParticle].Re = s*cos(phi);
  this->SpinorVCoordinates[nbrParticle].Im = s*sin(phi);  
}

double ParticleOnSphereCollectionGauged::GetRandomNumber()
{
  return Generator->GetRealRandomNumber();
}

// randomize particle positions
void ParticleOnSphereCollectionGauged::Randomize()
{
  double phi, theta,  c, s;
  for (int i=0; i<NbrParticles; ++i)
    {
      phi=Generator->GetRealRandomNumber()*2.0*M_PI;
      theta=acos(1.0-2.0*Generator->GetRealRandomNumber());
      s=sin(theta/2.); c = cos(theta/2.);
      this->SpinorUCoordinates[i].Re =c;
      this->SpinorUCoordinates[i].Im =0.0;      
      this->SpinorVCoordinates[i].Re = s*cos(phi);
      this->SpinorVCoordinates[i].Im = s*sin(phi);
      this->ThetaPhi[i<<1] = .0*acos(SpinorUCoordinates[i].Re);
      this->ThetaPhi[(i<<1)+1] = Arg(SpinorVCoordinates[i]);
    }
}

// get absolute values of all relative distances
// distances = matrix in which to return the distances
void ParticleOnSphereCollectionGauged::GetDistances(RealSymmetricMatrix &distances)
{
  if ((distances.GetNbrRow()!=NbrParticles)||(distances.GetNbrColumn()!=NbrParticles))
    distances.Resize(NbrParticles,NbrParticles);
  for (int i=0; i<NbrParticles; ++i)
    {
      distances(i,i)=0.0;
      for (int j=i+1; j<NbrParticles; ++j)
	{
	  distances(i,j)=(*this->Distances)(i,j)=Norm(this->SpinorUCoordinates[i].Re*this->SpinorVCoordinates[j]-this->SpinorUCoordinates[j].Re*this->SpinorVCoordinates[i]);
	}
    }
  this->DistancesUpToDate = false;
  return;
}


// get reference to internal matrix with particle distances
const RealSymmetricMatrix& ParticleOnSphereCollectionGauged::GetDistances()
{
  if (!DistancesUpToDate)
    {
      for (int i=0; i<NbrParticles; ++i)
	{
	  // this->Distances(i,i)=0.0;
	  for (int j=i+1; j<NbrParticles; ++j)
	    {
	      (*this->Distances)(i,j)=Norm(this->SpinorUCoordinates[i]*this->SpinorVCoordinates[j]-this->SpinorUCoordinates[j]*this->SpinorVCoordinates[i]);
	    }
	}
    }
  this->DistancesUpToDate = false;
  return *this->Distances;
}



// toggle positions of first N/2 particles with the remaining N/2 positions
//
void ParticleOnSphereCollectionGauged
::ToggleHalfHalf()
{
  if (NbrParticles&1) return;
  double TmpD;
  Complex TmpC;
  int NUp = this->NbrParticles/2;
  // exchange spin up and spin down
  for (int j = 0; j < NUp; ++j)
    {
      TmpD = ThetaPhi[j << 1];      
      ThetaPhi[j << 1] = ThetaPhi[(j+NUp) << 1];
      ThetaPhi[(j+NUp) << 1] = TmpD;
      TmpD = ThetaPhi[1 + (j << 1)];
      ThetaPhi[1+(j <<1)] = ThetaPhi[1+ ((j+NUp) << 1)];
      ThetaPhi[1+ ((j+NUp) << 1)] = TmpD;
      TmpD = SpinorUCoordinates[j].Re;
      SpinorUCoordinates[j].Re = SpinorUCoordinates[j+NUp].Re;
      SpinorUCoordinates[j+NUp].Re = TmpD;
      TmpC = SpinorVCoordinates[j];
      SpinorVCoordinates[j] = SpinorVCoordinates[j+NUp];
      SpinorVCoordinates[j+NUp] = TmpC;
    }
}
