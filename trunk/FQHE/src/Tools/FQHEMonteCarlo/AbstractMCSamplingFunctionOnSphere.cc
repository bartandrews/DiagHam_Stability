#include "AbstractMCSamplingFunctionOnSphere.h"
#include <iostream>
using std::cout;
using std::endl;

// virtual destructor
AbstractMCSamplingFunctionOnSphere::~AbstractMCSamplingFunctionOnSphere()
{
}

// register basic system of particles
AbstractParticleCollection * AbstractMCSamplingFunctionOnSphere::GetSystem()
{
  return ((AbstractParticleCollection *)this->System);
}

// register basic system of particles
void AbstractMCSamplingFunctionOnSphere::RegisterSystem(AbstractParticleCollection *system)
{
  cout << "Forwarding call in AbstractMCSamplingFunctionOnSphere::RegisterSystem(AbstractParticleCollection *system)"<<endl;
  this->RegisterSystem((AbstractParticleCollectionOnSphere *)system);
}
