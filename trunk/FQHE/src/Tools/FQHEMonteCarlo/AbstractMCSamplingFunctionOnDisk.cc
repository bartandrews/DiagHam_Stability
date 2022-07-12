#include "AbstractMCSamplingFunctionOnDisk.h"
#include <iostream>
using std::cout;
using std::endl;

// virtual destructor
AbstractMCSamplingFunctionOnDisk::~AbstractMCSamplingFunctionOnDisk()
{
}

// register basic system of particles
AbstractParticleCollection * AbstractMCSamplingFunctionOnDisk::GetSystem()
{
  return ((AbstractParticleCollection *)this->System);
}

// register basic system of particles
void AbstractMCSamplingFunctionOnDisk::RegisterSystem(AbstractParticleCollection *system)
{
  // cout << "Forwarding call in AbstractMCSamplingFunctionOnDisk::RegisterSystem(AbstractParticleCollection *system)"<<endl;
  if (system->GetCollectionType() != AbstractParticleCollection::OnDiskCollection)
    {
      std::cerr << "Error: wrong type of particle collection!"<<endl;
      exit(1);
    }
  this->RegisterSystem((AbstractParticleCollectionOnDisk *)system);
}
