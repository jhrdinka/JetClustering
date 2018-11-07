#include "GenParticleCollection.hh"
#include <iostream>

GenParticleCollection::GenParticleCollection() {}

GenParticleCollection::GenParticleCollection(
    const GenParticleCollection& coll) {
  _genparticles = std::vector<const GenParticle*>(coll.size());
  for (unsigned int i = 0; i < coll.size(); i++)
    _genparticles[i] = new GenParticle(*(coll._genparticles[i]));
}

GenParticleCollection::~GenParticleCollection() {
  size_t sz = _genparticles.size();
  for (size_t i = 0; i < sz; ++i) delete _genparticles[i];
}

const GenParticle* GenParticleCollection::AddGenParticle(
    const TLorentzVector& p4, int pdgid, int status, const TVector3& vertex,
    float charge) {
  const GenParticle* r = new GenParticle(p4, pdgid, status, vertex, charge);
  _genparticles.push_back(r);
  return r;
}

void GenParticleCollection::Add(const GenParticle* r) {
  _genparticles.push_back(r);
}

const GenParticle* GenParticleCollection::at(const unsigned int i) {
  return _genparticles[i];
}

unsigned int GenParticleCollection::size() const {
  return _genparticles.size();
}

void GenParticleCollection::Delete(const unsigned int i) {
  std::cout << i << std::endl;
  _genparticles.erase(_genparticles.begin() + i);
}
