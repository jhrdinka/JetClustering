#include <iostream>
#include "GenParticleCollection.hh"


GenParticleCollection::GenParticleCollection(){
}

GenParticleCollection::GenParticleCollection(GenParticleCollection& coll) {
    _genparticles = std::vector<GenParticle*>(coll.size());
    for (unsigned int i=0; i<coll.size(); i++)
        _genparticles[i] = new GenParticle(*(coll._genparticles[i]));
}

GenParticleCollection::~GenParticleCollection(){
     size_t sz = _genparticles.size();
     for (size_t i = 0; i < sz; ++i)
         delete _genparticles[i];
}

GenParticle* GenParticleCollection::AddGenParticle(TLorentzVector p4, int pdgid, int status){
    GenParticle* r = new GenParticle(p4, pdgid, status);
    _genparticles.push_back(r);
    return r;
}

void GenParticleCollection::Add(GenParticle* r){
    _genparticles.push_back(r);
}

GenParticle* GenParticleCollection::at(const unsigned int i){
    return _genparticles[i];
}

unsigned int GenParticleCollection::size(){
    return _genparticles.size();
}

void GenParticleCollection::Delete(const unsigned int i){
    std::cout<<i<<std::endl;
    _genparticles.erase(_genparticles.begin() + i);
}
