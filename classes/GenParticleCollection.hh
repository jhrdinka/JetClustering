#ifndef GENPARTICLECOLLECTION_H
#define GENPARTICLECOLLECTION_H

#include "GenParticle.hh"
#include "TLorentzVector.h"

class GenParticleCollection{

   private:
       std::vector<GenParticle*> _genparticles;

   public:

       GenParticleCollection();
       GenParticleCollection(GenParticleCollection& coll);

       ~GenParticleCollection();
       
       GenParticle* AddGenParticle(TLorentzVector p4, int pdgid, int status);
       void Add(GenParticle* r);
       
       GenParticle* at(const unsigned int i);       

       unsigned int size();
       void Delete(const unsigned int i);

};

#endif
