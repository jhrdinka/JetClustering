#ifndef GENPARTICLECOLLECTION_H
#define GENPARTICLECOLLECTION_H

#include "GenParticle.hh"
#include "TLorentzVector.h"

class GenParticleCollection {
 private:
  std::vector<const GenParticle*> _genparticles;

 public:
  GenParticleCollection();
  GenParticleCollection(const GenParticleCollection& coll);

  ~GenParticleCollection();

  const GenParticle* AddGenParticle(const TLorentzVector& p4, int pdgid,
                                    unsigned status, const TVector3& vertex,
                                    float charge, float deltaR, bool isPU);
  void Add(const GenParticle* r);

  const GenParticle* at(const unsigned int i);

  unsigned int size() const;
  void Delete(const unsigned int i);
};

#endif
