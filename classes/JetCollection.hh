#ifndef JETCOLLECTION_H
#define JETCOLLECTION_H

#include "Jet.hh"
#include "TLorentzVector.h"

class JetCollection {
 private:
  std::vector<Jet*> _jets;

 public:
  JetCollection();
  JetCollection(const JetCollection& coll);

  ~JetCollection();

  Jet* AddJet(const TLorentzVector& p4, float vertexZ = 0);
  void Add(Jet* r);

  const Jet* at(const unsigned int i) const;

  unsigned int size() const;
  void Delete(const unsigned int i);
};

#endif
