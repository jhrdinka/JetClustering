#include "JetCollection.hh"
#include <iostream>

JetCollection::JetCollection() {}

JetCollection::JetCollection(const JetCollection& coll) {
  _jets = std::vector<Jet*>(coll.size());
  for (unsigned int i = 0; i < coll.size(); i++)
    _jets[i] = new Jet(*(coll._jets[i]));
}

JetCollection::~JetCollection() {
  size_t sz = _jets.size();
  for (size_t i = 0; i < sz; ++i) delete _jets[i];
}

Jet* JetCollection::AddJet(const TLorentzVector& p4, float vertexZ) {
  Jet* r = new Jet(p4, vertexZ);
  _jets.push_back(r);
  return r;
}

void JetCollection::Add(Jet* r) { _jets.push_back(r); }

Jet* JetCollection::at(const unsigned int i) { return _jets[i]; }

unsigned int JetCollection::size() const { return _jets.size(); }

void JetCollection::Delete(const unsigned int i) {
  std::cout << i << std::endl;
  _jets.erase(_jets.begin() + i);
}
