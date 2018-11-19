#include "RecHitCollection.hh"
#include <iostream>

RecHitCollection::RecHitCollection() {}

RecHitCollection::RecHitCollection(RecHitCollection& coll) {
  _rechits = std::vector<RecHit*>(coll.size());
  for (unsigned int i = 0; i < coll.size(); i++)
    _rechits[i] = new RecHit(*(coll._rechits[i]));
}

RecHitCollection::~RecHitCollection() {
  size_t sz = _rechits.size();
  for (size_t i = 0; i < sz; ++i) delete _rechits[i];
}

RecHit* RecHitCollection::AddRecHit(const TLorentzVector p4,
                                    const TLorentzVector pos, int layer,
                                    int bits, float time) {
  RecHit* r = new RecHit(p4, pos, layer, bits, time);
  _rechits.push_back(r);
  return r;
}

void RecHitCollection::Add(RecHit* r) { _rechits.push_back(r); }

RecHit* RecHitCollection::at(const unsigned int i) { return _rechits[i]; }

unsigned int RecHitCollection::size() { return _rechits.size(); }

void RecHitCollection::Delete(const unsigned int i) {
  std::cout << i << std::endl;
  _rechits.erase(_rechits.begin() + i);
}
