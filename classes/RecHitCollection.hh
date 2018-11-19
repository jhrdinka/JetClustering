#ifndef RECHITCOLLECTION_H
#define RECHITCOLLECTION_H

#include "RecHit.hh"
#include "TLorentzVector.h"

class RecHitCollection {
 private:
  std::vector<RecHit*> _rechits;

 public:
  RecHitCollection();
  RecHitCollection(RecHitCollection& coll);

  ~RecHitCollection();

  RecHit* AddRecHit(const TLorentzVector p4, const TLorentzVector pos,
                    int layer, int bits, float time);
  void Add(RecHit* r);

  RecHit* at(const unsigned int i);

  unsigned int size() const;
  void Delete(const unsigned int i);
};

#endif
