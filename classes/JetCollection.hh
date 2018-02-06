#ifndef JETCOLLECTION_H
#define JETCOLLECTION_H

#include "Jet.hh"
#include "TLorentzVector.h"

class JetCollection{

   private:
       std::vector<Jet*> _jets;

   public:

       JetCollection();
       JetCollection(JetCollection& coll);

       ~JetCollection();
       
       Jet* AddJet(const TLorentzVector p4);
       void Add(Jet* r);
       
       Jet* at(const unsigned int i);       

       unsigned int size();
       void Delete(const unsigned int i);

};

#endif
