#ifndef CLUSTERCOLLECTION_H
#define CLUSTERCOLLECTION_H

#include "Cluster.hh"
#include "TLorentzVector.h"

class ClusterCollection{

   private:
       std::vector<Cluster*> _clusters;

   public:

       ClusterCollection();
       ClusterCollection(ClusterCollection& coll);

       ~ClusterCollection();
       
       Cluster* AddCluster(const TLorentzVector p4, const TLorentzVector pos);
       void Add(Cluster* r);
       
       Cluster* at(const unsigned int i);       

       unsigned int size();
       void Delete(const unsigned int i);

};

#endif
