#include <iostream>
#include "ClusterCollection.hh"


ClusterCollection::ClusterCollection(){
}

ClusterCollection::ClusterCollection(ClusterCollection& coll) {
    _clusters = std::vector<Cluster*>(coll.size());
    for (unsigned int i=0; i<coll.size(); i++)
        _clusters[i] = new Cluster(*(coll._clusters[i]));
}

ClusterCollection::~ClusterCollection(){
     size_t sz = _clusters.size();
     for (size_t i = 0; i < sz; ++i)
         delete _clusters[i];
}

Cluster* ClusterCollection::AddCluster(const TLorentzVector p4, const TLorentzVector pos){
    Cluster* r = new Cluster(p4, pos);
    _clusters.push_back(r);
    return r;
}

void ClusterCollection::Add(Cluster* r){
    _clusters.push_back(r);
}

Cluster* ClusterCollection::at(const unsigned int i){
    return _clusters[i];
}

unsigned int ClusterCollection::size(){
    return _clusters.size();
}

void ClusterCollection::Delete(const unsigned int i){
    std::cout<<i<<std::endl;
    _clusters.erase(_clusters.begin() + i);
}
