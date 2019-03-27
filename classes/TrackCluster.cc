#include "TrackCluster.hh"

TrackCluster::TrackCluster(const TVector3& position, unsigned nTracks)
    : _position(position), _nTracks(nTracks) {}

TrackCluster::TrackCluster(const TrackCluster& trackCluster) {
  _position = trackCluster._position;
  _nTracks = trackCluster._nTracks;
}
