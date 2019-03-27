#ifndef TRACKCLUSTER_H
#define TRACKCLUSTER_H

#include "TVector3.h"

class TrackCluster {
 public:
  /// constructor
  TrackCluster(const TVector3& position, unsigned nTracks);
  /// copy constructor
  TrackCluster(const TrackCluster& trackCluster);

  /// desctructor
  ~TrackCluster() = default;

  /// access position
  TVector3 position() const;

  /// access number of contributing tracks
  unsigned short nTracks() const;

 private:
  /// the cluster position
  TVector3 _position;
  /// the number of tracks contributing to this cluster
  unsigned short _nTracks;
};

inline TVector3 TrackCluster::position() const { return _position; }

inline unsigned short TrackCluster::nTracks() const { return _nTracks; }

#endif  // TRACKCLUSTER_H
