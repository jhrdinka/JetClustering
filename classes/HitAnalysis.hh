#ifndef HITANALYSIS_H
#define HITANALYSIS_H

#include <TDirectory.h>
#include <TString.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>

using namespace std;

class HitAnalysis {
 private:
  TString _treeName;
  /// layerID
  int _layerID;
  /// Optional parameter (e.g. for indicating position in r, eta)
  float _r;
  /// Optional parameter (e.g. for indicating position in z, phi)
  float _z;
  /// deltaR
  std::vector<float> _deltaR;
  /// mean
  float _mean;
  /// min
  float _min;
  /// max
  float _max;
  /// The output tree
  TTree* _outputTree;

 public:
  // constructors
  HitAnalysis(const TString treeName);
  void fill(const int& layerID, const float& r, const float& z,
            const std::vector<float>& deltaR, float mean, float min, float max);
  void write();
};

#endif
