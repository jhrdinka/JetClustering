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
  /// deltaR
  std::vector<float> _distancesRZ;
  /// mean
  float _meanDrz;
  /// min
  float _minDrz;
  /// max
  float _maxDrz;
  /// deltaPhi
  std::vector<float> _distancesS;
  /// mean
  float _meanDs;
  /// min
  float _minDs;
  /// max
  float _maxDs;
  /// The output tree
  TTree* _outputTree;

 public:
  // constructors
  HitAnalysis(const TString treeName);
  void fill(const int& layerID, const std::vector<float>& deltaRZ,
            const float& meanRZ, const float& minRZ, const float& maxRZ,
            const std::vector<float>& deltaS, const float& meanS,
            const float& minS, const float& maxS);
  void write();
};

#endif
