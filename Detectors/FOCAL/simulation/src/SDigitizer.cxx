#include "FOCALSimulation/SDigitizer.h"
#include "DataFormatsFOCAL/Digit.h"
#include "FOCALBase/Geometry.h"
#include "FOCALBase/Hit.h"
#include "MathUtils/Cartesian.h"
#include "SimulationDataFormat/MCCompLabel.h"

#include <climits>
#include <list>
#include <chrono>
#include <numeric>
#include <fairlogger/Logger.h>

ClassImp(o2::focal::SDigitizer);

using o2::focal::Digit;
using o2::focal::Hit;

using namespace o2::focal;
______________________________________________________________________
std::vector<o2::focal::LabeledDigit> SDigitizer::process(const std::vector<Hit>& hits)
{

  std::map<int, std::map<int, std::vector<o2::focal::Hit>>> hitsPerTowerPerParticleID;

  std::unordered_map<int, std::vector<LabeledDigit>> digitsPerTower;

  for (auto hit : hits) {
    hitsPerTowerPerParticleID[hit.GetDetectorID()][hit.GetTrackID()].push_back(hit);
  }

  std::vector<o2::focal::Hit> SHits;
  for (auto [towerID, hitsParticle] : hitsPerTowerPerParticleID) {
    for (auto [partID, Hits] : hitsParticle) {
      o2::focal::Hit SHit = std::accumulate(std::next(Hits.begin()), Hits.end(), Hits.front());
      SHits.push_back(SHit);
    }
  }

  for (auto hit : SHits) {
    try {

      int tower = hit.GetDetectorID();
      auto[inside, col, row, layer, segment] = mGeometry-> getVirtualInfo(hit.x(), hit.y(), hit.z());

      if (!inside) {
        LOG(warning) << "tower index out of range: " << tower;
        continue;
      }

      double energy = hit.GetEnergyLoss();

      Digit digit(col, row, layer, hit.GetTime(), energy);
      digit->setIndex(hit.DetID());

      MCLabel label(hit.GetTrackID(), mCurrEvID, mCurrSrcID, false, 1.0);
      if (digit.getAmplitude() < __DBL_EPSILON__) {
        label.setAmplitudeFraction(0);
      }
      LabeledDigit d(digit, label);

      digitsPerTower[tower].push_back(d);

    } catch (InvalidPositionException& e) {
      LOG(error) << "Error in creating the digit: " << e.what();
    }
  }

  std::vector<LabeledDigit> digitsVector;

  for (auto [towerID, labeledDigits] : digitsPerTower) {

    o2::focal::LabeledDigit Sdigit = std::accumulate(std::next(labeledDigits.begin()), labeledDigits.end(), labeledDigits.front());

    if (Sdigit.getAmplitude() < __DBL_EPSILON__) {
      continue;
    }

    digitsVector.push_back(Sdigit);
  }

  digitsPerTower.clear();

  return digitsVector;
}
_____________________________________________________________________
void SDigitizer::setCurrSrcID(int v)
{
  
  if (v > MCCompLabel::maxSourceID()) {
    LOG(fatal) << "MC source id " << v << " exceeds max storable in the label " << MCCompLabel::maxSourceID();
  }
  mCurrSrcID = v;
}
____________________________________________________________________
void SDigitizer::setCurrEvID(int v)
{
  
  if (v > MCCompLabel::maxEventID()) {
    LOG(fatal) << "MC event id " << v << " exceeds max storable in the label " << MCCompLabel::maxEventID();
  }
  mCurrEvID = v;
}
