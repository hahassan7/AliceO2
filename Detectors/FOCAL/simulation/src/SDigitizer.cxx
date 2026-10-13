#include "FOCALSimulation/SDigitizer.h"
#include "DataFormatsFOCAL/SDigit.h"
#include "FOCALBase/Geometry.h"
#include "FOCALBase/Hit.h"
#include "FOCALBase/VirtualCellIndex.h"
#include "MathUtils/Cartesian.h"
#include "SimulationDataFormat/MCCompLabel.h"

#include <climits>
#include <list>
#include <chrono>
#include <numeric>
#include <fairlogger/Logger.h>

ClassImp(o2::focal::SDigitizer);

// The digit subsystem is copied from the hit subsystem: both enums must use the same values
static_assert(static_cast<int>(o2::focal::Hit::Subsystem_t::EPADS) == o2::focal::SDigit::kEPads &&
                static_cast<int>(o2::focal::Hit::Subsystem_t::EPIXELS) == o2::focal::SDigit::kEPixels &&
                static_cast<int>(o2::focal::Hit::Subsystem_t::HCAL) == o2::focal::SDigit::kHCAL &&
                static_cast<int>(o2::focal::Hit::Subsystem_t::UNKNOWN) == o2::focal::SDigit::kUnknown,
              "o2::focal::SDigit::Subsystem_t and o2::focal::Hit::Subsystem_t must have the same values");

using o2::focal::SDigit;
using o2::focal::Hit;

using namespace o2::focal;

//_______________________________________________________________________
std::vector<o2::focal::LabeledSDigit> SDigitizer::process(const std::vector<Hit>& hits)
{
  // Reset the output of the previous call: every call handles ONE MC event, so digits of
  // different collisions / MC events must not be accumulated here (no double counting).
  mLabeledSDigits.clear();

  // Hits are grouped per cell. The cell key combines the subsystem and the detector ID of the hit, because
  // the detector IDs of the different subsystems are independent index spaces:
  //  - FOCAL-E pads   : detID = pad absID = layer*nRows*nCols + row*nCols + col on the grid of REAL pads
  //                     (FOCALBase/PadIndex.h, layer = geometric numbering 0..19, pads + pixels)
  //  - pixels and HCAL: detID = cell ID on the uniform virtual grid (FOCALBase/VirtualCellIndex.h)
  //                     TODO: to be replaced by a physical indexing together with the pixel/HCAL digitization
  auto cellKey = [](const o2::focal::Hit& hit) -> int64_t {
    return (static_cast<int64_t>(hit.getSubsystem()) << 32) | static_cast<uint32_t>(hit.GetDetectorID());
  };

  std::map<int64_t, std::map<int, std::vector<o2::focal::Hit>>> hitsPerTowerPerParticleID;

  std::unordered_map<int64_t, std::vector<LabeledSDigit>> digitsPerTower;

  for (auto hit : hits) {
    hitsPerTowerPerParticleID[cellKey(hit)][hit.GetTrackID()].push_back(hit);
  }

  const auto padIndex = mGeometry->getPadIndex();

  std::vector<o2::focal::Hit> SHits;
  for (auto [towerID, hitsParticle] : hitsPerTowerPerParticleID) {
    for (auto [partID, Hits] : hitsParticle) {
      o2::focal::Hit SHit = std::accumulate(std::next(Hits.begin()), Hits.end(), Hits.front());
      SHits.push_back(SHit);
    }
  }

  for (auto hit : SHits) {

    int tower = hit.GetDetectorID();
    // (layer,row,col): pads -> grid of real pads, pixels / HCAL -> virtual grid (unchanged behaviour)
    auto [layer, row, col] = hit.isPadHit() ? padIndex.getLayerRowColFromIndex(tower) : VirtualCellIndex::getLayerRowColFromIndex(tower);

    if (layer < 0) {
      LOG(warning) << "cell index out of range: " << tower << " (subsystem " << static_cast<int>(hit.getSubsystem()) << ")";
      continue;
    }

    double energy = hit.GetEnergyLoss();

    SDigit digit(col, row, layer, hit.GetTime(), energy); // energy in eV
    digit.setIndex(tower);                                 // unique cell index inside the subsystem (pads: absID, pixels/HCAL: cellID)
    digit.setSubsystem(static_cast<SDigit::Subsystem_t>(hit.getSubsystem()));

    MCLabel label(hit.GetTrackID(), mCurrEvID, mCurrSrcID, false, 1.0);
    if (digit.getEnergy() < __DBL_EPSILON__) {
      label.setAmplitudeFraction(0);
    }
    LabeledSDigit d(digit, label);

    digitsPerTower[cellKey(hit)].push_back(d);
  }

  for (auto [towerID, labeledDigits] : digitsPerTower) {

    o2::focal::LabeledSDigit Sdigit = std::accumulate(std::next(labeledDigits.begin()), labeledDigits.end(), labeledDigits.front());

    if (Sdigit.getEnergy() < __DBL_EPSILON__) {
      continue;
    }

    mLabeledSDigits.push_back(Sdigit);
  }

  digitsPerTower.clear();

  return mLabeledSDigits;
}

//_______________________________________________________________________
void SDigitizer::setCurrSrcID(int v)
{
  
  if (v > MCCompLabel::maxSourceID()) {
    LOG(fatal) << "MC source id " << v << " exceeds max storable in the label " << MCCompLabel::maxSourceID();
  }
  mCurrSrcID = v;
}

//_______________________________________________________________________
void SDigitizer::setCurrEvID(int v)
{
  
  if (v > MCCompLabel::maxEventID()) {
    LOG(fatal) << "MC event id " << v << " exceeds max storable in the label " << MCCompLabel::maxEventID();
  }
  mCurrEvID = v;
}

//_______________________________________________________________________
const std::vector<o2::focal::SDigit> SDigitizer::getSDigits()
{
  std::vector<o2::focal::SDigit> digits;
  for (const auto& labeledDigit : mLabeledSDigits) {
    digits.push_back(labeledDigit.getSDigit());
  }
  return digits;
}

//_______________________________________________________________________
const std::vector<o2::focal::MCLabel> SDigitizer::getMCLabels()
{
  std::vector<o2::focal::MCLabel> labels;
  for (const auto& ld : mLabeledSDigits) {
    const auto& l = ld.getLabels();
    labels.insert(labels.end(), l.begin(), l.end());
  }
  return labels;
}