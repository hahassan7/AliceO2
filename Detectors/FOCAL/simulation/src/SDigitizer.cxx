#include "FOCALSimulation/SDigitizer.h"
#include "DataFormatsFOCAL/SDigit.h"
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

  // Hits are grouped per cell. The cell key combines the subsystem and the full cell index, because the indices of
  // the different subsystems are independent index spaces:
  //  - FOCAL-E pads   : pad absID = layer*nRows*nCols + row*nCols + col on the grid of REAL pads
  //                     (FOCALBase/PadIndex.h, layer = geometric numbering 0..19, pads + pixels)
  //  - HCAL (sheets)  : SiPM index = row*nCols + col (FOCALBase/HCALIndex.h); layer = first HCAL layer (20), information only
  //  - HCAL (other designs): cell from the hit position on the uniform virtual grid (getVirtualInfo, old way),
  //                     index = (layer*nRows + row)*nCols + col
  //  - pixels         : SKIPPED. TODO: pixels need their own indexing (ITS approach: chip ID + position) and digitization.
  // The hit stores only (layer, index inside the layer), see Hit.h. This is the ONLY place where the full index
  // (pads: absID, HCAL: SiPM index) is rebuilt; everything after the SDigitizer uses the full index.
  struct Cell {
    int index = -1, layer = -1, row = -1, col = -1;
  };
  const bool hcalSheets = (mGeometry->getHCALDesign() == Geometry::HCALDesgin::Sheets);
  auto cellOf = [&](const o2::focal::Hit& hit) -> Cell {
    Cell c;
    if (hit.isHCALHit() && !hcalSheets) {
      auto [inside, col, row, layer, segment] = mGeometry->getVirtualInfo(hit.GetX(), hit.GetY(), hit.GetZ());
      auto [gridOK, nCols, nRows] = mGeometry->getVirtualNColRow(segment);
      if (inside && gridOK) {
        c = {(layer * nRows + row) * nCols + col, layer, row, col};
      }
      return c;
    }
    c.index = mGeometry->getIndexFromLayerAndIndexInLayer(hit.getSubsystem(), hit.getLayer(), hit.getIndexInLayer());
    std::tie(c.layer, c.row, c.col) = mGeometry->getLayerRowColFromIndex(hit.getSubsystem(), c.index);
    return c;
  };

  std::map<int64_t, std::map<int, std::vector<o2::focal::Hit>>> hitsPerTowerPerParticleID;
  std::unordered_map<int64_t, Cell> cells;
  std::unordered_map<int64_t, std::vector<LabeledSDigit>> digitsPerTower;

  for (const auto& hit : hits) {
    if (hit.isPixelHit()) {
      continue; // TODO: pixels (see above)
    }
    const Cell c = cellOf(hit);
    if (c.layer < 0) {
      LOG(warning) << "hit without valid cell: subsystem " << static_cast<int>(hit.getSubsystem()) << ", layer " << hit.getLayer() << ", index in layer " << hit.getIndexInLayer();
      continue;
    }
    const int64_t key = (static_cast<int64_t>(hit.getSubsystem()) << 32) | static_cast<uint32_t>(c.index);
    cells[key] = c;
    hitsPerTowerPerParticleID[key][hit.GetTrackID()].push_back(hit);
  }

  for (const auto& [key, hitsParticle] : hitsPerTowerPerParticleID) {
    const Cell& c = cells[key];
    for (const auto& [partID, Hits] : hitsParticle) {
      o2::focal::Hit hit = std::accumulate(std::next(Hits.begin()), Hits.end(), Hits.front()); // sum of the hits of one particle in this cell

      SDigit digit(c.col, c.row, c.layer, hit.GetTime(), hit.GetEnergyLoss()); // energy in eV
      digit.setIndex(c.index);                                                 // unique index inside the subsystem (pads: absID, HCAL: SiPM index)
      digit.setSubsystem(static_cast<SDigit::Subsystem_t>(hit.getSubsystem()));

      MCLabel label(hit.GetTrackID(), mCurrEvID, mCurrSrcID, false, 1.0);
      if (digit.getEnergy() < __DBL_EPSILON__) {
        label.setAmplitudeFraction(0);
      }
      digitsPerTower[key].emplace_back(digit, label);
    }
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