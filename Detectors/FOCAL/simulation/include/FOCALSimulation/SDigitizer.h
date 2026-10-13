// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#ifndef ALICEO2_FOCAL_SDIGITIZER_H
#define ALICEO2_FOCAL_SDIGITIZER_H

#include <memory>
#include <unordered_map>
#include <vector>
#include <list>

#include "Rtypes.h"
#include "TObject.h"

#include "DataFormatsFOCAL/SDigit.h"
#include "DataFormatsFOCAL/MCLabel.h"

#include "FOCALBase/Geometry.h"
#include "FOCALBase/Hit.h"
#include "FOCALSimulation/LabeledSDigit.h"

#include "SimulationDataFormat/MCTruthContainer.h"

namespace o2
{
namespace focal
{

/// \class SDigitizer
/// \brief FOCAL summed digitizer
/// \ingroup FOCALsimulation

class SDigitizer
{
 public:
  SDigitizer() = default;
  ~SDigitizer() = default;
  SDigitizer(const SDigitizer&) = delete;
  SDigitizer& operator=(const SDigitizer&) = delete;

  /// \brief Sum the hits of ONE MC event per cell
  /// \param hits Hits of the MC event
  /// \return Summed labeled digits of this MC event (also kept in mLabeledSDigits until the next call)
  ///
  /// The internal container is reset at the beginning of every call, so digits of different
  /// collisions / MC events are never double counted. The caller must collect the result
  /// (or getSDigits()/getMCLabels()) after each call.
  std::vector<o2::focal::LabeledSDigit> process(const std::vector<Hit>& hits);

  void setCurrSrcID(int v);
  int getCurrSrcID() const { return mCurrSrcID; }

  void setCurrEvID(int v);
  int getCurrEvID() const { return mCurrEvID; }

  void setGeometry(const o2::focal::Geometry* gm) { mGeometry = gm; }

  /// \brief Digits (without labels) of the last call to process()
  const std::vector<o2::focal::SDigit> getSDigits();
  /// \brief MC labels of the last call to process(), flattened over all digits
  const std::vector<o2::focal::MCLabel> getMCLabels();
  /// \brief Labeled digits of the last call to process()
  const std::vector<o2::focal::LabeledSDigit>& getLabeledSDigits() const { return mLabeledSDigits; }

 private:
  const Geometry* mGeometry = nullptr; 
  int mCurrSrcID = 0;                  
  int mCurrEvID = 0;
  std::vector<o2::focal::LabeledSDigit> mLabeledSDigits; ///< summed digits of the last call to process() (reset at every call)

  ClassDefNV(SDigitizer, 1);
};
} // namespace focal
} // namespace o2

#endif /* ALICEO2_FOCAL_SDIGITIZER_H */
