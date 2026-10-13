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
#ifndef ALICEO2_FOCAL_MCLABEL_H_
#define ALICEO2_FOCAL_MCLABEL_H_

#include "SimulationDataFormat/MCCompLabel.h"

namespace o2
{
namespace focal
{

/// \class MCLabel
/// \brief Monte-Carlo label for FOCAL clusters / digits
/// \ingroup FOCALDataFormat
class MCLabel : public o2::MCCompLabel
{
 private:
  double mAmplitudeFraction = 0.; ///< Fraction of the digit amplitude (energy) carried by this MC particle

 public:
  MCLabel() = default;
  MCLabel(int trackID, int eventID, int srcID, bool fake, double afraction) : o2::MCCompLabel(trackID, eventID, srcID, fake), mAmplitudeFraction(afraction) {}
  MCLabel(bool noise, double afraction) : o2::MCCompLabel(noise), mAmplitudeFraction(afraction) {}
  
  void setAmplitudeFraction(double afraction) { mAmplitudeFraction = afraction; }
  double getAmplitudeFraction() const { return mAmplitudeFraction; }

  ClassDefNV(MCLabel, 1);
};
} // namespace focal
} //namespace o2

#endif
