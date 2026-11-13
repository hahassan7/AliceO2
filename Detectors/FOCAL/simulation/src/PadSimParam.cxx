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

#include "FOCALSimulation/PadSimParam.h"
#include <iostream>

O2ParamImpl(o2::focal::PadSimParam);

namespace o2::focal
{

std::ostream& operator<<(std::ostream& stream, const PadSimParam& s)
{
  s.PrintStream(stream);
  return stream;
}

void PadSimParam::PrintStream(std::ostream& stream) const
{
  stream << "FOCAL pad simulation parameters:"
         << "\n FOCAL::PadSimParam::mEnergyPerPair = " << mEnergyPerPair
         << "\n FOCAL::PadSimParam::mCCE = " << mCCE
         << "\n FOCAL::PadSimParam::mCalibFileLowRange = " << mCalibFileLowRange
         << "\n FOCAL::PadSimParam::mCalibFileHighRange = " << mCalibFileHighRange
         << "\n FOCAL::PadSimParam::mADCTOTGap = " << mADCTOTGap
         << "\n FOCAL::PadSimParam::mTOTTransitionLow = " << mTOTTransitionLow
         << "\n FOCAL::PadSimParam::mTOTTransitionHigh = " << mTOTTransitionHigh
         << "\n FOCAL::PadSimParam::mNoiseADC = " << mNoiseADC
         << "\n FOCAL::PadSimParam::mCMNADC = " << mCMNADC
         << "\n FOCAL::PadSimParam::mTOTSmear = " << mTOTSmear
         << "\n FOCAL::PadSimParam::mTOAThreshold = " << mTOAThreshold
         << "\n FOCAL::PadSimParam::mReadoutAllSensors = " << (mReadoutAllSensors ? "yes" : "no")
         << "\n FOCAL::PadSimParam::mUseBadChannels = " << (mUseBadChannels ? "yes" : "no")
         << "\n FOCAL::PadSimParam::mDebugStream = " << (mDebugStream ? "yes" : "no") << "\n";
}

} // namespace o2::focal