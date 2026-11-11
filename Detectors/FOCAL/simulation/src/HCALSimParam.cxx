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

#include "FOCALSimulation/HCALSimParam.h"
#include <iostream>

O2ParamImpl(o2::focal::HCALSimParam);

using namespace o2::focal;

std::ostream& o2::focal::operator<<(std::ostream& stream, const o2::focal::HCALSimParam& s)
{
  s.PrintStream(stream);
  return stream;
}

void HCALSimParam::PrintStream(std::ostream& stream) const
{
  stream << "FOCALHSimParam.mApplyBirks = " << ((mApplyBirks) ? "true" : "false");
  stream << "\nFOCALHSimParam.mBirksConstant = " << mBirksConstant;
  stream << "\nFOCALHSimParam.mAttenuationLength = " << mAttenuationLength;
  stream << "\nFOCALHSimParam.mPhotoelectronsPerMeV = " << mPhotoelectronsPerMeV;
  stream << "\nFOCALHSimParam.mSiPMPixels = " << mSiPMPixels;
  stream << "\nFOCALHSimParam.mSiPMGain = " << mSiPMGain;
  stream << "\nFOCALHSimParam.mSiPMCrosstalk = " << mSiPMCrosstalk;
  stream << "\nFOCALHSimParam.mDarkCountRate = " << mDarkCountRate;
  stream << "\nFOCALHSimParam.mIntegrationWindow = " << mIntegrationWindow;
  stream << "\nFOCALHSimParam.mBoardCrosstalk = " << mBoardCrosstalk;
  stream << "\nFOCALHSimParam.mPulseRiseTime = " << mPulseRiseTime;
  stream << "\nFOCALHSimParam.mPulseFallTime = " << mPulseFallTime;
  stream << "\nFOCALHSimParam.mConveyorGain = " << mConveyorGain;
  stream << "\nFOCALHSimParam.mPedestal = " << mPedestal;
  stream << "\nFOCALHSimParam.mADCSlope = " << mADCSlope;
  stream << "\nFOCALHSimParam.mADCLinearEnd = " << mADCLinearEnd;
  stream << "\nFOCALHSimParam.mADCTOTGap = " << mADCTOTGap;
  stream << "\nFOCALHSimParam.mTOTConstLow = " << mTOTConstLow;
  stream << "\nFOCALHSimParam.mTOTSlopeLow = " << mTOTSlopeLow;
  stream << "\nFOCALHSimParam.mTOTConstHigh = " << mTOTConstHigh;
  stream << "\nFOCALHSimParam.mTOTSlopeHigh = " << mTOTSlopeHigh;
  stream << "\nFOCALHSimParam.mTOTTransitionLow = " << mTOTTransitionLow;
  stream << "\nFOCALHSimParam.mTOTTransitionHigh = " << mTOTTransitionHigh;
  stream << "\nFOCALHSimParam.mTOTSmear = " << mTOTSmear;
  stream << "\nFOCALHSimParam.mNoiseADC = " << mNoiseADC;
  stream << "\nFOCALHSimParam.mCMNADC = " << mCMNADC;
  stream << "\nFOCALHSimParam.mTOAThreshold = " << mTOAThreshold;
  stream << "\nFOCALHSimParam.mDebugStream = " << ((mDebugStream) ? "true" : "false");
}
