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
#ifndef ALICEO2_FOCAL_PADSIMPARAM_H_
#define ALICEO2_FOCAL_PADSIMPARAM_H_

#include <iosfwd>
#include <string>
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "Rtypes.h"

namespace o2
{
namespace focal
{
/// \class PadSimParam
/// \brief FoCal-E Pads simulation parameters
/// \ingroup FOCALsimulation

class PadSimParam : public o2::conf::ConfigurableParamHelper<PadSimParam>
{
 public:
  ~PadSimParam() override = default;

  Float_t getEnergyPerPair() const { return mEnergyPerPair; }
  Float_t getCCE() const { return mCCE; }
  const std::string& getCalibFileLow() const { return mCalibFileLow; }
  const std::string& getCalibFileHigh() const { return mCalibFileHigh; }
  Float_t getADCTOTGap() const { return mADCTOTGap; }
  Float_t getTOTTransitionLow() const { return mTOTTransitionLow; }
  Float_t getTOTTransitionHigh() const { return mTOTTransitionHigh; }
  Float_t getNoiseADC() const { return mNoiseADC; }
  Float_t getCMNADC() const { return mCMNADC; }
  Float_t getTOTSmear() const { return mTOTSmear; }
  Float_t getTOAThreshold() const { return mTOAThreshold; }
  Bool_t getReadoutAllSensors() const { return mReadoutAllSensors; }
  Bool_t getUseBadChannels() const { return mUseBadChannels; }
  Bool_t getDebugStream() const { return mDebugStream; }

  Float_t mEnergyPerPair{3.65};
  Float_t mCCE{1.0};
  std::string mCalibFileLow{"share/Detectors/FOC/simulation/data/FOCALEPads_Calibration_20230208.txt"};
  std::string mCalibFileHigh{"share/Detectors/FOC/simulation/data/calib_highinjection_pdfextracted.txt"};
  Float_t mADCTOTGap{0.1};
  Float_t mTOTTransitionLow{0.5};
  Float_t mTOTTransitionHigh{2.0};
  Float_t mNoiseADC{1.5};
  Float_t mCMNADC{4.5};
  Float_t mTOTSmear{2.0};
  Float_t mTOAThreshold{0.002};
  Bool_t mReadoutAllSensors{false};
  Bool_t mUseBadChannels{false};
  Bool_t mDebugStream{false};

  void PrintStream(std::ostream& stream) const;

  O2ParamDef(PadSimParam, "FOCALPadSimParam");
};

std::ostream& operator<<(std::ostream& stream, const PadSimParam& s);

} // namespace focal
} // namespace o2

#endif