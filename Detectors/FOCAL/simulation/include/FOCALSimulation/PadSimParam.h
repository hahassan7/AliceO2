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

  float getEnergyPerPair() const { return mEnergyPerPair; }
  float getCCE() const { return mCCE; }
  const std::string& getCalibFileLowRange() const { return mCalibFileLowRange; }
  const std::string& getCalibFileHighRange() const { return mCalibFileHighRange; }
  float getADCTOTGap() const { return mADCTOTGap; }
  float getTOTTransitionLow() const { return mTOTTransitionLow; }
  float getTOTTransitionHigh() const { return mTOTTransitionHigh; }
  float getNoiseADC() const { return mNoiseADC; }
  float getCMNADC() const { return mCMNADC; }
  float getTOTSmear() const { return mTOTSmear; }
  float getTOAThreshold() const { return mTOAThreshold; }
  bool getReadoutAllSensors() const { return mReadoutAllSensors; }
  bool getUseBadChannels() const { return mUseBadChannels; }
  bool getDebugStream() const { return mDebugStream; }

  int getNAsics() const { return mNAsics; }
  int getNChannels() const { return mNChannels; }
  int getADCMax() const { return mADCMax; }

  float getADCSlopeMin() const { return mADCSlopeMin; }
  float getADCSlopeMax() const { return mADCSlopeMax; }
  float getTOTSlopeMin() const { return mTOTSlopeMin; }
  float getTOTSlopeMax() const { return mTOTSlopeMax; }
  float getTOTSlopeHighRangeMin() const { return mTOTSlopeHighRangeMin; }
  float getTOTSlopeHighRangeMax() const { return mTOTSlopeHighRangeMax; }

  float mEnergyPerPair{3.65f};
  float mCCE{1.0f};
  std::string mCalibFileLowRange{"${O2_ROOT}/share/Detectors/FOC/simulation/data/FOCALEPads_Calibration_20230208.txt"};
  std::string mCalibFileHighRange{"${O2_ROOT}/share/Detectors/FOC/simulation/data/calib_highinjection_pdfextracted.txt"};
  float mADCTOTGap{0.1f};
  float mTOTTransitionLow{0.5f};
  float mTOTTransitionHigh{2.0f};
  float mNoiseADC{1.5f};
  float mCMNADC{4.5f};
  float mTOTSmear{2.0f};
  float mTOAThreshold{0.002f};
  bool mReadoutAllSensors{false};
  bool mUseBadChannels{false};
  bool mDebugStream{false};

  int mNAsics{18}; // test-beam 2023 setup
  int mNChannels{72};
  int mADCMax{1023}; // 10-bit ADC

  // sanity windows for the calibration fits (from the test-beam analysis)
  float mADCSlopeMin{2500.f};
  float mADCSlopeMax{4500.f};
  float mTOTSlopeMin{200.f};
  float mTOTSlopeMax{550.f};
  float mTOTSlopeHighRangeMin{150.f};
  float mTOTSlopeHighRangeMax{300.f};

  void PrintStream(std::ostream& stream) const;

  O2ParamDef(PadSimParam, "FOCALPadSimParam");
};

std::ostream& operator<<(std::ostream& stream, const PadSimParam& s);

} // namespace focal
} // namespace o2

#endif
