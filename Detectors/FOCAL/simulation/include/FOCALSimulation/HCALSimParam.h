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
#ifndef ALICEO2_FOCAL_HCALSIMPARAM_H_
#define ALICEO2_FOCAL_HCALSIMPARAM_H_

#include <iosfwd>
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "Rtypes.h"

namespace o2
{
namespace focal
{
/// \class HCALSimParam
/// \brief FOCAL-H (HCAL) simulation parameters: fibre, SiPM, H2GCROC
/// \ingroup FOCALsimulation

// Source of some defaults in the comments: "placeholder"s are not final, to be replaced by the HCAL settings/calibration (in progress)
// Many parameters are currently placeholders, put in to allow for testing of other classes as they are developed

class HCALSimParam : public o2::conf::ConfigurableParamHelper<HCALSimParam>
{
 public:
  ~HCALSimParam() override = default;

  void PrintStream(std::ostream& stream) const;

  // Fibre (applied at hit level in Detector::ProcessHitsHCAL)
  bool mApplyBirks{true};       // apply Birks' law at hit level
  float mBirksConstant{0.0126}; // Birks constant kB (cm/MeV) - polystyrene 0.126 mm/MeV
  float mAttenuationLength{0.}; // fibre attenuation length (cm), 0 = off - local sim uses 270 cm currently

  // SiPM (using datasheet for Hamamatsu S14160-6010PS)
  float mPhotoelectronsPerMeV{64.}; // mean p.e. per MeV of visible energy - placeholder
  int mSiPMPixels{359011};          // number of microcells - datasheet
  float mSiPMGain{1.8e5};           // electrons per fired pixel - datasheet at Vbr + 5 V, investigating
  float mSiPMCrosstalk{0.01};       // optical crosstalk probability - datasheet, < 1 % (upper bound), needs further investigation
  float mDarkCountRate{3.e6};       // dark count rate (Hz) - datasheet
  float mIntegrationWindow{25.};    // window in which dark counts add to the signal (ns)
  float mBoardCrosstalk{0.};        // fraction of mu leaking to each of the 4 neighbouring towers - currently placeholder 0, needs lab tests

  // Pulse shape, single collision
  float mPulseRiseTime{15.}; // SiPM pulse rise time (ns) - local sim, investiagion continuing
  float mPulseFallTime{80.}; // SiPM pulse fall time (ns) - local sim, investiagion continuing

  // H2GCROC
  float mConveyorGain{0.025};    // current-conveyor attenuation (0.025-0.375, placeholder from H2GCROC3b datasheet)
  float mPedestal{68.};          // pedestal (ADC) - placeholder from FOCAL-E
  float mADCSlope{4050.};        // ADC gain at the preamplifier input (ADC/pC) - placeholder from FOCAL-E
  float mADCLinearEnd{0.20};     // end of the linear ADC range Q_lin (pC) - placeholder from FOCAL-E
  float mADCTOTGap{0.100};       // gap between the end of the ADC range and the start of TOT (pC) - placeholder from FOCAL-E
  float mTOTConstLow{-16.};      // low-range TOT fit, constant - placeholder from FOCAL-E
  float mTOTSlopeLow{336.};      // low-range TOT fit, slope (per pC) - placeholder from FOCAL-E
  float mTOTConstHigh{152.};     // high-range TOT fit, constant - placeholder from FOCAL-E
  float mTOTSlopeHigh{257.};     // high-range TOT fit, slope (per pC) - placeholder from FOCAL-E
  float mTOTTransitionLow{0.5};  // end of the TOT low range (pC) - placeholder from FOCAL-E
  float mTOTTransitionHigh{2.0}; // start of the TOT high range (pC) - placeholder from FOCAL-E
  float mTOTSmear{2.};           // TOT resolution (TOT units) - placeholder
  float mNoiseADC{1.4};          // incoherent noise per channel (ADC) - placeholder, TDR (~1/730 of the ADC range)
  float mCMNADC{0.};             // common-mode noise per chip half (ADC) - placeholder
  float mTOAThreshold{0.002};    // TOA fires above this charge (pC) - placeholder

  // Output
  bool mDebugStream{false}; // write focalhdigitsDebug.root (one entry per digitized tower)

  O2ParamDef(HCALSimParam, "FOCALHSimParam");
};

std::ostream& operator<<(std::ostream& stream, const HCALSimParam& s);

} // namespace focal

namespace framework
{
template <typename T>
struct is_messageable;
template <>
struct is_messageable<o2::focal::HCALSimParam> : std::true_type {
};
} // namespace framework

} // namespace o2

#endif
