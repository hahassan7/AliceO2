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

#ifndef ALICEO2_FOCAL_HGCROCRESPONSE_H_
#define ALICEO2_FOCAL_HGCROCRESPONSE_H_

#include <array>
#include <string>
#include <cstdint>

namespace o2
{
namespace focal
{

class HgcrocResponse
{
 public:
  static constexpr int NASICS = 18; // test-beam 2023 setup
  static constexpr int NCHANNELS = 72;

  // container for the calibration parameters of each pad
  struct ChannelCalib {
    double adcConst = 0;          // ADC fit constant
    double adcSlope = 0;          // ADC fit slope
    double adcLowerFitRange = 0;  // lower bound (in pC) for ADC fit
    double adcHigherFitRange = 0; // upper bound (in pC) for ADC fit
    double totConst = 0;          // ToT fit constant for low range fit
    double totSlope = 0;          // ToT fit slope for low range fit
    double totLowerFitRange = 0;  // lower bound (in pC) for ToT fit
    double totHigherFitRange = 0; // upper bound (in pC) for ToT fit
    double totConstHigh = 0;      // ToT fit constant for high range fit
    double totSlopeHigh = 0;      // ToT fit slope for high range fit

    // derived parameters (computed from the above)
    double gain = 0;      // ADC counts per pC
    double qLin = 0;      // end of the linear ADC range
    double adcMargin = 0; // ADC counts left above the linear range
    double qSat = 0;      // charge at which ADC reaches ADCMAX
    double qTOT = 0;      // charge at which TOT starts

    bool seenLow = false, seenHigh = false; // calibration found in each file
    uint16_t badMask = 0;                   // which groups of fit parameters are bad (replaced by the half mean)
  };

  // bits used to set ChannelCalib::badMask, one per group of fit parameters
  enum BadCalibBits : uint16_t {
    BadADCFit = 1 << 0,
    BadTOTLowFit = 1 << 1,
    BadTOTHighFit = 1 << 2
  };

  void init();
  bool isInitialized() const { return mInitialized; }
  const std::array<std::array<ChannelCalib, NCHANNELS>, NASICS>& getCalibParams() const { return mCalibParams; }

 private:
  void loadCalibLowRange(const std::string& filepath);
  void loadCalibHighRange(const std::string& filepath);

  std::array<std::array<ChannelCalib, NCHANNELS>, NASICS> mCalibParams;
  bool mInitialized = false;
};

} // namespace focal
} // namespace o2

#endif