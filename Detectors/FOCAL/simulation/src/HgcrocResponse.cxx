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

#include "FOCALSimulation/HgcrocResponse.h"
#include "FOCALSimulation/PadSimParam.h"
#include <fairlogger/Logger.h>
#include <TSystem.h>
#include <TString.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>

namespace o2
{
namespace focal
{

void HgcrocResponse::init()
{
  mInitialized = false;
  mCalibParams = {};

  const PadSimParam& pars = PadSimParam::Instance();

  TString lowRangeCalibData = pars.getCalibFileLowRange();
  gSystem->ExpandPathName(lowRangeCalibData);
  TString highRangeCalibData = pars.getCalibFileHighRange();
  gSystem->ExpandPathName(highRangeCalibData);

  loadCalibLowRange(lowRangeCalibData.Data());
  loadCalibHighRange(highRangeCalibData.Data());
  validateCalibParams();
  regularizeBadCalibParams();

  mInitialized = true;
}

void HgcrocResponse::loadCalibLowRange(const std::string& filepath)
{
  std::ifstream in(filepath);
  std::string line;

  if (!in) {
    LOG(fatal) << "focal::HgcrocResponse: cannot read calibration file: " << filepath;
  }

  LOG(debug) << "focal::HgcrocResponse: reading calibration file: " << filepath;

  int asic, ch;
  double adcLowerFitRange, adcHigherFitRange, adcConst, adcSlope, totLowerFitRange, totHigherFitRange, totConst, totSlope;
  int lineNum = 0;

  while (std::getline(in, line)) {
    ++lineNum;
    if (line.empty()) {
      continue;
    }

    std::istringstream ss(line);

    if (!(ss >> asic >> ch >> adcLowerFitRange >> adcHigherFitRange >> adcConst >> adcSlope >> totLowerFitRange >> totHigherFitRange >> totConst >> totSlope)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":"
                 << lineNum << ": error parsing line";
    }
    if ((asic < 0) || (asic > NASICS - 1) || (ch < 0) || (ch > NCHANNELS - 1)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":"
                 << lineNum << ": out of range entry for ASIC " << asic << " channel " << ch;
    }
    if (mCalibParams[asic][ch].seenLow) {
      LOG(fatal) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch
                 << " has received more than one calibration entry from " << std::filesystem::path(filepath).filename().string();
    }

    mCalibParams[asic][ch].adcLowerFitRange = adcLowerFitRange;
    mCalibParams[asic][ch].adcHigherFitRange = adcHigherFitRange;
    mCalibParams[asic][ch].adcConst = adcConst;
    mCalibParams[asic][ch].adcSlope = adcSlope;
    mCalibParams[asic][ch].totLowerFitRange = totLowerFitRange;
    mCalibParams[asic][ch].totHigherFitRange = totHigherFitRange;
    mCalibParams[asic][ch].totConst = totConst;
    mCalibParams[asic][ch].totSlope = totSlope;
    mCalibParams[asic][ch].seenLow = true;
  }
}

void HgcrocResponse::loadCalibHighRange(const std::string& filepath)
{
  std::ifstream in(filepath);
  std::string line;

  if (!in) {
    LOG(fatal) << "focal::HgcrocResponse: cannot read calibration file: " << filepath;
  }

  LOG(debug) << "focal::HgcrocResponse: reading calibration file: " << filepath;

  int asic, ch;
  double totConstHigh, totSlopeHigh;
  int lineNum = 0;

  while (std::getline(in, line)) {
    ++lineNum;
    if (line.empty()) {
      continue;
    }

    std::istringstream ss(line);

    if (!(ss >> asic >> ch >> totSlopeHigh >> totConstHigh)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":"
                 << lineNum << ": error parsing line";
    }
    if ((asic < 0) || (asic > NASICS - 1) || (ch < 0) || (ch > NCHANNELS - 1)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":"
                 << lineNum << ": out of range entry for ASIC " << asic << " channel " << ch;
    }
    if (mCalibParams[asic][ch].seenHigh) {
      LOG(fatal) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch
                 << " has received more than one calibration entry from " << std::filesystem::path(filepath).filename().string();
    }

    mCalibParams[asic][ch].totConstHigh = totConstHigh;
    mCalibParams[asic][ch].totSlopeHigh = totSlopeHigh;
    mCalibParams[asic][ch].seenHigh = true;
  }
}

// validates whether calibration parameters are valid or nonsense and sets badMask accordingly
void HgcrocResponse::validateCalibParams()
{
  for (int asic = 0; asic < NASICS; ++asic) {
    for (int ch = 0; ch < NCHANNELS; ++ch) {
      auto& c = mCalibParams[asic][ch];
      c.badMask = 0;

      if (!c.seenLow) {
        c.badMask |= BadADCFit | BadTOTLowFit;
        LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch
                     << " did not receive calibration parameters for low charge injection fit, half mean will be used";
      } else {
        const double adcMargin = ADCMAX - c.adcConst - c.adcSlope * c.adcHigherFitRange;
        if (c.adcSlope < ADCSLOPE_MIN || c.adcSlope > ADCSLOPE_MAX || adcMargin <= 0) {
          c.badMask |= BadADCFit;
          LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch
                       << ": bad ADC fit (slope = " << c.adcSlope << ", adcMargin = " << adcMargin << "), will be replaced by the half mean";
        }
        if (c.totSlope < TOTSLOPE_MIN || c.totSlope > TOTSLOPE_MAX) {
          c.badMask |= BadTOTLowFit;
          LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch
                       << ": bad ToT low charge injection fit parameters (slope = " << c.totSlope << "), will be replaced by the half mean";
        }
      }

      if (!c.seenHigh) {
        c.badMask |= BadTOTHighFit;
        LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch
                     << " did not receive calibration parameters for high charge injection fit, half mean will be used";
      } else {
        if (c.totSlopeHigh < TOTSLOPEHIGH_MIN || c.totSlopeHigh > TOTSLOPEHIGH_MAX) {
          c.badMask |= BadTOTHighFit;
          LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch
                       << ": bad ToT high charge injection fit parameters (slope = " << c.totSlopeHigh << "), will be replaced by the half mean";
        }
      }
    }
  }
}

// replace the flagged groups of bad fit parameters by the mean of the good channels in the same ASIC half
void HgcrocResponse::regularizeBadCalibParams()
{
  using Param = double ChannelCalib::*;

  const std::vector<Param> adcParams{&ChannelCalib::adcLowerFitRange, &ChannelCalib::adcHigherFitRange, &ChannelCalib::adcConst, &ChannelCalib::adcSlope};
  const std::vector<Param> totLowParams{&ChannelCalib::totLowerFitRange, &ChannelCalib::totHigherFitRange, &ChannelCalib::totConst, &ChannelCalib::totSlope};
  const std::vector<Param> totHighParams{&ChannelCalib::totConstHigh, &ChannelCalib::totSlopeHigh};

  auto regularizeGroup = [this](int asic, int half, uint16_t bit, const std::vector<Param>& params) {
    const int firstChannel = half * NCHANNELS / 2;
    const int lastChannel = firstChannel + NCHANNELS / 2;

    int nGood = 0;
    std::vector<double> meanParams(params.size(), 0.);
    std::vector<int> badChannels;

    // calculating means over the good channels
    for (int ch = firstChannel; ch < lastChannel; ++ch) {
      const auto& c = mCalibParams[asic][ch];

      if (c.badMask & bit) {
        badChannels.push_back(ch);
        continue;
      }
      for (size_t i = 0; i < params.size(); ++i) {
        meanParams[i] += c.*params[i];
      }

      nGood++;
    }

    if (badChannels.empty()) {
      return; // nothing to replace
    }
    if (nGood == 0) {
      std::string group = "unknown fit group";
      if (bit == BadADCFit) {
        group = "ADC fit";
      } else if (bit == BadTOTLowFit) {
        group = "ToT low charge injection fit";
      } else if (bit == BadTOTHighFit) {
        group = "ToT high charge injection fit";
      }
      LOG(fatal) << "focal::HgcrocResponse: ASIC " << asic << " half " << half << " has no good channels to average for the " << group;
    }

    for (double& mean : meanParams) {
      mean /= nGood;
    }

    // replacing bad parameter with the mean for the ASIC half
    for (int ch : badChannels) {
      auto& c = mCalibParams[asic][ch];
      for (size_t i = 0; i < params.size(); ++i) {
        c.*params[i] = meanParams[i];
      }
      if (bit == BadADCFit) { // the half mean does not guarantee adcMargin > 0, so check it again after the substitution
        const double adcMargin = ADCMAX - c.adcConst - c.adcSlope * c.adcHigherFitRange;
        if (adcMargin <= 0) {
          LOG(fatal) << "focal::HgcrocResponse: the ADC fit parameters of ASIC " << asic << " half " << half
                     << " leave no ADC counts above the linear range (adcMargin = "
                     << adcMargin << " <= 0) after substitution of invalid fit parameters by half mean";
        }
      }
    }
  };

  for (int asic = 0; asic < NASICS; ++asic) {
    for (int half = 0; half < 2; ++half) {
      regularizeGroup(asic, half, BadADCFit, adcParams);
      regularizeGroup(asic, half, BadTOTLowFit, totLowParams);
      regularizeGroup(asic, half, BadTOTHighFit, totHighParams);
    }
  }
}

} // namespace focal
} // namespace o2