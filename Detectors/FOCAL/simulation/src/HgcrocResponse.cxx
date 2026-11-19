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
#include <filesystem>

namespace o2
{
namespace focal
{

template <typename SimParam>
void HgcrocResponse::init(const SimParam& pars)
{
  mInitialized = false;

  // loading config from SimParam class for pads or HCAL
  mConfig.mNAsics = pars.mNAsics;
  mConfig.mNChannels = pars.mNChannels;
  mConfig.mADCMax = pars.mADCMax;
  mConfig.mADCSlopeMin = pars.mADCSlopeMin;
  mConfig.mADCSlopeMax = pars.mADCSlopeMax;
  mConfig.mTOTSlopeMin = pars.mTOTSlopeMin;
  mConfig.mTOTSlopeMax = pars.mTOTSlopeMax;
  mConfig.mTOTSlopeHighRangeMin = pars.mTOTSlopeHighRangeMin;
  mConfig.mTOTSlopeHighRangeMax = pars.mTOTSlopeHighRangeMax;
  mConfig.mCalibFileLowRange = pars.mCalibFileLowRange;
  mConfig.mCalibFileHighRange = pars.mCalibFileHighRange;
  mConfig.mADCTOTGap = pars.mADCTOTGap;

  if (mConfig.mNAsics <= 0 || mConfig.mNChannels <= 0) {
    LOG(fatal) << "focal::HgcrocResponse: number of ASICs and channels must be positive";
  }
  if (mConfig.mNChannels % 2 != 0) {
    LOG(fatal) << "focal::HgcrocResponse: number of channels per ASIC must be even (the ASIC is split in two halves)";
  }

  mCalibParams.assign(mConfig.mNAsics, std::vector<ChannelCalib>(mConfig.mNChannels));

  TString lowRangeCalibData = mConfig.mCalibFileLowRange;
  gSystem->ExpandPathName(lowRangeCalibData);
  TString highRangeCalibData = mConfig.mCalibFileHighRange;
  gSystem->ExpandPathName(highRangeCalibData);

  loadCalibLowRange(lowRangeCalibData.Data());
  loadCalibHighRange(highRangeCalibData.Data());
  validateCalibParams();
  regularizeBadCalibParams();
  computeDerivedParams();

  mInitialized = true;
}

template void HgcrocResponse::init<PadSimParam>(const PadSimParam&);
// add HCAL here once it provides the same getters

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
    if ((asic < 0) || (asic > mConfig.mNAsics - 1) || (ch < 0) || (ch > mConfig.mNChannels - 1)) {
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
    if ((asic < 0) || (asic > mConfig.mNAsics - 1) || (ch < 0) || (ch > mConfig.mNChannels - 1)) {
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

// validates whether calibration parameters are valid or nonsense and sets badFitMask accordingly
void HgcrocResponse::validateCalibParams()
{
  int nMissingLow = 0, nMissingHigh = 0, nBadADC = 0, nBadTOTLow = 0, nBadTOTHigh = 0;
  std::ostringstream logMissingLow, logMissingHigh, logBadADC, logBadTOTLow, logBadTOTHigh;

  for (int asic = 0; asic < mConfig.mNAsics; ++asic) {
    for (int ch = 0; ch < mConfig.mNChannels; ++ch) {
      auto& c = mCalibParams[asic][ch];
      c.badFitMask = 0;

      if (!c.seenLow) {
        c.badFitMask |= BadADCFit | BadTOTLowFit;
        nMissingLow++;
        logMissingLow << " (" << asic << ", " << ch << ")";
      } else {
        const double adcMargin = mConfig.mADCMax - c.adcConst - c.adcSlope * c.adcHigherFitRange;
        if (c.adcSlope < mConfig.mADCSlopeMin || c.adcSlope > mConfig.mADCSlopeMax || adcMargin <= 0) {
          c.badFitMask |= BadADCFit;
          nBadADC++;
          logBadADC << " (" << asic << ", " << ch << ")";
        }
        if (c.totSlope < mConfig.mTOTSlopeMin || c.totSlope > mConfig.mTOTSlopeMax) {
          c.badFitMask |= BadTOTLowFit;
          nBadTOTLow++;
          logBadTOTLow << " (" << asic << ", " << ch << ")";
        }
      }

      if (!c.seenHigh) {
        c.badFitMask |= BadTOTHighFit;
        nMissingHigh++;
        logMissingHigh << " (" << asic << ", " << ch << ")";
      } else {
        if (c.totSlopeHigh < mConfig.mTOTSlopeHighRangeMin || c.totSlopeHigh > mConfig.mTOTSlopeHighRangeMax) {
          c.badFitMask |= BadTOTHighFit;
          nBadTOTHigh++;
          logBadTOTHigh << " (" << asic << ", " << ch << ")";
        }
      }
    }
  }

  if (nMissingLow > 0) {
    LOG(warning) << "focal::HgcrocResponse: " << nMissingLow << " channels with no low charge injection calibration, "
                 << "replaced by the half mean (asic, ch):" << logMissingLow.str();
  }
  if (nMissingHigh > 0) {
    LOG(warning) << "focal::HgcrocResponse: " << nMissingHigh << " channels with no high charge injection calibration, "
                 << "replaced by the half mean (asic, ch):" << logMissingHigh.str();
  }
  if (nBadADC > 0) {
    LOG(warning) << "focal::HgcrocResponse: " << nBadADC << " channels with a bad ADC fit, "
                 << "replaced by the half mean (asic, ch):" << logBadADC.str();
  }
  if (nBadTOTLow > 0) {
    LOG(warning) << "focal::HgcrocResponse: " << nBadTOTLow << " channels with a bad low charge injection ToT fit, "
                 << "replaced by the half mean (asic, ch):" << logBadTOTLow.str();
  }
  if (nBadTOTHigh > 0) {
    LOG(warning) << "focal::HgcrocResponse: " << nBadTOTHigh << " channels with a bad high charge injection ToT fit, "
                 << "replaced by the half mean (asic, ch):" << logBadTOTHigh.str();
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
    const int firstChannel = half * mConfig.mNChannels / 2;
    const int lastChannel = firstChannel + mConfig.mNChannels / 2;

    int nGood = 0;
    std::vector<double> meanParams(params.size(), 0.);
    std::vector<int> badFitChannels;

    // calculating means over the good channels
    for (int ch = firstChannel; ch < lastChannel; ++ch) {
      const auto& c = mCalibParams[asic][ch];

      if (c.badFitMask & bit) {
        badFitChannels.push_back(ch);
        continue;
      }
      for (size_t i = 0; i < params.size(); ++i) {
        meanParams[i] += c.*params[i];
      }

      nGood++;
    }

    if (badFitChannels.empty()) {
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
    for (int ch : badFitChannels) {
      auto& c = mCalibParams[asic][ch];
      for (size_t i = 0; i < params.size(); ++i) {
        c.*params[i] = meanParams[i];
      }
      if (bit == BadADCFit) { // the half mean does not guarantee adcMargin > 0, so check it again after the substitution
        const double adcMargin = mConfig.mADCMax - c.adcConst - c.adcSlope * c.adcHigherFitRange;
        if (adcMargin <= 0) {
          LOG(fatal) << "focal::HgcrocResponse: the ADC fit parameters of ASIC " << asic << " half " << half
                     << " leave no ADC counts above the linear range (adcMargin = "
                     << adcMargin << " <= 0) after substitution of invalid fit parameters by half mean";
        }
      }
    }
  };

  for (int asic = 0; asic < mConfig.mNAsics; ++asic) {
    for (int half = 0; half < 2; ++half) {
      regularizeGroup(asic, half, BadADCFit, adcParams);
      regularizeGroup(asic, half, BadTOTLowFit, totLowParams);
      regularizeGroup(asic, half, BadTOTHighFit, totHighParams);
    }
  }
}

void HgcrocResponse::computeDerivedParams()
{
  for (auto& row : mCalibParams) {
    for (auto& c : row) {
      c.gain = c.adcSlope;
      c.qLin = c.adcHigherFitRange;
      c.adcMargin = mConfig.mADCMax - c.adcConst - c.gain * c.qLin;
      c.qSat = c.qLin + 2 * c.adcMargin / c.gain;
      c.qTOT = c.totLowerFitRange;
    }
  }
}

} // namespace focal
} // namespace o2
