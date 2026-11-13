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

} // namespace focal
} // namespace o2