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

#include "FOCALSimulation/SiPMResponse.h"
#include "FOCALSimulation/HCALSimParam.h"

#include <TRandom.h>
#include <cmath>

ClassImp(o2::focal::SiPMResponse);

using namespace o2::focal;

void SiPMResponse::init()
{
  // Read the parameters from HCALSimParam and precompute the constants
  auto& param = HCALSimParam::Instance();
  mMeanDarkCounts = param.mDarkCountRate * param.mIntegrationWindow * 1.e-9;
  mChargePerFiredPixel = param.mSiPMGain * kElectronCharge;

  mPhotoelectronsPerMeV = param.mPhotoelectronsPerMeV;
  mSiPMCrosstalk = param.mSiPMCrosstalk;
  mSiPMPixels = param.mSiPMPixels;
  mConveyorGain = param.mConveyorGain;
}

double SiPMResponse::getMeanPhotoelectrons(double visibleEnergy) const
{
  // TDR: (light yield × capture × fibre transparency × interface) x SiPM PDE
  double meanPe = mPhotoelectronsPerMeV * visibleEnergy * 1.e-6; // MeV
  return meanPe;
}

SiPMSignal SiPMResponse::getResponse(double meanPe, double time) const
{
  SiPMSignal signal;

  // Use gRandom if no generator was set
  TRandom* randomGenerator = mRandomGenerator ? mRandomGenerator : gRandom;


  // Photoelectron statistics and dark counts
  double Npe = randomGenerator->Poisson(meanPe) + randomGenerator->Poisson(mMeanDarkCounts);

  // Optical crosstalk inside SiPM
  // Each avalanche triggers ~mSiPMCrosstalk extra pixels
  Npe = Npe + randomGenerator->Poisson(Npe * mSiPMCrosstalk);

  // Pixel saturation
  // Mean number of fired pixels when Npe photoelectrons land at random on Npix pixels
  // TODO: deal with total vs effective Npix issue, investigate with local sim
  double firedPixels = mSiPMPixels * (1 - std::exp(-Npe / mSiPMPixels));

  // Charge at preamp in pC
  double chargeSiPM = firedPixels * mChargePerFiredPixel; // charge at SiPM output
  double charge = mConveyorGain * chargeSiPM;             // charge at preamplifier input

  // Set signal values
  signal.mMeanPhotoelectrons = meanPe;
  signal.mTime = time;
  signal.mPhotoelectrons = Npe;
  signal.mFiredPixels = firedPixels;
  signal.mCharge = charge;

  return signal;
}

double SiPMResponse::getVisibleEnergyMeV(double charge) const
{
  double firedPixels = charge / (mChargePerFiredPixel * mConveyorGain); // charge at SiPM output
  // lower limit to avoid log(0) and negative Evis
  if (firedPixels >= mSiPMPixels) {
    firedPixels = mSiPMPixels - 1.e-6;
  }

  double Npe = -mSiPMPixels * std::log(1 - firedPixels / mSiPMPixels); // mean number of photoelectrons - lower limit if firedPixels >= Npix
  Npe = Npe / (1 + mSiPMCrosstalk);                                    // mean number of photoelectrons before crosstalk
  double Evis = Npe / mPhotoelectronsPerMeV;                           // visible energy (MeV)

  return Evis;
}
