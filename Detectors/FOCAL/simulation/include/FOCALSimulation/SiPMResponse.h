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
#ifndef ALICEO2_FOCAL_SIPMRESPONSE_H_
#define ALICEO2_FOCAL_SIPMRESPONSE_H_

#include "Rtypes.h"

class TRandom;

namespace o2
{
namespace focal
{

/// \struct SiPMSignal
/// \brief Response of the SiPM of one tower for one collision (all quantities kept for the debug stream)
struct SiPMSignal {
  double mMeanPhotoelectrons = 0.; // mu, after board crosstalk
  double mPhotoelectrons = 0.;     // N, after Poisson, dark counts and SiPM crosstalk
  double mFiredPixels = 0.;        // N_fired, after pixel saturation (mean value, no fluctuations in v1)
  double mCharge = 0.;             // Q at the preamplifier input, after the current conveyor (pC)
  double mTime = 0.;               // arrival time of the signal (ns), passed to the H2GCROC response
};

/// \class SiPMResponse
/// \brief FOCAL-H SiPM response of one tower: visible energy -> charge
/// \ingroup FOCALsimulation
///
/// Split in two steps because the board crosstalk acts on the mu of all towers of a collision at the same
/// time, and is therefore applied in HCALDigitizer between them:
///   1. getMeanPhotoelectrons(): E_vis -> mu
///      (board crosstalk on mu, in HCALDigitizer)
///   2. getResponse(): mu -> N (Poisson + dark counts) -> SiPM crosstalk -> pixel saturation -> charge
/// getVisibleEnergy() is the inverse (charge -> E_vis), for the analysis.
/// All parameters from HCALSimParam.
class SiPMResponse
{
 public:
  SiPMResponse() = default;
  ~SiPMResponse() = default;

  /// \brief Read the parameters from HCALSimParam and precompute the constants - tp be called before first collision.
  void init();

  /// \brief Random generator, owned by the caller (HCALDigitizer), so that one seed reproduces a collision
  void setRandomGenerator(TRandom* generator) { mRandomGenerator = generator; }

  /// \brief Return the mean number of photoelectrons (mu) of the tower
  /// \param visibleEnergy Visible energy of the tower (eV, as in the SDigit)
  /// \return mu
  double getMeanPhotoelectrons(double visibleEnergy) const;

  /// \brief Steps 3-6: photoelectron statistics, dark counts, SiPM crosstalk, pixel saturation, charge
  /// \param mu mu of the tower, after board crosstalk
  /// \param time Arrival time of the tower (ns); v1: passed through unchanged
  /// \return SiPM signal of the tower
  SiPMSignal getResponse(double mu, double time) const;

  /// \brief Charge -> visible energy
  /// \param Q Charge at the preamplifier input (pC)
  /// \return Visible energy (MeV)
  double getVisibleEnergyMeV(double Q) const;

  static constexpr double kElectronCharge = 1.602176634e-7; // one electron (pC)

 private:
  TRandom* mRandomGenerator = nullptr; //! random generator (not owned)
  double mMeanDarkCounts = 0.;         // mu_dark = mDarkCountRate * mIntegrationWindow
  double mChargePerFiredPixel = 0.;    // mSiPMGain * kElectronCharge (pC)

  // Parameters from HCALSimParam (read in init())
  double mPhotoelectronsPerMeV = 0.; // mean p.e. per MeV of visible energy
  double mSiPMCrosstalk = 0.;        // optical crosstalk probability
  int mSiPMPixels = 0;               // number of microcells
  double mConveyorGain = 0.;         // current-conveyor attenuation (0.

  ClassDefNV(SiPMResponse, 1);
};

} // namespace focal
} // namespace o2

#endif
