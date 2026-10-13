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

#ifndef STEER_DIGITIZERWORKFLOW_FOCALDIGITIZER_H_
#define STEER_DIGITIZERWORKFLOW_FOCALDIGITIZER_H_

#include <vector>

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "DataFormatsFOCAL/SDigit.h"
#include "FOCALBase/Hit.h"
#include "FOCALSimulation/SDigitizer.h"
#include <DetectorsBase/BaseDPLDigitizer.h>

class TChain;

namespace o2
{

namespace focal
{

/// \class DigitizerSpec
/// \brief DPL task running the FOCAL (S)digitization
///
/// Current status (v0): for every collision the hits of all MC events (parts) of that collision are
/// summed per cell by the SDigitizer and written out as summed digits (energy per cell), with one
/// TriggerRecord per collision and the MC labels. Every collision is treated as triggered.
///
/// TODO:
/// - Electronics response for the pads (charge, ADC/TOA/TOT, noise) -> PadDigitizer / HGCROCResponse
/// - Trigger handling (e.g. FIT minimum-bias) and pile-up / readout window
/// - Pixel and HCAL digitization
class DigitizerSpec final : public o2::base::BaseDPLDigitizer, public o2::framework::Task
{
 public:
  using o2::base::BaseDPLDigitizer::init;
  /// \brief Constructor
  DigitizerSpec() : o2::base::BaseDPLDigitizer(o2::base::InitServices::GEOM), o2::framework::Task() {}

  /// \brief Destructor
  ~DigitizerSpec() final = default;

  /// \brief init digitizer
  /// \param ctx Init context
  void initDigitizerTask(framework::InitContext& ctx) final;

  /// \brief run digitizer
  /// \param ctx Processing context
  ///
  /// One loop over collisions: for each collision, sum the hits of all its MC events per cell
  /// and append the summed digits, their MC labels and a TriggerRecord to the output.
  void run(framework::ProcessingContext& ctx) override;

 private:
  bool mFinished = false;              ///< Flag for digitization finished
  o2::focal::SDigitizer mSumDigitizer; ///< Summed digitizer
  std::vector<Hit> mHits;              ///< Vector with input hits
  std::vector<TChain*> mSimChains;     ///< Chains with the simulated hits
};

/// \brief Create the DPL processor spec of the FOCAL digitizer
/// \param channel Channel of the collision context
/// \param mctruth Whether to produce the MC labels
o2::framework::DataProcessorSpec getFOCALDigitizerSpec(int channel, bool mctruth = true);

} // namespace focal
} // end namespace o2

#endif /* STEER_DIGITIZERWORKFLOW_FOCALDIGITIZER_H_ */
