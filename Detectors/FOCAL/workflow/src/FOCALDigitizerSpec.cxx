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

#include "FOCALWorkflow/FOCALDigitizerSpec.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DataRefUtils.h"
#include "Framework/Lifetime.h"
#include "Headers/DataHeader.h"
#include "TStopwatch.h"
#include "Steer/HitProcessingManager.h" // for DigitizationContext
#include "TChain.h"
#include <TGeoManager.h>
#include <TString.h>
#include <TSystem.h>

#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsFOCAL/TriggerRecord.h"
#include "DataFormatsFOCAL/MCLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include <map>
#include <tuple>

using namespace o2::framework;
using SubSpecificationType = o2::framework::DataAllocator::SubSpecificationType;

namespace o2
{
namespace focal
{

void DigitizerSpec::initDigitizerTask(framework::InitContext& ctx)
{
  if (!gGeoManager) {
    LOG(error) << "Geometry needs to be loaded before";
  }

  // Use the same geometry description as the simulation (macro/build_geometry.C). Without a file,
  // Geometry::getInstance() only sets default parameters (FOCAL_Z0) and the pad dimensions are all zero.
  // TODO: make the geometry file configurable (same option as in the simulation)
  TString geoFile = "$O2_ROOT/share/Detectors/Geometry/FOC/geometryFiles/geometry_Sheets.txt";
  gSystem->ExpandPathName(geoFile);
  auto geom = o2::focal::Geometry::getInstance(geoFile.Data());

  mSumDigitizer.setGeometry(geom);

  mFinished = false;
}

void DigitizerSpec::run(framework::ProcessingContext& ctx)
{
  if (mFinished) {
    return;
  }

  auto context = ctx.inputs().get<o2::steer::DigitizationContext*>("collisioncontext");
  context->initSimChains(o2::detectors::DetID::FOC, mSimChains);

  auto& timesview = context->getEventRecords();  // one interaction record per collision
  auto& eventParts = context->getEventParts();   // MC events (signal, background, ...) of each collision
  LOG(info) << "FOCAL digitizer: got " << timesview.size() << " collisions";

  // Output containers for the whole time frame
  std::vector<o2::focal::SDigit> digitsOut;
  std::vector<o2::focal::TriggerRecord> triggerRecordsOut;
  o2::dataformats::MCTruthContainer<o2::focal::MCLabel> labelsOut;

  TStopwatch timer;
  timer.Start();

  // Single loop over collisions.
  // v0: every collision is treated as triggered (TODO: trigger selection, e.g. FIT minimum bias, and pile-up).
  for (int collID = 0; collID < timesview.size(); ++collID) {

    // Sum the digits of all MC events (parts) of this collision per cell.
    // Key = (subsystem, cell index): the cell index (SDigit::getIndex) is unique inside each subsystem.
    std::map<std::tuple<int, int>, o2::focal::LabeledSDigit> digitsInCollision;

    for (auto& part : eventParts[collID]) {
      mSumDigitizer.setCurrEvID(part.entryID);
      mSumDigitizer.setCurrSrcID(part.sourceID);

      mHits.clear();
      context->retrieveHits(mSimChains, "FOCHit", part.sourceID, part.entryID, &mHits);
      LOG(debug) << "Collision " << collID << ", source " << part.sourceID << ", event " << part.entryID << ": " << mHits.size() << " hits";

      for (const auto& sdigit : mSumDigitizer.process(mHits)) {
        auto key = std::make_tuple(int(sdigit.getSubsystem()), sdigit.getIndex());
        auto found = digitsInCollision.find(key);
        if (found == digitsInCollision.end()) {
          digitsInCollision.emplace(key, sdigit);
        } else {
          found->second += sdigit; // same cell hit by several MC events of this collision: add energy and labels
        }
      }
    }

    // Append the digits of this collision to the output and create its trigger record
    const int firstDigit = digitsOut.size();
    for (const auto& [key, labeledDigit] : digitsInCollision) {
      const int digitIndex = digitsOut.size();
      digitsOut.push_back(labeledDigit.getSDigit());
      for (const auto& label : labeledDigit.getLabels()) {
        labelsOut.addElement(digitIndex, label);
      }
    }
    const int nDigits = digitsOut.size() - firstDigit;
    triggerRecordsOut.emplace_back(timesview[collID], firstDigit, nDigits); // all subsystems (range mAllData)
  }

  LOG(info) << "FOCAL digitizer: " << digitsOut.size() << " digits in " << triggerRecordsOut.size() << " trigger records";

  ctx.outputs().snapshot(Output{"FOC", "DIGITS", 0}, digitsOut);
  ctx.outputs().snapshot(Output{"FOC", "TRGRDIG", 0}, triggerRecordsOut);
  if (ctx.outputs().isAllowed({"FOC", "DIGITSMCTR", 0})) {
    ctx.outputs().snapshot(Output{"FOC", "DIGITSMCTR", 0}, labelsOut);
  }

  const o2::parameters::GRPObject::ROMode roMode = o2::parameters::GRPObject::TRIGGERING;
  LOG(info) << "FOCAL: Sending ROMode= " << roMode << " to GRPUpdater";
  ctx.outputs().snapshot(Output{"FOC", "ROMode", 0}, roMode);

  timer.Stop();
  LOG(info) << "Digitization took " << timer.CpuTime() << "s";
  ctx.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  mFinished = true;
}

o2::framework::DataProcessorSpec getFOCALDigitizerSpec(int channel, bool mctruth)
{
  std::vector<OutputSpec> outputs;
  outputs.emplace_back("FOC", "DIGITS", 0, Lifetime::Timeframe);
  outputs.emplace_back("FOC", "TRGRDIG", 0, Lifetime::Timeframe);
  if (mctruth) {
    outputs.emplace_back("FOC", "DIGITSMCTR", 0, Lifetime::Timeframe);
  }
  outputs.emplace_back("FOC", "ROMode", 0, Lifetime::Timeframe);
  // Note: no "TRIGGERINPUT" output. In EMCAL it carries the EMCAL L0 trigger decisions to the CTP
  // digitizer (that is what EMCAL's first collision loop computes); FOCAL is not a CTP input.

  std::vector<o2::framework::InputSpec> inputs;
  inputs.emplace_back("collisioncontext", "SIM", "COLLISIONCONTEXT", static_cast<SubSpecificationType>(channel), Lifetime::Timeframe);

  return DataProcessorSpec{
    "FOCALDigitizer",
    inputs,
    outputs,
    AlgorithmSpec{o2::framework::adaptFromTask<DigitizerSpec>()},
    Options{}};
}
} // end namespace focal
} // end namespace o2
