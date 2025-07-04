
#include "FOCALWorkflow/FOCALDigitizerSpec.h"
#include "CommonConstants/Triggers.h"
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

#include "CommonDataFormat/EvIndex.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsFOCAL/TriggerRecord.h"
#include <Steer/MCKinematicsReader.h>

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

  auto geom = o2::focal::Geometry::GetInstance(); 

  mSumDigitizer.setGeometry(geom);
  mDigitizer.setGeometry(geom);

  if (ctx.options().get<bool>("debug-stream")) {
    mDigitizer.setDebugStreaming(true);
  }

  if (ctx.options().get<bool>("disable-dig")) {
    mRunDigitizer = false;
  }

  mFinished = false;
}

void DigitizerSpec::run(framework::ProcessingContext& ctx)
{
  if (mFinished) {
    return;
  }

  if (!mIsConfigured) {
    configure();
    mIsConfigured = true;
  }

  o2::focal::SimParam::Instance().printKeyValues(true, true);

  mDigitizer.flush();

  auto context = ctx.inputs().get<o2::steer::DigitizationContext*>("collisioncontext");

  auto intRate = context->getDigitizerInteractionRate();
  context->initSimChains(o2::detectors::DetID::FOC, mSimChains);

  if (mcReader == nullptr) {
    mcReader = new o2::steer::MCKinematicsReader(context.get());
  }

  auto& timesview = context->getEventRecords();
  LOG(debug) << "GOT " << timesview.size() << " COLLISSION TIMES";

  if (timesview.size() == 0) {
    return;
  }

  TStopwatch timer;
  timer.Start();

  auto& eventParts = context->getEventParts();

  int collisionN = 0;
  for (int collID = 0; collID < timesview.size(); ++collID) {

    if (mRunDigitizer == false) {
      break;
    }

    if (intRate < 10000.) {
      break;
    }

    collisionN++;

    mDigitizer.setEventTime(timesview[collID]);

    for (auto& part : eventParts[collID]) {

      mSumDigitizer.setCurrEvID(part.entryID);
      mSumDigitizer.setCurrSrcID(part.sourceID);

      mHits.clear();
      context->retrieveHits(mSimChains, "FOCHit", part.sourceID, part.entryID, &mHits);

      std::vector<o2::focal::LabeledDigit> summedLabeledDigits;
      std::vector<o2::focal::Digit> summedDigits;
      if (mRunSDitizer) {
        summedLabeledDigits = mSumDigitizer.process(mHits);
        for (auto labeledsummeddigit : summedLabeledDigits) {
          summedDigits.push_back(labeledsummeddigit.getDigit());
        }
      } else {
        for (auto& hit : mHits) {
          summedDigits.emplace_back(hit.GetDetectorID(), hit.GetEnergyLoss(), hit.GetTime());
        }
      }

      mDigitizer.process(summedDigits);
    }
  }
  mDigitizer.finish();

  LOG(info) << " CALLING FOCAL DIGITIZATION ";
  o2::dataformats::MCTruthContainer<o2::focal::MCLabel> labelAccum;

  for (int collID = 0; collID < timesview.size(); ++collID) {

    std::bitset<5> trigger{0x1}; 

    if (!trigger.any()) {
      continue;
    }

    mDigitizer.setEventTime(timesview[collID], trigger.any());
    if (!mDigitizer.isCurrentEventTriggered()) {
      LOG(debug) << "reject collision";
      continue;
    }
    LOG(debug) << "accept collision";

    for (auto& part : eventParts[collID]) {

      mSumDigitizer.setCurrEvID(part.entryID);
      mSumDigitizer.setCurrSrcID(part.sourceID);

      auto& mcEventHeader = mcReader->getMCEventHeader(part.sourceID, part.entryID);
      mcEventHeader.print();

      mHits.clear();
      context->retrieveHits(mSimChains, "FOCHit", part.sourceID, part.entryID, &mHits);

      LOG(info) << "For collision " << collID << " eventID " << part.entryID << " found " << mHits.size() << " hits ";

      std::vector<o2::focal::LabeledDigit> summedLabeledDigits;
      if (mRunSDitizer) {
        summedLabeledDigits = mSumDigitizer.process(mHits);
      } else {
        for (auto& hit : mHits) {
          o2::focal::MCLabel digitlabel(hit.GetTrackID(), part.entryID, part.sourceID, false, 1.);
          if (hit.GetEnergyLoss() < __DBL_EPSILON__) {
            digitlabel.setAmplitudeFraction(0);
          }
          summedLabeledDigits.emplace_back(hit.GetDetectorID(), hit.GetEnergyLoss(), hit.GetTime(), digitlabel);
        }
      }

      mDigitizer.process(summedLabeledDigits);
    }
  }
  mDigitizer.finish();

  ctx.outputs().snapshot(Output{"FOC", "DIGITS", 0}, mDigitizer.getDigits());
  if (ctx.outputs().isAllowed({"FOC", "DIGITSMCTR", 0})) {
    ctx.outputs().snapshot(Output{"FOC", "DIGITSMCTR", 0}, mDigitizer.getMCLabels());
  }
  const o2::parameters::GRPObject::ROMode roMode = o2::parameters::GRPObject::TRIGGERING;
  LOG(info) << "FOCAL: Sending ROMode= " << roMode << " to GRPUpdater";
  ctx.outputs().snapshot(Output{"FOC", "ROMode", 0}, roMode);

  timer.Stop();
  LOG(info) << "Digitization took " << timer.CpuTime() << "s";
  ctx.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  mFinished = true;
}

void DigitizerSpec::configure()
{
  mDigitizer.init();
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
  outputs.emplace_back("FOC", "TRIGGERINPUT", 0, Lifetime::Timeframe);

  std::vector<o2::framework::InputSpec> inputs;
  inputs.emplace_back("collisioncontext", "SIM", "COLLISIONCONTEXT", static_cast<SubSpecificationType>(channel), Lifetime::Timeframe);


  return DataProcessorSpec{
    "FOCALDigitizer",
    inputs,
    outputs,
    Options{
      {"pileup", VariantType::Int, 1, {"whether to run in continuous time mode"}},
      {"disable-dig", VariantType::Bool, false, {"Disable digitisation"}},
      {"debug-stream", VariantType::Bool, false, {"Enable debug streaming"}}}
  };
}
} // end namespace focal
} // end namespace o2
