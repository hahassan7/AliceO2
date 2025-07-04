#ifndef STEER_DIGITIZERWORKFLOW_FOCALDIGITIZER_H_
#define STEER_DIGITIZERWORKFLOW_FOCALDIGITIZER_H_

#include <memory>
#include <vector>

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "DataFormatsFOCAL/Digit.h"
#include "FOCALBase/Hit.h"
#include "FOCALSimulation/Digitizer.h" //we dont have this one yet
#include "FOCALSimulation/SDigitizer.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include <DetectorsBase/BaseDPLDigitizer.h>

class TChain;

namespace o2
{

namespace steer
{
class MCKinematicsReader;
}

namespace focal
{

class DigitizerSpec final : public o2::base::BaseDPLDigitizer, public o2::framework::Task
{
 public:
  using o2::base::BaseDPLDigitizer::init;
  /// \brief Constructor
  DigitizerSpec() : o2::base::BaseDPLDigitizer(o2::base::InitServices::GEOM), o2::framework::Task(){}

  /// \brief Destructor
  ~DigitizerSpec() final = default;

  /// \brief init digitizer
  /// \param ctx Init context
  void initDigitizerTask(framework::InitContext& ctx) final;

  void configure();

  /// \brief run digitizer
  /// \param ctx Processing context
  ///
  /// Handling of pileup events:
  /// - Open readout window when the event sets a trigger
  /// - Accumulate digits sampled via the time response from different bunch crossings
  /// - Retrieve digits when the readout window closes
  void run(framework::ProcessingContext& ctx) override;

 private:
  Bool_t mFinished = false;                   ///< Flag for digitization finished
  bool mIsConfigured = false;                 ///< Initialization status of the digitizer
  bool mRunSDitizer = false;                  ///< Run SDigitization
  Digitizer mDigitizer;                       ///< Digitizer object
  o2::focal::SDigitizer mSumDigitizer;        ///< Summed digitizer
  std::vector<Hit> mHits;                     ///< Vector with input hits
  std::vector<TChain*> mSimChains;
  o2::steer::MCKinematicsReader* mcReader; ///< reader to access MC collision information
  bool mRunDigitizer = true;
};

/// \brief Create new digitizer spec
/// \return Digitizer spec
o2::framework::DataProcessorSpec getFOCALDigitizerSpec(int channel, bool mctruth = true);

} // namespace focal
} // end namespace o2

#endif /* STEER_DIGITIZERWORKFLOW_FOCALDIGITIZER_H_ */
