/// @brief  Processor spec for a ROOT file writer for FOCAL digits

#include "FOCALWorkflow/FOCALDigitWriterSpec.h"
#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include <SimulationDataFormat/MCTruthContainer.h>
#include "DataFormatsFOCAL/Digit.h"
#include <DataFormatsFOCAL/MCLabel.h>
#include "DataFormatsFOCAL/TriggerRecord.h"

using namespace o2::framework;

namespace o2
{
namespace focal
{

template <typename T>
using BranchDefinition = framework::MakeRootTreeWriterSpec::BranchDefinition<T>;

/// create the processor spec
/// describing a processor receiving digits for FOCal writing them to file
DataProcessorSpec getFOCALDigitWriterSpec(bool mctruth)
{
  using InputSpec = framework::InputSpec;
  using MakeRootTreeWriterSpec = framework::MakeRootTreeWriterSpec;
  return MakeRootTreeWriterSpec("FOCALDigitWriter",
                                "focaldigits.root",
                                "o2sim",
                                1,
                                BranchDefinition<std::vector<o2::focal::Digit>>{InputSpec{"focaldigits", "FOC", "DIGITS"}, "FOCALDigit"},
                                BranchDefinition<std::vector<o2::focal::TriggerRecord>>{InputSpec{"trgrecorddigits", "FOC", "TRGRDIG"}, "FOCALDigitTRGR"},
                                BranchDefinition<o2::dataformats::MCTruthContainer<o2::focal::MCLabel>>{InputSpec{"focaldigitlabels", "FOC", "DIGITSMCTR"}, "FOCALDigitMCTruth", mctruth ? 1 : 0})();
}
} // end namespace focal
} // end namespace o2
