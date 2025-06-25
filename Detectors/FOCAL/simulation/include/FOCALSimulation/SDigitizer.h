#ifndef ALICEO2_FOCAL_SDIGITIZER_H
#define ALICEO2_FOCAL_SDIGITIZER_H

#include <memory>
#include <unordered_map>
#include <vector>
#include <list>

#include "Rtypes.h"
#include "TObject.h"

#include "DataFormatsFOCAL/Digit.h"
#include "FOCALBase/Geometry.h"
#include "FOCALBase/Hit.h"
#include "FOCALSimulation/LabeledDigit.h"

#include "SimulationDataFormat/MCTruthContainer.h"

namespace o2
{
namespace focal
{

/// \class SDigitizer
/// \brief FOCAL summed digitizer
/// \ingroup FOCALsimulation

class SDigitizer
{
 public:
  SDigitizer() = default;
  ~SDigitizer() = default;
  SDigitizer(const SDigitizer&) = delete;
  SDigitizer& operator=(const SDigitizer&) = delete;

  std::vector<o2::focal::LabeledDigit> process(const std::vector<Hit>& hits);

  void setCurrSrcID(int v);
  int getCurrSrcID() const { return mCurrSrcID; }

  void setCurrEvID(int v);
  int getCurrEvID() const { return mCurrEvID; }

  void setGeometry(const o2::focal::Geometry* gm) { mGeometry = gm; }

 private:
  const Geometry* mGeometry = nullptr; 
  int mCurrSrcID = 0;                  
  int mCurrEvID = 0;                   

  ClassDefNV(SDigitizer, 1);
};
} // namespace focal
} // namespace o2

#endif /* ALICEO2_FOCAL_SDIGITIZER_H */