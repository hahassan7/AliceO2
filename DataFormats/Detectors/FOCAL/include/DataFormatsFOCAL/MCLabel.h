#ifndef ALICEO2_FOCAL_MCLABEL_H_
#define ALICEO2_FOCAL_MCLABEL_H_

#include "SimulationDataFormat/MCCompLabel.h"

namespace o2
{
namespace focal
{

/// \class MCLabel
/// \brief Monte-Carlo label for FOCAL clusters / digits
/// \ingroup FOCALDataFormat
class MCLabel : public o2::MCCompLabel
{
 private:
  double mAmplitudeFraction;

 public:
  MCLabel() = default;
  MCLabel(int trackID, int eventID, int srcID, bool fake, double afraction) : o2::MCCompLabel(trackID, eventID, srcID, fake), mAmplitudeFraction(afraction) {}
  MCLabel(bool noise, double afraction) : o2::MCCompLabel(noise), mAmplitudeFraction(afraction) {}
  
  void setAmplitudeFraction(double afraction) { mAmplitudeFraction = afraction; }
  double getAmplitudeFraction() const { return mAmplitudeFraction; }

  ClassDefNV(MCLabel, 1);
};
} // namespace focal
} //namespace o2

#endif
