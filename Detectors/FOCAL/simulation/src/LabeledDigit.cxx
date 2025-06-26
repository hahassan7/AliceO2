#include "FOCALSimulation/LabeledDigit.h"
#include <iostream>

using namespace o2::focal;

LabeledDigit::LabeledDigit(Digit digit, o2::focal::MCLabel label)
  : mDigit(digit)
{
  mLabels.push_back(label);
}

LabeledDigit::LabeledDigit(int col, int row, int layer, double time, double amplitudeGeV, o2::focal::MCLabel label)
  : mDigit(col, row, layer, time, amplitudeGeV)
{
  mLabels.push_back(label);
}

LabeledDigit::LabeledDigit(int index, int layer, double time, double amplitudeGeV, o2::focal::MCLabel label)
  : mDigit(index, layer, time, amplitudeGeV)
{
  mLabels.push_back(label);
}

LabeledDigit& LabeledDigit::operator+=(const LabeledDigit& other)
{
  if (canAdd(other)) {
    double a1 = getEnergy();
    double a2 = other.getEnergy();
    double r = ((a1 + a2) != 0) ? 1.0 / (a1 + a2) : 0.0;
    mDigit += other.getDigit();

    for (int j = 0; j < mLabels.size(); j++) {
      mLabels.at(j).setAmplitudeFraction(mLabels.at(j).getAmplitudeFraction() * a1 * r);
    }

    for (auto label : other.getLabels()) {
      label.setAmplitudeFraction(label.getAmplitudeFraction() * a2 * r);
      mLabels.push_back(label);
    }
  }
  return *this;
}

void LabeledDigit::PrintStream(std::ostream& stream) const
{
  stream << "FOCAL LabeledDigit: Index " << getIndex() << ", Time " << getTimeStamp() << ", Energy " << getEnergy() << " GeV, Labels ( ";
  for (auto label : mLabels) {
    stream << label.getRawValue() << " ";
  }
  stream << ")";
}

std::ostream& operator<<(std::ostream& stream, const LabeledDigit& digi)
{
  digi.PrintStream(stream);
  return stream;
}