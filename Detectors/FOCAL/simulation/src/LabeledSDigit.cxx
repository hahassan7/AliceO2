#include "FOCALSimulation/LabeledSDigit.h"
#include <iostream>

using namespace o2::focal;

LabeledSDigit::LabeledSDigit(SDigit digit, o2::focal::MCLabel label)
  : mSDigit(digit)
{
  mLabels.push_back(label);
}

LabeledSDigit::LabeledSDigit(int col, int row, int layer, double time, double energy, o2::focal::MCLabel label)
  : mSDigit(col, row, layer, time, energy)
{
  mLabels.push_back(label);
}

LabeledSDigit::LabeledSDigit(int index, int layer, double time, double energy, o2::focal::MCLabel label)
  : mSDigit(index, layer, time, energy)
{
  mLabels.push_back(label);
}

LabeledSDigit& LabeledSDigit::operator+=(const LabeledSDigit& other)
{
  if (canAdd(other)) {
    double a1 = getEnergy();
    double a2 = other.getEnergy();
    double r = ((a1 + a2) != 0) ? 1.0 / (a1 + a2) : 0.0;
    mSDigit += other.getSDigit();

    for (size_t j = 0; j < mLabels.size(); j++) {
      mLabels.at(j).setAmplitudeFraction(mLabels.at(j).getAmplitudeFraction() * a1 * r);
    }

    for (auto label : other.getLabels()) {
      label.setAmplitudeFraction(label.getAmplitudeFraction() * a2 * r);
      mLabels.push_back(label);
    }
  }
  return *this;
}

void LabeledSDigit::PrintStream(std::ostream& stream) const
{
  stream << "FOCAL LabeledSDigit: Subsystem " << static_cast<int>(getSubsystem()) << ", Index " << getIndex() << ", Time " << getTimeStamp() << " ns, Energy " << getEnergy() << " eV, Labels ( ";
  for (auto label : mLabels) {
    stream << label.getRawValue() << " ";
  }
  stream << ")";
}

std::ostream& o2::focal::operator<<(std::ostream& stream, const LabeledSDigit& digi)
{
  digi.PrintStream(stream);
  return stream;
}