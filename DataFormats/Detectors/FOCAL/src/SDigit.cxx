#include "DataFormatsFOCAL/SDigit.h"
#include <iostream>

using namespace o2::focal;

SDigit::SDigit(int col, int row, int layer, double time, double energy)
  : SDigitBase(time), mEnergy(energy), mCol(col), mRow(row), mLayer(layer)
{
}

SDigit::SDigit(int index, int layer, double time, double energy)
  : SDigitBase(time), mEnergy(energy), mLayer(layer), mIndex(index)
{
}

SDigit& SDigit::operator+=(const SDigit& other)
{
  if (canAdd(other)) {
    mEnergy += other.mEnergy;
  }
  return *this;
}

void SDigit::PrintStream(std::ostream& stream) const
{
  stream << "FOCAL SDigit: Subsystem " << static_cast<int>(mSubsystem) << ", Index " << mIndex << ", Layer " << mLayer << ", Row " << mRow << ", Col " << mCol
         << ", Time " << getTimeStamp() << " ns, Energy " << getEnergy() << " eV";
}

std::ostream& o2::focal::operator<<(std::ostream& stream, const SDigit& digi)
{
  digi.PrintStream(stream);
  return stream;
}
