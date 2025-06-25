#include "DataFormatsFOCAL/Digit.h"
#include <iostream>

using namespace o2::focal;

Digit::Digit(int col, int row, int layer, double time, double amplitudeGeV)
  : DigitBase(time), mCol(col), mRow(row), mLayer(layer), mAmplitudeGeV(amplitudeGeV)
{
}

Digit::Digit(int index, int layer, double time, double amplitudeGeV) 
  : DigitBase(time), mIndex(index), mLayer(layer), mAmplitudeGeV(amplitudeGeV)
{
}

Digit& Digit::operator+=(const Digit& other)
{
  if (canAdd(other)) {
    mAmplitudeGeV += other.mAmplitudeGeV;
  }
  return *this;
}

void Digit::PrintStream(std::ostream& stream) const
{
  stream << "FOCAL Digit: Index " << mIndex << ", Energy " << getEnergy();
}

std::ostream& operator<<(std::ostream& stream, const Digit& digi)
{
  digi.PrintStream(stream);
  return stream;
}
