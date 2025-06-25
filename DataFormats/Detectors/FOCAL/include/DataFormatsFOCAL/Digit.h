#ifndef ALICEO2_FOCAL_DIGIT_H_
#define ALICEO2_FOCAL_DIGIT_H_

#include <iosfwd>
#include <cmath>
#include "Rtypes.h"
#include "CommonDataFormat/TimeStamp.h"
#include "DataFormatsFOCAL/Constants.h"

#include <boost/serialization/base_object.hpp>

namespace o2
{

namespace focal
{
using DigitBase = o2::dataformats::TimeStamp<double>;

/// \class Digit
/// \brief FOCAL digit implementation
/// \ingroup FOCALDataFormat
class Digit : public DigitBase
{
 public:
  Digit() = default;

  Digit(int col, int row, int layer, double time, double amplitudeGeV);
  Digit(int index, int layer, double time, double amplitudeGeV); 
  ~Digit() = default; 

  bool operator<(const Digit& other) const { return getTimeStamp() < other.getTimeStamp(); }
  bool operator>(const Digit& other) const { return getTimeStamp() > other.getTimeStamp(); }
  bool operator==(const Digit& other) const { return getTimeStamp() == other.getTimeStamp(); }

  bool canAdd(const Digit other)
  {
    return (mCol == other.getCol() && mRow == other.getRow() && mLayer == other.getLayer());
  }

  Digit& operator+=(const Digit& other);
  friend Digit operator+(Digit lhs, const Digit& rhs)
  {
    lhs += rhs;
    return lhs;
  }

  void setCol(int8_t col) { mCol = col; }
  int8_t getCol() const { return mCol; }

  void setRow(int8_t row) { mRow = row; }
  int8_t getRow() const { return mRow; }

  void setLayer(int8_t layer) { mLayer = layer; }
  int8_t getLayer() const { return mLayer; }

  void setIndex(int index) { mIndex = index; }
  int getIndex() const { return mIndex; }

  void setEnergy(double energy) { mAmplitudeGeV = energy; }
  double getEnergy() const { return mAmplitudeGeV; }

  void PrintStream(std::ostream& stream) const;

 private:
  friend class boost::serialization::access;

  double mAmplitudeGeV = 0.;
  int8_t mCol = -1;
  int8_t mRow = -1;
  int8_t mLayer = -1;
  int mIndex = -1;

  ClassDefNV(Digit, 3);
};

std::ostream& operator<<(std::ostream& stream, const Digit& dig);
} // namespace focal
} // namespace o2
#endif
