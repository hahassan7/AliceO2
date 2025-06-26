#ifndef ALICEO2_FOCAL_LABELEDDIGIT_H_
#define ALICEO2_FOCAL_LABELEDDIGIT_H_

#include <iosfwd>
#include <cmath>
#include "Rtypes.h"
#include "CommonDataFormat/TimeStamp.h"
#include "DataFormatsFOCAL/Constants.h"
#include "DataFormatsFOCAL/Digit.h"
#include "DataFormatsFOCAL/MCLabel.h"

#include <boost/serialization/base_object.hpp> 

namespace o2
{

namespace focal
{
/// \class LabeledDigit
/// \brief FOCAL labeled digit implementation
/// \ingroup FOCALsimulation

class LabeledDigit
{
 public:
  LabeledDigit() = default;

  LabeledDigit(Digit digit, o2::focal::MCLabel label);
  LabeledDigit(int col, int row, int layer,double time, double amplitudeGeV, o2::focal::MCLabel label);
  LabeledDigit(int index, int layer, double time, double amplitudeGeV, o2::focal::MCLabel label);

  ~LabeledDigit() = default;

  void setDigit(Digit d) { mDigit = d; }
  Digit getDigit() const { return mDigit; }

  void addLabel(o2::focal::MCLabel l) { mLabels.push_back(l); }
  int getNumberOfLabels() const { return mLabels.size(); }
  std::vector<o2::focal::MCLabel> getLabels() const { return mLabels; }

  bool operator<(const LabeledDigit& other) const { return getTimeStamp() < other.getTimeStamp(); }
  bool operator>(const LabeledDigit& other) const { return getTimeStamp() > other.getTimeStamp(); }
  bool operator==(const LabeledDigit& other) const { return (getTimeStamp() == other.getTimeStamp()); }

  bool canAdd(const LabeledDigit other)
  {
    return (getCol() == other.getCol() && getRow() == other.getRow() && getLayer() == other.getLayer());
  }

  LabeledDigit& operator+=(const LabeledDigit& other);                    
  friend LabeledDigit operator+(LabeledDigit lhs, const LabeledDigit& rhs) 
  {
    lhs += rhs;
    return lhs;
  }

  void setCol(int8_t col) { mDigit.setCol(col); }
  int8_t getCol() const { return mDigit.getCol(); }
  
  void setRow(int8_t row) { mDigit.setRow(row); }
  int8_t getRow() const { return mDigit.getRow(); }

  void setLayer(int8_t layer) { mDigit.setLayer(layer); }
  int8_t getLayer() const { return mDigit.getLayer(); }

  void setIndex(int index) { mDigit.setIndex(index); }
  int getIndex() const { return mDigit.getIndex(); }

  void setEnergy(double energy) { mDigit.setEnergy(energy); }
  double getEnergy() const { return mDigit.getEnergy(); }

  void setTimeStamp(double time) { mDigit.setTimeStamp(time); }
  double getTimeStamp() const { return mDigit.getTimeStamp(); }


  void PrintStream(std::ostream& stream) const;

 private:
  friend class boost::serialization::access;

  Digit mDigit;                            ///< Digit
  std::vector<o2::focal::MCLabel> mLabels; ///< Labels

  ClassDefNV(LabeledDigit, 2);
};

std::ostream& operator<<(std::ostream& stream, const LabeledDigit& dig);
} // namespace focal
} // namespace o2
#endif
