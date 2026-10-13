#ifndef ALICEO2_FOCAL_LABELEDSDIGIT_H_
#define ALICEO2_FOCAL_LABELEDSDIGIT_H_

#include <iosfwd>
#include <cmath>
#include <vector>
#include "Rtypes.h"
#include "CommonDataFormat/TimeStamp.h"
#include "DataFormatsFOCAL/Constants.h"
#include "DataFormatsFOCAL/SDigit.h"
#include "DataFormatsFOCAL/MCLabel.h"

#include <boost/serialization/base_object.hpp> 

namespace o2
{

namespace focal
{
/// \class LabeledSDigit
/// \brief FOCAL summed digit (SDigit) with its MC labels
/// \ingroup FOCALsimulation

class LabeledSDigit
{
 public:
  LabeledSDigit() = default;

  LabeledSDigit(SDigit digit, o2::focal::MCLabel label);
  LabeledSDigit(int col, int row, int layer, double time, double energy, o2::focal::MCLabel label);
  LabeledSDigit(int index, int layer, double time, double energy, o2::focal::MCLabel label);

  ~LabeledSDigit() = default;

  void setSDigit(SDigit d) { mSDigit = d; }
  SDigit getSDigit() const { return mSDigit; }

  void addLabel(o2::focal::MCLabel l) { mLabels.push_back(l); }
  int getNumberOfLabels() const { return mLabels.size(); }
  std::vector<o2::focal::MCLabel> getLabels() const { return mLabels; }

  bool operator<(const LabeledSDigit& other) const { return getTimeStamp() < other.getTimeStamp(); }
  bool operator>(const LabeledSDigit& other) const { return getTimeStamp() > other.getTimeStamp(); }
  bool operator==(const LabeledSDigit& other) const { return (getTimeStamp() == other.getTimeStamp()); }

  bool canAdd(const LabeledSDigit other)
  {
    return (getSubsystem() == other.getSubsystem() && getCol() == other.getCol() && getRow() == other.getRow() && getLayer() == other.getLayer());
  }

  LabeledSDigit& operator+=(const LabeledSDigit& other);                    
  friend LabeledSDigit operator+(LabeledSDigit lhs, const LabeledSDigit& rhs) 
  {
    lhs += rhs;
    return lhs;
  }

  void setCol(int16_t col) { mSDigit.setCol(col); }
  int16_t getCol() const { return mSDigit.getCol(); }
  
  void setRow(int16_t row) { mSDigit.setRow(row); }
  int16_t getRow() const { return mSDigit.getRow(); }

  void setLayer(int16_t layer) { mSDigit.setLayer(layer); }
  int16_t getLayer() const { return mSDigit.getLayer(); }

  void setIndex(int index) { mSDigit.setIndex(index); }
  int getIndex() const { return mSDigit.getIndex(); }

  void setSubsystem(SDigit::Subsystem_t subsystem) { mSDigit.setSubsystem(subsystem); }
  SDigit::Subsystem_t getSubsystem() const { return mSDigit.getSubsystem(); }

  /// \brief Energy in eV
  void setEnergy(double energy) { mSDigit.setEnergy(energy); }
  double getEnergy() const { return mSDigit.getEnergy(); }

  void setTimeStamp(double time) { mSDigit.setTimeStamp(time); }
  double getTimeStamp() const { return mSDigit.getTimeStamp(); }


  void PrintStream(std::ostream& stream) const;

 private:
  friend class boost::serialization::access;

  SDigit mSDigit;                            ///< SDigit
  std::vector<o2::focal::MCLabel> mLabels; ///< Labels

  ClassDefNV(LabeledSDigit, 3);
};

std::ostream& operator<<(std::ostream& stream, const LabeledSDigit& dig);
} // namespace focal
} // namespace o2
#endif
