#ifndef ALICEO2_FOCAL_SDIGIT_H_
#define ALICEO2_FOCAL_SDIGIT_H_

#include <iosfwd>
#include <cmath>
#include <cstdint>
#include "Rtypes.h"
#include "CommonDataFormat/TimeStamp.h"
#include "DataFormatsFOCAL/Constants.h"

#include <boost/serialization/base_object.hpp>

namespace o2
{

namespace focal
{
using SDigitBase = o2::dataformats::TimeStamp<double>;

/// \class SDigit
/// \brief FOCAL summed digit (SDigit): energy-level input of the digitizer
/// \ingroup FOCALDataFormat
///
/// Energy deposited in one cell by one collision (sum of the hits). The electronics-level
/// digit (ADC, TOA, TOT) is o2::focal::Digit, produced from the SDigits by the digitizer.
/// - energy  : in eV (same unit as the FOCAL hits, Edep*1e9)
/// - layer   : geometric layer numbering (0..19 ECAL pads + pixels, >= 20 HCAL)
/// - row/col : pads -> grid of real pads (row 0..87, col 0..89, see FOCALBase/PadIndex.h)
///             pixels/HCAL -> uniform virtual grid (see FOCALBase/VirtualCellIndex.h)
/// - index   : unique cell index inside its subsystem
///             pads -> PadIndex absID, pixels/HCAL -> VirtualCellIndex cellID
/// - subsystem: which index space the digit belongs to (same values as o2::focal::Hit::Subsystem_t)
class SDigit : public SDigitBase
{
 public:
  /// \enum Subsystem_t
  /// \brief Subsystem of the digit. Values MUST match o2::focal::Hit::Subsystem_t (checked in the SDigitizer)
  enum Subsystem_t : uint8_t {
    kEPads = 0,   ///< ECAL pads
    kEPixels = 1, ///< ECAL pixels
    kHCAL = 2,    ///< HCAL
    kUnknown = 3  ///< Undefined
  };

  SDigit() = default;

  /// \param col Column of the cell
  /// \param row Row of the cell
  /// \param layer Layer of the cell
  /// \param time Time of the digit (ns)
  /// \param energy Deposited energy (eV)
  SDigit(int col, int row, int layer, double time, double energy);
  SDigit(int index, int layer, double time, double energy);
  ~SDigit() = default;

  bool operator<(const SDigit& other) const { return getTimeStamp() < other.getTimeStamp(); }
  bool operator>(const SDigit& other) const { return getTimeStamp() > other.getTimeStamp(); }
  bool operator==(const SDigit& other) const { return getTimeStamp() == other.getTimeStamp(); }

  /// \brief Two digits can be added if they belong to the same cell of the same subsystem
  bool canAdd(const SDigit other)
  {
    return (mSubsystem == other.getSubsystem() && mCol == other.getCol() && mRow == other.getRow() && mLayer == other.getLayer());
  }

  SDigit& operator+=(const SDigit& other);
  friend SDigit operator+(SDigit lhs, const SDigit& rhs)
  {
    lhs += rhs;
    return lhs;
  }

  void setCol(int16_t col) { mCol = col; }
  int16_t getCol() const { return mCol; }

  void setRow(int16_t row) { mRow = row; }
  int16_t getRow() const { return mRow; }

  void setLayer(int16_t layer) { mLayer = layer; }
  int16_t getLayer() const { return mLayer; }

  /// \brief Unique cell index inside the subsystem (pads: PadIndex absID, pixels/HCAL: VirtualCellIndex cellID)
  void setIndex(int index) { mIndex = index; }
  int getIndex() const { return mIndex; }

  void setSubsystem(Subsystem_t subsystem) { mSubsystem = subsystem; }
  Subsystem_t getSubsystem() const { return static_cast<Subsystem_t>(mSubsystem); }
  bool isPadDigit() const { return mSubsystem == kEPads; }
  bool isPixelDigit() const { return mSubsystem == kEPixels; }
  bool isHCALDigit() const { return mSubsystem == kHCAL; }

  /// \brief Energy in eV
  void setEnergy(double energy) { mEnergy = energy; }
  /// \brief Energy in eV
  double getEnergy() const { return mEnergy; }

  void PrintStream(std::ostream& stream) const;

 private:
  friend class boost::serialization::access;

  double mEnergy = 0.;            ///< Deposited energy (eV)
  int16_t mCol = -1;              ///< Column of the cell
  int16_t mRow = -1;              ///< Row of the cell
  int16_t mLayer = -1;            ///< Layer of the cell (geometric numbering)
  int mIndex = -1;                ///< Unique cell index inside the subsystem
  uint8_t mSubsystem = kUnknown;  ///< Subsystem (Subsystem_t)

  ClassDefNV(SDigit, 1);
};

std::ostream& operator<<(std::ostream& stream, const SDigit& dig);
} // namespace focal
} // namespace o2
#endif
