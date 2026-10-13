// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#ifndef ALICEO2_FOCAL_HIT_H
#define ALICEO2_FOCAL_HIT_H

#include <iosfwd>
#include "SimulationDataFormat/BaseHits.h"
#include "CommonUtils/ShmAllocator.h"
#include <boost/container_hash/hash.hpp>

namespace o2::focal
{

/// \class Hit
/// \brief Common FOCAL hit class for the detector simulation
/// \ingroup FOCALbase
/// \author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
/// \since May 29, 2024
///
/// Cell index: the detector ID of the O2 base hit (GetDetectorID) is only 16 bits (max 65535), too small for the
/// full FOCAL-E pad index (absID up to 158400). A hit therefore stores the cell as (layer, index inside the layer):
///   - getIndexInLayer() = the 16-bit base detector ID:
///       pads -> row * nCols + col (PadIndex::getIndexInLayer, max 7919),
///       HCAL (sheets) -> SiPM index (HCALIndex, one layer, max 8000 with 3x3 bundles),
///       pixels -> TRUNCATED to 16 bits (TODO: the pixel digitization must follow the ITS approach: chip ID + position);
///   - getLayer() = layer (one byte): pads 0..19 (geometric numbering), HCAL 20.
/// The full index (e.g. the pad absID) is rebuilt ONLY in the SDigitizer, with Geometry::getIndexFromLayerAndIndexInLayer.
/// Digits always carry the full index.
class Hit : public o2::BasicXYZEHit<float>
{
 public:
  /// \enum Subsystem_t
  /// \brief Subsystem index of the Hit
  enum class Subsystem_t {
    EPADS,   ///< ECAL pads
    EPIXELS, ///< ECAL pixels
    HCAL,    ///< HCAL
    UNKNOWN  ///< Undefined
  };

  /// \struct HitID
  /// \brief Mapped information of a channel
  struct HitID {
    int mParentID;      ///< parentID of the track creating this hit
    uint16_t mRow;      ///< Row of the hit in the calorimeter (uint16: pixel/HCAL grids can exceed 255)
    uint16_t mColumn;   ///< Column of the hit in the calorimeter
    uint8_t mLayer;     ///< Layer the was hit
    uint8_t mSubsystem; ///< Subsystem (EPADS/EPIXELS/HCAL), avoids collisions between the index spaces of different subsystems

    bool operator==(const HitID& other) const
    {
      return mParentID == other.mParentID && mRow == other.mRow && mColumn == other.mColumn && mLayer == other.mLayer && mSubsystem == other.mSubsystem;
    }
    friend std::ostream& operator<<(std::ostream& stream, const Hit::HitID& channel);
  };

  /// \struct HitIDHasher
  /// \brief Hash functor for hit ID
  struct HitIDHasher {

    /// \brief Functor implementation
    /// \param s Hit for which to determine a hash value
    /// \return hash value for channel ID
    size_t operator()(const HitID& s) const
    {
      std::size_t seed = 0;
      boost::hash_combine(seed, s.mParentID);
      boost::hash_combine(seed, s.mRow);
      boost::hash_combine(seed, s.mColumn);
      boost::hash_combine(seed, s.mLayer);
      boost::hash_combine(seed, s.mSubsystem);
      return seed;
    }
  };

  /// \brief Dummy constructor
  Hit() = default;

  /// \brief Construction of the FOCAL hit with full information
  /// \param primary Index of the incoming primary particle
  /// \param trackID Index of the MC particle (also shower particle) responsible for the energy loss
  /// \param indexInLayer Index of the cell inside its layer (16 bits, stored as the detector ID of the base hit)
  /// \param layer Layer of the cell
  /// \param subsystem FOCAL Subdetector (E-Pads, E-Pixels, HCAL)
  /// \param pos Geometric hit position (global coordintate system)
  /// \param tof Time-of-flight of the particle to the FOCAL
  /// \param eLoss Energy loss
  Hit(int primary, int trackID, int indexInLayer, int layer, Subsystem_t subsystem, double initialEnergy, const math_utils::Point3D<float>& pos,
      double tof, double eLoss);

  /// \brief Destructor
  ~Hit() = default;

  /// \brief Comparison operator for equalness
  /// \param other Hit to compare to
  /// \return True if subsytem, module and MC particle ID match, false otherwise
  bool operator==(const Hit& other) const;

  /// \brief Comparison operator for smaller
  /// \param other Hit to compare to
  /// \return True if other hit is smaller (first track ID, then, subsystem ID, then module ID), false otherwise
  bool operator<(const Hit& other) const;

  /// \brief Operator for incremental sum, adding energy loss of the other hit to this energy loss
  /// \param other Hit to add to this one
  /// \return This hit containing the sum of the two energy losses
  Hit& operator+=(const Hit& other);

  /// \brief Get the type of the subsystem for which the hit was created
  /// \return Subsystem type
  Subsystem_t getSubsystem() const { return static_cast<Subsystem_t>(mSubSystem); }

  /// \brief Get the index of the cell INSIDE ITS LAYER (= GetDetectorID, 16 bits). NOT the full index (e.g. pad absID),
  ///        use Geometry::getIndexFromLayerAndIndexInLayer for that.
  /// \return Index inside the layer (pads: row*nCols+col, HCAL: SiPM index)
  int getIndexInLayer() const { return GetDetectorID(); }

  /// \brief Get the layer of the cell
  /// \return Layer (pads 0..19, HCAL 20)
  int getLayer() const { return mLayer; }

  /// \brief Check if the hit is a FOCAL-E pixel hit
  /// \return True if the hit is a FOCAL-E pixel hit, false otherwise
  bool isPixelHit() const { return getSubsystem() == Subsystem_t::EPIXELS; }

  /// \brief Check if the hit is a FOCAL-E pad hit
  /// \return True if the hit is a FOCAL-E pad hit, false otherwise
  bool isPadHit() const { return getSubsystem() == Subsystem_t::EPADS; }

  /// \brief Check if the hit is a FOCAL-H hit
  /// \return True if the hit is a FOCAL-H, false otherwise
  bool isHCALHit() const { return getSubsystem() == Subsystem_t::HCAL; }

  /// \brief Get index of the incomimg primary particle associated with the hit
  /// \return Associated primary particle
  int getPrimary() const { return mPrimary; }

  /// \brief Get energy of the incoming primary particle at the entrance of FOCAL
  /// \return Initial energy
  double getInitialEnergy() const { return mInitialEnergy; }

  /// \brief Set energy of the incoming primary particle at the entrance of FOCAL
  /// \param energy Initial energy
  void setInitialEnergy(double energy) { mInitialEnergy = energy; }

  /// \brief Set index of the incomimg primary particle associated with the hit
  /// \param primary Associated primary particle
  void setPrimary(int primary) { mPrimary = primary; }

  /// \brief Print information of this hit on the output stream
  /// \param stream Stream to print on
  void printStream(std::ostream& stream) const;

 private:
  uint8_t mLayer = 0;                                                     ///< Layer of the cell (see class description)
  uint8_t mSubSystem = static_cast<uint8_t>(Subsystem_t::UNKNOWN);        ///< FOCAL subdetector (Subsystem_t), one byte
  int mPrimary = -1;                                                      ///< Primary particles at the origin of the hit
  double mInitialEnergy = 0.;                                             ///< Energy of the parent particle that entered the EMCAL

  ClassDefNV(Hit, 2);
};

/// @brief Sum operator, creating a new hit with the sum of the two energy losses
/// @param lhs Left-hand side of the sum
/// @param rhs Right-hand side of the sum
/// @return New hit with the properties of the lhs hit and the summed energy loss of the two hits
Hit operator+(const Hit& lhs, const Hit& rhs);

/// \brief Output stream operator for FOCAL hits
/// \param stream Stream to write on
/// \param point Hit to be printed
/// \return Stream after printing
std::ostream& operator<<(std::ostream& stream, const Hit& point);
} // namespace o2::focal

#ifdef USESHM
namespace std
{
template <>
class allocator<o2::focal::Hit> : public o2::utils::ShmAllocator<o2::focal::Hit>
{
};
} // namespace std
#endif
#endif
