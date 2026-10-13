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
#ifndef ALICEO2_FOCAL_VIRTUALCELLINDEX_H
#define ALICEO2_FOCAL_VIRTUALCELLINDEX_H

#include <tuple>

namespace o2::focal
{

/// \class VirtualCellIndex
/// \brief Unique cell index for the FOCAL-E pixel layers and FOCAL-H hits on the "virtual" grid
/// \ingroup FOCALbase
///
/// Pixels and HCAL still use the uniform virtual grid of Geometry::getVirtualInfo (col, row, layer).
/// This class packs that triplet into ONE unique, invertible integer (replacing the old
/// detID = row*col + col, which was not unique and did not contain the layer):
///
///   cellID = (layer << 24) | ((row + OFFSET) << 12) | (col + OFFSET)
///
/// - 12 bits for row and col with an OFFSET of 128: rows/cols from -128 to 3967 are allowed
///   (pixel grid 0.05 cm: ~1880 x 2024 cells; getVirtualInfo accepts |x| up to half the FOCAL size
///   plus the 5 cm middle-tower offset, i.e. columns from -100 to ~1980).
/// - 7 bits for the layer (0..127; the geometry has up to 54 layers).
///
/// TODO: to be replaced by a proper physical indexing for pixels and HCAL (as done for the pads,
/// see PadIndex.h) by the people implementing the pixel and HCAL digitization.
class VirtualCellIndex
{
 public:
  static constexpr int OFFSET = 128;
  static constexpr int BITS = 12;
  static constexpr int MASK = (1 << BITS) - 1;
  static constexpr int MAXLAYER = 127;

  static bool isValid(int layer, int row, int col)
  {
    return layer >= 0 && layer <= MAXLAYER && row + OFFSET >= 0 && row + OFFSET <= MASK && col + OFFSET >= 0 && col + OFFSET <= MASK;
  }

  /// \brief (layer, row, col) -> cell index, -1 if out of range
  static int getIndexFromLayerRowCol(int layer, int row, int col)
  {
    if (!isValid(layer, row, col)) {
      return -1;
    }
    return (layer << (2 * BITS)) | ((row + OFFSET) << BITS) | (col + OFFSET);
  }

  /// \brief cell index -> {layer, row, col}; all -1 if the index is negative
  static std::tuple<int /*layer*/, int /*row*/, int /*col*/> getLayerRowColFromIndex(int index)
  {
    if (index < 0) {
      return {-1, -1, -1};
    }
    return {index >> (2 * BITS), ((index >> BITS) & MASK) - OFFSET, (index & MASK) - OFFSET};
  }
};

} // namespace o2::focal

#endif // ALICEO2_FOCAL_VIRTUALCELLINDEX_H
