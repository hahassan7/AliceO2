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
#ifndef ALICEO2_FOCAL_HCALINDEX_H
#define ALICEO2_FOCAL_HCALINDEX_H

#include <tuple>

namespace o2::focal
{

/// \class HCALIndex
/// \brief Index of the FOCAL-H SiPMs (sheet design)
/// \ingroup FOCALbase
///
/// One SiPM reads one bundle of nFibresX x nFibresY scintillating fibres (GLOBAL_HCAL_FIBRES_PER_SIPM_X/Y in the
/// geometry file, default 4 x 4). The SiPMs are numbered on a (row, col) grid:
///
///   row 0 .. nRows-1   bottom -> top
///   col 0 .. nCols-1   left -> right (increasing global x)
///   (row 0, col 0) = the bottom-left SiPM
///
///   index = row * nCols + col
///
/// How the fibres are grouped (see Geometry::getHCALRowColFromCopyNumbers):
///   - fibre rows are grouped from the bottom: row = fibreRow / nFibresY;
///   - in each half, fibre columns are grouped starting from the beam-pipe end of a full-width plate: full bundles
///     next to the beam pipe, a partial bundle (if any) at the outer edge.
/// The columns are the same for all pitches. In the fibre rows beside the beam-pipe opening the plates are shorter:
/// the SiPM columns in the middle (in the hole) are empty there, and the column next to the hole may be partial
/// (TODO: to be confirmed with the engineers, see Geometry::getHCALRowColFromCopyNumbers).
///
/// TODO: a readout mapper (as PadMapper for the pads) from (row, col) to the electronics:
///       SiPM board, H2GCROC, half, channel, lpGBT, link, CRU. Needs the cabling documents from the engineers.
class HCALIndex
{
 public:
  HCALIndex() = default;
  HCALIndex(int nRows, int nCols) : mNRows(nRows), mNCols(nCols) {}

  int getNRows() const { return mNRows; }
  int getNCols() const { return mNCols; }
  int getNChannels() const { return mNRows * mNCols; }

  bool isValid(int row, int col) const { return row >= 0 && row < mNRows && col >= 0 && col < mNCols; }
  bool isValid(int index) const { return index >= 0 && index < getNChannels(); }

  /// \brief (row, col) -> SiPM index. Returns -1 if out of range.
  int getIndexFromRowCol(int row, int col) const { return isValid(row, col) ? row * mNCols + col : -1; }

  /// \brief SiPM index -> {row, col}. Returns {-1, -1} if out of range.
  std::tuple<int /*row*/, int /*col*/> getRowColFromIndex(int index) const
  {
    if (!isValid(index)) {
      return {-1, -1};
    }
    return {index / mNCols, index % mNCols};
  }

 private:
  int mNRows = 0;
  int mNCols = 0;
};

} // namespace o2::focal

#endif // ALICEO2_FOCAL_HCALINDEX_H
