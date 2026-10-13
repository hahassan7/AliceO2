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
#ifndef ALICEO2_FOCAL_PADINDEX_H
#define ALICEO2_FOCAL_PADINDEX_H

#include <tuple>

namespace o2::focal
{

/// \class PadIndex
/// \brief Channel indexing of the FOCAL-E pad layers
/// \ingroup FOCALbase
///
/// The pads are indexed on a global grid that counts ONLY real pads
/// (gaps between sensors/modules, guard rings and the beam-pipe hole take no index):
///
///   col   0 .. nCols-1   left -> right (increasing global x)
///   row   0 .. nRows-1   bottom -> top (increasing global y)
///   layer 0 .. nLayers-1 ALL ECAL layers (pads + pixels), front -> back, i.e. the same layer
///                        numbering as the geometry (0..19). Pixel layers (4 and 9) simply never
///                        contain pad digits; they are kept so that there is only ONE layer numbering.
///
///   index (absID) = layer * (nRows * nCols) + row * nCols + col
///
/// with nCols = nModulesX * nSensorsX * nPadsX (2*5*9 = 90) and
///      nRows = nModulesY * nSensorsY * nPadsY (11*1*8 = 88), i.e. 20*88*90 = 158400 indices
///      (of which 18*88*90 = 142560 are real pads).
///
/// The location of a pad in the detector hardware (side, module, sensor, pad inside the sensor)
/// is obtained with getPadLocation(); it is a different VIEW of the same pad, never a second index.
///
/// =====================================================================================
/// TODO (IMPORTANT): the PadLocation conventions below are PROVISIONAL working assumptions
/// (2026-09-29) and MUST be updated from the engineers' mechanical drawings, PCB schematics
/// and cabling documents before use with real data:
///   - side   : 0 = left (x < 0), 1 = right (x > 0)
///   - module : 0 .. nModulesY-1, bottom -> top
///   - sensor : 0 = closest to the beam pipe, nSensorsX-1 = furthest (mirror-symmetric)
///   - padX   : 0 .. nPadsX-1 counted outwards from the beam pipe
///   - padY   : 0 .. nPadsY-1 bottom -> top
/// All of these conventions live ONLY in getPadLocation()/getIndexFromPadLocation().
/// =====================================================================================
class PadIndex
{
 public:
  /// \struct PadLocation
  /// \brief Location of a pad in the detector hardware
  struct PadLocation {
    int layer = -1;  ///< ECAL layer (0 .. nLayers-1, same numbering as the geometry)
    int side = -1;   ///< 0 = left (x<0), 1 = right (x>0)
    int module = -1; ///< module (tower) row, bottom -> top
    int sensor = -1; ///< sensor inside the module-layer, 0 = closest to the beam pipe
    int padX = -1;   ///< pad column inside the sensor, counted outwards from the beam pipe
    int padY = -1;   ///< pad row inside the sensor, bottom -> top
  };

  PadIndex() = default;

  /// \param nLayers   number of ECAL layers, pads + pixels (= geometric layers, 20)
  /// \param nModulesX number of modules (towers) in x (= 2 halves)
  /// \param nModulesY number of modules (towers) in y
  /// \param nSensorsX number of sensors (wafers) per module in x
  /// \param nSensorsY number of sensors (wafers) per module in y
  /// \param nPadsX    number of pads per sensor in x
  /// \param nPadsY    number of pads per sensor in y
  PadIndex(int nLayers, int nModulesX, int nModulesY, int nSensorsX, int nSensorsY, int nPadsX, int nPadsY)
    : mNLayers(nLayers), mNModulesX(nModulesX), mNModulesY(nModulesY), mNSensorsX(nSensorsX), mNSensorsY(nSensorsY), mNPadsX(nPadsX), mNPadsY(nPadsY)
  {
  }

  int getNLayers() const { return mNLayers; }
  int getNCols() const { return mNModulesX * mNSensorsX * mNPadsX; }
  int getNRows() const { return mNModulesY * mNSensorsY * mNPadsY; }
  int getNChannelsPerLayer() const { return getNRows() * getNCols(); }
  int getNChannels() const { return mNLayers * getNChannelsPerLayer(); }
  int getNPadsXInSensor() const { return mNPadsX; }
  int getNPadsYInSensor() const { return mNPadsY; }

  bool isValid(int layer, int row, int col) const
  {
    return layer >= 0 && layer < mNLayers && row >= 0 && row < getNRows() && col >= 0 && col < getNCols();
  }
  bool isValid(int index) const { return index >= 0 && index < getNChannels(); }

  /// \brief (layer, row, col) -> pad index. Returns -1 if out of range.
  int getIndexFromLayerRowCol(int layer, int row, int col) const
  {
    if (!isValid(layer, row, col)) {
      return -1;
    }
    return layer * getNChannelsPerLayer() + row * getNCols() + col;
  }

  /// \brief (row, col) -> index of the pad INSIDE its layer = row * nCols + col (< 65536, fits the 16-bit hit detector ID).
  ///        Returns -1 if out of range. absID = layer * getNChannelsPerLayer() + indexInLayer.
  int getIndexInLayer(int row, int col) const
  {
    if (row < 0 || row >= getNRows() || col < 0 || col >= getNCols()) {
      return -1;
    }
    return row * getNCols() + col;
  }

  /// \brief (layer, index inside the layer) -> pad index (absID). Returns -1 if out of range.
  int getIndexFromLayerAndIndexInLayer(int layer, int indexInLayer) const
  {
    if (layer < 0 || layer >= mNLayers || indexInLayer < 0 || indexInLayer >= getNChannelsPerLayer()) {
      return -1;
    }
    return layer * getNChannelsPerLayer() + indexInLayer;
  }

  /// \brief pad index -> {layer, row, col}. Returns {-1,-1,-1} if out of range.
  std::tuple<int /*layer*/, int /*row*/, int /*col*/> getLayerRowColFromIndex(int index) const
  {
    if (!isValid(index)) {
      return {-1, -1, -1};
    }
    const int layer = index / getNChannelsPerLayer();
    const int inLayer = index % getNChannelsPerLayer();
    return {layer, inLayer / getNCols(), inLayer % getNCols()};
  }

  /// \brief Pad index from the indices used to build the geometry in the simulation
  ///        (see Geometry::getPadRowColLayerFromCopyNumbers)
  /// \param layer   ECAL layer (geometric numbering, pads + pixels)
  /// \param moduleX module (tower) index in x, 0 = left (x<0)
  /// \param moduleY module (tower) index in y, 0 = bottom
  /// \param waferX  sensor index inside the module along +x (geometry convention)
  /// \param waferY  sensor index inside the module counted from the TOP (geometry convention)
  /// \param ix      pad index inside the sensor along +x (geometry convention)
  /// \param iy      pad index inside the sensor counted from the TOP (geometry convention)
  int getIndexFromVolumeIndices(int layer, int moduleX, int moduleY, int waferX, int waferY, int ix, int iy) const
  {
    int col = (moduleX * mNSensorsX + waferX) * mNPadsX + ix;
    int row = (moduleY * mNSensorsY + (mNSensorsY - 1 - waferY)) * mNPadsY + (mNPadsY - 1 - iy);
    return getIndexFromLayerRowCol(layer, row, col);
  }

  /// \brief (layer, row, col) -> location in the hardware (PROVISIONAL conventions, see class TODO)
  PadLocation getPadLocation(int layer, int row, int col) const
  {
    PadLocation loc;
    if (!isValid(layer, row, col)) {
      return loc;
    }
    const int colsPerSide = mNSensorsX * mNPadsX;
    loc.layer = layer;
    loc.side = col / colsPerSide;
    int u = col % colsPerSide; // 0 .. colsPerSide-1, left -> right inside this side
    if (loc.side == 0) {
      u = colsPerSide - 1 - u; // left side: count outwards = towards -x
    }
    loc.sensor = u / mNPadsX;
    loc.padX = u % mNPadsX;
    const int rowsPerModule = mNSensorsY * mNPadsY;
    loc.module = row / rowsPerModule;
    loc.padY = row % rowsPerModule; // with nSensorsY = 1 this is the pad row inside the sensor
    return loc;
  }
  /// \brief pad index -> location in the hardware
  PadLocation getPadLocation(int index) const
  {
    auto [layer, row, col] = getLayerRowColFromIndex(index);
    return getPadLocation(layer, row, col);
  }

  /// \brief Location in the hardware -> pad index (inverse of getPadLocation). Returns -1 if out of range.
  int getIndexFromPadLocation(const PadLocation& loc) const
  {
    const int colsPerSide = mNSensorsX * mNPadsX;
    const int rowsPerModule = mNSensorsY * mNPadsY;
    if (loc.side < 0 || loc.side >= mNModulesX || loc.sensor < 0 || loc.sensor >= mNSensorsX || loc.padX < 0 || loc.padX >= mNPadsX ||
        loc.module < 0 || loc.module >= mNModulesY || loc.padY < 0 || loc.padY >= rowsPerModule) {
      return -1;
    }
    int u = loc.sensor * mNPadsX + loc.padX;
    if (loc.side == 0) {
      u = colsPerSide - 1 - u;
    }
    const int col = loc.side * colsPerSide + u;
    const int row = loc.module * rowsPerModule + loc.padY;
    return getIndexFromLayerRowCol(loc.layer, row, col);
  }

 private:
  int mNLayers = 0;
  int mNModulesX = 0;
  int mNModulesY = 0;
  int mNSensorsX = 0;
  int mNSensorsY = 0;
  int mNPadsX = 0;
  int mNPadsY = 0;
};

} // namespace o2::focal

#endif // ALICEO2_FOCAL_PADINDEX_H
