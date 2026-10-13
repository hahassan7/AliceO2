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

#ifndef ALICEO2_FOCAL_GEOMETRY_H_
#define ALICEO2_FOCAL_GEOMETRY_H_

#include <string>
#include <vector>
#include <array>
#include <list>
#include <tuple>

#include "FOCALBase/Composition.h"
#include "FOCALBase/PadIndex.h"
#include "FOCALBase/HCALIndex.h"
#include "FOCALBase/Hit.h"

namespace o2
{
namespace focal
{

class VirtualSegment
{

 public:
  VirtualSegment()
  {
  }

  VirtualSegment& operator=(const VirtualSegment& segment)
  {
    if (this != &segment) {
      mMinLayer = segment.mMinLayer;
      mMaxLayer = segment.mMaxLayer;
      mPadSize = segment.mPadSize;
      mRelativeSensitiveThickness = segment.mRelativeSensitiveThickness;
      mIsPixel = segment.mIsPixel;
      mPixelTreshold = segment.mPixelTreshold;
    }
    return *this;
  }

  ~VirtualSegment() {}

  int mMinLayer = -1;
  int mMaxLayer = -1;
  float mPadSize = -1.0;
  float mRelativeSensitiveThickness = -1.0;
  int mIsPixel = 0; // 0: pad or strip; 1: pixel; 2: HCAL
  float mPixelTreshold = -1.0;
}; // class VirtualSegment

class Geometry
{
 public:
  enum HCALDesgin {
    Sandwich = 0,  // Sandwich HCAL design
    Spaghetti = 1, // Spaghetti design
    Sheets = 2     // Sheets design
  };

  Geometry() = default;
  Geometry(Geometry* geo);
  Geometry(const Geometry& geo) = default;
  Geometry& operator=(const Geometry& geo) = default;
  ~Geometry() = default;

  static Geometry* getInstance();
  static Geometry* getInstance(const std::string name);

  void init(const std::string geoFile);
  void init();
  void buildComposition();
  void setParameters(const std::string geoFile);
  void setParameters();

  // geometry helper functions
  // TODO: return tuples instead of using references
  std::tuple<double /*x*/, double /*y*/, double /*z*/> getGeoPadCenter(int tower, int layer, int stack, int row, int col) const;
  std::tuple<double /*x*/, double /*y*/> getGeoPadCenterLocal(int towerX, int towerY, int row, int col) const;
  std::tuple<double /*x*/, double /*y*/> getGeoPixCenterLocal(int row, int col) const;
  std::tuple<double /*x*/, double /*y*/, double /*z*/> getGeoPixelCenter(int pixel_id, int tower, int layer, int stack, int row, int col) const;
  std::tuple<double /*x*/, double /*y*/, double /*z*/> getGeoCompositionCenter(int tower, int layer, int stack) const;
  std::tuple<int /*row*/, int /*col*/, int /*stack*/, int /*layer*/, int /*seg*/, int /*waferx*/, int /*wafery*/> getPadPositionId2RowColStackLayer(int id) const;
  int getPixelNumber(int vol0, int vol1, int vol2, double x, double y) const;
  std::tuple<double /*x*/, double /*y*/, double /*z*/> getGeoTowerCenter(int tower, int segment = -1) const;
  bool disabledTower(int tower);
  std::tuple<bool /*statusCode*/, int /*col*/, int /*row*/, int /*layer*/, int /*segment*/> getVirtualInfo(double x, double y, double z) const;
  std::tuple<bool /*statusCode*/, double /*x*/, double /*y*/, double /*z*/> getXYZFromColRowSeg(int col, int row, int segment) const;
  std::tuple<bool /*statusCode*/, int /*nCol*/, int /*nRow*/> getVirtualNColRow(int segment) const;
  std::tuple<bool /*statusCode*/, int /*layer*/, int /*segment*/> getVirtualLayerSegment(float z) const;
  std::tuple<bool /*statusCode*/, int /*segment*/> getVirtualSegmentFromLayer(int layer) const;
  int getVirtualSegment(float z) const;

  // getters
  int getNumberOfPads() const { return mGlobal_PAD_NX * mGlobal_PAD_NY; }
  int getNumberOfPADsInX() const { return mGlobal_PAD_NX_Tower; }
  int getNumberOfPADsInY() const { return mGlobal_PAD_NY_Tower; }
  int getNumberOfPIXsInX() const { return mGlobal_PIX_NX_Tower; }
  int getNumberOfPIXsInY() const { return mGlobal_PIX_NY_Tower; }
  float getHCALTowerSize() const { return mGlobal_HCAL_Tower_Size; }
  int getHCALTowersInX() const { return mGlobal_HCAL_Tower_NX; }
  int getHCALTowersInY() const { return mGlobal_HCAL_Tower_NY; }
  int getNumberOfSegments() const { return mNumberOfSegments; } // NOTE: These are not the virtual segments, but just total number of layers as read from the geometry file. Need to disambiguate
  int getNumberOfPadLayers() const { return mNPadLayers; }
  int getNumberOfPixelLayers() const { return mNPixelLayers; }
  int getNumberOfHCalLayers() const { return mNHCalLayers; }
  int getNumberOfLayers() const { return mNPadLayers + mNPixelLayers + mNHCalLayers; }
  int getNumberOfLayerSeg() const { return mLayerSeg; }

  // ---- FOCAL-E pad indexing (global grid of real pads, see FOCALBase/PadIndex.h) ----
  int getNumberOfPadsInSensorX() const { return mGlobal_PAD_NX; }         ///< pads per sensor (wafer) in x (9)
  int getNumberOfPadsInSensorY() const { return mGlobal_PAD_NY; }         ///< pads per sensor (wafer) in y (8)
  int getNumberOfSensorsInModuleX() const { return mGlobal_PAD_NX_Tower; } ///< sensors per module (tower) in x (5)
  int getNumberOfSensorsInModuleY() const { return mGlobal_PAD_NY_Tower; } ///< sensors per module (tower) in y (1)
  /// \brief Pad index helper configured with the dimensions of this geometry
  ///
  /// The layer numbering is the geometric one (pads + pixels, 0..19); pixel layers never contain pad digits.
  PadIndex getPadIndex() const
  {
    return PadIndex(mNPadLayers + mNPixelLayers, mGlobal_Tower_NX, mGlobal_Tower_NY, mGlobal_PAD_NX_Tower, mGlobal_PAD_NY_Tower, mGlobal_PAD_NX, mGlobal_PAD_NY);
  }
  /// \brief Check whether an ECAL layer (geometric numbering, pads + pixels) is a pixel layer
  bool isPixelLayer(int layer) const;
  /// \brief Global pad (layer,row,col) from the volume copy numbers of the pad (EWSIPAD1) and of its module (EMSC1/EMSC2)
  /// \param padCopy   copy number of the pad volume = 1 + id + (stack << 12) + (layer << 16)
  /// \param moduleCopy copy number of the module volume = 1 + moduleY * nModulesX + moduleX
  /// \return {valid, layer (geometric numbering 0..19), row, col}
  std::tuple<bool /*valid*/, int /*layer*/, int /*row*/, int /*col*/> getPadRowColLayerFromCopyNumbers(int padCopy, int moduleCopy) const;

  double getFOCALSizeX() const;
  double getFOCALSizeY() const;
  double getTowerSize() const { return mTowerSizeX; }
  double getTowerSizeX() const;
  double getTowerSizeY() const;
  double getFOCALSizeZ() const;
  double getECALSizeZ() const;
  double getECALCenterZ() const;
  double getHCALSizeZ() const;
  double getHCALCenterZ() const;
  double getFOCALSegmentZ(int seg) const;
  double getFOCALZ0() const { return mGlobal_FOCAL_Z0; }
  int getNumberOfTowersInX() const { return mGlobal_Tower_NX; }
  int getNumberOfTowersInY() const { return mGlobal_Tower_NY; }
  double getTowerGapSize() const { return mGlobal_TOWER_TOLX; }
  double getTowerGapSizeX() const { return mGlobal_TOWER_TOLX; }
  double getTowerGapSizeY() const { return mGlobal_TOWER_TOLY; }
  double getGlobalPixelSize() const { return mGlobal_Pixel_Size; }
  double getGlobalPixelWaferSizeX() const { return mGlobal_PIX_SizeX; }
  double getGlobalPixelWaferSizeY() const { return mGlobal_PIX_SizeY; }
  double getGlobalPixelSkin() const { return mGlobal_PIX_SKIN; }
  double getGlobalPixelOffsetX() const { return mGlobal_PIX_OffsetX; }
  double getGlobalPadSize() const { return mGlobal_Pad_Size; }
  float getMiddleTowerOffset() const { return mGlobal_Middle_Tower_Offset; }
  bool getInsertFrontPadLayers() const { return mInsertFrontPadLayers; }
  bool getInsertHCalReadoutMaterial() const { return mInsertFrontHCalReadoutMaterial; }
  float getHCALPitchSize() const { return mGlobal_HCAL_Pitch_Size; }
  float getHCALBeamPipeHoleSize() const { return mGlobal_HCAL_BeamPipeHole_Size; }

  // ---- FOCAL-H, sheet design: dimensions (COMPOSITION_HCAL_S0/S1/S2 and GLOBAL_HCAL_* lines of the geometry file) ----
  /// One pitch = sheet with grooves + full sheet + sheet with grooves + full sheet, with one fibre row in each grooved sheet
  static constexpr int HCAL_FIBRE_ROWS_PER_PITCH = 2;
  const Composition* getHCALComposition(int stack) const; ///< HCAL composition COMPOSITION_HCAL_S<stack> of the geometry file, nullptr if not defined
  double getHCALPlateWidth() const;          ///< width of a full-width plate in x (COMPOSITION_HCAL_S0 X, 49.81 cm)
  double getHCALPlateLength() const;         ///< length of the plates and fibres in z (COMPOSITION_HCAL_S0 Z, 110 cm)
  double getHCALGroovedSheetThickness() const; ///< thickness of the sheet with grooves (COMPOSITION_HCAL_S0 Y, 0.20 cm)
  double getHCALFullSheetThickness() const;  ///< thickness of the full sheet (COMPOSITION_HCAL_S2 Y, 0.15 cm)
  double getHCALFibreRadius() const;         ///< fibre radius (COMPOSITION_HCAL_S1 X / 2, 0.05 cm)
  double getHCALPitchHeight() const { return HCAL_FIBRE_ROWS_PER_PITCH * (getHCALGroovedSheetThickness() + getHCALFullSheetThickness()); } ///< 0.7 cm
  float getHCALSizeY() const { return mGlobal_HCAL_SizeY; } ///< total height of FOCAL-H (GLOBAL_HCAL_SIZE_Y, default 100.8 cm)
  /// \brief Number of pitches over the full height = GLOBAL_HCAL_SIZE_Y / pitch height (144). The HCAL height is a parameter
  ///        of its own (not taken from the ECAL), so that ECAL descoping scenarios do not change FOCAL-H.
  int getHCALNumberOfPitches() const;
  int getHCALNumberOfPitchesBeamPipe() const;    ///< pitches beside the beam-pipe opening (GLOBAL_HCAL_BEAMPIPE / pitch height, 12)
  int getHCALNumberOfFibresPerRow() const;       ///< fibre positions in one row of a full-width plate = width / fibre pitch (124)
  int getHCALNumberOfFibreRows() const { return HCAL_FIBRE_ROWS_PER_PITCH * getHCALNumberOfPitches(); } ///< fibre rows over the full height (288)

  // ---- FOCAL-H SiPM indexing, sheet design (see FOCALBase/HCALIndex.h) ----
  int getHCALFibresPerSiPMX() const { return mGlobal_HCAL_FibresPerSiPM_X; } ///< fibres per SiPM bundle in x (default 4)
  int getHCALFibresPerSiPMY() const { return mGlobal_HCAL_FibresPerSiPM_Y; } ///< fibre rows per SiPM bundle in y (default 4)
  /// \brief SiPM index helper configured with the dimensions of this geometry
  HCALIndex getHCALIndex() const;
  /// \brief (layer, row, col) of the full cell index of a subsystem:
  ///        pads -> PadIndex; HCAL (sheets) -> HCALIndex, with layer = first HCAL layer (20, information only).
  ///        Pixels and the other HCAL designs have no index (use getVirtualInfo with the position): {-1, -1, -1}.
  std::tuple<int /*layer*/, int /*row*/, int /*col*/> getLayerRowColFromIndex(Hit::Subsystem_t subsystem, int index) const;
  /// \brief Full cell index of a subsystem from the (layer, index inside the layer) stored in a hit (see Hit.h):
  ///        pads -> PadIndex absID; HCAL (sheets) -> SiPM index (= index inside the layer).
  ///        Pixels and the other HCAL designs: -1 (no index, see getLayerRowColFromIndex). -1 if not valid.
  int getIndexFromLayerAndIndexInLayer(Hit::Subsystem_t subsystem, int layer, int indexInLayer) const;
  /// \brief SiPM (row, col) from the volume copy numbers of a fibre (sheet design, see Detector::CreatePitchAssembly)
  /// \param fibreCopy copy number of the fibre volume inside its pitch:
  ///                  fibreCopy = rowInPitch * nFibresPerRow + (fibre position counted from the beam-pipe end of a full-width plate)
  ///                  with rowInPitch = 0 (lower row) or 1 (upper row) and nFibresPerRow = getHCALNumberOfFibresPerRow().
  ///                  The shorter plates beside the beam pipe use the same positions (their fibres closest to the beam pipe
  ///                  are simply missing), so every row has nFibresPerRow positions and the copy number is unique in the pitch.
  /// \param pitchCopy copy number of its pitch assembly = pitch number counted from the bottom (0 .. getHCALNumberOfPitches()-1)
  /// \param halfCopy  copy number of its half: 0 = left (x < 0), 1 = right (x > 0)
  std::tuple<bool /*valid*/, int /*row*/, int /*col*/> getHCALRowColFromCopyNumbers(int fibreCopy, int pitchCopy, int halfCopy) const;

  float getDetectorOpeningRight() const { return mGlobal_DetectorOpening_Right; }
  float getDetectorOpeningLeft() const { return mGlobal_DetectorOpening_Left; }

  std::vector<const Composition*> getFOCALMicroModule(int layer) const;
  const Composition* getComposition(int layer, int stack) const;
  std::string_view getTowerGapMaterial() const { return mGlobal_Gap_Material; }

  int getVirtualNSegments() const;

  float getVirtualPadSize(int segment) const;
  float getVirtualRelativeSensitiveThickness(int segment) const;
  float getVirtualPixelTreshold(int segment) const;
  float getVirtualSegmentSizeZ(int segment) const;
  float getVirtualSegmentZ(int segment) const;
  bool getVirtualIsPixel(int segment) const;
  bool getVirtualIsHCal(int segment) const;
  int getVirtualNLayersInSegment(int segment) const;
  int getVirtualMinLayerInSegment(int segment) const;
  int getVirtualMaxLayerInSegment(int segment) const;

  void setUpLayerSegmentMap();
  void setUpTowerWaferSize();

  HCALDesgin getHCALDesign() const { return mHCALDesign; }

 protected:
  std::vector<Composition> mGeometryComposition;
  std::vector<Composition> mFrontMatterCompositionBase;
  std::vector<Composition> mPadCompositionBase;
  std::vector<Composition> mPixelCompositionBase;
  std::vector<Composition> mHCalCompositionBase;

  // PAD setup
  float mGlobal_Pad_Size = 0.0; // pad size
  int mGlobal_PAD_NX = 0;       // number of X pads in wafer
  int mGlobal_PAD_NY = 0;       // number of Y pads in wafer
  int mGlobal_PAD_NX_Tower = 0; // number of X wafers in tower
  int mGlobal_PAD_NY_Tower = 0; // number of Y wafers in tower
  float mGlobal_PPTOL = 0.0;    // tolerance between the wafers
  float mGlobal_PAD_SKIN = 0.0; // dead area (guard ring) on the wafer
  float mWaferSizeX = 0.0;      // Wafer X size
  float mWaferSizeY = 0.0;      // Wafer Y size

  float mGlobal_DetectorOpening_Right = 0.0; // detector opening in X
  float mGlobal_DetectorOpening_Left = 0.0;  // detector opening in Y

  // PIX setup
  float mGlobal_Pixel_Size = 0.0;  // pixel size
  float mGlobal_PIX_SizeX = 0.0;   // sensor size X
  float mGlobal_PIX_SizeY = 0.0;   // sensor size Y
  float mGlobal_PIX_OffsetX = 0.0; // offset for pixel layers in X
  float mGlobal_PIX_OffsetY = 0.0; // offset for pixel layers in Y
  float mGlobal_PIX_SKIN = 0.0;
  int mGlobal_PIX_NX_Tower = 0;       // number of sensors in X
  int mGlobal_PIX_NY_Tower = 0;       // number of sensors in Y
  bool mGlobal_Pixel_Readout = false; // readout on

  // Tower setup
  int mNPadLayers = 0;                      // total number of pad layers
  int mNPixelLayers = 0;                    // number of pixel layers
  std::array<int, 20> mPixelLayerLocations; // location of the pixel layers
  int mGlobal_Tower_NX = 0;                 // How many towers in X
  int mGlobal_Tower_NY = 0;                 // How many towers in Y
  float mTowerSizeX = 0.0;                  // X size of tower
  float mTowerSizeY = 0.0;                  // Y size of tower
  float mGlobal_TOWER_TOLX = 0.0;           // X - tolarance around tower
  float mGlobal_TOWER_TOLY = 0.0;           // Y - tolarance around tower
  float mGlobal_Middle_Tower_Offset = 0.0;  // if odd layers, the middle tower is offset due to the beampipe
  std::string mGlobal_Gap_Material;         // gap filling material NOTE: currently not used

  float mGlobal_HCAL_Tower_Size = 0.0;
  int mGlobal_HCAL_Tower_NX = 0;              // Number of HCAL towers on X
  int mGlobal_HCAL_Tower_NY = 0;              // Number of HCAL towers on Y
  float mGlobal_HCAL_Pitch_Size = 0.0;        // Distance between two fibers
  float mGlobal_HCAL_BeamPipeHole_Size = 0.0; // beam pipe hole size in HCAL
  float mGlobal_HCAL_SizeY = 100.8;           // total height of FOCAL-H (sheet design), independent of the ECAL size
  int mGlobal_HCAL_FibresPerSiPM_X = 4;       // fibres per SiPM bundle in x
  int mGlobal_HCAL_FibresPerSiPM_Y = 4;       // fibre rows per SiPM bundle in y
  HCALDesgin mHCALDesign = Sandwich;          // HCAL design type

  float mGlobal_FOCAL_Z0 = 0.0;

  bool mInsertFrontPadLayers = false;           // Have 2 pad layers in front of ECAL for charged particle veto
  bool mInsertFrontHCalReadoutMaterial = false; // if true, insert an 1cm thick aluminium layer at 2cm behind HCal to simulate the material introduced by the readout

  int mLayerSeg = 0;
  int mNHCalLayers = 0;                           // number of HCalLayers
  std::array<int, 100> mSegments;                 //  which layer belongs to which segment
  std::array<int, 100> mNumberOfLayersInSegments; // nymber of layers in each segment
  int mNumberOfSegments = 0;                      // number of long. segements
  int mNFrontMatterCompositionBase = 0;
  std::array<float, 100> mLocalLayerZ;    //// layer location in z
  std::array<float, 100> mLocalSegmentsZ; /// segment location in z
  float mFrontMatterLayerThickness = 0.0;
  float mPadLayerThickness = 0.0;
  float mPixelLayerThickness = 0.0;
  float mHCalLayerThickness = 0.0;
  std::array<float, 100> mLayerThickness; // thickenss of the layers
  std::list<int> mDisableTowers;

  int mVirtualNSegments = 0;
  std::vector<VirtualSegment> mVirtualSegmentComposition;

 private:
  static Geometry* sGeom;
  static bool sInit;
}; // Geometry

} // namespace focal
} // namespace o2
#endif
