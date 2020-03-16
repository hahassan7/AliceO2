// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file AnalysisCluster.cxx

#include "FairLogger.h"
#include <gsl/span>
#include <array>
#include <TLorentzVector.h>
#include "DataFormatsEMCAL/AnalysisCluster.h"

using namespace o2::emcal;

//_______________________________________________________________________
void AnalysisCluster::clear()
{
  //if(mTracksMatched) delete mTracksMatched;
  //mTracksMatched = 0;
  mCellsAmpFraction = 0x0;
  mCellsAbsId = 0x0;
}

//_______________________________________________________________________
TLorentzVector AnalysisCluster::getMomentum(std::array<const double, 3> vertex) const
{

  TLorentzVector p;

  Double32_t pos[3] = {mGlobalPos.X(), mGlobalPos.Y(), mGlobalPos.Z()};
  pos[0] -= vertex[0];
  pos[1] -= vertex[1];
  pos[2] -= vertex[2];

  double r = TMath::Sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);

  if (r > 0)
    p.SetPxPyPzE(mEnergy * pos[0] / r, mEnergy * pos[1] / r, mEnergy * pos[2] / r, mEnergy);
  else
    LOG(INFO) << "Null cluster radius, momentum calculation not possible";

  return p;
}

//_______________________________________________________________________
void AnalysisCluster::setCellsAbsId(gsl::span<unsigned short> array)
{
  mCellsAbsId = gsl::span<unsigned short>(&array[0], mNCells);
}

//_______________________________________________________________________
void AnalysisCluster::setCellsAmplitudeFraction(gsl::span<Double32_t> array)
{
  mCellsAmpFraction = gsl::span<Double32_t>(&array[0], mNCells);
}

//______________________________________________________________________________
void AnalysisCluster::setGlobalPosition(Point3D<double> x)
{
  mGlobalPos.SetX(x.X());
  mGlobalPos.SetY(x.Y());
  mGlobalPos.SetZ(x.Z());
}

//______________________________________________________________________________
void AnalysisCluster::setLocalPosition(Point3D<double> x)
{
  mLocalPos.SetX(x.X());
  mLocalPos.SetY(x.Y());
  mLocalPos.SetZ(x.Z());
}