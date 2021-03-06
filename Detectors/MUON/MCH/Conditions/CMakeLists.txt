# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

o2_add_library(
  MCHConditions
  SOURCES src/DCSNamer.cxx
  PUBLIC_LINK_LIBRARIES fmt::fmt)

if(BUILD_TESTING)
  o2_add_test(
    dcs-namer
    SOURCES test/testDCSNamer.cxx test/HVAliases.cxx test/LVAliases.cxx
    COMPONENT_NAME mch
    LABELS "muon;mch;dcs"
    PUBLIC_LINK_LIBRARIES O2::MCHConditions)
endif()

o2_add_executable(dcs-check-ccdb
                  COMPONENT_NAME mch
                  SOURCES src/dcs-check-ccdb.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsDCS O2::MCHConditions O2::CCDB)

o2_add_executable(mch-dcs-processor-workflow
                  COMPONENT_NAME calibration
                  SOURCES workflow/dcs-processor-workflow.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsDCS O2::MCHConditions)


o2_add_executable(mch-dcs-sim-workflow
                  COMPONENT_NAME calibration
                  SOURCES workflow/dcs-sim-workflow.cxx
                  PUBLIC_LINK_LIBRARIES O2::DCStestWorkflow O2::MCHConditions)

