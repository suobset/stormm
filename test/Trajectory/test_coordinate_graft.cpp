#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hpc_config.h"
#include "Accelerator/hybrid.h"
#include "Constants/fixed_precision.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "FileManagement/file_listing.h"
#include "Math/rounding.h"
#include "Math/series_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Random/random.h"
#include "Reporting/summary_file.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Trajectory/coordinate_copy.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/coordinate_graft.h"
#include "Trajectory/trajectory_enumerators.h"
#include "UnitTesting/approx.h"
#include "UnitTesting/test_environment.h"
#include "UnitTesting/test_system_manager.h"
#include "UnitTesting/unit_test.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::numerics;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Prepare the index map for a coordinate graft between two objects.
//
// Arguments:
//   xrs:               Source of random numbers to generate mapping assignments
//   orig_natom:        The number of atoms in the original object
//   make_dest_larger:  Indicate whether the new object will be larger or smaller
//   orig_tier:         Indicate whether to take data from the original object on the CPU host or
//                      GPU device
//   dest_tier:         Indicate whether to take data from the destination object on the CPU host
//                      or GPU device
//   dest_offset:       Optional offset to be applied to all destination indices, for mapping
//                      into objects with multiple systems
//   orig_offset:       Optional offset to be applied to all origin indices, for mapping out of an
//                      object with multipl esystems
//-------------------------------------------------------------------------------------------------
Hybrid<int2> prepCoordinateMap(Xoshiro256ppGenerator *xrs, const int orig_natom,
                               const bool make_dest_larger, const HybridTargetLevel orig_tier,
                               const HybridTargetLevel dest_tier, const int dest_offset = 0,
                               const int orig_offset = 0) {
  Hybrid<int2> result(orig_natom);
  result.resize(0);
  int dest_idx = 0;
  for (int j = 0; j < orig_natom; j++) {
    if (make_dest_larger) {

      // The destination object will be larger than the original.  Some elements of the
      // destination will be skipped.
      dest_idx += static_cast<int>(2.0 * xrs->uniformRandomNumber());
      result.pushBack({ dest_idx + dest_offset, j + orig_offset });
      dest_idx++;
    }
    else {

      // Every element of the destination will be covered by one of the transfers.
      if (xrs->uniformRandomNumber() < 0.5) {
        result.pushBack({ dest_idx + dest_offset, j + orig_offset });
        dest_idx++;
      }
    }
  }
  if (make_dest_larger == false) {
    result.resize(dest_idx);
  }
#ifdef STORMM_USE_HPC
  if (orig_tier == HybridTargetLevel::DEVICE || dest_tier == HybridTargetLevel::DEVICE) {
    result.upload();
  }
#endif
  return result;
}

//-------------------------------------------------------------------------------------------------
// Test grafts into CoordinateFrame objects.  Begin by making the transfer map, then design       
// CoordinateFrame objects of the appropriate size to receive the transfers.
//
// Arguments:
//   tsm:               The collection of test systems to study
//   xrs:               Source of random numbers needed to design the transfer map
//   dest_tier:         Indicate whether the destination object will receive the coordinates in
//                      memory on the CPU host or on the GPU device
//   orig_tier:         Indicate whether the original object will transmit coordinates from memory
//                      on the CPU host or the GPU device
//   make_dest_larger:  The destination object will be created larger than the original object,
//                      meaning that not all atoms of the destination object may be affected by
//                      the coordinate graft
//   gpos_scale_bits:   Scaling factor to place on global positions in the PhaseSpaceSynthesis
//                      tests
//   vel_scale_bits:    Scaling factor to place on velocities in the PhaseSpaceSynthesis tests
//   frc_scale_bits:    Scaling factor to place on force values in the PhaseSpaceSynthesis tests
//-------------------------------------------------------------------------------------------------
void transferToCoordinateFrame(const TestSystemManager &tsm, Xoshiro256ppGenerator *xrs,
                               const HybridTargetLevel dest_tier = HybridTargetLevel::HOST,
                               const HybridTargetLevel orig_tier = HybridTargetLevel::HOST,
                               const bool make_dest_larger = false, const int gpos_scale_bits = 40,
                               const int vel_scale_bits = default_velocity_scale_bits,
                               const int frc_scale_bits = default_force_scale_bits,
                               const GpuDetails gpu = null_gpu) {

  // Test CoordinateFrame ==> CoordinateFrame grafting
  for (size_t i = 0; i < tsm.getSystemCount(); i++) {
    CoordinateFrame orig = tsm.exportCoordinateFrame(i);
    const CoordinateFrameReader origr = orig.data();

    // Create the map, to be used for all transfers
    const Hybrid<int2> crd_map = prepCoordinateMap(xrs, origr.natom, make_dest_larger, orig_tier,
                                                   dest_tier);
#ifdef STORMM_USE_HPC
    if (orig_tier == HybridTargetLevel::DEVICE) {
      orig.upload();
    }
#endif
    // Create a CoordinateFrame destination object and check the transfer
    const int dest_natom =  (make_dest_larger) ? crd_map.readHost(crd_map.size() - 1).x + 1 :
                                                 crd_map.size();
    CoordinateFrame dest(dest_natom);
    coordGraft(&dest, orig, crd_map, dest_tier, orig_tier, gpu);
    const std::vector<double> dest_holdings = dest.getInterlacedCoordinates(dest_tier);
    const std::vector<double> orig_holdings = orig.getInterlacedCoordinates(orig_tier);
    const int ngraft = crd_map.size();
    std::vector<double> orig_x(ngraft), orig_y(ngraft), orig_z(ngraft);
    std::vector<double> dest_x(ngraft), dest_y(ngraft), dest_z(ngraft);
    for (int j = 0; j < ngraft; j++) {
      const size_t dest_idx = crd_map.readHost(j).x;
      const size_t orig_idx = crd_map.readHost(j).y;
      orig_x[j] = orig_holdings[(3 * orig_idx)    ];
      orig_y[j] = orig_holdings[(3 * orig_idx) + 1];
      orig_z[j] = orig_holdings[(3 * orig_idx) + 2];
      dest_x[j] = dest_holdings[(3 * dest_idx)    ];
      dest_y[j] = dest_holdings[(3 * dest_idx) + 1];
      dest_z[j] = dest_holdings[(3 * dest_idx) + 2];
    }
    check(dest_x, RelationalOperator::EQUAL, orig_x, "Cartesian X coordinates were transferred "
          "incorrectly in a CoordinateFrame ==> CoordinateFrame graft.", tsm.getTestingStatus());
    check(dest_y, RelationalOperator::EQUAL, orig_y, "Cartesian Y coordinates were transferred "
          "incorrectly in a CoordinateFrame ==> CoordinateFrame graft.", tsm.getTestingStatus());
    check(dest_z, RelationalOperator::EQUAL, orig_z, "Cartesian Z coordinates were transferred "
          "incorrectly in a CoordinateFrame ==> CoordinateFrame graft.", tsm.getTestingStatus());
  }

  // Test PhaseSpace ==> CoordinateFrame grafting
  const std::vector<TrajectoryKind> tkinds = { TrajectoryKind::POSITIONS,
                                               TrajectoryKind::VELOCITIES,
                                               TrajectoryKind::FORCES };
  const std::vector<CoordinateCycle> tcyc = { CoordinateCycle::WHITE, CoordinateCycle::BLACK };
  for (size_t i = 0; i < tsm.getSystemCount(); i++) {
    PhaseSpace orig = tsm.exportPhaseSpace(i);
    PhaseSpaceWriter origw = orig.data(CoordinateCycle::WHITE);
    for (size_t j = 0; j < tkinds.size(); j++) {
      for (size_t k = 0; k < tcyc.size(); k++) {

        // Populate the force and velocities with random non-trivial numbers, if necessary.
        switch (tkinds[j]) {
        case TrajectoryKind::POSITIONS:
          break;
        case TrajectoryKind::VELOCITIES:
          switch (tcyc[k]) {
          case CoordinateCycle::WHITE:
            for (int m = 0; m < origw.natom; m++) {
              origw.xvel[m] = 0.25 * xrs->gaussianRandomNumber();
              origw.yvel[m] = 0.25 * xrs->gaussianRandomNumber();
              origw.zvel[m] = 0.25 * xrs->gaussianRandomNumber();
            }
            break;
          case CoordinateCycle::BLACK:
            for (int m = 0; m < origw.natom; m++) {
              origw.vxalt[m] = 0.25 * xrs->gaussianRandomNumber();
              origw.vyalt[m] = 0.25 * xrs->gaussianRandomNumber();
              origw.vzalt[m] = 0.25 * xrs->gaussianRandomNumber();
            }
            break;
          }
          break;
        case TrajectoryKind::FORCES:
          switch (tcyc[k]) {
          case CoordinateCycle::WHITE:
            for (int m = 0; m < origw.natom; m++) {
              origw.xfrc[m] = 3.0 * xrs->gaussianRandomNumber();
              origw.yfrc[m] = 3.0 * xrs->gaussianRandomNumber();
              origw.zfrc[m] = 3.0 * xrs->gaussianRandomNumber();
            }
            break;
          case CoordinateCycle::BLACK:
            for (int m = 0; m < origw.natom; m++) {
              origw.fxalt[m] = 3.0 * xrs->gaussianRandomNumber();
              origw.fyalt[m] = 3.0 * xrs->gaussianRandomNumber();
              origw.fzalt[m] = 3.0 * xrs->gaussianRandomNumber();
            }
            break;
          }
          break;
        }

        // Create the map, to be used for all transfers
        const Hybrid<int2> crd_map = prepCoordinateMap(xrs, origw.natom, make_dest_larger,
                                                       orig_tier, dest_tier);
#ifdef STORMM_USE_HPC
        if (orig_tier == HybridTargetLevel::DEVICE) {
          orig.upload();
        }
#endif
        // Create a CoordinateFrame destination object and check the transfer
        const int dest_natom =  (make_dest_larger) ? crd_map.readHost(crd_map.size() - 1).x + 1 :
                                                     crd_map.size();
        CoordinateFrame dest(dest_natom);
        coordGraft(&dest, orig, crd_map, tkinds[j], tcyc[k], dest_tier, orig_tier, gpu);
        const std::vector<double> dest_holdings = dest.getInterlacedCoordinates(dest_tier);
        const std::vector<double> orig_holdings = orig.getInterlacedCoordinates(tcyc[k], tkinds[j],
                                                                                orig_tier);
        const int ngraft = crd_map.size();
        std::vector<double> orig_x(ngraft), orig_y(ngraft), orig_z(ngraft);
        std::vector<double> dest_x(ngraft), dest_y(ngraft), dest_z(ngraft);
        for (int m = 0; m < ngraft; m++) {
          const size_t dest_idx = crd_map.readHost(m).x;
          const size_t orig_idx = crd_map.readHost(m).y;
          orig_x[m] = orig_holdings[(3 * orig_idx)    ];
          orig_y[m] = orig_holdings[(3 * orig_idx) + 1];
          orig_z[m] = orig_holdings[(3 * orig_idx) + 2];
          dest_x[m] = dest_holdings[(3 * dest_idx)    ];
          dest_y[m] = dest_holdings[(3 * dest_idx) + 1];
          dest_z[m] = dest_holdings[(3 * dest_idx) + 2];
        }
        check(dest_x, RelationalOperator::EQUAL, orig_x, "Cartesian X coordinates were "
              "transferred incorrectly in a PhaseSpace (" +
              getEnumerationName(tkinds[j]) + ", " + getEnumerationName(tcyc[k]) + " ) ==> "
              "CoordinateFrame graft.", tsm.getTestingStatus());
        check(dest_y, RelationalOperator::EQUAL, orig_y, "Cartesian Y coordinates were "
              "transferred incorrectly in a PhaseSpace (" +
              getEnumerationName(tkinds[k]) + ", " + getEnumerationName(tcyc[k]) + " ) ==> "
              "CoordinateFrame graft.", tsm.getTestingStatus());
        check(dest_z, RelationalOperator::EQUAL, orig_z, "Cartesian Z coordinates were "
              "transferred incorrectly in a PhaseSpace (" +
              getEnumerationName(tkinds[j]) + ", " + getEnumerationName(tcyc[k]) + " ) ==> "
              "CoordinateFrame graft.", tsm.getTestingStatus());
      }
    }
  }

  // Test CoordinateSeries ==> CoordinateFrame grafting
  for (size_t i = 0; i < tsm.getSystemCount(); i++) {
    CoordinateFrame tmp = tsm.exportCoordinateFrame(i);
    CoordinateSeries<float> orig(tmp, 12);
    CoordinateSeriesWriter<float> origw = orig.data();
    const int padded_frame_size = roundUp(orig.getAtomCount(), warp_size_int);
    for (int j = 0; j < orig.getFrameCount(); j++) {
      for (int k = 0; k < orig.getAtomCount(); k++) {
        origw.xcrd[(j * padded_frame_size) + k] += (0.25 * xrs->gaussianRandomNumber());
        origw.ycrd[(j * padded_frame_size) + k] += (0.25 * xrs->gaussianRandomNumber());
        origw.zcrd[(j * padded_frame_size) + k] += (0.25 * xrs->gaussianRandomNumber());
      }
    }
    
    // Create the map, to be used for all transfers
    const int choice_frm = 12.0 * xrs->uniformRandomNumber();
    Hybrid<int2> crd_map = prepCoordinateMap(xrs, orig.getAtomCount(), make_dest_larger,
                                             orig_tier, dest_tier, 0,
                                             choice_frm * padded_frame_size);
#ifdef STORMM_USE_HPC
    if (orig_tier == HybridTargetLevel::DEVICE) {
      orig.upload();
    }
#endif
    // Create a CoordinateFrame destination object and check the transfer
    const int dest_natom =  (make_dest_larger) ? crd_map.readHost(crd_map.size() - 1).x + 1 :
                                                 crd_map.size();
    CoordinateFrame dest(dest_natom);
    coordGraft(&dest, orig, crd_map, dest_tier, orig_tier, gpu);
    const std::vector<double> dest_holdings = dest.getInterlacedCoordinates(dest_tier);
    const std::vector<double> orig_holdings = orig.getInterlacedCoordinates<double>(choice_frm, -1,
                                                                                    orig_tier);
    const int ngraft = crd_map.size();
    std::vector<double> orig_x(ngraft), orig_y(ngraft), orig_z(ngraft);
    std::vector<double> dest_x(ngraft), dest_y(ngraft), dest_z(ngraft);
    for (int j = 0; j < ngraft; j++) {
      const size_t dest_idx = crd_map.readHost(j).x;
      const size_t orig_idx = crd_map.readHost(j).y - (choice_frm * padded_frame_size);
      orig_x[j] = orig_holdings[(3 * orig_idx)    ];
      orig_y[j] = orig_holdings[(3 * orig_idx) + 1];
      orig_z[j] = orig_holdings[(3 * orig_idx) + 2];
      dest_x[j] = dest_holdings[(3 * dest_idx)    ];
      dest_y[j] = dest_holdings[(3 * dest_idx) + 1];
      dest_z[j] = dest_holdings[(3 * dest_idx) + 2];
    }
    check(dest_x, RelationalOperator::EQUAL, Approx(orig_x).margin(1.0e-5), "Cartesian X "
          "coordinates were transferred incorrectly in a CoordinateSeries ==> CoordinateFrame "
          "graft.", tsm.getTestingStatus());
    check(dest_y, RelationalOperator::EQUAL, Approx(orig_y).margin(1.0e-5), "Cartesian Y "
          "coordinates were transferred incorrectly in a CoordinateSeries ==> CoordinateFrame "
          "graft.", tsm.getTestingStatus());
    check(dest_z, RelationalOperator::EQUAL, Approx(orig_z).margin(1.0e-5), "Cartesian Z "
          "coordinates were transferred incorrectly in a CoordinateSeries ==> CoordinateFrame "
          "graft.", tsm.getTestingStatus());
  }

  // Test PhaseSpaceSynthesis ==> CoordinateFrame grafting
  PhaseSpaceSynthesis poly_ps =
    tsm.exportPhaseSpaceSynthesis(incrementingSeries(0, tsm.getSystemCount()), 0.05, 716389,
                                  gpos_scale_bits, vel_scale_bits, frc_scale_bits);
  PsSynthesisWriter poly_psw = poly_ps.data(CoordinateCycle::WHITE);
  const double vscl = poly_psw.vel_scale;
  const double fscl = poly_psw.frc_scale;
  for (size_t i = 0; i < tsm.getSystemCount(); i++) {
    for (size_t j = 0; j < tkinds.size(); j++) {
      for (size_t k = 0; k < tcyc.size(); k++) {

        // Populate the force and velocities with random non-trivial numbers, if necessary.
        switch (tkinds[j]) {
        case TrajectoryKind::POSITIONS:
          break;
        case TrajectoryKind::VELOCITIES:
          {
            const int95_t tmp_xv = hostDoubleToInt95(xrs->gaussianRandomNumber() * vscl);
            const int95_t tmp_yv = hostDoubleToInt95(xrs->gaussianRandomNumber() * vscl);
            const int95_t tmp_zv = hostDoubleToInt95(xrs->gaussianRandomNumber() * vscl);
            switch (tcyc[k]) {
            case CoordinateCycle::WHITE:
              for (int m = 0; m < poly_psw.atom_counts[i]; m++) {
                poly_psw.xvel[m] = tmp_xv.x;
                poly_psw.yvel[m] = tmp_yv.x;
                poly_psw.zvel[m] = tmp_zv.x;
                poly_psw.xvel_ovrf[m] = tmp_xv.y;
                poly_psw.yvel_ovrf[m] = tmp_yv.y;
                poly_psw.zvel_ovrf[m] = tmp_zv.y;
              }
              break;
            case CoordinateCycle::BLACK:
              for (int m = 0; m < poly_psw.atom_counts[i]; m++) {
                poly_psw.vxalt[m] = tmp_xv.x;
                poly_psw.vyalt[m] = tmp_yv.x;
                poly_psw.vzalt[m] = tmp_zv.x;
                poly_psw.vxalt_ovrf[m] = tmp_xv.y;
                poly_psw.vyalt_ovrf[m] = tmp_yv.y;
                poly_psw.vzalt_ovrf[m] = tmp_zv.y;
              }
              break;
            }
          }
          break;
        case TrajectoryKind::FORCES:
          {
            const int95_t tmp_xf = hostDoubleToInt95(xrs->gaussianRandomNumber() * fscl);
            const int95_t tmp_yf = hostDoubleToInt95(xrs->gaussianRandomNumber() * fscl);
            const int95_t tmp_zf = hostDoubleToInt95(xrs->gaussianRandomNumber() * fscl);
            switch (tcyc[k]) {
            case CoordinateCycle::WHITE:
              for (int m = 0; m < poly_psw.atom_counts[i]; m++) {
                poly_psw.xfrc[m] = tmp_xf.x;
                poly_psw.yfrc[m] = tmp_yf.x;
                poly_psw.zfrc[m] = tmp_zf.x;
                poly_psw.xfrc_ovrf[m] = tmp_xf.y;
                poly_psw.yfrc_ovrf[m] = tmp_yf.y;
                poly_psw.zfrc_ovrf[m] = tmp_zf.y;
              }
              break;
            case CoordinateCycle::BLACK:
              for (int m = 0; m < poly_psw.atom_counts[i]; m++) {
                poly_psw.fxalt[m] = tmp_xf.x;
                poly_psw.fyalt[m] = tmp_yf.x;
                poly_psw.fzalt[m] = tmp_zf.x;
                poly_psw.fxalt_ovrf[m] = tmp_xf.y;
                poly_psw.fyalt_ovrf[m] = tmp_yf.y;
                poly_psw.fzalt_ovrf[m] = tmp_zf.y;
              }
              break;
            }
          }
          break;
        }

        // Create the map, to be used for all transfers
        const Hybrid<int2> crd_map = prepCoordinateMap(xrs, poly_psw.atom_counts[i],
                                                       make_dest_larger, orig_tier, dest_tier, 0,
                                                       poly_psw.atom_starts[i]);
#ifdef STORMM_USE_HPC
        if (orig_tier == HybridTargetLevel::DEVICE) {
          poly_ps.upload();
        }
#endif
        // Create a CoordinateFrame destination object and check the transfer
        const int dest_natom =  (make_dest_larger) ? crd_map.readHost(crd_map.size() - 1).x + 1 :
                                                     crd_map.size();
        CoordinateFrame dest(dest_natom);
        coordGraft(&dest, poly_ps, crd_map, tkinds[j], tcyc[k], dest_tier, orig_tier, gpu);
        const std::vector<double> dest_holdings = dest.getInterlacedCoordinates(dest_tier);
        const std::vector<double> orig_holdings = poly_ps.getInterlacedCoordinates(i, tcyc[k],
                                                                                   tkinds[j],
                                                                                   orig_tier);
        const int ngraft = crd_map.size();
        std::vector<double> orig_x(ngraft), orig_y(ngraft), orig_z(ngraft);
        std::vector<double> dest_x(ngraft), dest_y(ngraft), dest_z(ngraft);
        for (int m = 0; m < ngraft; m++) {
          const size_t dest_idx = crd_map.readHost(m).x;
          const size_t orig_idx = crd_map.readHost(m).y - poly_ps.getAtomOffset(i);
          orig_x[m] = orig_holdings[(3 * orig_idx)    ];
          orig_y[m] = orig_holdings[(3 * orig_idx) + 1];
          orig_z[m] = orig_holdings[(3 * orig_idx) + 2];
          dest_x[m] = dest_holdings[(3 * dest_idx)    ];
          dest_y[m] = dest_holdings[(3 * dest_idx) + 1];
          dest_z[m] = dest_holdings[(3 * dest_idx) + 2];
        }
        check(dest_x, RelationalOperator::EQUAL, orig_x, "Cartesian X coordinates were "
              "transferred incorrectly in a PhaseSpaceSynthesis (" +
              getEnumerationName(tkinds[j]) + ", " + getEnumerationName(tcyc[k]) + " ) ==> "
              "CoordinateFrame graft.", tsm.getTestingStatus());
        check(dest_y, RelationalOperator::EQUAL, orig_y, "Cartesian Y coordinates were "
              "transferred incorrectly in a PhaseSpaceSynthesis (" +
              getEnumerationName(tkinds[k]) + ", " + getEnumerationName(tcyc[k]) + " ) ==> "
              "CoordinateFrame graft.", tsm.getTestingStatus());
        check(dest_z, RelationalOperator::EQUAL, orig_z, "Cartesian Z coordinates were "
              "transferred incorrectly in a PhaseSpaceSynthesis (" +
              getEnumerationName(tkinds[j]) + ", " + getEnumerationName(tcyc[k]) + " ) ==> "
              "CoordinateFrame graft.", tsm.getTestingStatus());
      }
    }
  }

  // Test Condensate ==> CoordinateFrame grafting
  const std::vector<PrecisionModel> all_prec = { PrecisionModel::DOUBLE, PrecisionModel::SINGLE };
  for (size_t i = 0; i < all_prec.size(); i++) {
    Condensate cdns(poly_ps, PrecisionModel::SINGLE);
#ifdef STORMM_USE_HPC
    cdns.upload();
#endif
    const CondensateReader cdnsr = cdns.data();
    for (int j = 0; j < cdnsr.system_count; j++) {

      // Create the map, to be used for all transfers
      const Hybrid<int2> crd_map = prepCoordinateMap(xrs, cdnsr.atom_counts[j], make_dest_larger,
                                                     orig_tier, dest_tier, 0,
                                                     cdnsr.atom_starts[j]);
#ifdef STORMM_USE_HPC
      if (orig_tier == HybridTargetLevel::DEVICE) {
        cdns.upload();
      }
#endif
      // Create a CoordinateFrame destination object and check the transfer
      const int dest_natom =  (make_dest_larger) ? crd_map.readHost(crd_map.size() - 1).x + 1 :
                                                   crd_map.size();
      CoordinateFrame dest(dest_natom);
      coordGraft(&dest, cdns, crd_map, dest_tier, orig_tier, gpu);
      const std::vector<double> dest_holdings = dest.getInterlacedCoordinates(dest_tier);
      const std::vector<double> orig_holdings = cdns.getInterlacedCoordinates(j, orig_tier);
      const int ngraft = crd_map.size();
      std::vector<double> orig_x(ngraft), orig_y(ngraft), orig_z(ngraft);
      std::vector<double> dest_x(ngraft), dest_y(ngraft), dest_z(ngraft);
      for (int k = 0; k < ngraft; k++) {
        const size_t dest_idx = crd_map.readHost(k).x;
        const size_t orig_idx = crd_map.readHost(k).y - cdns.getAtomOffset(j);
        orig_x[k] = orig_holdings[(3 * orig_idx)    ];
        orig_y[k] = orig_holdings[(3 * orig_idx) + 1];
        orig_z[k] = orig_holdings[(3 * orig_idx) + 2];
        dest_x[k] = dest_holdings[(3 * dest_idx)    ];
        dest_y[k] = dest_holdings[(3 * dest_idx) + 1];
        dest_z[k] = dest_holdings[(3 * dest_idx) + 2];
      }
      check(dest_x, RelationalOperator::EQUAL, orig_x, "Cartesian X coordinates were transferred "
            "incorrectly in a Condensate (" + getEnumerationName(cdns.getMode()) +
            ") ==> CoordinateFrame graft.", tsm.getTestingStatus());
      check(dest_y, RelationalOperator::EQUAL, orig_y, "Cartesian Y coordinates were transferred "
            "incorrectly in a Condensate (" + getEnumerationName(cdns.getMode()) +
            ") ==> CoordinateFrame graft.", tsm.getTestingStatus());
      check(dest_z, RelationalOperator::EQUAL, orig_z, "Cartesian Z coordinates were transferred "
            "incorrectly in a Condensate (" + getEnumerationName(cdns.getMode()) +
            ") ==> CoordinateFrame graft.", tsm.getTestingStatus());
    }
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
#else
  const GpuDetails gpu = null_gpu;
#endif

  // Section 1
  section("Grafts into CoordinateFrame objects");

  // Section 2
  section("Grafts into PhaseSpace objects");

  // Section 3
  section("Grafts into CoordinateSeries objects");

  // Section 4
  section("Grafts into PhaseSpaceSynthesis objects");

  // Section 5
  section("Grafts into Condensate objects");

  // Read in systems.  Make perturbations and, if HPC mode is engaged, upload the coordinates.
  const std::vector<std::string> iso_systems = { "drug_example_iso", "med_1", "med_2", "med_3",
                                                 "med_4", "med_5", "dhfr_cmap", "trpcage" };
  const std::vector<std::string> pbc_systems = { "bromobenzene", "drug_example_vs", "tamavidin",
                                                 "symmetry_C1_in_water", "symmetry_C6_in_water",
                                                 "tip3p", "tip4p" };
  const char osc = osSeparator();
  const std::string top_dir = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string crd_dir = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  TestSystemManager tsm_iso(top_dir, "top", iso_systems, crd_dir, "inpcrd", iso_systems);
  TestSystemManager tsm_pbc(top_dir, "top", pbc_systems, crd_dir, "inpcrd", pbc_systems);
  Xoshiro256ppGenerator xrs;

  // Test grafts into CoordinateFrame objects.  Begin by making the transfer map, then design
  // CoordinateFrame objects of the appropriate size to receive the transfers.
  const std::vector<HybridTargetLevel> tiers = { HybridTargetLevel::HOST /*, };
                                                 HybridTargetLevel::DEVICE */ };
  for (size_t i = 0; i < tiers.size(); i++) {
    for (size_t j = 0; j < tiers.size(); j++) {
      transferToCoordinateFrame(tsm_iso, &xrs, tiers[i], tiers[j], true,
                                globalpos_scale_nonoverflow_bits - 1,
                                velocity_scale_nonoverflow_bits - 1,
                                force_scale_nonoverflow_bits - 1, gpu);
      transferToCoordinateFrame(tsm_iso, &xrs, tiers[i], tiers[j], false,
                                globalpos_scale_nonoverflow_bits + 1,
                                velocity_scale_nonoverflow_bits + 1,
                                force_scale_nonoverflow_bits + 1, gpu);
    }
  }

  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
