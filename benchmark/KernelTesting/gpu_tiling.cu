#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <nvml.h>
#include <string>
#include <vector>
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/Constants/scaling.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Namelists/command_line_parser.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/hpc_valence_potential.h"
#include "../../src/Potential/hpc_nonbonded_potential.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Synthesis/implicit_solvent_workspace.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/static_mask_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/nonbonded_workunit.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "src/Accelerator/hybrid.h"
#include "src/Namelists/namelist_emulator.h"
#include "src/UnitTesting/test_environment.h"
#include "src/UnitTesting/unit_test_enumerators.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::diskutil;
using namespace stormm::errors;
using namespace stormm::energy;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::numerics;
using namespace stormm::parse;
using namespace stormm::restraints;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// GPU Kernel to compute interactions within a tile
//-------------------------------------------------------------------------------------------------
__global__ void computeTileInteractions (PsSynthesisWriter poly_psw,
                                         const SyNonbondedKit<float, float2> poly_nbk,
                                         const int4 *t_assgn, const int num_tiles,
                                         const LocalExclusionMask lem) {
  __shared__ volatile int total_warps;

  // Get the total number of warps (blocks * threads)
  if(threadIdx.x == 0){
    total_warps = (blockDim.x >> warp_bits) * gridDim.x;
  }
  __syncthreads();
  
  const int warp_idx = (blockDim.x >> warp_bits) * blockIdx.x + (threadIdx.x >> warp_bits);  // Warp index
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);  // Lane within the warp

  // Iterate over the tiles, each block handles multiple tiles
  for(int i = warp_idx; i < num_tiles; i += total_warps) {
    int4 instr = t_assgn[i];  // Instruction to fetch atom counts for current tile
    const int atom_start = instr.x;
    const int atom_count = instr.y;
    
    // Interactions between atoms within the tile
    for (int j = atom_start + lane_idx; j < atom_start + atom_count; j++) {
      // Reading atom coordinates and charges for atom j
      const double atom_jx = hostInt95ToDouble(poly_psw.xcrd[j], poly_psw.xcrd_ovrf[j]);
      const double atom_jy = hostInt95ToDouble(poly_psw.ycrd[j], poly_psw.ycrd_ovrf[j]);
      const double atom_jz = hostInt95ToDouble(poly_psw.zcrd[j], poly_psw.zcrd_ovrf[j]);
      const double atom_jq = poly_nbk.charge[j] * poly_nbk.coulomb;

      const int ljidx_j = poly_nbk.lj_idx[j];
    }
  }
}

//-------------------------------------------------------------------------------------------------
// C++ Function to compute interactions within a system
//-------------------------------------------------------------------------------------------------
void computeTileInteractions (PsSynthesisWriter *poly_psw, 
                              const SyNonbondedKit<float, float2> &poly_nbk,
                              const LocalExclusionMaskReader &lemr) {
  for(int i = 0; i < poly_psw->system_count; ++i) {
    const int atom_count = poly_psw->atom_counts[i];
    const int atom_start = poly_psw->atom_starts[i];
    const int ljabc_start = poly_nbk.ljabc_offsets[i];
    const int nlj_types = poly_nbk.n_lj_types[i];
    
    for(int j = atom_start; j < atom_start + atom_count; ++j){
      const double atom_jx = hostInt95ToDouble(poly_psw->xcrd[j], poly_psw->xcrd_ovrf[j]);
      const double atom_jy = hostInt95ToDouble(poly_psw->ycrd[j], poly_psw->ycrd_ovrf[j]);
      const double atom_jz = hostInt95ToDouble(poly_psw->zcrd[j], poly_psw->zcrd_ovrf[j]);
      const double atom_jq = poly_nbk.charge[j] * poly_nbk.coulomb;

      const int ljidx_j = poly_nbk.lj_idx[j];

      for(int k = atom_start; k < j; ++k){
        double fmag;
        if(testExclusion(lemr, j, k)) {
          fmag = 0.0;
        }
        else {
          const double atom_kx = hostInt95ToDouble(poly_psw->xcrd[k], poly_psw->xcrd_ovrf[k]);
          const double atom_ky = hostInt95ToDouble(poly_psw->ycrd[k], poly_psw->ycrd_ovrf[k]);
          const double atom_kz = hostInt95ToDouble(poly_psw->zcrd[k], poly_psw->zcrd_ovrf[k]);
          const double atom_kq = poly_nbk.charge[k] * poly_nbk.coulomb;

          const int ljidx_k = poly_nbk.lj_idx[k];
          const int ljparm_jk = nlj_types * ljidx_j + ljidx_k + ljabc_start;
          const float2 ljab_jk = poly_nbk.ljab_coeff[ljparm_jk];

          const double dx = atom_kx - atom_jx;
          const double dy = atom_ky - atom_jy;
          const double dz = atom_kz - atom_jz;
          double distance_sq = dx*dx + dy*dy + dz*dz;
          double distance = sqrt(distance_sq);

          double fmag = -(atom_jq * atom_kq) / (distance * distance_sq);
          const double inv_distance = 1.0 / distance;
          const double inv_distance_sq = inv_distance * inv_distance;
          const double inv_distance_4 = inv_distance_sq * inv_distance_sq;
          const double inv_distance_6 = inv_distance_sq * inv_distance_4;
          const double inv_distance_8 = inv_distance_4  * inv_distance_4;

          fmag += ((6.0 * ljab_jk.y) - (12.0 * ljab_jk.x * inv_distance_6)) * inv_distance_8;
        }

        fmag *= poly_psw->frc_scale;
        const int95_t force_x = (fmag * dx);
        const int95_t force_y = (fmag * dy);
        const int95_t force_z = (fmag * dz);
        

        // Update the forces on atom j and atom k
        hostSplitFPSum(poly_psw->xfrc[j], poly_psw->xfrc_ovrf[j], force_x);
        hostSplitFPSum(poly_psw->yfrc[j], poly_psw->yfrc_ovrf[j], force_y);
        hostSplitFPSum(poly_psw->zfrc[j], poly_psw->zfrc_ovrf[j], force_z);

        hostSplitFPSubtract(poly_psw->xfrc[k], poly_psw->xfrc_ovrf[k], force_x);
        hostSplitFPSubtract(poly_psw->yfrc[k], poly_psw->yfrc_ovrf[k], force_y);
        hostSplitFPSubtract(poly_psw->zfrc[k], poly_psw->zfrc_ovrf[k], force_z);
        
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Main function
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // baseline
  StopWatch timer;
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> engage_gpu(1);

  // Program Specific Command Line Inputs
  CommandLineParser clip("gpu_tiling", "A benchmarking program to test GPU tile strategies.");
  clip.addStandardAmberInputs("-p", "-c", "-ig_seed");
  NamelistEmulator *t_nml = clip.getNamelistPointer();

  // Initialize the testing environment
  TestEnvironment oe(argc, argv, &clip, TmpdirStatus::NOT_REQUIRED, ExceptionResponse::SILENT);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Read command line instructions
  clip.parseUserInput(argc, argv);
  const std::string inpcrd_file = t_nml->getStringValue("-c");
  const std::string topology_file = t_nml->getStringValue("-p");
  
  // Check that a system has been provided and load the basic topology and input coordinates
#ifdef STORMM_USE_HPC
  TestSystemManager tsm(std::vector<std::string>(1, topology_file),
                        std::vector<std::string>(1, inpcrd_file), ExceptionResponse::DIE);
  timer.assignTime(time_load);

  // Create basic objects and load system data
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(std::vector<int>(1, 0),
                                                              0.0, ig_seed, 40, 44, fp_bits);

  int total_tile_count = 0;
  for(int i = 0; i < poly_ps.getSystemCount(); ++i){
    PhaseSpace ps = poly_ps.exportSystem(i);
    PhaseSpaceWriter psw = ps.data();
    imageCoordinates(psw.xcrd, psw.ycrd, psw.zcrd, psw.natom, psw.umat, psw.invu, 
                     ImagingMethod::PRIMARY_UNIT_CELL);
    poly_ps.importSystem(ps, i);
    int atom_counts = psw.natom;
    int tiles = (atom_counts + 31) / 32;
    total_tile_count += (tiles * (tiles + 1)) / 2;
  }

  // x: Atom Start for the receiving atoms
  // y: Atom Start for the sending atoms
  // z: Number of receiving atoms
  // w: Number of sending atoms
  Hybrid<int4> tile_assignments(total_tile_count, "tile_assignments");

  int counter = 0;
  for(int i = 0; i < poly_ps.getSystemCount(); ++i){
    PhaseSpace ps = poly_ps.exportSystem(i);
    PhaseSpaceWriter psw = ps.data();

    // Range over all tiles in the current system
    int atom_start = psw.atom_starts[i];
    int atom_count = psw.natom;
    int tiles = (atom_count + 31) / 32;

    // Assign tile assignments for each tile in the system
    for (int j = 0; j < tiles; ++j) {
      int start_recv = atom_start + (j * 32);  // Starting index for receiving atoms
      int end_recv = std::min(atom_start + (j + 1) * 32, atom_start + atom_count);  
      int num_recv = end_recv - start_recv;  // Number of receiving atoms

      for (int k = j; k < tiles; ++k) {
        int start_send = atom_start + (k * 32);  // Starting index for sending atoms
        int end_send = std::min(atom_start + (k + 1) * 32, atom_start + atom_count);  
        int num_send = end_send - start_send;  // Number of sending atoms

        // Fill the tile assignments
        tile_assignments[counter] = Hybrid<int4>(start_recv, start_send, num_recv, num_send);
        counter++;
      }
    }
  }
  PsSynthesisWriter host_psw = poly_ps.data();
  Xoshiro256ppGenerator xrs(ig_seed);
  for (int i = 1; i < host_psw.system_count; i++) {
    const size_t aoff = host_psw.atom_starts[i];
    addRandomNoise(&xrs, &host_psw.xcrd[aoff], &host_psw.xcrd_ovrf[aoff], &host_psw.ycrd[aoff],
                   &host_psw.ycrd_ovrf[aoff], &host_psw.zcrd[aoff], &host_psw.zcrd_ovrf[aoff],
                   host_psw.atom_counts[i], 0.01, host_psw.gpos_scale_f);
  }
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(std::vector<int>(1, 0));
  
  LocalExclusionMask lem(poly_ag);
  
  // Upload basic objects to the GPU
  poly_ps.upload();
  poly_ag.upload();
  tile_assignments.upload();
  lem.upload();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps.data(devc);
  SyNonbondedKit<float, float2> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc);
  const int4 *t_assgn = tile_assignments.data(devc);
  timer.assignTime(time_build);

}
