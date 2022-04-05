#include "Math/vector_ops.h"
#include "Topology/atomgraph_enumerators.h"
#include "eval_valence_workunit.h"
#include "valence_potential.h"

namespace omni {
namespace energy {

using math::crossProduct;
using math::readBitFromMask;
using topology::TorsionKind;
using trajectory::PhaseSpaceReader;
using trajectory::PhaseSpaceWriter;
  
//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                          const NonbondedKit<double> nbk, const double* xcrd, const double* ycrd,
                          const double* zcrd, const double* umat, const double* invu,
                          const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                          ScoreCard *ecard, const int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          const EvaluateForce eval_force, const VwuTask activity) {  

  // Initialize energy accumulators
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();
  llint bond_acc = 0LL;
  llint angl_acc = 0LL;
  llint dihe_acc = 0LL;
  llint impr_acc = 0LL;
  llint ubrd_acc = 0LL;
  llint cimp_acc = 0LL;
  llint cmap_acc = 0LL;
  llint qq14_acc = 0LL;
  llint lj14_acc = 0LL;

  // Loop over each work unit
  for (size_t vidx = 0LLU; vidx < vwu_list.size(); vidx++) {

    // Get all tasks from the work unit
    const std::vector<int> task_counts = vwu_list[vidx].getTaskCounts();

    // Make local arrays of the coordinates and forces.  Import coordinates and extant forces
    // from the global arrays.
    const int max_atoms = vwu_list[vidx].getMaxAtoms();
    std::vector<double> sh_xcrd(max_atoms), sh_ycrd(max_atoms), sh_zcrd(max_atoms);
    std::vector<double> sh_xfrc(max_atoms), sh_yfrc(max_atoms), sh_zfrc(max_atoms);
    std::vector<double> sh_charges(max_atoms);
    std::vector<int> sh_lj_idx(max_atoms);
    const int n_imp_atoms = vwu_list[vidx].getImportedAtomCount();
    for (int i = 0; i < n_imp_atoms; i++) {
      const size_t atom_idx = vwu_list[vidx].getImportedAtomIndex(i);
      sh_xcrd[i] = xcrd[atom_idx];
      sh_ycrd[i] = ycrd[atom_idx];
      sh_zcrd[i] = zcrd[atom_idx];
      sh_xfrc[i] = xfrc[atom_idx];
      sh_yfrc[i] = yfrc[atom_idx];
      sh_zfrc[i] = zfrc[atom_idx];
      sh_charges[i] = nbk.charge[atom_idx];
      sh_lj_idx[i] = nbk.lj_idx[atom_idx];
    }
    
    // Evaluate bonds and Urey-Bradley terms
    if (activity == VwuTask::BOND || activity == VwuTask::UBRD || activity == VwuTask::ALL_TASKS) {
      const int ncbnd = task_counts[static_cast<int>(VwuTask::CBND)];
      const std::vector<uint> cbnd_acc_mask = vwu_list[vidx].getAccumulationFlags(VwuTask::CBND);
      for (int pos = 0; pos < ncbnd; pos++) { 
        const uint2 tinsr = vwu_list[vidx].getCompositeBondInstruction(pos);
        const bool is_urey_bradley = ((tinsr.x >> 20) & 0x1);
        
        // Skip Urey-Bradley interactions if only bonds are desired, and skip bonds if only
        // Urey-Bradley interactions are of interest.
        if ((activity == VwuTask::BOND && is_urey_bradley) ||
            (activity == VwuTask::UBRD && (! is_urey_bradley))) {
          continue;
        }
        const int i_atom = (tinsr.x & 0x3ff);
        const int j_atom = ((tinsr.x >> 10) & 0x3ff);
        const int param_idx = tinsr.y;
        const double keq = (is_urey_bradley) ? vk.ubrd_keq[param_idx] : vk.bond_keq[param_idx];
        const double leq = (is_urey_bradley) ? vk.ubrd_leq[param_idx] :
                                               fabs(vk.bond_leq[param_idx]);
        const double dx = sh_xcrd[j_atom] - sh_xcrd[i_atom];
        const double dy = sh_ycrd[j_atom] - sh_ycrd[i_atom];
        const double dz = sh_zcrd[j_atom] - sh_zcrd[i_atom];
        const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
        const double dl = dr - leq;
        const double du = keq * dl * dl;
        if (readBitFromMask(cbnd_acc_mask, pos) == 1) {
          if (is_urey_bradley) {
            ubrd_acc += static_cast<llint>(round(du * nrg_scale_factor));
          }
          else {
            bond_acc += static_cast<llint>(round(du * nrg_scale_factor));
          }
        }

        // Compute forces
        if (eval_force == EvaluateForce::YES) {
          const double fmag = 2.0 * keq * dl / dr;
          const double fmag_dx = fmag * dx;
          const double fmag_dy = fmag * dy;
          const double fmag_dz = fmag * dz;
          sh_xfrc[i_atom] += fmag_dx;
          sh_yfrc[i_atom] += fmag_dy;
          sh_zfrc[i_atom] += fmag_dz;
          sh_xfrc[j_atom] -= fmag_dx;
          sh_yfrc[j_atom] -= fmag_dy;
          sh_zfrc[j_atom] -= fmag_dz;
        }
      }
    }
    
    // Evaluate harmonic bond angles
    if (activity == VwuTask::ANGL || activity == VwuTask::ALL_TASKS) {
      const int nangl = task_counts[static_cast<int>(VwuTask::ANGL)];
      const std::vector<uint> angl_acc_mask = vwu_list[vidx].getAccumulationFlags(VwuTask::ANGL);
      for (int pos = 0; pos < nangl; pos++) {
        const uint2 tinsr = vwu_list[vidx].getAngleInstruction(pos);
        const int i_atom = (tinsr.x & 0x3ff);
        const int j_atom = ((tinsr.x >> 10) & 0x3ff);
        const int k_atom = ((tinsr.x >> 20) & 0x3ff);
        const int param_idx = tinsr.y;
        const double keq = vk.angl_keq[param_idx];
        const double theta0 = vk.angl_theta[param_idx];

        // Compute displacements
        double ba[3], bc[3];
        ba[0] = sh_xcrd[i_atom] - sh_xcrd[j_atom];
        ba[1] = sh_ycrd[i_atom] - sh_ycrd[j_atom];
        ba[2] = sh_zcrd[i_atom] - sh_zcrd[j_atom];
        bc[0] = sh_xcrd[k_atom] - sh_xcrd[j_atom];
        bc[1] = sh_ycrd[k_atom] - sh_ycrd[j_atom];
        bc[2] = sh_zcrd[k_atom] - sh_zcrd[j_atom];

        // On to the angle force computation
        const double mgba = (ba[0] * ba[0]) + (ba[1] * ba[1]) + (ba[2] * ba[2]);
        const double mgbc = (bc[0] * bc[0]) + (bc[1] * bc[1]) + (bc[2] * bc[2]);
        const double invbabc = 1.0 / sqrt(mgba * mgbc);
        double costheta = ((ba[0] * bc[0]) + (ba[1] * bc[1]) + (ba[2] * bc[2])) * invbabc;
        costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
        const double theta = acos(costheta);
        const double dtheta = theta - theta0;
        const double du = keq * dtheta * dtheta;
        if (readBitFromMask(angl_acc_mask, pos) == 1) {
          angl_acc += static_cast<llint>(round(du * nrg_scale_factor));
        }
        
        // Compute forces
        if (eval_force == EvaluateForce::YES) {
          const double dA = -2.0 * keq * dtheta / sqrt(1.0 - costheta * costheta);
          const double sqba = dA / mgba;
          const double sqbc = dA / mgbc;
          const double mbabc = dA * invbabc;
          double adf[3], cdf[3];
          for (int i = 0; i < 3; i++) {
            adf[i] = (bc[i] * mbabc) - (costheta * ba[i] * sqba);
            cdf[i] = (ba[i] * mbabc) - (costheta * bc[i] * sqbc);
          }
          sh_xfrc[i_atom] -= adf[0];
          sh_yfrc[i_atom] -= adf[1];
          sh_zfrc[i_atom] -= adf[2];
          sh_xfrc[j_atom] += adf[0] + cdf[0];
          sh_yfrc[j_atom] += adf[1] + cdf[1];
          sh_zfrc[j_atom] += adf[2] + cdf[2];
          sh_xfrc[k_atom] -= cdf[0];
          sh_yfrc[k_atom] -= cdf[1];
          sh_zfrc[k_atom] -= cdf[2];
        }
      }
    }

    // Evaluate cosine-based dihedrals and CHARMM improper dihedrals
    if (activity == VwuTask::DIHE || activity == VwuTask::CIMP || activity == VwuTask::INFR14 ||
        activity == VwuTask::ALL_TASKS) {
      const int ncdhe = task_counts[static_cast<int>(VwuTask::CDHE)];
      const std::vector<uint> cdhe_acc_mask = vwu_list[vidx].getAccumulationFlags(VwuTask::CDHE);
      for (int pos = 0; pos < ncdhe; pos++) {
        const uint3 tinsr = vwu_list[vidx].getCompositeDihedralInstruction(pos);

        // Compute 1:4 interactions and go on if that is all that is required.  Here, the
        // 1:4 interactions are being *inferred* from the dihedral that controls them, as this
        // routine can compute all van-der Waals and electrostatic 1:4 interactions as specific
        // energy quantities.
        const int i_atom = (tinsr.x & 0x3ff);
        const int l_atom = (tinsr.y & 0x3ff);
        if (activity == VwuTask::INFR14 || activity == VwuTask::ALL_TASKS) {
          const int attn_idx = ((tinsr.y >> 10) & 0x1f);
          if (attn_idx > 0) {
            const double2 uc = evaluateAttenuated14Pair(i_atom, l_atom, attn_idx,
                                                        nbk.coulomb_constant, sh_charges.data(),
                                                        sh_lj_idx.data(), vk.attn14_elec,
                                                        vk.attn14_vdw, nbk.lja_14_coeff,
                                                        nbk.ljb_14_coeff, nbk.n_lj_types,
                                                        sh_xcrd.data(), sh_ycrd.data(),
                                                        sh_zcrd.data(), umat, invu, unit_cell,
                                                        sh_xfrc.data(), sh_yfrc.data(),
                                                        sh_zfrc.data(), eval_force, eval_force);
            if (readBitFromMask(cdhe_acc_mask, pos) == 1) {
              qq14_acc += static_cast<llint>(round(uc.x * nrg_scale_factor));
              lj14_acc += static_cast<llint>(round(uc.y * nrg_scale_factor));
            }
          }
          if (activity == VwuTask::INFR14) {
            continue;
          }
        }
        
        const bool is_charmm_improper = ((tinsr.x >> 30) & 0x1);

        // Skip CHARMM improper interactions if only standard dihedrals are of interest, and
        // skip standard dihedrals if only CHARMM impropers are of interest.
        if ((activity == VwuTask::DIHE && is_charmm_improper) ||
            (activity == VwuTask::CIMP && (! is_charmm_improper))) {
          continue;
        }
        const int j_atom = ((tinsr.x >> 10) & 0x3ff);
        const int k_atom = ((tinsr.x >> 20) & 0x3ff);
        const int param_idx = ((tinsr.y >> 16) & 0xffff);
        const bool second_term = ((tinsr.y >> 15) & 0x1);
        const TorsionKind kind = ((tinsr.x >> 31) & 0x1) ? TorsionKind::IMPROPER :
                                                           TorsionKind::PROPER;
        
        // Compute displacements
        double ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
        ab[0] = sh_xcrd[j_atom] - sh_xcrd[i_atom];
        ab[1] = sh_ycrd[j_atom] - sh_ycrd[i_atom];
        ab[2] = sh_zcrd[j_atom] - sh_zcrd[i_atom];
        bc[0] = sh_xcrd[k_atom] - sh_xcrd[j_atom];
        bc[1] = sh_ycrd[k_atom] - sh_ycrd[j_atom];
        bc[2] = sh_zcrd[k_atom] - sh_zcrd[j_atom];
        cd[0] = sh_xcrd[l_atom] - sh_xcrd[k_atom];
        cd[1] = sh_ycrd[l_atom] - sh_ycrd[k_atom];
        cd[2] = sh_zcrd[l_atom] - sh_zcrd[k_atom];

        // Compute cross products and then the angle between the planes
        crossProduct(ab, bc, crabbc);
        crossProduct(bc, cd, crbccd);
        double costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
        costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                         (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
        crossProduct(crabbc, crbccd, scr);
        costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
        const double theta = (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) ?  acos(costheta) :
                                                                                  -acos(costheta);
        double sangle, sangle_ii, stiffness, stiffness_ii;
        if (is_charmm_improper) {
          stiffness = vk.cimp_keq[param_idx];
          sangle = theta - vk.cimp_phi[param_idx];
          double contrib = stiffness * sangle * sangle;
          if (readBitFromMask(cdhe_acc_mask, pos) == 1) {
            cimp_acc += static_cast<llint>(round(contrib * nrg_scale_factor));          
          }
        }
        else {
          const double ampl = (param_idx < 65535) ? vk.dihe_amp[param_idx] : 0.0;
          const double freq = vk.dihe_freq[param_idx];
          const double phi  = vk.dihe_phi[param_idx];
          stiffness = ampl * freq;
          sangle = (freq * theta) - phi;
          double contrib = ampl * (1.0 + cos(sangle));
          if (second_term) {
            const int param_ii_idx = (tinsr.z & 0xfffff);
            const double ampl_ii = vk.dihe_amp[param_ii_idx];
            const double freq_ii = vk.dihe_freq[param_ii_idx];
            const double phi_ii  = vk.dihe_phi[param_ii_idx];
            stiffness_ii = ampl_ii * freq_ii;
            sangle_ii = (freq_ii * theta) - phi_ii;
            contrib += ampl_ii * (1.0 + cos(sangle_ii));
          }
          if (readBitFromMask(cdhe_acc_mask, pos) == 1) {
            if (kind == TorsionKind::PROPER) {
              dihe_acc += static_cast<llint>(round(contrib * nrg_scale_factor));          
            }
            else {
              impr_acc += static_cast<llint>(round(contrib * nrg_scale_factor));
            }
          }
        }

        // Compute forces, if requested
        if (eval_force == EvaluateForce::YES) {
          double fr;
          if (is_charmm_improper) {
            fr = -2.0 * stiffness * sangle;
          }
          else {
            fr = stiffness * sin(sangle);
            if (second_term) {
              fr += stiffness_ii * sin(sangle_ii);
            }
          }
          const double mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
          const double invab = 1.0 / mgab;
          const double mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
          const double invbc = 1.0 / mgbc;
          const double mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
          const double invcd = 1.0 / mgcd;
          const double cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) * invab * invbc;
          const double isinb2 = (cosb * cosb < asymptotic_to_one_lf) ?
                                fr / (1.0 - (cosb * cosb)) : fr * inverse_one_minus_asymptote_lf;
          const double cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
          const double isinc2 = (cosc * cosc < asymptotic_to_one_lf) ?
                                fr / (1.0 - (cosc * cosc)) : fr * inverse_one_minus_asymptote_lf;
          const double invabc = invab * invbc;
          const double invbcd = invbc * invcd;
          for (int i = 0; i < 3; i++) {
            crabbc[i] *= invabc;
            crbccd[i] *= invbcd;
          }

          // Transform the rotational derivatives to cartesian coordinates
          const double fa = -invab * isinb2;
          const double fb1 = (mgbc - (mgab * cosb)) * invabc * isinb2;
          const double fb2 = cosc * invbc * isinc2;
          const double fc1 = (mgbc - (mgcd * cosc)) * invbcd * isinc2;
          const double fc2 = cosb * invbc * isinb2;
          const double fd = -invcd * isinc2;
          sh_xfrc[i_atom] += crabbc[0] * fa;
          sh_xfrc[j_atom] += (fb1 * crabbc[0]) - (fb2 * crbccd[0]);
          sh_xfrc[k_atom] += (fc2 * crabbc[0]) - (fc1 * crbccd[0]);
          sh_xfrc[l_atom] -= fd * crbccd[0];
          sh_yfrc[i_atom] += crabbc[1] * fa;
          sh_yfrc[j_atom] += (fb1 * crabbc[1]) - (fb2 * crbccd[1]);
          sh_yfrc[k_atom] += (fc2 * crabbc[1]) - (fc1 * crbccd[1]);
          sh_yfrc[l_atom] -= fd * crbccd[1];
          sh_zfrc[i_atom] += crabbc[2] * fa;
          sh_zfrc[j_atom] += (fb1 * crabbc[2]) - (fb2 * crbccd[2]);
          sh_zfrc[k_atom] += (fc2 * crabbc[2]) - (fc1 * crbccd[2]);
          sh_zfrc[l_atom] -= fd * crbccd[2];
        }
      }
    }

    // Evaluate CMAP-based terms 
    if (activity == VwuTask::CMAP || activity == VwuTask::ALL_TASKS) {
      const int ncmap = task_counts[static_cast<int>(VwuTask::CMAP)];
      const std::vector<uint> cmap_acc_mask = vwu_list[vidx].getAccumulationFlags(VwuTask::CMAP);
      for (int pos = 0; pos < ncmap; pos++) {
        const uint2 tinsr = vwu_list[vidx].getCmapInstruction(pos);
        const int i_atom = (tinsr.x & 0x3ff);
        const int j_atom = ((tinsr.x >> 10) & 0x3ff);
        const int k_atom = ((tinsr.x >> 20) & 0x3ff);
        const int l_atom = (tinsr.y & 0x3ff);
        const int m_atom = ((tinsr.y >> 10) & 0x3ff);
        const int surf_idx = (tinsr.y >> 20);

        // Assume that, by this point, the imported atoms in the ValenceWorkUnit have been
        // arranged into an image that can be trusted for all interactions.  Avoid further
        // re-imaging of displacements inside the CMAP calculation.
        const double contrib = evalCmap(vk.cmap_patches, vk.cmap_patch_bounds, surf_idx,
                                        vk.cmap_dim[surf_idx], i_atom, j_atom, k_atom, l_atom,
                                        m_atom, sh_xcrd.data(), sh_ycrd.data(), sh_zcrd.data(),
                                        nullptr, nullptr, UnitCellType::NONE, sh_xfrc.data(),
                                        sh_yfrc.data(), sh_zfrc.data(), eval_force);
        if (readBitFromMask(cmap_acc_mask, pos) == 1) {
          cmap_acc += static_cast<llint>(round(contrib * nrg_scale_factor));
        }
      }
    }
    
    // Add accumulated forces back to the global arrays (this is not done by all GPU kernels, as
    // in some cases the ValenceWorkUnits also move atoms and then leave the global force arrays
    // zero'ed for atoms they are responsible for moving).
    for (int i = 0; i < n_imp_atoms; i++) {
      const size_t atom_idx = vwu_list[vidx].getImportedAtomIndex(i);
      xfrc[atom_idx] += sh_xfrc[i];
      yfrc[atom_idx] += sh_yfrc[i];
      zfrc[atom_idx] += sh_zfrc[i];
    }

    // Evaluate remaining 1:4 interactions.  Calling this routine by the "INFR14" enumeration will
    // compute all 1:4 interactions, a redefinition of "inferred 1:4" relative to other situations.
    // In calling this function, "inferred 1:4" means *inferring* 1:4 interactions from dihedrals
    // that control them as well as extra 1:4 interactions that cannot be attached to a dihedral
    // interaction.  Elsewhere, including the following conditional statement, "inferred 1:4"
    // refers to only the extra interactions (which had to be *inferred* from virtual site
    // exclusion rules).  This minor bit of confusion was accepted for convenience, to reuse
    // VwuTask and avoid making a totally separate enumerator.  The StateVariable enumerator
    // doesn't serve quite the right definitions, either.
    if (activity == VwuTask::INFR14 || activity == VwuTask::ALL_TASKS) {
      const int ninfr14 = task_counts[static_cast<int>(VwuTask::INFR14)];
      const std::vector<uint> infr_acc_mask = vwu_list[vidx].getAccumulationFlags(VwuTask::INFR14);
      for (int pos = 0; pos < ninfr14; pos++) {
        const uint tinsr = vwu_list[vidx].getInferred14Instruction(pos);
        const int i_atom = (tinsr & 0x3ff);
        const int l_atom = ((tinsr >> 10) & 0x3ff);
        const int attn_idx = (tinsr >> 20);
        if (attn_idx == 0) {
          continue;
        }
        const double2 uc = evaluateAttenuated14Pair(i_atom, l_atom, attn_idx, nbk.coulomb_constant,
                                                    sh_charges.data(), sh_lj_idx.data(),
                                                    vk.attn14_elec, vk.attn14_vdw,
                                                    nbk.lja_14_coeff, nbk.ljb_14_coeff,
                                                    nbk.n_lj_types, sh_xcrd.data(), sh_ycrd.data(),
                                                    sh_zcrd.data(), umat, invu, unit_cell,
                                                    sh_xfrc.data(), sh_yfrc.data(), sh_zfrc.data(),
                                                    eval_force, eval_force);
        if (readBitFromMask(infr_acc_mask, pos) == 1) {
          qq14_acc += static_cast<llint>(round(uc.x * nrg_scale_factor));
          lj14_acc += static_cast<llint>(round(uc.y * nrg_scale_factor));
        }
      }
    }
  }

  // Contribute results as the instantaneous states
  switch (activity) {
  case VwuTask::BOND:
    ecard->contribute(StateVariable::BOND, bond_acc, sysid);
    break;
  case VwuTask::ANGL:
    ecard->contribute(StateVariable::ANGLE, angl_acc, sysid);
    break;
  case VwuTask::DIHE:
    ecard->contribute(StateVariable::PROPER_DIHEDRAL, dihe_acc, sysid);
    ecard->contribute(StateVariable::IMPROPER_DIHEDRAL, impr_acc, sysid);
    break;
  case VwuTask::UBRD:
    ecard->contribute(StateVariable::UREY_BRADLEY, ubrd_acc, sysid);
    break;
  case VwuTask::CIMP:
    ecard->contribute(StateVariable::CHARMM_IMPROPER, cimp_acc, sysid);
    break;
  case VwuTask::CMAP:
    ecard->contribute(StateVariable::CMAP, cmap_acc, sysid);
    break;
  case VwuTask::INFR14:
    ecard->contribute(StateVariable::ELECTROSTATIC_ONE_FOUR, qq14_acc, sysid);
    ecard->contribute(StateVariable::VDW_ONE_FOUR, lj14_acc, sysid);
    break;
  case VwuTask::RPOSN:
  case VwuTask::RBOND:
  case VwuTask::RANGL:
  case VwuTask::RDIHE:
    break;
  case VwuTask::SETTLE:
  case VwuTask::CGROUP:
  case VwuTask::VSITE:
  case VwuTask::CDHE:
  case VwuTask::CBND:
    break;
  case VwuTask::ALL_TASKS:
    ecard->contribute(StateVariable::BOND, bond_acc, sysid);
    ecard->contribute(StateVariable::ANGLE, angl_acc, sysid);
    ecard->contribute(StateVariable::PROPER_DIHEDRAL, dihe_acc, sysid);
    ecard->contribute(StateVariable::IMPROPER_DIHEDRAL, impr_acc, sysid);
    ecard->contribute(StateVariable::UREY_BRADLEY, ubrd_acc, sysid);
    ecard->contribute(StateVariable::CHARMM_IMPROPER, cimp_acc, sysid);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard, const int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          const EvaluateForce eval_force, const VwuTask activity) {
  PhaseSpaceWriter psw = ps->data();
  evalValenceWorkUnits(ag->getDoublePrecisionValenceKit(), ag->getDoublePrecisionVirtualSiteKit(),
                       ag->getDoublePrecisionNonbondedKit(), psw.xcrd, psw.ycrd, psw.zcrd,
                       psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, ecard,
                       sysid, vwu_list, eval_force, activity);
}

//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const AtomGraph &ag, const PhaseSpace &ps, ScoreCard *ecard,
                          const int sysid, const std::vector<ValenceWorkUnit> &vwu_list,
                          const VwuTask activity) {
  PhaseSpaceReader psr = ps.data();
  evalValenceWorkUnits(ag.getDoublePrecisionValenceKit(), ag.getDoublePrecisionVirtualSiteKit(),
                       ag.getDoublePrecisionNonbondedKit(), psr.xcrd, psr.ycrd, psr.zcrd, psr.umat,
                       psr.invu, psr.unit_cell, nullptr, nullptr, nullptr, ecard, sysid, vwu_list,
                       EvaluateForce::NO, activity);
}

} // namespace energy
} // namespace omni
