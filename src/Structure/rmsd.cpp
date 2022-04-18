#include "rmsd.h"
#include "structure_utils.h"

namespace omni {
namespace structure {
  
//-------------------------------------------------------------------------------------------------
double rmsd(const PhaseSpaceReader &psr_a, const PhaseSpaceReader &psr_b,
            const ChemicalDetailsKit &cdk, const RmsdMethod method, const int lower_limit,
            const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, cdk.natom, "rmsd");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cdk.natom : upper_limit;
  return rmsd<double, double>(psr_a.xcrd, psr_a.ycrd, psr_a.zcrd, psr_b.xcrd, psr_b.ycrd,
                              psr_b.zcrd, cdk.masses, method, lower_limit, actual_upper_limit);
}

//-------------------------------------------------------------------------------------------------
double rmsd(const PhaseSpaceWriter &psw_a, const PhaseSpaceWriter &psw_b,
            const ChemicalDetailsKit &cdk, const RmsdMethod method, const int lower_limit,
            const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, cdk.natom, "rmsd");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cdk.natom : upper_limit;
  return rmsd<double, double>(psw_a.xcrd, psw_a.ycrd, psw_a.zcrd, psw_b.xcrd, psw_b.ycrd,
                              psw_b.zcrd, cdk.masses, method, lower_limit, actual_upper_limit);
}

//-------------------------------------------------------------------------------------------------
double rmsd(const PhaseSpace &ps_a, const PhaseSpace &ps_b, const AtomGraph &ag,
            const RmsdMethod method, const int lower_limit, const int upper_limit) {
  const PhaseSpaceReader psr_a = ps_a.data();
  const PhaseSpaceReader psr_b = ps_b.data();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  coordinateBoundsCheck(lower_limit, upper_limit, cdk.natom, "rmsd");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cdk.natom : upper_limit;
  return rmsd<double, double>(psr_a.xcrd, psr_a.ycrd, psr_a.zcrd, psr_b.xcrd, psr_b.ycrd,
                              psr_b.zcrd, cdk.masses, method, lower_limit, actual_upper_limit);
}

//-------------------------------------------------------------------------------------------------
double rmsd(const CoordinateFrameReader &cfr_a, const CoordinateFrameReader &cfr_b,
            const ChemicalDetailsKit &cdk, const RmsdMethod method, const int lower_limit,
            const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, cdk.natom, "rmsd");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cdk.natom : upper_limit;
  return rmsd<double, double>(cfr_a.xcrd, cfr_a.ycrd, cfr_a.zcrd, cfr_b.xcrd, cfr_b.ycrd,
                              cfr_b.zcrd, cdk.masses, method, lower_limit, actual_upper_limit);
}

//-------------------------------------------------------------------------------------------------
double rmsd(const CoordinateFrameWriter &cfw_a, const CoordinateFrameWriter &cfw_b,
            const ChemicalDetailsKit &cdk, const RmsdMethod method, const int lower_limit,
            const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, cdk.natom, "rmsd");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cdk.natom : upper_limit;
  return rmsd<double, double>(cfw_a.xcrd, cfw_a.ycrd, cfw_a.zcrd, cfw_b.xcrd, cfw_b.ycrd,
                              cfw_b.zcrd, cdk.masses, method, lower_limit, actual_upper_limit);
}

//-------------------------------------------------------------------------------------------------
double rmsd(const CoordinateFrame &cf_a, const CoordinateFrame &cf_b, const AtomGraph &ag,
            const RmsdMethod method, const int lower_limit, const int upper_limit) {
  const CoordinateFrameReader cfr_a = cf_a.data();
  const CoordinateFrameReader cfr_b = cf_b.data();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  coordinateBoundsCheck(lower_limit, upper_limit, cdk.natom, "rmsd");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cdk.natom : upper_limit;
  return rmsd<double, double>(cfr_a.xcrd, cfr_a.ycrd, cfr_a.zcrd, cfr_b.xcrd, cfr_b.ycrd,
                              cfr_b.zcrd, cdk.masses, method, lower_limit, actual_upper_limit);
}

} // namespace structure
} // namespace omni
