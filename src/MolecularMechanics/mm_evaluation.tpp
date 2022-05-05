// -*-c++-*-
namespace omni {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void evalValeMM(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, double* umat, double* invu,
                const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                ScoreCard *sc, const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                const EvaluateForce eval_force, const int system_index,
                const Tcalc inv_gpos_factor, const Tcalc force_factor) {
  evaluateBondTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc,
                                           zfrc, sc, eval_force, system_index, inv_gpos_factor,
                                           force_factor);
  evaluateAngleTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                            yfrc, zfrc, sc, eval_force, system_index,
                                            inv_gpos_factor, force_factor);
  evaluateDihedralTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                               yfrc, zfrc, sc, eval_force, system_index,
                                               inv_gpos_factor, force_factor);
  evaluateUreyBradleyTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell,
                                                  xfrc, yfrc, zfrc, sc, eval_force, system_index,
                                                  inv_gpos_factor, force_factor);
  evaluateCharmmImproperTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell,
                                                     xfrc, yfrc, zfrc, sc, eval_force,
                                                     system_index, inv_gpos_factor, force_factor);
  evaluateCmapTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc,
                                           zfrc, sc, eval_force, system_index, inv_gpos_factor,
                                           force_factor);
  evaluateAttenuated14Terms<Tcoord, Tforce, Tcalc>(vk, nbk, xcrd, ycrd, zcrd, umat, invu,
                                                   unit_cell, xfrc, yfrc, zfrc, sc, eval_force,
                                                   eval_force, system_index, inv_gpos_factor,
                                                   force_factor);
}

} // namespace mm
} // namespace omni
