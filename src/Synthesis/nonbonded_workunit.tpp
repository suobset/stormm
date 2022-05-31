//-*-c++-*-
namespace omni {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename Tmask> std::vector<NonbondedWorkUnit>
enumerateNonbondedWorkUnits(const Tmask &se, const int target_tile_count,
                            const int nbwu_estimate, const std::vector<int> &atom_counts) {
  std::vector<NonbondedWorkUnit> result;
  result.reserve(nbwu_estimate);
  const int system_count = atom_counts.size();
  std::vector<int3> tile_list(target_tile_count);
  const int max_imports = small_block_max_atoms / tile_length;
  std::vector<int> import_list(max_imports);
  int current_tile_count = 0;
  for (int sysid = 0; sysid < system_count; sysid++) {
    const int total_tile_lengths = (atom_counts[sysid] + tile_length - 1) / tile_length;
    std::vector<int> import_coverage(total_tile_lengths, 0);

    // The strategy is to list tiles in the typical order, keeping track of import requirements.
    // Get as close as possible to the target work unit size, keeping to multiples of eight if at
    // all possible, then commit work units to the list as they reach the target size or it
    // becomes impossible to grow the work unit without exceeding the import limits.
    const int ordi_stride = (target_tile_count ==  8) ? 2 :
                            (target_tile_count == 16) ? 4 :
                            (target_tile_count == 32) ? 6 :
                            (target_tile_count == 64) ? 8 : tile_lengths_per_supertile;
    int ordi_llim = 0;
    int ordi_hlim = std::min(ordi_llim + ordi_stride, total_tile_lengths);
    int absc_increment = 1;
    int import_count = 0;
    while (ordi_llim < total_tile_lengths) {
      int istart, ilimit;
      if (absc_increment == 1) {
        istart = 0;
        ilimit = total_tile_lengths;
      }
      else {
        istart = total_tile_lengths - 1;
        ilimit = -1;
      }
      for (int i = istart; i != ilimit; i += absc_increment) {
        for (int j = ordi_llim; j < ordi_hlim; j++) {
          if (j > i) {
            continue;
          }

          // Test whether this new tile can be added within the limits of the work unit.
          if (addTileToWorkUnitList(tile_list.data(), import_coverage.data(), &import_count,
                                    &current_tile_count, i, j, sysid)) {
            if (current_tile_count == target_tile_count) {

              // This list of tiles has reached its target size and should be converted into a
              // work unit.
              result.emplace_back(se, tile_list);
              for (int k = 0 ; k < current_tile_count; k++) {
                if (tile_list[k].z == sysid) {
                  import_coverage[tile_list[k].x] = 0;
                  import_coverage[tile_list[k].y] = 0;
                }
              }
              import_count = 0;
              current_tile_count = 0;
            }
          }
          else {

            // Backtrack to the most recent multiple of the target batch size.
            const int nbatch = current_tile_count / small_block_tile_width;
            const int fallback_size = nbatch * small_block_tile_width;
            std::vector<int3> tmp_tile_list(fallback_size);
            for (int k = 0; k < fallback_size; k++) {
              tmp_tile_list[k] = tile_list[k];
            }
            result.emplace_back(se, tmp_tile_list);
            for (int k = 0 ; k < fallback_size; k++) {
              if (tile_list[k].z == sysid) {
                import_coverage[tile_list[k].x] = 0;
                import_coverage[tile_list[k].y] = 0;
              }
            }
            for (int k = fallback_size; k < current_tile_count; k++) {
              tile_list[k - fallback_size] = tile_list[k];
            }
            current_tile_count -= fallback_size;
            import_count = 0;
            for (int k = 0; k < current_tile_count; k++) {
              const int kx_import = tile_list[k].x;
              const int ky_import = tile_list[k].y;
              const int kz_system = tile_list[k].z;
              bool x_unique = true;
              bool y_unique = true;
              for (int m = 0; m < k; m++) {
                if (tile_list[m].z == kz_system) {
                  x_unique = (x_unique && (kx_import != tile_list[m].x) &&
                              (kx_import != tile_list[m].y));
                  y_unique = (y_unique && (ky_import != tile_list[m].x) &&
                              (ky_import != tile_list[m].y));
                }
              }
              if (kx_import == ky_import) {
                import_count += x_unique;
              }
              else {
                import_count += x_unique + y_unique;
              }
            }

            // Add the tile now.  There is guaranteed to be sufficient room in the work unit
            // based on the comparison of small_block_max_imports and small_block_tile_width.
            addTileToWorkUnitList(tile_list.data(), import_coverage.data(), &import_count,
                                  &current_tile_count, i, j, sysid);
          }
        }
      }
      absc_increment *= -1;
      ordi_llim = ordi_hlim;
      ordi_hlim = std::min(ordi_llim + ordi_stride, total_tile_lengths);
    }
  }
  
  // Commit the final tile, then return the result
  tile_list.resize(current_tile_count);
  result.emplace_back(se, tile_list);
  return result;
}

} // namespace synthesis
} // namespace omni
