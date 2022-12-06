#!/bin/bash

echo "&files" > cgen.in
for SYS in gly_lys gly_gly ala arg gly lys phe pro trp tyr gly_ala gly_arg gly_gly gly_lys \
           gly_phe gly_pro gly_trp gly_tyr ; do
  echo "  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/${SYS}.top" >> cgen.in
  echo "         -c ${STORMM_SOURCE}/test/Namelists/coord/${SYS}.inpcrd }" >> cgen.in
done
cat >> cgen.in << EOF
&end

&conformer
  rotation_samples 6,
  max_system_trials 24, final_states 24,
  core_mask { atoms "@N,CA,C & !(:ACE,NME)" }
&end

&solvent
  igb 5
&end

&minimize
  ncyc 2, maxcyc 5
&end
EOF

cat >> cgen.in << EOF

EOF

valgrind ${STORMM_HOME}/apps/bin/conformer.stormm -i cgen.in -warn
