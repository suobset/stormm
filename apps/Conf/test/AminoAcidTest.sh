#!/bin/bash

echo "&files" > cgen.in
for SYS in gly_lys gly_gly ala arg gly lys phe pro trp tyr gly_ala gly_arg gly_gly gly_lys gly_phe \
           gly_pro gly_trp gly_tyr ; do
  echo "  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/${SYS}.top" >> cgen.in
  echo "         -c ${STORMM_SOURCE}/test/Namelists/coord/${SYS}.inpcrd }" >> cgen.in
done
cat >> cgen.in << EOF
&end

&conformer
  rotation_samples 6,
  max_system_trials 1024,
&end

&solvent
  igb 5
&end

&minimize
  ncyc 50, maxcyc 500
&end
EOF

cat >> cgen.in << EOF

EOF

${STORMM_HOME}/apps/bin/conformer.stormm.cuda -i cgen.in -warn
