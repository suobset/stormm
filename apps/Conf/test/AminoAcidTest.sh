#!/bin/bash

echo "&files" > cgen.in
for SYS in gly_lys ; do
  echo "  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/${SYS}.top" >> cgen.in
  echo "         -c ${STORMM_SOURCE}/test/Namelists/coord/${SYS}.inpcrd }" >> cgen.in
done
cat >> cgen.in << EOF
&end

&conformer
  rotation_samples 3,
  max_system_trials 1024,
&end
EOF

cat >> cgen.in << EOF

EOF

valgrind ${STORMM_HOME}/apps/bin/conformer.stormm -i cgen.in -warn
