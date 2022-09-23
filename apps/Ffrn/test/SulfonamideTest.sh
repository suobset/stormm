#!/bin/bash

cat > minimize.in << EOF
&files
  -sys { -p ${STORMM_SOURCE}/test/MoleculeFormat/sulfonamide.top
         -c ${STORMM_SOURCE}/test/MoleculeFormat/sulfonamide_rots.sdf }
&end

&ffrefine
  rotation_samples 600,
&end
EOF

valgrind ${STORMM_HOME}/apps/bin/ffrefine.stormm -i minimize.in -warn
