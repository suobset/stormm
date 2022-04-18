#!/bin/bash

cat > cgen.in << EOF
&files
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
&end

&conformer
  rotation_samples 600,
&end
EOF

cat >> cgen.in << EOF

EOF

valgrind ${OMNI_HOME}/apps/bin/conformer.omni -i cgen.in -warn
