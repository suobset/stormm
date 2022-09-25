#!/bin/bash

cat > cgen.in << EOF
&files
  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${STORMM_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${STORMM_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${STORMM_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${STORMM_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${STORMM_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${STORMM_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${STORMM_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${STORMM_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
&end

&conformer
  rotation_samples 60,
&end
EOF

cat >> cgen.in << EOF

EOF

valgrind ${STORMM_HOME}/apps/bin/conformer.stormm -i cgen.in -warn
