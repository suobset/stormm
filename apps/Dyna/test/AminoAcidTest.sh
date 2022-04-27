#!/bin/bash

cat > cgen.in << EOF
&files
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/phe.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/phe.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/phe.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/phe.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/phe.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/phe.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/trp.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/trp.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/trp.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/trp.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/trp.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/trp.inpcrd }
  -sys { -p ${OMNI_SOURCE}/test/Namelists/topol/gly_arg.top
         -c ${OMNI_SOURCE}/test/Namelists/coord/gly_arg.inpcrd }
&end

&minimize
  ncyc 100,
  maxcyc 500,
&end
EOF

cat >> cgen.in << EOF

EOF

${OMNI_HOME}/apps/bin/dynamics.omni -i cgen.in -warn
