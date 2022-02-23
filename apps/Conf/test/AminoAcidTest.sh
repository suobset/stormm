#!/bin/bash

cat > cgen.in << EOF
&files
  -p ${OMNI_SOURCE}/test/Namelists/topol/.*.top
  -c ${OMNI_SOURCE}/test/Namelists/coord/.*.inpcrd
&end
EOF

cat >> cgen.in << EOF

EOF

valgrind ${OMNI_HOME}/apps/bin/conformer.omni -i cgen.in -warn
