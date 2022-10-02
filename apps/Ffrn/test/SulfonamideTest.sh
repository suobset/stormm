#!/bin/bash

cat > ffld.in << EOF
&files
  -sys { -p ${STORMM_SOURCE}/test/Topology/sulfonamide.top
         -c ${STORMM_SOURCE}/test/MoleculeFormat/sulfonamide_rots.sdf
         -label sulfonamide frame_end -1 }
&end

&ffrefine
  
&end

&minimize
  ncyc 50,  maxcyc 1000,
&end

&restraint
  label sulfonamide
  mask1 @O1, mask2 @S1, mask3 @C2, mask4 @C3
  
&end
EOF

valgrind ${STORMM_HOME}/apps/bin/ffrefine.stormm -O -i ffld.in -warn
