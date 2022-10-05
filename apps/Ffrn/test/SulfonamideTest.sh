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
  ncyc 50,  maxcyc 500,
&end

&restraint
  system sulfonamide,
  ensemble heavy_dihedrals,
  mask '@O1,S1,C2,C3',
  penalty 50.0, fbhw 0.0,
&end

%&restraint
%  system sulfonamide
%  ensemble heavy_dihedrals
%  penalty 20.0, fbhw 0.1,
%&end
EOF

${STORMM_HOME}/apps/bin/ffrefine.stormm -O -i ffld.in -warn
