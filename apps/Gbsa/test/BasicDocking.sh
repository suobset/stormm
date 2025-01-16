#!/bin/bash

cat > nrg_est.in << EOF
&files
#  -sys { -p ${STORMM_SOURCE}/test/Topology/trpcage.top
#         -c ${STORMM_SOURCE}/test/Trajectory/trpcage.inpcrd -label ligand}
  -sys { -p ${STORMM_SOURCE}/test/Topology/symmetry_C1.top
         -c ${STORMM_SOURCE}/test/Trajectory/symmetry_C1.inpcrd -label ligand}
  -sys { -p ${STORMM_SOURCE}/test/Topology/symmetry_C2.top
         -c ${STORMM_SOURCE}/test/Trajectory/symmetry_C2.inpcrd -label ligand}
&end
EOF

valgrind ${STORMM_BUILD}/apps/Gbsa/mmgbsa.stormm.cuda -O -i nrg_est.in \
                                        -rec_p ${STORMM_SOURCE}/test/Topology/symmetry_C3.top \
                                        -rec_c ${STORMM_SOURCE}/test/Trajectory/symmetry_C3.inpcrd

              
