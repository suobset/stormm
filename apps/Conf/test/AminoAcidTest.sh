#!/bin/bash

echo "&files" > cgen.in
for SYS in  gly_lys gly_gly ala arg gly lys phe pro trp tyr gly_ala gly_arg gly_gly gly_lys \
            gly_phe gly_pro gly_trp gly_tyr ; do
  echo "  -sys { -p ${STORMM_SOURCE}/test/Namelists/topol/${SYS}.top" >> cgen.in
  echo "         -c ${STORMM_SOURCE}/test/Namelists/coord/${SYS}.inpcrd " >> cgen.in
  echo "         -x conf_${SYS}.crd x_kind AMBER_CRD }" >> cgen.in
done
cat >> cgen.in << EOF
  -x conf.crd
  x_kind amber_crd
&end

&conformer
  rotation_samples 6,
  max_system_trials 8, final_states 2,
  core_mask { atoms "@N,CA,C & !(:ACE,NME)" }
&end

&solvent
  igb 8
&end

&minimize
  ncyc 50, cdcyc 50, maxcyc 200, ntpr = 1,
  clash_vdw_ratio 0.65,
&end

&random
  igseed 9183025
&end
EOF

cat >> cgen.in << EOF

EOF

if [ -e ${STORMM_BUILD}/apps/Conf/conformer.stormm.cuda ] ; then
  ${STORMM_BUILD}/apps/Conf/conformer.stormm.cuda -i cgen.in -warn
elif [ -e ${STORMM_BUILD}/apps/Conf/conformer.stormm ] ; then
  ${STORMM_BUILD}/apps/Conf/conformer.stormm -i cgen.in -warn
fi
