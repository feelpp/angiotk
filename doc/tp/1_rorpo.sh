mkdir -p /home/feelpp/tp/1_rorpo
(time RORPO_multiscale_usage /home/feelpp/tp/0_iso/pcp_iso.nii /home/feelpp/tp/1_rorpo/pcp_25_1.34_7_4.nii 25 1.34 7 --core 4 --verbose) 2>&1 | tee /home/feelpp/tp/log/1_rorpo.log
