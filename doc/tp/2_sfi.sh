(time angiotk_meshing_surfacefromimage --input.image.filename /home/feelpp/tp/1_rorpo/pcp_25_1.34_7_4.mha \
      --output.path /home/feelpp/tp/2_surfacefromimage/model.stl \
      --config-file /home/feelpp/tp/cfg/surfacefromimage.cfg
) 2>&1 | tee /home/feelpp/tp/log/2_sfi.log
