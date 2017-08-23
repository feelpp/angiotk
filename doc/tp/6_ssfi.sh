(time angiotk_meshing_surfacefromimage --config-file /home/feelpp/tp/cfg/surfacefromimage2.cfg \
      --input.image.filename /home/feelpp/tp/5_imagefromcenterlines/model_centerlines_up_spacing0.25.mha \
      --output.directory /home/feelpp/tp/6_surfacefromimage2 \
) 2>&1 | tee /home/feelpp/tp/log/6_ssfi.log
