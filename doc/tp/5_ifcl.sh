(time angiotk_meshing_imagefromcenterlines --config-file /home/feelpp/tp/cfg/imagefromcenterlines.cfg \
      --input.centerlines.filename /home/feelpp/tp/4_centerlinesmanager/model_centerlines_up.vtk \
      --output.directory /home/feelpp/tp/5_imagefromcenterlines \
) 2>&1 | tee /home/feelpp/tp/log/5_ifcl.log
