(time angiotk_meshing_volumefromstlandcenterlines --config-file /home/feelpp/tp/cfg/volumefromstlandcenterlines.cfg \
      --input.surface.filename /home/feelpp/tp/7_surfaceremeshing/model_centerlines_up_spacing0.25_remeshGMSHpt15.stl \
      --input.centerlines.filename /home/feelpp/tp/4_centerlinesmanager/model_centerlines_up.vtk \
      --output.directory /home/feelpp/tp/8_volumemeshing \
) 2>&1 | tee /home/feelpp/tp/log/8_vp.log
