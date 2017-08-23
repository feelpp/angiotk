(time angiotk_meshing_remeshstl --config-file /home/feelpp/tp/cfg/remeshstlgmsh.cfg \
      --input.surface.filename /home/feelpp/tp/6_surfacefromimage2/model_centerlines_up_spacing0.25.stl \
      --gmsh.centerlines.filename /home/feelpp/tp/4_centerlinesmanager/model_centerlines_up.vtk \
      --output.directory /home/feelpp/tp/7_surfaceremeshing \
) 2>&1 | tee /home/feelpp/tp/log/7_sp.log
