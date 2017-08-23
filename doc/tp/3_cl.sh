(time angiotk_meshing_centerlines --config-file /home/feelpp/tp/cfg/centerlines.cfg \
      --input.surface.filename /home/feelpp/tp/2_surfacefromimage/model.stl \
      --input.pointpair.filename /home/feelpp/tp/cfg/pointpair1.data \
      --output.directory /home/feelpp/tp/3_centerlines/part1 \
      --delaunay-tessellation.output.directory /home/feelpp/tp/3_centerlines/part1 \
      --source-ids 0 --target-ids 1 \
) 2>&1 | tee /home/feelpp/tp/log/3_cl.log
