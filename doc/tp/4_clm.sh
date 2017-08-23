(time angiotk_meshing_centerlinesmanager --config-file /home/feelpp/tp/cfg/centerlinesmanager.cfg \
				   --input.surface.filename /home/feelpp/tp/2_surfacefromimage/model.stl \
				   --output.directory /home/feelpp/tp/4_centerlinesmanager \
				   --input.centerlines.filename /home/feelpp/tp/3_centerlines/part1/model_centerlines.vtk
 ) 2>&1 | tee /home/feelpp/tp/log/4_clm.log
