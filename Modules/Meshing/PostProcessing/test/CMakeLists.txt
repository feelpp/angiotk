
add_test(NAME test_meshing_postprocessing COMMAND meshing_postprocessing 
  --extract-submesh.input.filename=@ExternalData_BINARY_ROOT@/Modules/Meshing/PostProcessing/simpleTbifurcationFSI.msh
  --extract-submesh.output.directory=angiotk/meshing/postprocessing/test/simpleTbifurcation/postprocessing
  --extract-submesh.force-rebuild=true
  --mesh-partitioner.npartitions=2
  --mesh-partitioner.npartitions=4
  --mesh-partitioner.npartitions=8
  --mesh-partitioner.force-rebuild=true )

