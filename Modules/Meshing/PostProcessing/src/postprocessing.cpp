
#include <postprocessing.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>

namespace Feel
{

MeshPartitioner::MeshPartitioner( std::string prefix )
  :
  M_prefix( prefix ),
  M_inputPath( Environment::expand( soption(_name="input.filename",_prefix=this->prefix()) ) ),
  M_outputPath( 1,Environment::expand( soption(_name="output.filename",_prefix=this->prefix()) ) ),
  M_outputDirectory( Environment::expand( soption(_name="output.directory",_prefix=this->prefix()) ) ),
  M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
  //M_nPartitions( ioption(_name="npartitions",_prefix=this->prefix()) ),
  M_nPartitions( Environment::vm()[prefixvm(prefix,"npartitions").c_str()].as<std::vector<int> >() ),
  M_partitioner( ioption(_name="partitioner",_prefix=this->prefix()) )
{
  if ( M_nPartitions.size()==0 )
    M_nPartitions.resize(1,this->worldComm().size());

  // todo M_inputPath is relative

  if ( !M_inputPath.empty() && M_outputPath[0].empty() )
  {
    if ( M_outputDirectory.empty() )
      this->updateOutputPathFromInputPath();
    else
      this->updateOutputPathFromInputFileName();
  }
  else if ( !M_outputPath[0].empty() )
  {
      this->updateOutputPathFromPath( M_outputPath[0] );
  }

}
MeshPartitioner::MeshPartitioner( MeshPartitioner const& m )
  :
  M_prefix( m.M_prefix ),
  M_inputPath( m.M_inputPath ),
  M_outputPath( m.M_outputPath ),
  M_outputDirectory( m.M_outputDirectory ),
  M_forceRebuild( m.M_forceRebuild ),
  M_nPartitions( m.M_nPartitions ),
  M_partitioner( m.M_partitioner )
{}

void
MeshPartitioner::updateOutputPathFromInputPath()
{
  CHECK( !M_inputPath.empty() ) << "input path is empty";

  fs::path gp = M_inputPath;
  std::string nameMeshFile = gp.stem().string();
  fs::path directory = gp.parent_path();
  M_outputPath.resize( M_nPartitions.size() );
  for ( int k=0;k<M_nPartitions.size();++k)
  {
    std::string newFileName = (boost::format("%1%_p%2%.msh")%nameMeshFile%this->nPartitions(k)).str();
    fs::path outputPath = directory / fs::path(newFileName);
    M_outputPath[k] = outputPath.string();
  }
}

void
MeshPartitioner::updateOutputPathFromInputFileName()
{
  CHECK( !M_inputPath.empty() ) << "input path is empty";

  // define output directory
  fs::path meshesdirectories = Environment::rootRepository();
  if ( M_outputDirectory.empty() )
  {
    meshesdirectories /= fs::path("/angiotk/meshing/meshfile-garbage/" + this->prefix() );
  }
  else meshesdirectories /= fs::path(M_outputDirectory);

  // get filename without extension
  fs::path gp = M_inputPath;
  std::string nameMeshFile = gp.stem().string();

  M_outputPath.resize( M_nPartitions.size() );
  for ( int k=0;k<M_nPartitions.size();++k)
  {
    std::string newFileName = (boost::format("%1%_p%2%.msh")%nameMeshFile%this->nPartitions(k)).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath[k] = outputPath.string();
  }
}

void
MeshPartitioner::updateOutputPathFromPath( std::string outputpathbase )
{
  //CHECK(false) << "TODO";
  fs::path meshesdirectories;
  if ( fs::path(outputpathbase).is_relative() )
  {
    // define output directory
    meshesdirectories = Environment::rootRepository();
    if ( M_outputDirectory.empty() )
    {
      meshesdirectories /= fs::path("/angiotk/meshing/meshfile-garbage/" + this->prefix() );
    }
    else meshesdirectories /= fs::path(M_outputDirectory);
  }
  else
  {
    meshesdirectories = fs::path(outputpathbase).parent_path();
  }

  // get filename without extension
  fs::path gp = outputpathbase;
  std::string nameOutputBase = gp.stem().string();

  M_outputPath.resize( M_nPartitions.size() );
  for ( int k=0;k<M_nPartitions.size();++k)
  {
    std::string newFileName = (boost::format("%1%_p%2%.msh")%nameOutputBase%this->nPartitions(k)).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath[k] = outputPath.string();
  }

}

void
MeshPartitioner::run()
{
  if ( !fs::exists( this->inputPath() ) )
  {
    if ( this->worldComm().isMasterRank() )
      std::cout << "WARNING : partitioning not done because this input path not exist :" << this->inputPath() << "\n";
    return;
  }

  for ( int k=0;k<M_nPartitions.size();++k)
  {
    if ( M_outputPath[k].empty() )
    {
      if ( this->worldComm().isMasterRank() )
	std::cout << "WARNING : partitioning not done because the output path is empty\n";
      return;  
    }


    if ( this->worldComm().isMasterRank() )
    {
      std::ostringstream coutStr;

      coutStr << "\n"
	      << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
	      << "---------------------------------------\n"
	      << "[MeshPartitioner] : run (start) \n"
	      << "---------------------------------------\n";
      coutStr << "input path                          : " << this->inputPath() << "\n"
	      << "number of partition                 : " << this->nPartitions(k) << "\n"
	      << "Gmsh partitioner (1=CHACO, 2=METIS) : " << this->partitioner() << "\n"
	      << "output path                         : " << this->outputPath(k) << "\n"
	      << "---------------------------------------\n"
	      << "---------------------------------------\n";
      std::cout << coutStr.str();
    }


    if ( !fs::exists( this->outputPath(k) ) || this->forceRebuild() )
    {
      // wait all process to be sure that previous test on file existance is idem for all proc
      this->worldComm().globalComm().barrier();
      CHECK( this->inputPath() != this->outputPath(k) ) << "not allow to use same name";

      // build directories if necessary
      if ( !this->outputPath(k).empty() && this->worldComm().isMasterRank() )
      {
	fs::path directory = fs::path(this->outputPath(k)).parent_path();
	if ( !fs::exists( directory ) )
	  fs::create_directories( directory );
      }
      // // wait all process
      this->worldComm().globalComm().barrier();

      // partioning mesh
      Gmsh gmsh( 3,//mesh_fluid_type::nDim,
		 1,//mesh_fluid_type::nOrder,
		 this->worldCommPtr() );
      gmsh.setNumberOfPartitions( this->nPartitions(k) );
      gmsh.setPartitioner( (GMSH_PARTITIONER)this->partitioner() );

      gmsh.rebuildPartitionMsh( this->inputPath(), this->outputPath(k) );
    }
    else
    {
      if ( this->worldComm().isMasterRank() )
	std::cout << "[MeshPartitioner] : partitioning already done : skip this step\n";
    }

    if ( this->worldComm().isMasterRank() )
    {
      std::ostringstream coutStr;
      coutStr << "---------------------------------------\n"
	      << "---------------------------------------\n"
	      << "[MeshPartitioner] : run (finish) \n"
	      << "---------------------------------------\n"
	      << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n";
      std::cout << coutStr.str();
    }

  }
}


po::options_description
MeshPartitioner::options( std::string const& prefix )
{
  po::options_description desc_options("mesh-partitioner options");
  desc_options.add_options()
    (prefixvm(prefix,"input.filename").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) path of input file")
    (prefixvm(prefix,"output.filename").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) path of output file")
    (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
    (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
    //(prefixvm(prefix,"npartitions").c_str(), Feel::po::value<int>()->default_value(-1), "(int) number of partition")
    (prefixvm(prefix,"npartitions").c_str(), po::value<std::vector<int> >()->multitoken(), "number of partition" )
    (prefixvm(prefix,"partitioner").c_str(), Feel::po::value<int>()->default_value(GMSH_PARTITIONER_DEFAULT), "(int) Gmsh partitioner (1=CHACO, 2=METIS)")
    ;

  return desc_options;
}




template< class MeshType >
ExtractSubMeshFromFSIMesh<MeshType>::ExtractSubMeshFromFSIMesh( std::string prefix )
  :
  M_prefix( prefix ),
  M_inputBloodFlowMeshFilename( Environment::expand( soption(_name="input.filename",_prefix=this->prefix()) ) ),
  M_outputDirectory( Environment::expand( soption(_name="output.directory",_prefix=this->prefix()) ) ),
  M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
  M_markerNameLumenVolume( soption(_name="marker.lumen",_prefix=this->prefix()) ),
  M_markerNameArterialWallVolume( soption(_name="marker.arterial-wall",_prefix=this->prefix()) )
{
  this->updateOutputPathFromInputFileName();
}

template< class MeshType >
ExtractSubMeshFromFSIMesh<MeshType>::ExtractSubMeshFromFSIMesh( ExtractSubMeshFromFSIMesh const& m )
  :
  M_prefix( m.M_prefix ),
  M_inputBloodFlowMeshFilename( m.M_inputBloodFlowMeshFilename ),
  M_outputDirectory( m.M_outputDirectory ),
  M_outputPathLumen( m.M_outputPathLumen ),
  M_outputPathArterialWall( m.M_outputPathArterialWall ),
  M_forceRebuild( m.M_forceRebuild ),
  M_markerNameLumenVolume( m.M_markerNameLumenVolume ),
  M_markerNameArterialWallVolume( m.M_markerNameArterialWallVolume ),
  M_meshLumen( m.M_meshLumen ),
  M_meshArterialWall( m.M_meshArterialWall )
{}



template< class MeshType >
void
ExtractSubMeshFromFSIMesh<MeshType>::updateOutputPathFromInputFileName()
{
  // define output directory
  fs::path meshesdirectories = Environment::rootRepository();
  if ( M_outputDirectory.empty() )
  {
    meshesdirectories /= fs::path("/angiotk/meshing/meshfile/" + this->prefix() );
  }
  else meshesdirectories /= fs::path(M_outputDirectory);

  // get filename without extension
  std::string mshfileFSIpathBase = M_inputBloodFlowMeshFilename;
  fs::path gp = mshfileFSIpathBase;
  std::string nameMeshFile = gp.stem().string();

  // define lumen and arterial wall filename
  //std::string lumenFileName = nameMeshFile + "_fluid_" + shape_type::name() + ".msh";
  //std::string arterialWallFileName = nameMeshFile + "_wall_" + shape_type::name() + ".msh";
  std::string lumenFileName = nameMeshFile + "_lumen.msh";
  std::string arterialWallFileName = nameMeshFile + "_wall.msh";
  M_outputPathLumen = ( meshesdirectories / fs::path(lumenFileName) ).string();
  M_outputPathArterialWall = ( meshesdirectories / fs::path(arterialWallFileName) ).string();
}


template< class MeshType >
void
ExtractSubMeshFromFSIMesh<MeshType>::run()
{

  if ( !fs::exists( this->inputPath() ) )
  {
    if ( this->worldComm().isMasterRank() )
      std::cout << "WARNING : sub mesh extraction not done because this input path not exist :" << this->inputPath() << "\n";
    return;
  }

  if ( this->outputPathLumen().empty() && this->outputPathArterialWall().empty() )
  {
    if ( this->worldComm().isMasterRank() )
      std::cout << "WARNING : sub mesh extraction not done because the output path of lumen and arterial-wall are empty\n";
    return;  
  }

  if ( this->worldComm().isMasterRank() )
  {
    std::ostringstream coutStr;
    coutStr << "\n"
	    << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
	    << "---------------------------------------\n"
	    << "[ExtractSubMeshFromFSIMesh] : run (start) \n"
	    << "---------------------------------------\n";
    coutStr << "input path                  : " << this->inputPath() << "\n"
	    << "marker name (Lumen)         : " << this->markerNameLumenVolume() << "\n"
	    << "marker name (Arterial Wall) : " << this->markerNameArterialWallVolume() << "\n"
	    << "output path (Lumen)         : " << this->outputPathLumen() << "\n"
	    << "output path (Arterial Wall) : " << this->outputPathArterialWall() << "\n"
	    << "---------------------------------------\n"
	    << "---------------------------------------\n";
    std::cout << coutStr.str();
  }


  if ( this->worldComm().isMasterRank() )
  {
    if ( !fs::exists( M_outputPathLumen ) || !fs::exists( M_outputPathArterialWall ) || this->forceRebuild() )
    {
      // build directories if necessary
      if ( !this->outputPathLumen().empty() )
      {
	fs::path directoryLumen = fs::path(this->outputPathLumen()).parent_path();
	if ( !fs::exists( directoryLumen ) )
	  fs::create_directories( directoryLumen );
      }
      if ( !this->outputPathArterialWall().empty() )
      {
	fs::path directoryArterialWall = fs::path(this->outputPathArterialWall()).parent_path();
	if ( !fs::exists( directoryArterialWall ) )
	  fs::create_directories( directoryArterialWall );
      }

      //-----------------------------------------------------------------------//
      std::cout << "[ExtractSubMeshFromFSIMesh] : load fsi mesh (start)\n";
      auto meshFSI = loadMesh(_mesh=new mesh_type(this->worldComm().subWorldCommSeqPtr()),
			      _filename=this->inputPath(),
			      _straighten=false,
			      _worldcomm=this->worldComm().subWorldCommSeqPtr(),
			      _rebuild_partitions=false );
      std::cout << "[ExtractSubMeshFromFSIMesh] : load fsi mesh (finish)\n";
      //-----------------------------------------------------------------------//
      //-----------------------------------------------------------------------//
      std::cout << "[ExtractSubMeshFromFSIMesh] : build lumen and arterial wall submesh (start)\n";
      // create lumen submesh
      if ( !this->outputPathLumen().empty() )
      {
	M_meshLumen = createSubmesh(meshFSI,markedelements(meshFSI,this->markerNameLumenVolume()));
	std::cout << "[ExtractSubMeshFromFSIMesh] : number of mesh element in lumen :" << M_meshLumen/*fluidMesh_temp*/->numGlobalElements()<<"\n";
	saveGMSHMesh(_mesh=M_meshLumen,_filename=this->outputPathLumen());
      }
      if ( !this->outputPathArterialWall().empty() )
      {
	// create arterial wall submesh
	M_meshArterialWall = createSubmesh(meshFSI,markedelements(meshFSI,this->markerNameArterialWallVolume()));
	std::cout << "[ExtractSubMeshFromFSIMesh] : number of mesh element in arterial wall :" << M_meshArterialWall/*wallMesh_temp*/->numGlobalElements()<<"\n";
	saveGMSHMesh(_mesh=M_meshArterialWall,_filename=this->outputPathArterialWall());
      }
      //-----------------------------------------------------------------------//
      std::cout << "[ExtractSubMeshFromFSIMesh] : build lumen and arterial wall submesh (finish)\n";
    }
    else
    {
      std::cout << "[ExtractSubMeshFromFSIMesh] : extraction already done : skip this step\n";
    }
  } // isMasterRank()

  // wait writing meshes done
  this->worldComm().globalComm().barrier();

  if ( this->worldComm().isMasterRank() )
  {
    std::ostringstream coutStr;
    coutStr << "---------------------------------------\n"
	    << "---------------------------------------\n"
	    << "[ExtractSubMeshFromFSIMesh] : run (finish) \n"
	    << "---------------------------------------\n"
	    << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n";
    std::cout << coutStr.str();
  }
}

template< class MeshType >
po::options_description
ExtractSubMeshFromFSIMesh<MeshType>::options( std::string const& prefix )
{
  po::options_description desc_options("extract-submesh options");
  desc_options.add_options()
    (prefixvm(prefix,"input.filename").c_str(), Feel::po::value<std::string>(), "(string) path of fsi msh file")
    (prefixvm(prefix,"output.directory").c_str(),Feel::po::value< std::string >()->default_value( "" ), "(string) output directory")
    //(prefixvm(prefix,"output.filename").c_str(),Feel::po::value< std::string >()->default_value( "" ), "(string) subdir export")
    (prefixvm(prefix,"marker.lumen").c_str(), Feel::po::value<std::string>()->default_value("lumenVolume"), "(string) marker fluid")
    (prefixvm(prefix,"marker.arterial-wall").c_str(), Feel::po::value<std::string>()->default_value("wallVolume"), "(string) marker wall")
    (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
    ;
  return desc_options;
}


template class ExtractSubMeshFromFSIMesh< Mesh<Simplex<3> > >;


} //Feel
