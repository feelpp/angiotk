#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions = VolumeMeshing::options("");

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_volumefromstlandcenterlines",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));


    VolumeMeshing myMeshingVolume("");
    if ( myMeshingVolume.inputInletOutletDescPath().empty() )
      {
	std::string inputSurfacePath = myMeshingVolume.inputSTLPath();
	std::string nameWithoutExt = fs::path(inputSurfacePath).stem().string();
	std::string outputDescPath = (fs::path(myMeshingVolume.outputPath()).parent_path()/fs::path(nameWithoutExt+".desc")).string();

	if ( !fs::exists( outputDescPath ) || myMeshingVolume.forceRebuild() )
	  {
	    InletOutletDesc ioDesc;
	    ioDesc.loadFromSTL( inputSurfacePath );
	    ioDesc.save(outputDescPath);
	  }
	myMeshingVolume.setInputInletOutletDescPath(outputDescPath);
      }
    myMeshingVolume.run();

    return 0;
}

