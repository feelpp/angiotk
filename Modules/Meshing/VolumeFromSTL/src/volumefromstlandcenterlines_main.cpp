#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace AngioTk;

    po::options_description myoptions = VolumeMeshing::options("");

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_volumefromstlandcenterlines",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));


    VolumeMeshing myMeshingVolume("");
    if ( myMeshingVolume.inputInletOutletDescPath().empty() )
      {
	std::string inputSurfacePath = myMeshingVolume.inputSurfacePath();
	std::string nameWithoutExt = fs::path(inputSurfacePath).stem().string();
	std::string outputDescPath = (fs::path(myMeshingVolume.outputPath()).parent_path()/fs::path(nameWithoutExt+".desc")).string();

    if( !fs::exists(myMeshingVolume.inputSurfacePath()))
    {
        std::cout << "WARNING: The input surface path file does not exists." << std::endl
                  << "Please set it with the \"input.surface.filename\" option." << std::endl;
    }

	if ( !fs::exists( outputDescPath ) || myMeshingVolume.forceRebuild() )
	  {
	    InletOutletDesc ioDesc;
        /* Handle non-zero returns (errors) */
	    if(ioDesc.loadFromSTL( inputSurfacePath ))
        {
            return 1;
        }
	    ioDesc.save(outputDescPath);
	  }
	myMeshingVolume.setInputInletOutletDescPath(outputDescPath);
      }

    /* Create json file path */
    std::string baseDir = fs::path(myMeshingVolume.inputInletOutletDescPath()).parent_path().string();
    std::string jsonFile = fs::path(myMeshingVolume.inputInletOutletDescPath()).stem().string() + ".json";
	std::string outputJSONPath = baseDir + "/" + jsonFile;

    /* write json file */
    InletOutletDesc ioDesc(myMeshingVolume.inputInletOutletDescPath());
    ioDesc.saveJSON(outputJSONPath);

    myMeshingVolume.run();

    return 0;
}

