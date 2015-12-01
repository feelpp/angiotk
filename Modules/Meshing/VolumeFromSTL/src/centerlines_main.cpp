#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace AngioTk;

    po::options_description myoptions;
    myoptions.add( CenterlinesFromSurface::options("") );
    myoptions.add_options()
      ("pre-process.subdivide-surface", Feel::po::value<bool>()->default_value(false), "use subdivide-surface ")
      ("post-process.convert-centerlines", Feel::po::value<bool>()->default_value(false), "convert-centerlines ");
    myoptions.add( SubdivideSurface::options("subdivide-surface") ).add( CenterlinesManager::options("convert-centerlines") );

    AngioTkEnvironment env( Feel::_argc=argc, Feel::_argv=argv,
			    Feel::_desc=myoptions,
			    Feel::_about=Feel::about(Feel::_name="meshing_centerlines",
						     Feel::_author="Feel++ Consortium",
						     Feel::_email="feelpp-devel@feelpp.org"));

    bool preProcessSubdivideSurface = boption(Feel::_name="pre-process.subdivide-surface");
    bool postProcessConvertCenterlines = boption(Feel::_name="post-process.convert-centerlines");

    CenterlinesFromSurface centerlines("");
    std::string finalOutputPath = centerlines.outputPath();
    std::string finalOutputFileName = fs::path(finalOutputPath).stem().string();

    if ( preProcessSubdivideSurface )
      {
	SubdivideSurface mySubdivideSurface("subdivide-surface");
	mySubdivideSurface.setInputSurfacePath( centerlines.inputSurfacePath() );
	mySubdivideSurface.setOutputDirectory( fs::path(centerlines.outputPath()).parent_path().string() );
	mySubdivideSurface.updateOutputPathFromInputFileName();
	mySubdivideSurface.run();

	centerlines.setInputSurfacePath( mySubdivideSurface.outputPath() );
	//centerlines.updateOutputPathFromInputFileName();
      }

    if ( postProcessConvertCenterlines )
      {
	centerlines.setOutputPath( (fs::path(finalOutputPath).parent_path()/ fs::path(finalOutputFileName+"_temporary.vtk")).string() );
      }
    centerlines.run();

    if ( postProcessConvertCenterlines )
      {
	CenterlinesManager myCM("convert-centerlines");
	myCM.setInputSurfacePath( centerlines.inputSurfacePath() );
	myCM.setInputCenterlinesPath( centerlines.outputPath() );
	myCM.setOutputPath( finalOutputPath );
	if ( centerlines.forceRebuild() )
	  myCM.setForceRebuild( true );
	myCM.run();
      }


    return 0;
}

