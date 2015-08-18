#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions;
    myoptions.add( RemeshSTL::options("") );
    myoptions.add_options()
      ("pre-process.open-surface", Feel::po::value<bool>()->default_value(false), "preprocess : open-surface ");
    myoptions.add( OpenSurface::options("open-surface") );

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_remeshstl",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));

    bool preProcessOpenSurface = boption(_name="pre-process.open-surface");

    RemeshSTL myRemeshSTL("");
    if ( preProcessOpenSurface )
      {
	if ( myRemeshSTL.inputCenterlinesPath().empty() )
	  {
	    std::cout << "no inputCenterlinesPath -> exit";
	    return 0;
	  }

	OpenSurface myOpenSurface("open-surface");
	myOpenSurface.setInputSurfacePath(myRemeshSTL.inputSurfacePath());
	myOpenSurface.setInputCenterlinesPath(myRemeshSTL.inputCenterlinesPath());

	myOpenSurface.setOutputDirectory( fs::path(myRemeshSTL.outputPath()).parent_path().string() );
	myOpenSurface.updateOutputPathFromInputFileName();
	myOpenSurface.run();

	// update RemeshSTL
	myRemeshSTL.setInputSurfacePath( myOpenSurface.outputPath() );
	if ( myOpenSurface.forceRebuild() )
	  myRemeshSTL.setForceRebuild( true );
      }

    myRemeshSTL.run();

    return 0;
}

