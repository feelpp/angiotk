#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace AngioTk;

    po::options_description myoptions;
    myoptions.add( RemeshSurface::options("") );
    myoptions.add_options()
      ("pre-process.open-surface", Feel::po::value<bool>()->default_value(false), "preprocess : open-surface ");
    myoptions.add( OpenSurface::options("open-surface") );

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_remeshstl",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));

    bool preProcessOpenSurface = boption(_name="pre-process.open-surface");

    RemeshSurface myRemeshSurface("");
    if ( preProcessOpenSurface )
      {
	if ( myRemeshSurface.inputCenterlinesPath().empty() )
	  {
	    std::cout << "no inputCenterlinesPath -> exit";
	    return 0;
	  }

	OpenSurface myOpenSurface("open-surface");
	myOpenSurface.setInputSurfacePath(myRemeshSurface.inputSurfacePath());
	myOpenSurface.setInputCenterlinesPath(myRemeshSurface.inputCenterlinesPath());

	myOpenSurface.setOutputDirectory( fs::path(myRemeshSurface.outputPath()).parent_path().string() );
	myOpenSurface.updateOutputPathFromInputFileName();
	myOpenSurface.run();

	// update RemeshSurface
	myRemeshSurface.setInputSurfacePath( myOpenSurface.outputPath() );
	if ( myOpenSurface.forceRebuild() )
	  myRemeshSurface.setForceRebuild( true );
      }

    myRemeshSurface.run();

    return 0;
}

