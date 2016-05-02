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
        if ( myRemeshSurface.inputSurfacePath().empty() )
        {
            std::cout << "inputSurfacePath is not specified. The program cannot continue." << std::endl
                << "Please set it with the \"input.surface.filename\" option." << std::endl;
            return 1;
        }
        if ( myRemeshSurface.inputCenterlinesPath().empty() )
        {
            std::cout << "inputCenterlinesPath is not specified. The program cannot continue." << std::endl
                << "Please set it with the \"gmsh.centerlines.filename\" option." << std::endl;
            return 1;
        }

        OpenSurface myOpenSurface("open-surface");
        myOpenSurface.setInputSurfacePath(myRemeshSurface.inputSurfacePath());
        myOpenSurface.setInputCenterlinesPath(myRemeshSurface.inputCenterlinesPath());

        myOpenSurface.setOutputDirectory( fs::path(myRemeshSurface.outputPath()).parent_path().string() );
        myOpenSurface.updateOutputPathFromInputFileName();
        /* Checking for errors, exiting if on is returned */
        if(myOpenSurface.run())
        {
            return 1;
        }

        // update RemeshSurface
        myRemeshSurface.setInputSurfacePath( myOpenSurface.outputPath() );
        if ( myOpenSurface.forceRebuild() )
            myRemeshSurface.setForceRebuild( true );
    }

    /* Checking for errors, exiting if on is returned */
    if(myRemeshSurface.run())
    {
        return 1;
    }

    return 0;
}

