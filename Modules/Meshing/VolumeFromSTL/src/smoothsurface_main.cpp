#include <volumefromstl.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions = SmoothSurface::options("");

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_smoothsurface",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));

    SmoothSurface mySmoothSurface("");
    mySmoothSurface.run();

    return 0;
}

