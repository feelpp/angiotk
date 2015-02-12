#include <volumefromstl.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions = OpenSurface::options("");

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_opensurface",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));

    OpenSurface myOpenSurface("");
    myOpenSurface.run();

    return 0;
}

