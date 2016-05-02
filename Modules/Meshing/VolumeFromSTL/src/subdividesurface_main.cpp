#include <volumefromstl.hpp>

int main( int argc, char** argv )
{
    using namespace AngioTk;

    po::options_description myoptions = SubdivideSurface::options("");

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_subdividesurface",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));

    SubdivideSurface mySubdivideSurface("");
    mySubdivideSurface.run();

    return 0;
}

