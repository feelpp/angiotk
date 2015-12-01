#include <volumefromstl.hpp>

int main( int argc, char** argv )
{
    using namespace AngioTk;

    po::options_description myoptions = TubularExtension::options("");

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_tubularextension",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));

    TubularExtension myTubularExtension("");
    myTubularExtension.run();

    return 0;
}

