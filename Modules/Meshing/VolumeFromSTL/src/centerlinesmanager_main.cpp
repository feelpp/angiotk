#include <volumefromstl.hpp>

int main( int argc, char** argv )
{
    using namespace AngioTk;
    using namespace Feel;

    po::options_description myoptions = CenterlinesManager::options("");

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_centerlinesmanager",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));

    CenterlinesManager myCM("");
    myCM.run();

    return 0;
}

