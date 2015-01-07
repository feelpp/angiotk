#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions = CenterlinesFromSTL::options("");

    Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
		     _about=about(_name="centerlines",
				  _author="Feel++ Consortium",
				  _email="feelpp-devel@feelpp.org"));


    CenterlinesFromSTL centerlines("");
    centerlines.run();

    return 0;
}

