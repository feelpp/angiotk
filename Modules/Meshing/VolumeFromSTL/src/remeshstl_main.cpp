#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions = RemeshSTL::options("");

    Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
		     _about=about(_name="meshing_remeshstl",
				  _author="Feel++ Consortium",
				  _email="feelpp-devel@feelpp.org"));

    RemeshSTL myRemeshSTL("");
    myRemeshSTL.run();

    return 0;
}

