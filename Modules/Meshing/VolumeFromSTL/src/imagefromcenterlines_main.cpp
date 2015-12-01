#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace AngioTk;

    po::options_description myoptions = ImageFromCenterlines::options("");

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_imagefromcenterlines",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));

    ImageFromCenterlines imageFromCenterline("");
    imageFromCenterline.run();

    return 0;
}

