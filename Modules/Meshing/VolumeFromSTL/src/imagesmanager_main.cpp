#include <volumefromstl.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions = ImagesManager::options("");

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_imagesmanager",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));

    ImagesManager myImagesManager("");
    myImagesManager.run();

    return 0;
}

