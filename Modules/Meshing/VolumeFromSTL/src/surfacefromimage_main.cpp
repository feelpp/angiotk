#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions = SurfaceFromImage::options("");



#if 0
    Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
		     _about=about(_name="meshing_surfacefromimage",
				  _author="Feel++ Consortium",
				  _email="feelpp-devel@feelpp.org"));
#endif
    AngioTkEnvironment env2( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_surfacefromimage",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));


    SurfaceFromImage surfaceFromImage("");
    surfaceFromImage.run();

    return 0;
}

