#include <volumefromstl.hpp>
#include <centerlinesmanagerwindowinteractor.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions;
    myoptions.add_options()
      ("input.centerlines.path", po::value<std::vector<std::string>>()->multitoken(), "(vector<string>) input centerline filename" )
      ("input.surface.path", po::value<std::string>()->default_value( "" ), "(string) input surface filename" )
      ("input.point-set.path", po::value<std::string>()->default_value( "" ), "(string) input point-set filename" )
      ("input.point-pair.path", po::value<std::string>()->default_value( "" ), "input.pointpair.filename" )
      ("window-width", Feel::po::value<int>()->default_value(1024), "(int) window width")
      ("window-height", Feel::po::value<int>()->default_value(768), "(int) window height")
      ;

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_centerlinesmanagergui",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));

    std::vector<std::string> inputCenterlinesPath;
    if ( AngioTkEnvironment::vm().count("input.centerlines.path" ) )
      inputCenterlinesPath = Environment::vm()["input.centerlines.path"].as<std::vector<std::string> >();

    std::string inputSurfacePath = AngioTkEnvironment::expand( soption(_name="input.surface.path") );
    std::string inputPointSetPath = AngioTkEnvironment::expand( soption(_name="input.point-set.path") );
    std::string inputPointPairPath = AngioTkEnvironment::expand( soption(_name="input.point-pair.path") );

    CenterlinesManagerWindowInteractor windowInteractor;
    windowInteractor.setWindowWidth( ioption(_name="window-width" ) );
    windowInteractor.setWindowHeight( ioption(_name="window-height" ) );

    if ( !inputSurfacePath.empty() && fs::exists( inputSurfacePath ) )
    {
      windowInteractor.setInputSurfacePath( inputSurfacePath );
    }
    else
    {
        std::cout << "WARNING: The surface data you specified is either empty (" << inputSurfacePath.empty() 
                  << ") or the path doesn't exist (" << (!fs::exists( inputSurfacePath )) << ")." << std::endl
                  << "Please specify it with the --input.surface.path option." << std::endl;
    }
    if ( !inputPointSetPath.empty() && fs::exists( inputPointSetPath ) )
      windowInteractor.setInputPointSetPath( inputPointSetPath );
    if ( !inputPointPairPath.empty() && fs::exists( inputPointPairPath ) )
      windowInteractor.setInputPointPairPath( inputPointPairPath );

    std::vector<std::string> centerlinesPath;
    for (int k = 0;k<inputCenterlinesPath.size();++k)
      if ( !inputCenterlinesPath[k].empty() && fs::exists( inputCenterlinesPath[k] ) )
	centerlinesPath.push_back( inputCenterlinesPath[k] );
    windowInteractor.setInputCenterlinesPath( centerlinesPath );

    windowInteractor.run();

    return 0;
}

