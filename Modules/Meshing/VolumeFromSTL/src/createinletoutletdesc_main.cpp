
#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    //std::cout << "fs::current_path() " << fs::current_path() << "\n";
    fs::path initialPath = fs::current_path();

    po::options_description myoptions;// = CenterlinesFromSTL::options("");
    myoptions.add_options()
      ( "input.filename", po::value<std::string>()->default_value( "" ), "stl filename" );

    Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
		     _about=about(_name="meshing_centerlines",
				  _author="Feel++ Consortium",
				  _email="feelpp-devel@feelpp.org"));

    //std::cout << "fs::current_path() " << fs::current_path() << "\n";
    std::string inputPath = soption(_name="input.filename");
    if ( inputPath.empty() )
    {
        std::cout <<  "not done because inputPath is empty\n";
	return 0;
    }
    else if ( fs::path( inputPath ).is_relative() )
    {
        inputPath = ( initialPath/fs::path(inputPath) ).string();
    }
    if ( !fs::exists(inputPath) )
    {
        std::cout <<  "not done because inputPath does not exist " << inputPath << "\n";
	return 0;	
    }

    InletOutletDesc ioDesc;
    ioDesc.loadFromSTL(inputPath);

    std::string nameWithoutExt = fs::path(inputPath).stem().string();
    std::string outputPath = (fs::path(inputPath).parent_path()/fs::path(nameWithoutExt+".desc")).string();
    ioDesc.save(outputPath);

    return 0;
}

