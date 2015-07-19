#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions;
    myoptions.add( CenterlinesFromSTL::options("") );
    myoptions.add_options()
      ("pre-process.subdivide-surface", Feel::po::value<bool>()->default_value(false), "use subdivide-surface ")
      ("post-process.convert-centerlines", Feel::po::value<bool>()->default_value(false), "convert-centerlines ");
      myoptions.add( SubdivideSurface::options("subdivide-surface") ).add( CenterlinesManager::options("convert-centerlines") );
    Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
		     _about=about(_name="meshing_centerlines",
				  _author="Feel++ Consortium",
				  _email="feelpp-devel@feelpp.org"));

    bool preProcessSubdivideSurface = boption(_name="pre-process.subdivide-surface");
    bool postProcessConvertCenterlines = boption(_name="post-process.convert-centerlines");

    CenterlinesFromSTL centerlines("");
    std::string finalOutputPath = centerlines.outputPath();
    std::string finalOutputFileName = fs::path(finalOutputPath).stem().string();

    if ( preProcessSubdivideSurface )
      {
	SubdivideSurface mySubdivideSurface("subdivide-surface");
	mySubdivideSurface.setInputPath( centerlines.inputPath() );
	mySubdivideSurface.setOutputDirectory( fs::path(centerlines.outputPath()).parent_path().string() );
	mySubdivideSurface.updateOutputPathFromInputFileName();
	mySubdivideSurface.run();

	centerlines.setStlFileName( mySubdivideSurface.outputPath() );
	//centerlines.updateOutputPathFromInputFileName();
      }

    if ( postProcessConvertCenterlines )
      {
	centerlines.setOutputPath( (fs::path(finalOutputPath).parent_path()/ fs::path(finalOutputFileName+"_vmtkformat.vtk")).string() );
      }
    centerlines.run();

    if ( postProcessConvertCenterlines )
      {
	CenterlinesManager myCM("convert-centerlines");
	myCM.setInputSurfacePath( centerlines.inputPath() );
	myCM.setInputCenterlinesPath( centerlines.outputPath() );
	myCM.setOutputPath( finalOutputPath );
	myCM.run();
      }


    return 0;
}

