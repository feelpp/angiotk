#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions = SurfaceFromImage::options("");
    myoptions.
      add( SubdivideSurface::options("subdivide-surface") ).
      add( SmoothSurface::options("smooth-surface") ).
      add( RemeshSTL::options("remesh-surface") );
    myoptions.add_options()
      ("post-process.subdivide-surface", Feel::po::value<bool>()->default_value(false), "subdivide-surface")
      ("post-process.smooth-surface", Feel::po::value<bool>()->default_value(false), "smooth-surface")
      ("post-process.remesh-surface", Feel::po::value<bool>()->default_value(false), "remesh-surface");

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_surfacefromimage",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));


    bool postProcessSubdivideSurface = boption(_name="post-process.subdivide-surface");
    bool postProcessSmoothSurface = boption(_name="post-process.smooth-surface");
    bool postProcessRemeshSurface = boption(_name="post-process.remesh-surface");

    SurfaceFromImage surfaceFromImage("");

    std::string finalOutputPath = surfaceFromImage.outputPath();
    std::string finalOutputFileName = fs::path(finalOutputPath).stem().string();
    std::string lastOutputPath = finalOutputPath;
    if ( postProcessSubdivideSurface || postProcessSmoothSurface || postProcessRemeshSurface )
      {
	surfaceFromImage.setOutputPath( (fs::path(finalOutputPath).parent_path()/ fs::path(finalOutputFileName+"_surfaceFromImage.stl")).string() );
      }
    surfaceFromImage.run();


    SubdivideSurface mySubdivideSurface("subdivide-surface");
    if ( postProcessSubdivideSurface )
      {
	mySubdivideSurface.setInputSurfacePath( surfaceFromImage.outputPath() );

	if ( !postProcessSmoothSurface && !postProcessRemeshSurface )
	  mySubdivideSurface.setOutputPath( finalOutputPath );
	else
	  mySubdivideSurface.setOutputPath( (fs::path(finalOutputPath).parent_path()/ fs::path(finalOutputFileName+"_subdivideSurface.stl")).string() );

	if ( surfaceFromImage.forceRebuild() )
	  mySubdivideSurface.setForceRebuild( true );
	mySubdivideSurface.run();
      }


    SmoothSurface mySmoothSurface("smooth-surface");
    if ( postProcessSmoothSurface )
      {
	if ( postProcessSubdivideSurface )
	  mySmoothSurface.setInputSurfacePath( mySubdivideSurface.outputPath() );
	else
	  mySmoothSurface.setInputSurfacePath( surfaceFromImage.outputPath() );

	if ( !postProcessRemeshSurface )
	  mySmoothSurface.setOutputPath( finalOutputPath );
	else
	  mySubdivideSurface.setOutputPath( (fs::path(finalOutputPath).parent_path()/ fs::path(finalOutputFileName+"_smoothSurface.stl")).string() );

	if ( surfaceFromImage.forceRebuild() || ( postProcessSubdivideSurface && mySubdivideSurface.forceRebuild() ) )
	  mySmoothSurface.setForceRebuild( true );
	mySmoothSurface.run();
      }

    RemeshSTL myRemeshSurface("remesh-surface");
    if ( postProcessRemeshSurface )
      {
	myRemeshSurface.setPackageType("vmtk");

	if ( postProcessSmoothSurface )
	  myRemeshSurface.setInputSurfacePath( mySmoothSurface.outputPath() );
	else if ( postProcessSubdivideSurface )
	  myRemeshSurface.setInputSurfacePath( mySubdivideSurface.outputPath() );
	else
	  myRemeshSurface.setInputSurfacePath( surfaceFromImage.outputPath() );

	myRemeshSurface.setOutputPath( finalOutputPath );

	if ( surfaceFromImage.forceRebuild() ||
	     ( postProcessSubdivideSurface && mySubdivideSurface.forceRebuild() ) ||
	     ( postProcessSmoothSurface && mySmoothSurface.forceRebuild() ) )
	  myRemeshSurface.setForceRebuild( true );

	myRemeshSurface.run();
      }
    return 0;
}

