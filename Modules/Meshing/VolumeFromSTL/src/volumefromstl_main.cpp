
#include <feel/feelcore/environment.hpp>

#include <volumefromstl.hpp>


int main( int argc, char** argv )
{
    using namespace AngioTk;

    po::options_description myoptions( "blood flow mesh from STL options" );
    myoptions.add_options()
      ( "compute-centerlines", po::value<bool>()->default_value( true ), "(bool) compute-centerlines" )
      ( "remesh-stl-for-centerlines", po::value<bool>()->default_value( false ), "(bool) remesh-stl-for-centerlines" )
      ( "remesh-surface", po::value<bool>()->default_value( true ), "(bool) remesh-surface" )
      ( "mesh-volume", po::value<bool>()->default_value( true ), "(bool) mesh-volume" )
      ( "force-rebuild", po::value<bool>()->default_value( false ), "(bool) force rebuild" )
      ;

    myoptions.add( CenterlinesFromSurface::options("centerlines") )
      .add( RemeshSurface::options("mesh-surface") )
      .add( VolumeMeshing::options("mesh-volume") );

    AngioTkEnvironment env( _argc=argc, _argv=argv,
			    _desc=myoptions,
			    _about=about(_name="meshing_volumefromstl",
					 _author="Feel++ Consortium",
					 _email="feelpp-devel@feelpp.org"));


    bool doForceRebuild = boption(_name="force-rebuild");
    bool doRemeshSurfaceForCenterlines=boption(_name="remesh-stl-for-centerlines") || doForceRebuild;
    bool doComputeCenterlines = boption(_name="compute-centerlines") || doForceRebuild;
    bool doRemeshSurface=boption(_name="remesh-surface") || doForceRebuild;
    bool doMeshVolume=boption(_name="mesh-volume") || doForceRebuild;

    CenterlinesFromSurface centerlines("centerlines");
    if ( doForceRebuild )
      centerlines.setForceRebuild(true);
    if ( doComputeCenterlines )
    {
        if ( doRemeshSurfaceForCenterlines )
        {
	    RemeshSurface remshVMTK("mesh-surface");
	    if ( doForceRebuild )
	      remshVMTK.setForceRebuild(true);
	    remshVMTK.setPackageType("vmtk");
            remshVMTK.setInputSurfacePath( centerlines.inputSurfacePath() );
	    remshVMTK.updateOutputPathFromInputFileName();
            remshVMTK.run();
            centerlines.setInputSurfacePath( remshVMTK.outputPath() );
	    centerlines.updateOutputPathFromInputFileName();
        }
        centerlines.run();
    }

    RemeshSurface remshGMSH("mesh-surface");
    if ( doForceRebuild )
      remshGMSH.setForceRebuild(true);
    remshGMSH.setPackageType("gmsh");
    if ( doRemeshSurface )
    {
        if ( remshGMSH.inputSurfacePath().empty() )
	{
            remshGMSH.setInputSurfacePath( centerlines.inputSurfacePath() );
	    remshGMSH.updateOutputPathFromInputFileName();
	}
        remshGMSH.setInputCenterlinesPath( centerlines.outputPath() );
        remshGMSH.run();
    }

    VolumeMeshing meshVolume("mesh-volume");
    if ( doForceRebuild )
      meshVolume.setForceRebuild(true);
    if ( doMeshVolume )
    {
        meshVolume.setInputSurfacePath( remshGMSH.outputPath() );
        meshVolume.setInputCenterlinesPath( centerlines.outputPath() );
	meshVolume.updateOutputPathFromInputFileName();
        meshVolume.run();
    }

    return 0;
}
