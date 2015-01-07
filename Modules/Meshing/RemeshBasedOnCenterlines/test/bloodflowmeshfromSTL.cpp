
#include <feel/feelcore/environment.hpp>

#include <toolboxbloodflowmesh.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions( "blood flow mesh from STL options" );
    myoptions.add_options()
      ( "compute-centerlines", po::value<bool>()->default_value( true ), "compute-centerlines" )
      ( "remesh-stl-for-centerlines", po::value<bool>()->default_value( false ), "remesh-stl-for-centerlines" )
      ( "remesh-surface", po::value<bool>()->default_value( true ), "remesh-surface" )
      ( "mesh-volume", po::value<bool>()->default_value( true ), "mesh-volume" )
      ;

    myoptions.add( ToolBoxCenterlines::options("centerlines") )
      .add( ToolBoxBloodFlowReMeshSTL::options("mesh-surface") )
      .add( ToolBoxBloodFlowMesh::options("mesh-volume") );

    Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
		     _about=about(_name="bloodflowmeshfromSTL",
				  _author="Feel++ Consortium",
				  _email="feelpp-devel@feelpp.org"));



    ToolBoxCenterlines centerlines("centerlines");

    bool remeshSTLForCenterlines=boption(_name="remesh-stl-for-centerlines");
    bool doComputeCenterlines = boption(_name="compute-centerlines");
    bool doRemeshSurface=boption(_name="remesh-surface");
    bool doMeshVolume=boption(_name="mesh-volume");

    if ( doComputeCenterlines )
    {
        if ( remeshSTLForCenterlines )
        {
	    ToolBoxBloodFlowReMeshSTL remshVMTK("mesh-surface");
	    remshVMTK.setPackageType("vmtk");
            remshVMTK.setInputPath( centerlines.inputPath() );
	    remshVMTK.updateOutputPathFromInputFileName();
            remshVMTK.run();
            centerlines.setStlFileName( remshVMTK.outputPath() );
	    centerlines.updateOutputPathFromInputFileName();
        }
        centerlines.run();
    }

    ToolBoxBloodFlowReMeshSTL remshGMSH("mesh-surface");
    remshGMSH.setPackageType("gmsh");

    if ( doRemeshSurface )
    {
        if ( remshGMSH.inputPath().empty() )
	{
            remshGMSH.setInputPath( centerlines.inputPath() );
	    remshGMSH.updateOutputPathFromInputFileName();
	}
        remshGMSH.setCenterlinesFileName( centerlines.outputPath() );
        remshGMSH.run();
    }

    ToolBoxBloodFlowMesh meshVolume("mesh-volume");
    if ( doMeshVolume )
    {
        meshVolume.setInputSTLPath( remshGMSH.outputPath() );
        meshVolume.setInputCenterlinesPath( centerlines.outputPath() );
	meshVolume.updateOutputPathFromInputFileName();
        meshVolume.run();
    }

    return 0;
}
