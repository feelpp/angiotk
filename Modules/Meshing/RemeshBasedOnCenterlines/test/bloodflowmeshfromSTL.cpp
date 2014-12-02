
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

    myoptions.add( ToolBoxCenterlines::options() )
      .add( ToolBoxBloodFlowReMeshSTL::options() )
      .add( ToolBoxBloodFlowMesh::options() );

    Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
		     _about=about(_name="bloodflowmeshfromSTL",
				  _author="Feel++ Consortium",
				  _email="feelpp-devel@feelpp.org"));


#if 1
    ToolBoxCenterlines centerlines;

    bool remeshSTLForCenterlines=boption(_name="remesh-stl-for-centerlines");
    bool doComputeCenterlines = boption(_name="compute-centerlines");
    bool doRemeshSurface=boption(_name="remesh-surface");
    bool doMeshVolume=boption(_name="mesh-volume");

    if ( doComputeCenterlines )
    {
        if ( remeshSTLForCenterlines )
        {
            ToolBoxBloodFlowReMeshSTL remshVMTK("vmtk");
            remshVMTK.setStlFileName( centerlines.stlFileName() );
            remshVMTK.run();
            centerlines.setStlFileName( remshVMTK.outputFileName() );
        }
        if ( false )
        {
            // for cut3
            centerlines.setTargetids( { 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12 } );
            centerlines.setSourceids( { 9 } );
        }
        else
        {
            centerlines.setTargetids( { 1, 2 } );
            centerlines.setSourceids( { 0 } );

        }
        centerlines.run();
    }

    ToolBoxBloodFlowReMeshSTL remshGMSH("gmsh");
    if ( doRemeshSurface )
    {
        if ( remshGMSH.stlFileName().empty() )
            remshGMSH.setStlFileName( centerlines.stlFileName() );
        remshGMSH.setCenterlinesFileName( centerlines.outputFileName() );
        remshGMSH.run();
    }

    ToolBoxBloodFlowMesh meshVolume;
    if ( doMeshVolume )
    {
        meshVolume.setStlFileName( remshGMSH.outputFileName() );
        //if ( doComputeCenterlines )
        meshVolume.setCenterlinesFileName( centerlines.outputFileName() );
        meshVolume.run();
    }
#endif
    return 0;
}
