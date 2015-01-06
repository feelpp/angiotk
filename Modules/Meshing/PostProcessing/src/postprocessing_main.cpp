
#include <feel/feelcore/environment.hpp>

#include <postprocessing.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions;
    myoptions.add( ExtractSubMeshFromFSIMesh<>::options( "extract-submesh" ) );
    myoptions.add( MeshPartitioner::options( "mesh-partitioner" ) );

    Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
		     _about=about(_name="postprocessing",
				  _author="Feel++ Consortium",
				  _email="feelpp-devel@feelpp.org"));

    // extraction
    ExtractSubMeshFromFSIMesh< Mesh<Simplex<3> > > myExtract( "extract-submesh");
    myExtract.run();

    // partitioning
    MeshPartitioner myPartitioner("mesh-partitioner");
    myPartitioner.setInputPath( myExtract.outputPathLumen() );
    myPartitioner.updateOutputPathFromInputPath();
    myPartitioner.run();
    myPartitioner.setInputPath( myExtract.outputPathArterialWall() );
    myPartitioner.updateOutputPathFromInputPath();
    myPartitioner.run();

    return 0;
}
