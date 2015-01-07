#include <feel/feelcore/environment.hpp>

#include <toolboxbloodflowmesh.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description myoptions = ToolBoxBloodFlowReMeshSTL::options("");

    Environment env( _argc=argc, _argv=argv,
                     _desc=myoptions,
		     _about=about(_name="remeshstl",
				  _author="Feel++ Consortium",
				  _email="feelpp-devel@feelpp.org"));

    ToolBoxBloodFlowReMeshSTL myRemeshSTL("");
    myRemeshSTL.run();

    return 0;
}

