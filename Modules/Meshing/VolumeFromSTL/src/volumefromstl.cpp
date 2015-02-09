/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <volumefromstl.hpp>

namespace Feel
{

InletOutletDesc::InletOutletDesc( std::string const& path )
{
    std::ifstream fileloaded(path.c_str(), std::ios::in);  // load file
    if( fileloaded ) // if open sucess
    {
        std::string markerLumen,markerArterialWall;
        double ptx,pty,ptz;
        while ( !fileloaded.eof() )
        {
            fileloaded >> markerLumen >> markerArterialWall >> ptx >> pty >> ptz;
            this->add( InletOutletData( markerLumen,markerArterialWall,ptx,pty,ptz ) );
	    }
        fileloaded.close();
    }
}

void
InletOutletDesc::add( InletOutletData const& data )
{
    this->push_back( data );
}

void
InletOutletDesc::loadFromSTL( std::string inputPath )
{

    std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
    std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );

    //if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )

    //fs::path gp = M_inputPath;
    std::string nameWithoutExt = fs::path(inputPath).stem().string();

    std::string outputPath = (fs::path(inputPath).parent_path()/fs::path(nameWithoutExt+".dat")).string();

    std::ostringstream __str;
    __str << pythonExecutable << " ";
    //vmtkboundaryreferencesystems -ifile myinput.stl -ofile myoutput.dat
    __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtkboundaryreferencesystems ";
    __str << "-ifile " << inputPath << " -ofile " << outputPath;

    std::cout << "---------------------------------------\n"
              << "run in system : \n" << __str.str() << "\n"
              << "---------------------------------------\n";
    auto err = ::system( __str.str().c_str() );

    if ( !fs::exists(inputPath) )
    {
        std::cout << "WARNING!, outputPath does not exist : " << outputPath <<"\n";
        return;
    }

    std::ifstream fileDat(outputPath.c_str(), std::ios::in); // load file .dat
    int ncol = 13;
    std::vector<std::string> detailCol(ncol);
    for (int k=0;k<ncol;++k)
        fileDat >> detailCol[k];


    while ( !fileDat.eof() )
    {
        std::vector<double> valueOnCurrentRow( ncol );
        fileDat >> valueOnCurrentRow[0];
        if ( fileDat.eof() ) break;
        for (int k=1;k<ncol;++k)
            fileDat >> valueOnCurrentRow[k];

        int descId = this->size();
        InletOutletData thedat( (boost::format("markerLumen%1%")%descId).str(),
                                (boost::format("markerArterialWall%1%")%descId).str(),
                                valueOnCurrentRow[0],valueOnCurrentRow[1],valueOnCurrentRow[2] );
        this->add( thedat );
    }

    fileDat.close();
    //std::cout << "number of inlet-outlet " << this->size() << "\n";
}
void
InletOutletDesc::save( std::string outputPath )
{
    std::cout << "save desc in : " << outputPath << "\n";
    std::ofstream fileWrited( outputPath, std::ios::out | std::ios::trunc);
    if( fileWrited )
    {
        for ( auto const& desc : *this )
            fileWrited << desc.markerLumen() << " "
                       << desc.markerArterialWall() << " "
                       << desc.nodeX() << " "
                       << desc.nodeY() << " "
                       << desc.nodeZ() << "\n";
        fileWrited.close();
    }
    else
        std::cout << "writing save desc fail\n";

}


CenterlinesFromSTL::CenterlinesFromSTL( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputPath( soption(_name="input.filename",_prefix=this->prefix()) ),
    M_inputInletOutletDescPath( soption(_name="input.desc.filename",_prefix=this->prefix()) ),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
    M_viewResults( boption(_name="view-results",_prefix=this->prefix() ) ),
    M_viewResultsWithSurface( boption(_name="view-results.with-surface",_prefix=this->prefix() ) )
{
    std::vector<int> sids,tids;
    if ( Environment::vm().count(prefixvm(this->prefix(),"source-ids").c_str()) )
        sids = Environment::vm()[prefixvm(this->prefix(),"source-ids").c_str()].as<std::vector<int> >();
    if ( Environment::vm().count(prefixvm(this->prefix(),"target-ids").c_str() ) )
        tids = Environment::vm()[prefixvm(this->prefix(),"target-ids").c_str()].as<std::vector<int> >();
    for ( int id : sids )
        M_sourceids.insert( id );
    for ( int id : tids )
        M_targetids.insert( id );

    if ( !M_inputPath.empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}

CenterlinesFromSTL::CenterlinesFromSTL( CenterlinesFromSTL const& e )
    :
    M_prefix( e.M_prefix ),
    M_inputPath( e.M_inputPath ),
    M_inputInletOutletDescPath( e.M_inputInletOutletDescPath ),
    M_outputPath( e.M_outputPath ),
    M_outputDirectory( e.M_outputDirectory ),
    M_targetids( e.M_targetids ),
    M_sourceids( e.M_sourceids ),
    M_forceRebuild( e.M_forceRebuild ),
    M_viewResults( e.M_viewResults )
{}

void
CenterlinesFromSTL::updateOutputPathFromInputFileName()
{
    CHECK( !M_inputPath.empty() ) << "input path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( M_outputDirectory.empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(M_outputDirectory).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(M_outputDirectory);
    else
        meshesdirectories = fs::path(M_outputDirectory);

    // get filename without extension
    fs::path gp = M_inputPath;
    std::string nameMeshFile = gp.stem().string();

    std::string newFileName = (boost::format("%1%_centerlines.vtk")%nameMeshFile ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath = outputPath.string();
}


void
CenterlinesFromSTL::run()
{
    if ( !fs::exists( this->inputPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : centerlines computation not done because this input path not exist :" << this->inputPath() << "\n";
        return;
    }
    if ( this->sourceids().empty() && this->targetids().empty() )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : centerlines computation not done because this sourceids and targetids are empty\n";
        return;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run CenterlinesFromSTL \n"
            << "---------------------------------------\n";
    coutStr << "inputPath          : " << this->inputPath() << "\n";
    coutStr << "targetids : ";
    for ( int id : this->targetids() )
        coutStr << id << " ";
    coutStr << "\n";
    coutStr << "sourceids : ";
    for ( int id : this->sourceids() )
        coutStr << id << " ";
    coutStr << "\n"
            << "output path       : " << this->outputPath() << "\n"
            << "---------------------------------------\n"
            << "---------------------------------------\n";
    std::cout << coutStr.str();
#if 0
    if ( fs::exists( this->outputPath() ) && !this->forceRebuild() )
    {
        std::cout << "already file exist, ignore centerline\n"
                  << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
        return;
    }
#endif
    fs::path directory;
    // build directories if necessary
    if ( !this->outputPath().empty() && this->worldComm().isMasterRank() )
    {
        directory = fs::path(this->outputPath()).parent_path();
        if ( !fs::exists( directory ) )
            fs::create_directories( directory );
    }
    // // wait all process
    this->worldComm().globalComm().barrier();


    //fs::path stlNamePath = fs::path(this->inputPath());
    fs::path outputFileNamePath = fs::path(this->outputPath());
    std::string name = outputFileNamePath.stem().string();
    std::string pathVTP = (directory/fs::path(name+".vtp")).string();

    // source ~/packages/vmtk/vmtk.build2/Install/vmtk_env.sh
    std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
    std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );

    if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )
    {
        std::ostringstream __str;
        __str << pythonExecutable << " ";
        __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtkcenterlines ";

        if ( !this->inputInletOutletDescPath().empty() )
        {
            InletOutletDesc iodesc( this->inputInletOutletDescPath() );

            __str <<"-seedselector pointlist ";
            __str << "-sourcepoints ";
            for ( int id : this->sourceids() )
            {
                CHECK( id < iodesc.size() ) << "id : " << id << "not valid! must be < " << iodesc.size();
                auto iodata = iodesc[id];
                __str << iodata.nodeX() << " " << iodata.nodeY() << " " << iodata.nodeZ() << " ";;
            }
            __str << "-targetpoints ";
            for ( int id : this->targetids() )
            {
                CHECK( id < iodesc.size() ) << "id : " << id << "not valid! must be < " << iodesc.size();
                auto iodata = iodesc[id];
                __str << iodata.nodeX() << " " << iodata.nodeY() << " " << iodata.nodeZ() << " ";;
            }
        }
        else
        {
            __str <<"-seedselector profileidlist ";
            __str << "-targetids ";
            for ( int id : this->targetids() )
                __str << id << " ";
            __str << "-sourceids ";
            for ( int id : this->sourceids() )
                __str << id << " ";
        }

        __str << " -ifile " << this->inputPath() << " ";
        __str << " -ofile " << pathVTP << " " //name << ".vtp "
              << " --pipe " << dirBaseVmtk << "/vmtksurfacewriter "
              << " -ifile " << pathVTP << " " //name << ".vtp "
              << " -ofile " << this->outputPath() << " " //name << ".vtk "
              << " -mode ascii ";

        std::cout << "---------------------------------------\n"
                  << "run in system : \n" << __str.str() << "\n"
                  << "---------------------------------------\n";
        auto err = ::system( __str.str().c_str() );
    }
    else
    {
        std::cout << "already file exist, ignore centerline\n";
    }


    if ( M_viewResults )
    {
        std::ostringstream ostrView;
        ostrView << pythonExecutable << " ";// << dirBaseVmtk << "/vmtk "
        if ( M_viewResultsWithSurface )
        {
            ostrView << dirBaseVmtk << "/vmtksurfacereader -ifile " << this->inputPath() << " --pipe "
                     << dirBaseVmtk << "/vmtkrenderer" << " --pipe "
                     << dirBaseVmtk << "/vmtksurfaceviewer -opacity 0.25 " << " --pipe "
                     << dirBaseVmtk << "/vmtksurfaceviewer -ifile " << this->outputPath() << " -array MaximumInscribedSphereRadius ";
        }
        else
        {
            ostrView << dirBaseVmtk << "/vmtkcenterlineviewer -ifile " << this->outputPath() << " -pointarray MaximumInscribedSphereRadius";
        }


        std::cout << "---------------------------------------\n"
                  << "run in system : \n" << ostrView.str() << "\n"
                  << "---------------------------------------\n";
        auto errView = ::system( ostrView.str().c_str() );
    }
    std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n";
}


po::options_description
CenterlinesFromSTL::options( std::string const& prefix )
{
    po::options_description myCenterlinesOptions( "Centerlines from STL for blood flow mesh options" );
    myCenterlinesOptions.add_options()
        (prefixvm(prefix,"input.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input filename" )
        (prefixvm(prefix,"input.desc.filename").c_str(), po::value<std::string>()->default_value( "" ), "inletoutlet-desc.filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"source-ids").c_str(), po::value<std::vector<int> >()->multitoken(), "(vector of int) source ids" )
        (prefixvm(prefix,"target-ids").c_str(), po::value<std::vector<int> >()->multitoken(), "(vector of int) target ids" )
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        (prefixvm(prefix,"view-results").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) view-results")
        (prefixvm(prefix,"view-results.with-surface").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) view-results with surface")

        ;
    return myCenterlinesOptions;
}




RemeshSTL::RemeshSTL( std::string prefix )
    :
    M_prefix( prefix ),
    M_packageType(soption(_name="package-type",_prefix=this->prefix())),
    M_inputPath(soption(_name="input.filename",_prefix=this->prefix())),
    M_centerlinesFileName(soption(_name="centerlines.filename",_prefix=this->prefix())),
    M_remeshNbPointsInCircle( ioption(_name="nb-points-in-circle",_prefix=this->prefix()) ),
    M_area( doption(_name="area",_prefix=this->prefix()) ),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) )
{
    CHECK( M_packageType == "gmsh" || M_packageType == "vmtk" ) << "error on packageType : " << M_packageType;
    if ( !M_inputPath.empty() && M_outputPathGMSH.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}

RemeshSTL::RemeshSTL( RemeshSTL const& e )
    :
    M_prefix( e.M_prefix ),
    M_packageType( e.M_packageType ),
    M_inputPath( e.M_inputPath ),
    M_centerlinesFileName( e.M_centerlinesFileName ),
    M_remeshNbPointsInCircle( e.M_remeshNbPointsInCircle ),
    M_area( e.M_area ),
    M_outputPathGMSH( e.M_outputPathGMSH ), M_outputPathVMTK( e.M_outputPathVMTK),
    M_outputDirectory( e.M_outputDirectory ),
    M_forceRebuild( e.M_forceRebuild )
{}


void
RemeshSTL::updateOutputPathFromInputFileName()
{
    CHECK( !M_inputPath.empty() ) << "input path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( M_outputDirectory.empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(M_outputDirectory).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(M_outputDirectory);
    else
        meshesdirectories = fs::path(M_outputDirectory);

    // get filename without extension
    fs::path gp = M_inputPath;
    std::string nameMeshFile = gp.stem().string();

    std::string specGMSH = (boost::format( "remeshGMSHpt%1%") %this->remeshNbPointsInCircle() ).str();
    std::string specVMTK = (boost::format( "remeshVMTKarea%1%") %this->area() ).str();
    std::string newFileNameGMSH = (boost::format("%1%_%2%.stl")%nameMeshFile %specGMSH ).str();
    std::string newFileNameVMTK = (boost::format("%1%_%2%.stl")%nameMeshFile %specVMTK ).str();

    M_outputPathGMSH = ( meshesdirectories / fs::path(newFileNameGMSH) ).string();
    M_outputPathVMTK = ( meshesdirectories / fs::path(newFileNameVMTK) ).string();
}

void
RemeshSTL::run()
{
    std::cout << "\n"
              << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
              << "---------------------------------------\n"
              << "run RemeshSTL \n"
              << "---------------------------------------\n";
    std::cout << "type                 : " << M_packageType << "\n"
              << "inputPath          : " << this->inputPath() << "\n";
    if ( this->packageType() == "gmsh" )
        std::cout << "centerlinesFileName  : " << this->centerlinesFileName() << "\n";
    else if ( this->packageType() == "vmtk" )
        std::cout << "area  : " << this->area() << "\n";
    std::cout << "output path       : " << this->outputPath() << "\n"
              << "---------------------------------------\n"
              << "---------------------------------------\n";

    if ( !this->forceRebuild() && fs::exists( this->outputPath() ) )
    {
        std::cout << "already file exist, ignore remeshSTL\n"
                  << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
        return;
    }


    fs::path directory;
    // build directories if necessary
    if ( !this->outputPath().empty() && this->worldComm().isMasterRank() )
    {
        directory = fs::path(this->outputPath()).parent_path();
        if ( !fs::exists( directory ) )
            fs::create_directories( directory );
    }
    // // wait all process
    this->worldComm().globalComm().barrier();


    if ( this->packageType() == "gmsh" )
        this->runGMSH();
    else if ( this->packageType() == "vmtk" )
        this->runVMTK();

    std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n";
}

void
RemeshSTL::runVMTK()
{
    CHECK( !this->inputPath().empty() ) << "inputPath is empty";

    std::ostringstream __str;
    // source ~/packages/vmtk/vmtk.build2/Install/vmtk_env.sh
    std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
    __str << pythonExecutable << " ";
    //std::string dirBaseVmtk = "/Users/vincentchabannes/packages/vmtk/vmtk.build2/Install/bin/";
    std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );
    __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtksurfaceremeshing ";
    __str << "-ifile " << this->inputPath() << " ";
    __str << "-ofile " << this->outputPath() << " ";
    __str << "-area " << this->area() << " ";
    auto err = ::system( __str.str().c_str() );

    //std::cout << "hola\n"<< __str.str() <<"\n";
}

void
RemeshSTL::runGMSH()
{
    CHECK( !this->inputPath().empty() ) << "inputPath is empty";
    CHECK( !this->centerlinesFileName().empty() ) << "centerlinesFileName is empty";

    std::ostringstream geodesc;
    geodesc << "Mesh.Algorithm = 6; //(1=MeshAdapt, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) \n"
            << "Mesh.Algorithm3D = 1; //(1=tetgen, 4=netgen, 7=MMG3D, 9=R-tree) \n";

    //geodesc << "Merge \"stl_remesh_vmtk/fluidskin3.stl\"\n";
    geodesc << "Merge \""<< this->inputPath() <<"\";\n";

    geodesc << "Field[1] = Centerline;\n";


    //geodesc << "Field[1].FileName = \"../centerlines/fluidskin3.vtk\"\n";
    geodesc << "Field[1].FileName = \"" << this->centerlinesFileName() << "\";\n";
    geodesc << "Field[1].nbPoints = "<< this->remeshNbPointsInCircle() << ";//15//25 //number of mesh elements in a circle\n";
    //geodesc << "Field[1].nbPoints = 15;//25 //number of mesh elements in a circle\n";

    //Remesh the initial stl
    geodesc << "Field[1].reMesh =1;\n"
            << "Field[1].run;\n"
            << "Background Field = 1;\n";

    fs::path outputMeshNamePath = fs::path(this->outputPath());
    std::string _name = outputMeshNamePath.stem().string();
    std::string geoname=_name+".geo";

    std::ofstream geofile( geoname.c_str() );
    geofile << geodesc.str();
    geofile.close();

    //std::cout << " this->outputFileName() " << this->outputFileName() << "\n";
    detail::generateMeshFromGeo(geoname/*remeshGeoFileName*/,this->outputPath()/*remeshMshFileName*/,2);
}

po::options_description
RemeshSTL::options( std::string const& prefix )
{
    po::options_description myMeshSurfaceOptions( "Mesh surface blood flow from STL and Centerlines options" );
    myMeshSurfaceOptions.add_options()

        ( prefixvm(prefix,"package-type").c_str(), po::value<std::string>()->default_value( "vmtk" ), "force-remesh" )
        //( prefixvm(prefix,"force-remesh").c_str(), po::value<bool>()->default_value( false ), "force-remesh" )
        ( prefixvm(prefix,"nb-points-in-circle").c_str(), po::value<int>()->default_value( 15 ), "nb-points-in-circle" )
        ( prefixvm(prefix,"area").c_str(), po::value<double>()->default_value( 0.5 ), "area" )
        ( prefixvm(prefix,"input.filename").c_str(), po::value<std::string>()->default_value( "" ), "stl.filename" )
        ( prefixvm(prefix,"centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "centerlines.filename" )

        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")

        ;
    return myMeshSurfaceOptions;
}

namespace detail
{
void
generateMeshFromGeo( std::string inputGeoName,std::string outputMeshName,int dim )
{
    //fs::path outputMeshNamePath=fs::path(Environment::findFile(outputMeshName));
    fs::path outputMeshNamePath=fs::path(outputMeshName);
    std::string _name = outputMeshNamePath.stem().string();
    //std::cout << "\n _name " << _name << "\n";
    //std::cout << " _ext " << outputMeshNamePath.extension() << "\n";

    GmshInitialize();

    GModel * M_gmodel = new GModel();

    M_gmodel->setName( _name );
#if !defined( __APPLE__ )
    M_gmodel->setFileName( _name );
#endif
    M_gmodel->readGEO( inputGeoName );
    M_gmodel->mesh( dim );

    CHECK(M_gmodel->getMeshStatus() > 0)  << "Invalid Gmsh Mesh, Gmsh status : " << M_gmodel->getMeshStatus() << " should be > 0. Gmsh mesh cannot be written to disk\n";

    //M_gmodel->writeMSH( _name+"."+formatExportMesh, 2.2, CTX::instance()->mesh.binary );
    if ( outputMeshNamePath.extension() == ".msh" )
        M_gmodel->writeMSH( outputMeshName, 2.2, CTX::instance()->mesh.binary );
    else if ( outputMeshNamePath.extension() == ".stl" )
        M_gmodel->writeSTL( outputMeshName, CTX::instance()->mesh.binary );
    else
        CHECK( false ) << "error \n";


    delete M_gmodel;
}

} //  namespace detail



VolumeMeshing::VolumeMeshing( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputSTLPath(soption(_name="input.stl.filename",_prefix=this->prefix()) ),
    M_inputCenterlinesPath(soption(_name="input.centerlines.filename",_prefix=this->prefix())),
    M_inputInletOutletDescPath(soption(_name="input.desc.filename",_prefix=this->prefix())),
    M_remeshNbPointsInCircle( ioption(_name="nb-points-in-circle",_prefix=this->prefix() ) ),
    M_extrudeWall( boption(_name="extrude-wall",_prefix=this->prefix() ) ),
    M_extrudeWallNbElemLayer( ioption(_name="extrude-wall.nb-elt-layer",_prefix=this->prefix() ) ),
    M_extrudeWallhLayer( doption(_name="extrude-wall.h-layer",_prefix=this->prefix()) ),
                                                M_outputPath(),
                                                M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
                                                M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) )
{
    if ( !M_inputSTLPath.empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}

VolumeMeshing::VolumeMeshing( VolumeMeshing const& e )
    :
    M_prefix( e.M_prefix ),
    M_inputSTLPath( e.M_inputSTLPath ),
    M_inputCenterlinesPath( e.M_inputCenterlinesPath ),
    M_inputInletOutletDescPath( e.M_inputInletOutletDescPath ),
    M_remeshNbPointsInCircle( e.M_remeshNbPointsInCircle ),
    M_extrudeWall( e.M_extrudeWall ),
    M_extrudeWallNbElemLayer( e.M_extrudeWallNbElemLayer ),
    M_extrudeWallhLayer( e.M_extrudeWallhLayer ),
    M_outputPath( e.M_outputPath ),
    M_outputDirectory( e.M_outputDirectory ),
    M_forceRebuild( e.M_forceRebuild )
{}


void
VolumeMeshing::updateOutputPathFromInputFileName()
{
    CHECK( !this->inputSTLPath().empty() ) << "input path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( M_outputDirectory.empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(M_outputDirectory).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(M_outputDirectory);
    else
        meshesdirectories = fs::path(M_outputDirectory);

    // get filename without extension
    fs::path gp = M_inputSTLPath;
    std::string nameMeshFile = gp.stem().string();


    std::string arterialWallSpecStr = (this->extrudeWall())? "WithArterialWall" : "WithoutArterialWall";
    if ( this->extrudeWall() )
        arterialWallSpecStr += (boost::format("hLayer%1%nEltLayer%2%")%this->extrudeWall_hLayer() %this->extrudeWall_nbElemLayer() ).str();

    std::string newFileName = (boost::format( "%1%_%2%pt%3%%4%.msh")
                               %nameMeshFile
                               %"volumeMesh"
                               %this->remeshNbPointsInCircle()
                               %arterialWallSpecStr
                               ).str();

    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath = outputPath.string();
}

void
VolumeMeshing::run()
{
    std::cout << "\n"
              << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
              << "---------------------------------------\n"
              << "run VolumeMeshing \n"
              << "---------------------------------------\n";
    std::cout << "input path             : " << this->inputSTLPath() << "\n"
              << "centerlinesFileName     : " << this->inputCenterlinesPath() << "\n"
              << "inletOutletDescFileName : " << this->inputInletOutletDescPath() << "\n"
              << "output path          : " << this->outputPath() << "\n"
              << "---------------------------------------\n"
              << "---------------------------------------\n";

    CHECK( !this->inputSTLPath().empty() ) << "inputSTLPath is empty";
    CHECK( !this->inputCenterlinesPath().empty() ) << "centerlinesFileName is empty";
    CHECK( !this->inputInletOutletDescPath().empty() ) << "inletOutletDescFileName is empty";

    if ( fs::exists( this->outputPath() ) && !this->forceRebuild() )
    {
        std::cout << "already file exist, ignore meshVolume\n"
                  << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
        return;
    }

        fs::path directory;
        // build directories if necessary
        if ( !this->outputPath().empty() && this->worldComm().isMasterRank() )
        {
            directory = fs::path(this->outputPath()).parent_path();
            if ( !fs::exists( directory ) )
                fs::create_directories( directory );
        }
        // // wait all process
        this->worldComm().globalComm().barrier();


        std::string volumeGeoFileName = "secondStep.geo";
        std::string volumeMshFileName = "secondStep.msh";
        this->generateGeoFor3dVolumeFromSTLAndCenterlines(volumeGeoFileName);
        detail::generateMeshFromGeo(volumeGeoFileName,this->outputPath()/*volumeMshFileName*/,3);

        std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";

}

void
VolumeMeshing::generateGeoFor3dVolumeFromSTLAndCenterlines(std::string geoname)
{
    std::ostringstream geodesc;
    geodesc << "General.ExpertMode=1;\n";
    geodesc << "Mesh.Algorithm = 6; //(1=MeshAdapt, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) \n"
            << "Mesh.Algorithm3D = 1; //(1=tetgen, 4=netgen, 7=MMG3D, 9=R-tree) \n";

    //Mesh.CharacteristicLengthFactor=0.015;

    //Merge "stl_remesh_gmsh/fluidskin_P15.stl";
    geodesc << "Merge \""<< this->inputSTLPath() <<"\";\n";

    geodesc << "Field[1] = Centerline;\n";
    // centerline file in vtk format
    geodesc << "Field[1].FileName = \"" << this->inputCenterlinesPath() << "\";\n";

    geodesc << "Field[1].nbPoints = "<< this->remeshNbPointsInCircle() << ";//15//25 //number of mesh elements in a circle\n";

    // close inlets and outlets with planar faces
    geodesc << "Field[1].closeVolume =1;\n";

    // extrude in the outward direction a vessel wall
    geodesc << "Field[1].extrudeWall ="<< ( this->extrudeWall() ? 1 : 0 ) <<";\n";
    if ( this->extrudeWall() )
        {
            geodesc << "Field[1].nbElemLayer = "<< this->extrudeWall_nbElemLayer() <<"; //number of layers\n";
            geodesc << "Field[1].hLayer = "<< this->extrudeWall_hLayer() <<";// extrusion thickness given as percent of vessel radius\n";
        }

    // identify inlets/outlets boundaries
    geodesc << "Field[1].descInletOutlet = \"" << this->inputInletOutletDescPath() << "\";\n";

    // apply field
    geodesc << "Field[1].run;\n"
            << "Background Field = 1;\n";


    std::ofstream geofile( geoname.c_str() );
    geofile << geodesc.str();
    geofile.close();
}

po::options_description
VolumeMeshing::options( std::string const& prefix )
{
    po::options_description myMeshVolumeOptions( "Mesh volume blood flow from STL and Centerlines options" );
    myMeshVolumeOptions.add_options()
        ( prefixvm(prefix,"input.stl.filename").c_str(), po::value<std::string>()->default_value( "" ), "stl.filename" )
        ( prefixvm(prefix,"input.centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "centerlines.filename" )
        ( prefixvm(prefix,"input.desc.filename").c_str(), po::value<std::string>()->default_value( "" ), "inletoutlet-desc.filename" )
        ( prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        ( prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        ( prefixvm(prefix,"nb-points-in-circle").c_str(), po::value<int>()->default_value( 15 ), "nb-points-in-circle" )
        ( prefixvm(prefix,"extrude-wall").c_str(),po::value<bool>()->default_value( true ), "extrude-wall" )
        ( prefixvm(prefix,"extrude-wall.nb-elt-layer").c_str(), po::value<int>()->default_value( 2 ), "nb-elt-layer" )
        ( prefixvm(prefix,"extrude-wall.h-layer").c_str(), po::value<double>()->default_value( 0.2 ), "h-layer" )
        ;
    return myMeshVolumeOptions;
}





} // namespace Feel
