/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <volumefromstl.hpp>
#include <AngioTkCenterlineField.h>

namespace detail
{
Feel::fs::path AngioTkEnvironment::S_pathInitial;
boost::shared_ptr<Feel::Environment> AngioTkEnvironment::S_feelEnvironment;
}

namespace Feel
{

InletOutletDesc::InletOutletDesc( std::string const& path )
{
    std::ifstream fileloaded(path.c_str(), std::ios::in);  // load file
    if( fileloaded ) // if open sucess
    {
        std::string bctype,markerLumen,markerArterialWall;
        double ptx,pty,ptz;
        while ( !fileloaded.eof() )
        {
            fileloaded >> bctype;
            if ( fileloaded.eof() ) break;
            fileloaded >> markerLumen >> markerArterialWall >> ptx >> pty >> ptz;
            CHECK( bctype == "INLET" || bctype == "OUTLET" ) << "invalid type " << bctype;
            BCType mybctype = ( bctype == "INLET" )? BCType::BC_INLET : BCType::BC_OUTLET;
            this->add( InletOutletData( mybctype,markerLumen,markerArterialWall,ptx,pty,ptz ) );
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
        InletOutletData thedat( BCType::BC_INLET,
                                (boost::format("markerLumen%1%")%descId).str(),
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
            fileWrited << desc.bcTypeStr() << " "
                       << desc.markerLumen() << " "
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
    M_inputCenterlinesPointSetPath( soption(_name="input.pointset.filename",_prefix=this->prefix()) ),
    M_inputInletOutletDescPath( soption(_name="input.desc.filename",_prefix=this->prefix()) ),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
    M_useInteractiveSelection( boption(_name="use-interactive-selection",_prefix=this->prefix()) ),
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

std::tuple< std::vector<std::vector<double> >, std::vector<std::vector<double> > >
CenterlinesFromSTL::loadFromCenterlinesPointSetFile()
{
    std::vector<std::vector<double> > sourcePts, targetPts;

    std::ifstream fileLoaded( this->inputCenterlinesPointSetPath(), std::ios::in);
    while ( !fileLoaded.eof() )
      {
          std::vector<double> pt(3);
          double radius;
          int typePt = -1;
          fileLoaded >> typePt;
          if ( fileLoaded.eof() ) break;
          fileLoaded >> pt[0] >> pt[1] >> pt[2] >> radius;
#if 1
          if ( typePt == 0 )
              sourcePts.push_back(pt);
          else if ( typePt == 1 )
              targetPts.push_back(pt);
#else
          if ( typePt == 0 )
              {
              sourcePts.push_back(pt);
              sourcePts.push_back(pt);
              sourcePts.push_back(pt);
              }
          //else if ( typePt == 1 )
              {
              targetPts.push_back(pt);
              targetPts.push_back(pt);
              targetPts.push_back(pt);
              }

#endif
      }
    fileLoaded.close();  // on ferme le fichier

    auto res = std::make_tuple( sourcePts, targetPts );
    return res;
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
    if ( this->inputCenterlinesPointSetPath().empty() && this->sourceids().empty() && this->targetids().empty() )
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
    if ( !M_useInteractiveSelection )
    {
        coutStr << "targetids : ";
        for ( int id : this->targetids() )
            coutStr << id << " ";
        coutStr << "\n";
        coutStr << "sourceids : ";
        for ( int id : this->sourceids() )
            coutStr << id << " ";
        coutStr << "\n";
    }
    coutStr << "output path       : " << this->outputPath() << "\n"
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


        if ( M_useInteractiveSelection )
        {
            if ( false )
            __str <<"-seedselector openprofiles ";
            else
            __str <<"-seedselector pickpoint ";
        }
        else if ( !this->inputCenterlinesPointSetPath().empty() )
        {
            auto pointset = loadFromCenterlinesPointSetFile();
            __str <<"-seedselector pointlist ";
            __str << "-sourcepoints ";
            for ( std::vector<double> const& pt : std::get<0>(pointset) ) // source
            {
                __str << pt[0] << " " << pt[1] << " " << pt[2] << " ";
            }
            __str << "-targetpoints ";
            for ( std::vector<double> const& pt : std::get<1>(pointset) ) // target
            {
                __str << pt[0] << " " << pt[1] << " " << pt[2] << " ";
            }
        }
        else if ( !this->inputInletOutletDescPath().empty() )
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
        (prefixvm(prefix,"input.pointset.filename").c_str(), po::value<std::string>()->default_value( "" ), "input.pointset.filename" )
        (prefixvm(prefix,"input.desc.filename").c_str(), po::value<std::string>()->default_value( "" ), "inletoutlet-desc.filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"use-interactive-selection").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) use-interactive-selection")

        (prefixvm(prefix,"source-ids").c_str(), po::value<std::vector<int> >()->multitoken(), "(vector of int) source ids" )
        (prefixvm(prefix,"target-ids").c_str(), po::value<std::vector<int> >()->multitoken(), "(vector of int) target ids" )
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")

        (prefixvm(prefix,"view-results").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) view-results")
        (prefixvm(prefix,"view-results.with-surface").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) view-results with surface")
        ;
    return myCenterlinesOptions;
}


CenterlinesManager::CenterlinesManager( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputPath( soption(_name="input.filename",_prefix=this->prefix()) ),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) )
{
    std::vector<int> removeids;
    if ( Environment::vm().count(prefixvm(this->prefix(),"remove-branch-ids").c_str()) )
        removeids = Environment::vm()[prefixvm(this->prefix(),"remove-branch-ids").c_str()].as<std::vector<int> >();
    for ( int id : removeids )
        M_removeBranchIds.insert( id );

    if ( !M_inputPath.empty() && fs::path(M_inputPath).is_relative() )
        M_inputPath = (AngioTkEnvironment::pathInitial()/fs::path(M_inputPath) ).string();

    if ( !M_inputPath.empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}
CenterlinesManager::CenterlinesManager( CenterlinesManager const& e )
    :
    M_prefix( e.M_prefix ),
    M_inputPath( e.M_inputPath ),
    M_outputDirectory( e.M_outputDirectory ), M_outputPath( e.M_outputPath ),
    M_forceRebuild( e.M_forceRebuild )
{}

void
CenterlinesManager::updateOutputPathFromInputFileName()
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

    std::string newFileName = (boost::format("%1%_up.vtk")%nameMeshFile ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath = outputPath.string();
}

void
CenterlinesManager::run()
{
    if ( !fs::exists( this->inputPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : Centerlines Manager not run because this input path not exist :" << this->inputPath() << "\n";
        return;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run CenterlinesManager \n"
            << "---------------------------------------\n";
    coutStr << "inputPath          : " << this->inputPath() << "\n";
    if ( M_removeBranchIds.size() > 0 )
    {
        coutStr << "remove branch ids :";
        for ( int id : M_removeBranchIds )
            coutStr << " " << id;
        coutStr << "\n";
    }
    coutStr << "output path       : " << this->outputPath() << "\n"
            << "---------------------------------------\n"
            << "---------------------------------------\n";
    std::cout << coutStr.str();


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


    if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )
    {
        GmshInitialize();
        // if(!Msg::GetGmshClient())
        CTX::instance()->terminal = 1;
        //GmshBatch();

        int verbosityLevel = 5;
        Msg::SetVerbosity( verbosityLevel );

        AngioTkCenterline centerlinesTool;
        centerlinesTool.updateCenterlinesFromFile( this->inputPath() );
        centerlinesTool.removeBranchIds( M_removeBranchIds );
        centerlinesTool.addFieldBranchIds();
        centerlinesTool.writeCenterlinesVTK( this->outputPath() );
    }
}
po::options_description
CenterlinesManager::options( std::string const& prefix )
{
    po::options_description myCenterlinesManagerOptions( "Centerlines Manager from Image options" );

    myCenterlinesManagerOptions.add_options()
        (prefixvm(prefix,"input.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"remove-branch-ids").c_str(), po::value<std::vector<int> >()->multitoken(), "(vector of int) remove branch ids" )
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        ;
    return myCenterlinesManagerOptions;
}



ImageFromCenterlines::ImageFromCenterlines( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputPath( soption(_name="input.filename",_prefix=this->prefix()) ),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_dimX( doption(_name="dim.x",_prefix=this->prefix() ) ),
    M_dimY( doption(_name="dim.y",_prefix=this->prefix() ) ),
    M_dimZ( doption(_name="dim.z",_prefix=this->prefix() ) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) )
{
    if ( !M_inputPath.empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}
ImageFromCenterlines::ImageFromCenterlines( ImageFromCenterlines const& e )
    :
    M_prefix( e.M_prefix ),
    M_inputPath( e.M_inputPath ),
    M_outputDirectory( e.M_outputDirectory ), M_outputPath( e.M_outputPath ),
    M_dimX( e.M_dimX ), M_dimY( e.M_dimY ), M_dimZ( e.M_dimZ ),
    M_forceRebuild( e.M_forceRebuild )
{}

void
ImageFromCenterlines::updateOutputPathFromInputFileName()
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

    std::string newFileName = (boost::format("%1%_%2%-%3%-%4%.mha")%nameMeshFile %M_dimX %M_dimY %M_dimZ ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath = outputPath.string();
}

void
ImageFromCenterlines::run()
{
    if ( !fs::exists( this->inputPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : image building not done because this input path not exist :" << this->inputPath() << "\n";
        return;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run ImageFromCenterlines \n"
            << "---------------------------------------\n";
    coutStr << "inputPath          : " << this->inputPath() << "\n";
    coutStr << "dimensions : [" << M_dimX << "," << M_dimY << "," << M_dimZ << "]\n";
    coutStr << "output path       : " << this->outputPath() << "\n"
            << "---------------------------------------\n"
            << "---------------------------------------\n";
    std::cout << coutStr.str();


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


    if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )
    {
        std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
        std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );

        std::ostringstream __str;
        __str << pythonExecutable << " ";
        __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtkcenterlinemodeller ";
        __str << "-ifile " << this->inputPath() << " -radiusarray MaximumInscribedSphereRadius "
              << "-negate 1 -dimensions " << M_dimX << " " << M_dimY << " " << M_dimZ << " ";
        __str << "-ofile " << this->outputPath();

        std::cout << "---------------------------------------\n"
                  << "run in system : \n" << __str.str() << "\n"
                  << "---------------------------------------\n";
        auto err = ::system( __str.str().c_str() );
    }
}

po::options_description
ImageFromCenterlines::options( std::string const& prefix )
{
    po::options_description myImageFromCenterlinesOptions( "Create Image from Centerlines options" );

    myImageFromCenterlinesOptions.add_options()
        (prefixvm(prefix,"input.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"dim.x").c_str(), Feel::po::value<double>()->default_value(64), "(bool) force-rebuild")
        (prefixvm(prefix,"dim.y").c_str(), Feel::po::value<double>()->default_value(64), "(bool) force-rebuild")
        (prefixvm(prefix,"dim.z").c_str(), Feel::po::value<double>()->default_value(64), "(bool) force-rebuild")
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        ;
    return myImageFromCenterlinesOptions;
}

SurfaceFromImage::SurfaceFromImage( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputPath( soption(_name="input.filename",_prefix=this->prefix()) ),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_hasThresholdLower(false), M_hasThresholdUpper(false),
    M_thresholdLower(0.0),M_thresholdUpper(0.0),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) )
{

    if ( !M_inputPath.empty() && fs::path(M_inputPath).is_relative() )
        M_inputPath = (AngioTkEnvironment::pathInitial()/fs::path(M_inputPath) ).string();

    if ( Environment::vm().count(prefixvm(this->prefix(),"threshold.lower").c_str()) )
    {
        M_hasThresholdLower = true;
        M_thresholdLower = doption(_name="threshold.lower",_prefix=this->prefix());
    }
    if ( Environment::vm().count(prefixvm(this->prefix(),"threshold.upper").c_str()) )
    {
        M_hasThresholdUpper = true;
        M_thresholdUpper = doption(_name="threshold.upper",_prefix=this->prefix());
    }


    if ( !M_inputPath.empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}
SurfaceFromImage::SurfaceFromImage( SurfaceFromImage const& e )
    :
    M_prefix( e.M_prefix ),
    M_inputPath( e.M_inputPath ),
    M_outputDirectory( e.M_outputDirectory ), M_outputPath( e.M_outputPath ),
    M_thresholdLower( e.M_thresholdLower ), M_thresholdUpper( e.M_thresholdUpper ),
    M_hasThresholdLower( e.M_hasThresholdLower), M_hasThresholdUpper( e.M_hasThresholdUpper),
    M_forceRebuild( e.M_forceRebuild )
{}

void
SurfaceFromImage::updateOutputPathFromInputFileName()
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

    std::string newFileName = (boost::format("%1%.stl")%nameMeshFile ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath = outputPath.string();
}

void
SurfaceFromImage::run()
{
    if ( !fs::exists( this->inputPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : surface segmentation not done because this input path not exist :" << this->inputPath() << "\n";
        return;
    }

    if ( !M_hasThresholdLower && !M_hasThresholdUpper )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : surface segmentation not done because no threshold has been given\n";
        return;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run ImageFromCenterlines \n"
            << "---------------------------------------\n";
    coutStr << "inputPath         : " << this->inputPath() << "\n";
    if ( M_hasThresholdLower )
        coutStr << "threshold lower   : " << M_thresholdLower << "\n";
    if ( M_hasThresholdUpper )
        coutStr << "threshold upper   : " << M_thresholdUpper << "\n";
    coutStr << "output path       : " << this->outputPath() << "\n"
            << "---------------------------------------\n"
            << "---------------------------------------\n";
    std::cout << coutStr.str();


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

    if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )
    {
        std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
        std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );

        std::string nameImageInit = fs::path(this->outputPath()).stem().string()+"_levelsetInit.vti";
        std::string outputPathImageInit = (fs::path(this->outputPath()).parent_path()/fs::path(nameImageInit)).string();

        std::ostringstream __str;
        __str << pythonExecutable << " ";
        __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtkimageinitialization ";
        __str << "-ifile " << this->inputPath() << " -interactive 0 -method threshold ";
        if ( M_hasThresholdLower )
            __str << " -lowerthreshold " << M_thresholdLower << " ";
        if ( M_hasThresholdUpper )
            __str << " -upperthreshold " << M_thresholdUpper << " ";

        __str << "-olevelsetsfile " << outputPathImageInit;

        std::cout << "---------------------------------------\n"
                  << "run in system : \n" << __str.str() << "\n"
                  << "---------------------------------------\n";
        auto err = ::system( __str.str().c_str() );

        std::ostringstream __str2;
        __str2 << pythonExecutable << " ";
        __str2 << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtklevelsetsegmentation ";
        __str2 << "-ifile " << this->inputPath() << " -initiallevelsetsfile " << outputPathImageInit << " -iterations 1 ";
        //__str2 << "--pipe vmtkmarchingcubes -i @.o -ofile " << this->outputPath();
        __str2 << "--pipe vmtkmarchingcubes -i @.o ";
        if ( true )
            __str2 << "--pipe " << dirBaseVmtk << "/vmtksurfaceconnectivity ";
        __str2 << "-ofile " << this->outputPath();

        std::cout << "---------------------------------------\n"
                  << "run in system : \n" << __str2.str() << "\n"
                  << "---------------------------------------\n";
        auto err2 = ::system( __str2.str().c_str() );
    }

    //vmtkimageinitialization -ifile hh32.mha -method threshold -lowerthreshold -1 -interactive 0 -osurfacefile imginit.stl -olevelsetsfile levelsetinit.vti
    //vmtklevelsetsegmentation -ifile hh32.mha -initiallevelsetsfile levelsetinit2.vti -iterations 1  --pipe vmtkmarchingcubes -i @.o -ofile bid.stl
}

po::options_description
SurfaceFromImage::options( std::string const& prefix )
{
    po::options_description mySurfaceFromImageOptions( "Create Surface from Image options" );

    mySurfaceFromImageOptions.add_options()
        (prefixvm(prefix,"input.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"threshold.lower").c_str(), Feel::po::value<double>(), "(double) threshold lower")
        (prefixvm(prefix,"threshold.upper").c_str(), Feel::po::value<double>(), "(double) threshold upper")
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        ;
    return mySurfaceFromImageOptions;
}



SubdivideSurface::SubdivideSurface( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputPath( soption(_name="input.filename",_prefix=this->prefix()) ),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
    M_method( soption(_name="method",_prefix=this->prefix()) ),
    M_nSubdivisions( ioption(_name="subdivisions",_prefix=this->prefix()) )
{
    CHECK( M_method == "linear" || M_method == "butterfly" || M_method == "loop" ) << "invalid method " << M_method << "\n";
    if ( !M_inputPath.empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}


SubdivideSurface::SubdivideSurface( SubdivideSurface const& e )
    :
    M_prefix( e.M_prefix ),
    M_inputPath( e.M_inputPath ),
    M_outputDirectory( e.M_outputDirectory ), M_outputPath( e.M_outputPath ),
    M_forceRebuild( e.M_forceRebuild ),
    M_method( e.M_method ),
    M_nSubdivisions( e.M_nSubdivisions )
{}

void
SubdivideSurface::updateOutputPathFromInputFileName()
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
    //fs::path gp = M_inputPath;
    std::string nameMeshFile = fs::path(this->inputPath()).stem().string();

    std::string newFileName = (boost::format("%1%_subdivide%2%%3%.stl")%nameMeshFile %M_nSubdivisions %M_method ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath = outputPath.string();
}

void
SubdivideSurface::run()
{
    if ( !fs::exists( this->inputPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : smoothing surface not done because this input surface path not exist :" << this->inputPath() << "\n";
        return;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run SubdivideSurface \n"
            << "---------------------------------------\n";
    coutStr << "inputPath     : " << this->inputPath() << "\n";
    coutStr << "method        : " << M_method << "\n"
            << "nSubdivisions : " << M_nSubdivisions << "\n";
    coutStr << "output path   : " << this->outputPath() << "\n"
            << "---------------------------------------\n"
            << "---------------------------------------\n";
    std::cout << coutStr.str();


    fs::path directory;
    // build directories if necessary
    if ( !this->outputPath().empty() && this->worldComm().isMasterRank() )
    {
        directory = fs::path(this->outputPath()).parent_path();
        if ( !fs::exists( directory ) )
            fs::create_directories( directory );
    }
    // wait all process
    this->worldComm().globalComm().barrier();

    if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )
    {
        std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
        std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );

        std::ostringstream __str;
        __str << pythonExecutable << " ";
        __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtksurfacesubdivision ";
        __str << "-ifile " << this->inputPath() << " "
              << "-method " << M_method << " "
              << "-subdivisions " << M_nSubdivisions << " ";
        __str << "-ofile " << this->outputPath();
        std::cout << "---------------------------------------\n"
                  << "run in system : \n" << __str.str() << "\n"
                  << "---------------------------------------\n";
        auto err = ::system( __str.str().c_str() );
    }

}

po::options_description
SubdivideSurface::options( std::string const& prefix )
{
    po::options_description mySubdivideSurfaceOptions( "Subdivide the surface options" );

    mySubdivideSurfaceOptions.add_options()
        (prefixvm(prefix,"input.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        (prefixvm(prefix,"method").c_str(), Feel::po::value<std::string>()->default_value("butterfly"), "(string) linear, butterfly, loop")
        (prefixvm(prefix,"subdivisions").c_str(), Feel::po::value<int>()->default_value(1), "(int) number of subdivisions")
        ;
    return mySubdivideSurfaceOptions;

}

SmoothSurface::SmoothSurface( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputPath( soption(_name="input.filename",_prefix=this->prefix()) ),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
    M_method( soption(_name="method",_prefix=this->prefix()) ),
    M_nIterations( ioption(_name="iterations",_prefix=this->prefix()) ),
    M_taubinPassBand( doption(_name="taubin.passband",_prefix=this->prefix()) ),
    M_laplaceRelaxationFactor( doption(_name="laplace.relaxation",_prefix=this->prefix()) )
{
    CHECK( M_method == "taubin" || M_method == "laplace" ) << "invalid method " << M_method << "\n";
    if ( !M_inputPath.empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}

SmoothSurface::SmoothSurface( SmoothSurface const& e )
    :
    M_prefix( e.M_prefix ),
    M_inputPath( e.M_inputPath ),
    M_outputDirectory( e.M_outputDirectory ), M_outputPath( e.M_outputPath ),
    M_forceRebuild( e.M_forceRebuild ),
    M_method( e.M_method ),
    M_nIterations( e.M_nIterations ),
    M_taubinPassBand( e.M_taubinPassBand ),
    M_laplaceRelaxationFactor( e.M_laplaceRelaxationFactor )
{}

void
SmoothSurface::updateOutputPathFromInputFileName()
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

    std::string newFileName;
    if ( M_method == "taubin")
        newFileName = (boost::format("%1%_smooth%2%Taubin%3%.stl")%nameMeshFile %M_nIterations %M_taubinPassBand ).str();
    else
        newFileName = (boost::format("%1%_smooth%2%Laplace%3%.stl")%nameMeshFile %M_nIterations %M_laplaceRelaxationFactor ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath = outputPath.string();
}

void
SmoothSurface::run()
{
    if ( !fs::exists( this->inputPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : smoothing surface not done because this input surface path not exist :" << this->inputPath() << "\n";
        return;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run SmoothSurface \n"
            << "---------------------------------------\n";
    coutStr << "inputPath     : " << this->inputPath() << "\n";
    coutStr << "method        : " << M_method << "\n"
            << "nIterations   : " << M_nIterations << "\n";
    if ( M_method == "taubin" )
        coutStr << "passband      : " << M_taubinPassBand << "\n";
    if ( M_method == "laplace" )
        coutStr << "relaxation factor : " << M_laplaceRelaxationFactor << "\n";
    coutStr << "output path          : " << this->outputPath() << "\n"
            << "---------------------------------------\n"
            << "---------------------------------------\n";
    std::cout << coutStr.str();


    fs::path directory;
    // build directories if necessary
    if ( !this->outputPath().empty() && this->worldComm().isMasterRank() )
    {
        directory = fs::path(this->outputPath()).parent_path();
        if ( !fs::exists( directory ) )
            fs::create_directories( directory );
    }
    // wait all process
    this->worldComm().globalComm().barrier();

    if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )
    {
        std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
        std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );

        std::ostringstream __str;
        __str << pythonExecutable << " ";
        __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtksurfacesmoothing ";
        __str << "-ifile " << this->inputPath() << " "
              << "-iterations " << M_nIterations << " ";
        if ( M_method == "taubin" )
            __str << "-method taubin -passband " << M_taubinPassBand << " ";
        if ( M_method == "laplace" )
            __str << "-method laplace -relaxation " << M_laplaceRelaxationFactor << " ";
        __str << "-ofile " << this->outputPath();
        std::cout << "---------------------------------------\n"
                  << "run in system : \n" << __str.str() << "\n"
                  << "---------------------------------------\n";
        auto err = ::system( __str.str().c_str() );
    }

}

po::options_description
SmoothSurface::options( std::string const& prefix )
{
    po::options_description mySmoothSurfaceOptions( "Smooth the surface options" );

    mySmoothSurfaceOptions.add_options()
        (prefixvm(prefix,"input.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        (prefixvm(prefix,"method").c_str(), Feel::po::value<std::string>()->default_value("taubin"), "(string) taubin or laplace")
        (prefixvm(prefix,"iterations").c_str(), Feel::po::value<int>()->default_value(30), "(int) number of iterations")
        (prefixvm(prefix,"taubin.passband").c_str(), Feel::po::value<double>()->default_value(0.1), "(double) taubin passband")
        (prefixvm(prefix,"laplace.relaxation").c_str(), Feel::po::value<double>()->default_value(0.01), "(double) laplace.relaxation")
        ;
    return mySmoothSurfaceOptions;

}



OpenSurface::OpenSurface( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputSurfacePath( soption(_name="input.surface.filename",_prefix=this->prefix()) ),
    M_inputCenterlinesPath( soption(_name="input.centerlines.filename",_prefix=this->prefix()) ),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) )
{
    if ( !M_inputSurfacePath.empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}

OpenSurface::OpenSurface( OpenSurface const& e )
    :
    M_prefix( e.M_prefix ),
    M_inputSurfacePath( e.M_inputSurfacePath ),
    M_inputCenterlinesPath( e.M_inputCenterlinesPath ),
    M_outputDirectory( e.M_outputDirectory ), M_outputPath( e.M_outputPath ),
    M_forceRebuild( e.M_forceRebuild )
{}

void
OpenSurface::updateOutputPathFromInputFileName()
{
    CHECK( !M_inputSurfacePath.empty() ) << "input path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( M_outputDirectory.empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(M_outputDirectory).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(M_outputDirectory);
    else
        meshesdirectories = fs::path(M_outputDirectory);

    // get filename without extension
    fs::path gp = M_inputSurfacePath;
    std::string nameMeshFile = gp.stem().string();

    std::string newFileName = (boost::format("%1%_open.stl")%nameMeshFile ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath = outputPath.string();
}

void
OpenSurface::run()
{
    if ( !fs::exists( this->inputSurfacePath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : opening surface not done because this input surface path for centerlines not exist :" << this->inputSurfacePath() << "\n";
        return;
    }
    if ( !fs::exists( this->inputCenterlinesPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : opening surface not done because this input centerlines path for centerlines not exist :" << this->inputCenterlinesPath() << "\n";
        return;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run OpenSurface \n"
            << "---------------------------------------\n";
    coutStr << "inputSurfacePath     : " << this->inputSurfacePath() << "\n";
    coutStr << "inputCenterlinesPath : " << this->inputCenterlinesPath() << "\n";
    coutStr << "output path          : " << this->outputPath() << "\n"
            << "---------------------------------------\n"
            << "---------------------------------------\n";
    std::cout << coutStr.str();


    fs::path directory;
    // build directories if necessary
    if ( !this->outputPath().empty() && this->worldComm().isMasterRank() )
    {
        directory = fs::path(this->outputPath()).parent_path();
        if ( !fs::exists( directory ) )
            fs::create_directories( directory );
    }
    // wait all process
    this->worldComm().globalComm().barrier();

    if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )
    {
        if ( true )
        {
            this->runGMSH();
        }
        else
        {
            std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
            std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );

            std::ostringstream __str;
            __str << pythonExecutable << " ";
            __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtkendpointextractor ";
            __str << "-ifile " << this->inputCenterlinesPath() << " "
                  << "-radiusarray MaximumInscribedSphereRadius -numberofendpointspheres 1 ";
            __str << "--pipe ";
            __str << dirBaseVmtk << "/vmtkbranchclipper "
                  << "-ifile " << this->inputSurfacePath() << " "
                  << "-groupidsarray TractIds -blankingarray Blanking -radiusarray MaximumInscribedSphereRadius -insideout 0 -interactive 0 ";
            __str << "--pipe ";
            __str << dirBaseVmtk << "/vmtksurfaceconnectivity "
                  << " -cleanoutput 1 ";
            __str << "--pipe "
                  <<  dirBaseVmtk << "/vmtksurfacewriter "
                  << "-ofile " << this->outputPath();
            std::cout << "---------------------------------------\n"
                      << "run in system : \n" << __str.str() << "\n"
                      << "---------------------------------------\n";
            auto err = ::system( __str.str().c_str() );
        }
    }

    // vmtkendpointextractor -ifile ~/Desktop/MaillagesOdysée/meshing2/centerlines/cut6_centerlines.vtk -radiusarray MaximumInscribedSphereRadius -numberofendpointspheres 1 --pipe vmtkbranchclipper -interactive 0 -ifile ~/feel/blabll/hh128.stl -groupidsarray TractIds -blankingarray Blanking -radiusarray MaximumInscribedSphereRadius   -insideout 0 --pipe vmtksurfaceconnectivity -cleanoutput 1 --pipe vmtksurfacewriter -ofile hola3.stl
}

void
OpenSurface::runGMSH()
{
    std::ostringstream geodesc;
    geodesc << "Mesh.Algorithm = 6; //(1=MeshAdapt, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) \n"
            << "Mesh.Algorithm3D = 1; //(1=tetgen, 4=netgen, 7=MMG3D, 9=R-tree) \n";

    geodesc << "Merge \""<< this->inputSurfacePath() <<"\";\n";
    geodesc << "Field[1] = AngioTkCenterline;\n";
    geodesc << "Field[1].FileName = \"" << this->inputCenterlinesPath() << "\";\n";
    geodesc << "Field[1].clipMesh =1;\n"
            << "Field[1].run;\n"
            << "Background Field = 1;\n";

    fs::path outputMeshNamePath = fs::path(this->outputPath());
    std::string _name = outputMeshNamePath.stem().string();
    std::string geoname=_name+".geo";

    std::ofstream geofile( geoname.c_str() );
    geofile << geodesc.str();
    geofile.close();

    detail::generateMeshFromGeo(geoname,this->outputPath(),2);
}

po::options_description
OpenSurface::options( std::string const& prefix )
{
    po::options_description myOpenSurfaceOptions( "Open the surface options" );

    myOpenSurfaceOptions.add_options()
        (prefixvm(prefix,"input.centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"input.surface.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        ;
    return myOpenSurfaceOptions;
}




RemeshSTL::RemeshSTL( std::string prefix )
    :
    M_prefix( prefix ),
    M_packageType(soption(_name="package-type",_prefix=this->prefix())),
    M_inputPath(soption(_name="input.filename",_prefix=this->prefix())),
    M_centerlinesFileName(soption(_name="centerlines.filename",_prefix=this->prefix())),
    M_remeshNbPointsInCircle( ioption(_name="nb-points-in-circle",_prefix=this->prefix()) ),
    M_area( doption(_name="area",_prefix=this->prefix()) ),
    M_nIterationVMTK( ioption(_name="vmtk.n-iteration",_prefix=this->prefix()) ),
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
    M_area( e.M_area ),M_nIterationVMTK( e.M_nIterationVMTK ),
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
    __str << "-iterations " << M_nIterationVMTK << " ";
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

    //geodesc << "Field[1] = Centerline;\n";
    geodesc << "Field[1] = AngioTkCenterline;\n";


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
        ( prefixvm(prefix,"vmtk.n-iteration").c_str(), po::value<int>()->default_value( 10 ), "maxit" )

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

#if 0
    char * cmd2[3];// = new char("-v 2");
    cmd2[0] = new char[4];cmd2[0][0]='h';cmd2[0][1]='o';cmd2[0][2]='l';cmd2[0][3]='a';
    cmd2[1] = new char[2];cmd2[1][0]='-';cmd2[1][1]='v';
    cmd2[2] = new char[1];cmd2[2][0]='5';
    for (int i = 0; i < 3; ++i) std::cout << cmd2[i] << std::endl;
    GmshInitialize(3,cmd2 );
#else
    GmshInitialize();
#endif

    // if(!Msg::GetGmshClient())
    CTX::instance()->terminal = 1;
    //GmshBatch();

    int verbosityLevel = 5;
    Msg::SetVerbosity( verbosityLevel );

    GModel * M_gmodel = new GModel();

    M_gmodel->setName( _name );
#if !defined( __APPLE__ )
    M_gmodel->setFileName( _name );
#endif

    GModel::current()->getFields()->map_type_name["AngioTkCenterline"] = new FieldFactoryT<AngioTkCenterline>();

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

    //geodesc << "Field[1] = Centerline;\n";
    geodesc << "Field[1] = AngioTkCenterline;\n";

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
