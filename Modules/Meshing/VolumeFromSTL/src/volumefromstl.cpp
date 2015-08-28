/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <volumefromstl.hpp>
#include <AngioTkCenterlineField.h>
#include <centerlinesmanagerwindowinteractor.hpp>

#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkMetaImageReader.h> // mha
#include <vtkXMLImageDataReader.h> // vti
#include <vtkMarchingCubes.h>
#include <vtkSTLWriter.h>
#include <vtkImageTranslateExtent.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPointData.h>
#include <vtkCleanPolyData.h>

//#include <vtkvmtkPolyBallModeller.h>
#include <angiotkPolyBallModeller.h>

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
    fs::path dir = fs::path(outputPath).parent_path();
    if ( !fs::exists( dir ) )
        fs::create_directories( dir );

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
    M_inputPath( AngioTkEnvironment::expand( soption(_name="input.filename",_prefix=this->prefix()) ) ),
    M_inputCenterlinesPointSetPath( AngioTkEnvironment::expand( soption(_name="input.pointset.filename",_prefix=this->prefix()) ) ),
    M_inputInletOutletDescPath( AngioTkEnvironment::expand( soption(_name="input.desc.filename",_prefix=this->prefix()) ) ),
    M_inputGeoCenterlinesPath( AngioTkEnvironment::expand( soption(_name="input.geo-centerlines.filename",_prefix=this->prefix()) ) ),
    M_outputDirectory( AngioTkEnvironment::expand( soption(_name="output.directory",_prefix=this->prefix()) ) ),
    M_costFunctionExpr( soption(_name="cost-function.expression",_prefix=this->prefix()) ),
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
    if ( this->inputCenterlinesPointSetPath().empty() && this->sourceids().empty() && this->targetids().empty() && this->inputGeoCenterlinesPath().empty() )
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
    coutStr << "input surface path   : " << this->inputPath() << "\n";
    if ( !M_useInteractiveSelection )
    {
        if ( !this->inputGeoCenterlinesPath().empty() )
            coutStr << "input GeoCenterlines : " << this->inputGeoCenterlinesPath() << "\n";
        else if ( !this->inputCenterlinesPointSetPath().empty() )
            coutStr << "input GeoCenterlines : " << this->inputCenterlinesPointSetPath() << "\n";
        else
        {
            coutStr << "targetids            : ";
            for ( int id : this->targetids() )
                coutStr << id << " ";
            coutStr << "\n";
            coutStr << "sourceids            : ";
            for ( int id : this->sourceids() )
                coutStr << id << " ";
            coutStr << "\n";
        }
    }

    coutStr << "output path          : " << this->outputPath() << "\n"
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
        if ( !this->inputGeoCenterlinesPath().empty() && fs::exists( this->inputGeoCenterlinesPath() ) )
        {
            GmshInitialize();
            CTX::instance()->terminal = 1;
            int verbosityLevel = 5;
            Msg::SetVerbosity( verbosityLevel );

            AngioTkCenterline centerlinesTool;
            centerlinesTool.createFromGeoCenterlinesFile( this->inputGeoCenterlinesPath(),this->inputPath() );
            centerlinesTool.addFieldBranchIds();
            centerlinesTool.addFieldRadiusMin("RadiusMin");
            centerlinesTool.addFieldRadiusMin("MaximumInscribedSphereRadius");
            centerlinesTool.writeCenterlinesVTK( this->outputPath() );
        }
        else
        {
            std::ostringstream __str;
            __str << pythonExecutable << " ";
            __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtkcenterlines ";
            //__str << "-usetetgen 1 -simplifyvoronoi 1 ";

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

            __str << " -costfunction " << M_costFunctionExpr << " ";
#if 0
            //std::string costFunction = "1/R";
            //std::string costFunction = "(R-0.5)*(R-0.5)";
            //std::string costFunction = "R*R-R+0.25";
            //std::string costFunction = "R*R-2*0.8*R+0.8*0.8";
            std::string costFunction = "R*R-2*1.1*R+1.1*1.1";
            __str << " -costfunction " << costFunction << " ";
            //__str << " -delaunaytessellationfile " << costFunction << " ";
#endif

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
        (prefixvm(prefix,"input.geo-centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "geo-centerlines.filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"use-interactive-selection").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) use-interactive-selection")
        (prefixvm(prefix,"cost-function.expression").c_str(), Feel::po::value<std::string>()->default_value("1/R"), "(string) cost-function")
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
    M_inputCenterlinesPath(),// 1, soption(_name="input.centerlines.filename",_prefix=this->prefix()) ),
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.surface.filename",_prefix=this->prefix()) ) ),
    M_inputPointSetPath( AngioTkEnvironment::expand( soption(_name="input.point-set.filename",_prefix=this->prefix()) ) ),
    M_outputDirectory( AngioTkEnvironment::expand( soption(_name="output.directory",_prefix=this->prefix()) ) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
    M_useWindowInteractor( boption(_name="use-window-interactor",_prefix=this->prefix() ) ),
    M_applyThresholdMinRadius( doption(_name="threshold-radius.min",_prefix=this->prefix() ) ),
    M_applyThresholdMaxRadius( doption(_name="threshold-radius.max",_prefix=this->prefix() ) ),
    M_applyThresholdZonePointSetPath( AngioTkEnvironment::expand( soption(_name="thresholdzone-radius.point-set.filename",_prefix=this->prefix()) ) ),
    M_applyThresholdZoneMinRadius( doption(_name="thresholdzone-radius.min",_prefix=this->prefix() ) ),
    M_applyThresholdZoneMaxRadius( doption(_name="thresholdzone-radius.max",_prefix=this->prefix() ) )
{
    if ( Environment::vm().count(prefixvm(this->prefix(),"input.centerlines.filename").c_str()) )
        M_inputCenterlinesPath = Environment::vm()[prefixvm(this->prefix(),"input.centerlines.filename").c_str()].as<std::vector<std::string> >();

    std::vector<int> removeids;
    if ( Environment::vm().count(prefixvm(this->prefix(),"remove-branch-ids").c_str()) )
        removeids = Environment::vm()[prefixvm(this->prefix(),"remove-branch-ids").c_str()].as<std::vector<int> >();
    for ( int id : removeids )
        M_removeBranchIds.insert( id );

    for (int k = 0;k<this->inputCenterlinesPath().size();++k)
    {
        M_inputCenterlinesPath[k] = AngioTkEnvironment::expand( M_inputCenterlinesPath[k] );
        if ( !this->inputCenterlinesPath(k).empty() && fs::path(this->inputCenterlinesPath(k)).is_relative() )
            M_inputCenterlinesPath[k] = (AngioTkEnvironment::pathInitial()/fs::path(this->inputCenterlinesPath(k)) ).string();
    }

    if ( !this->inputCenterlinesPath().empty() && !this->inputCenterlinesPath(0).empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}

void
CenterlinesManager::updateOutputPathFromInputFileName()
{
    CHECK( !M_inputCenterlinesPath.empty() ) << "input path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( M_outputDirectory.empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(M_outputDirectory).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(M_outputDirectory);
    else
        meshesdirectories = fs::path(M_outputDirectory);

    // get filename without extension
    fs::path gp = this->inputCenterlinesPath(0);
    std::string nameMeshFile = gp.stem().string();

    std::string newFileName = (boost::format("%1%_up.vtk")%nameMeshFile ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath = outputPath.string();
}

std::map<int,std::vector<std::tuple<double,double,double> > >
CenterlinesManager::loadPointSetFile( std::string const& filepath )
{
    CHECK( fs::exists(filepath) ) << "inputTubularColisionPointSetPath not exists : " << filepath;
    std::ifstream fileLoaded( filepath.c_str(), std::ios::in);
    std::map<int,std::vector<std::tuple<double,double,double> > > pointPair;
    while ( !fileLoaded.eof() )
      {
          std::vector<double> pt(3);
          double radius;
          int typePt = -1;
          fileLoaded >> typePt;
          if ( fileLoaded.eof() ) break;
          fileLoaded >> pt[0] >> pt[1] >> pt[2] >> radius;
          pointPair[typePt].push_back( std::make_tuple(pt[0],pt[1],pt[2]) );
      }
    fileLoaded.close();
    return pointPair;
}

void
CenterlinesManager::run()
{
    if ( M_useWindowInteractor )
    {
        if ( !this->inputSurfacePath().empty() && fs::exists( this->inputSurfacePath() ) )
        {
            CenterlinesManagerWindowInteractor windowInteractor;
            windowInteractor.setInputSurfacePath( this->inputSurfacePath() );
            if ( !this->inputPointSetPath().empty() && fs::exists( this->inputPointSetPath() ) )
                windowInteractor.setInputPointSetPath( this->inputPointSetPath() );

            std::vector<std::string> centerlinesPath;
            for (int k = 0;k<this->inputCenterlinesPath().size();++k)
                if ( !this->inputCenterlinesPath(k).empty() && fs::exists( this->inputCenterlinesPath(k) ) )
                    centerlinesPath.push_back( this->inputCenterlinesPath(k) );
            windowInteractor.setInputCenterlinesPath( centerlinesPath );

            windowInteractor.run();
        }
        else if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : Centerlines Manager not run because this input surface path not exist :" << this->inputSurfacePath() << "\n";

        return;
    }

    if ( this->inputCenterlinesPath().empty() )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : Centerlines Manager not run because this input centerlines is empty\n";
        return;
    }
    if ( !fs::exists( this->inputCenterlinesPath(0) ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : Centerlines Manager not run because this input centerlines path not exist :" << this->inputCenterlinesPath(0) << "\n";
        return;
    }
    if ( this->inputCenterlinesPath().size() > 1 && ( this->inputSurfacePath().empty() || !fs::exists( this->inputSurfacePath() ) ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : Centerlines Manager (fusion case) not run because this input centerlines path not exist :" << this->inputCenterlinesPath(0) << "\n";
        return;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run CenterlinesManager \n"
            << "---------------------------------------\n";
    for ( int k=0;k<this->inputCenterlinesPath().size();++k)
        coutStr << "inputCenterlinesPath[" << k << "]  : " << this->inputCenterlinesPath(k) << "\n";
    if ( M_removeBranchIds.size() > 0 )
    {
        coutStr << "remove branch ids :";
        for ( int id : M_removeBranchIds )
            coutStr << " " << id;
        coutStr << "\n";
    }
    if ( M_applyThresholdMinRadius > 0 )
        coutStr << "threshold min-radius     : " << M_applyThresholdMinRadius << "\n";
    if ( M_applyThresholdMaxRadius > 0 )
        coutStr << "threshold max-radius     : " << M_applyThresholdMaxRadius << "\n";

    bool hasThresholdZonePointSetPath = !M_applyThresholdZonePointSetPath.empty() && fs::exists( M_applyThresholdZonePointSetPath );
    if ( hasThresholdZonePointSetPath )
    {
        coutStr << "thresholdzonezone point set path : " << M_applyThresholdZonePointSetPath << "\n";
        if ( M_applyThresholdZoneMinRadius > 0 )
            coutStr << "thresholdzone min-radius         : " << M_applyThresholdZoneMinRadius << "\n";
        if ( M_applyThresholdZoneMaxRadius > 0 )
            coutStr << "thresholdzone max-radius         : " << M_applyThresholdZoneMaxRadius << "\n";
    }

    coutStr << "output path              : " << this->outputPath() << "\n"
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


    bool surfaceMeshExist = !this->inputSurfacePath().empty() && fs::exists( this->inputSurfacePath() );

    if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )
    {
        std::shared_ptr<AngioTkCenterline> centerlinesTool;
        for ( int k=0;k<this->inputCenterlinesPath().size();++k)
        {
            if ( !centerlinesTool )
            {
                if ( true )
                {
                    CTX::instance()->terminal = 1;
                    int verbosityLevel = 5;
                    Msg::SetVerbosity( verbosityLevel );
                }
                centerlinesTool.reset( new AngioTkCenterline );
                if ( surfaceMeshExist )
                    centerlinesTool->importSurfaceFromFile( this->inputSurfacePath() );
            }
            centerlinesTool->importFile( this->inputCenterlinesPath(k) );
        }
        centerlinesTool->removeBranchIds( M_removeBranchIds );

        centerlinesTool->addFieldBranchIds();
        if ( surfaceMeshExist )
            centerlinesTool->addFieldRadiusMin();

        if ( M_applyThresholdMinRadius > 0 || M_applyThresholdMaxRadius > 0 )
        {
            std::vector<std::string> fieldnames = { "RadiusMin","MaximumInscribedSphereRadius" };
            if ( M_applyThresholdMinRadius > 0 )
                centerlinesTool->applyFieldThresholdMin( fieldnames,M_applyThresholdMinRadius );
            if ( M_applyThresholdMaxRadius > 0 )
                centerlinesTool->applyFieldThresholdMax( fieldnames,M_applyThresholdMaxRadius );
        }
        if ( hasThresholdZonePointSetPath )
        {
            auto pointSetData = this->loadPointSetFile( M_applyThresholdZonePointSetPath );
            std::vector<std::string> fieldnames = { "RadiusMin","MaximumInscribedSphereRadius" };
            if ( M_applyThresholdZoneMinRadius > 0 )
                centerlinesTool->applyFieldThresholdZoneMin( fieldnames,M_applyThresholdZoneMinRadius,pointSetData );
            if ( M_applyThresholdZoneMaxRadius > 0 )
                centerlinesTool->applyFieldThresholdZoneMax( fieldnames,M_applyThresholdZoneMaxRadius,pointSetData );
        }

        centerlinesTool->writeCenterlinesVTK( this->outputPath() );
    }

}

po::options_description
CenterlinesManager::options( std::string const& prefix )
{
    po::options_description myCenterlinesManagerOptions( "Centerlines Manager from Image options" );

    myCenterlinesManagerOptions.add_options()
        //(prefixvm(prefix,"input.centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"input.centerlines.filename").c_str(), po::value<std::vector<std::string>>()->multitoken(), "(vector<string>) input centerline filename" )
        (prefixvm(prefix,"input.surface.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input surface filename" )
        (prefixvm(prefix,"input.point-set.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input point-set filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"remove-branch-ids").c_str(), po::value<std::vector<int> >()->multitoken(), "(vector of int) remove branch ids" )
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        (prefixvm(prefix,"use-window-interactor").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) use-window-interactor")
        (prefixvm(prefix,"threshold-radius.min").c_str(), Feel::po::value<double>()->default_value(-1), "(double) threshold-radius.min")
        (prefixvm(prefix,"threshold-radius.max").c_str(), Feel::po::value<double>()->default_value(-1), "(double) threshold-radius.max")
        (prefixvm(prefix,"thresholdzone-radius.point-set.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input point-set filename" )
        (prefixvm(prefix,"thresholdzone-radius.max").c_str(), Feel::po::value<double>()->default_value(-1), "(double) threshold-radius.max")
        (prefixvm(prefix,"thresholdzone-radius.min").c_str(), Feel::po::value<double>()->default_value(-1), "(double) threshold-radius.max")
        ;
    return myCenterlinesManagerOptions;
}



ImageFromCenterlines::ImageFromCenterlines( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputCenterlinesPath( AngioTkEnvironment::expand( soption(_name="input.filename",_prefix=this->prefix()) ) ),
    M_outputDirectory( AngioTkEnvironment::expand( soption(_name="output.directory",_prefix=this->prefix()) ) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
    M_dimX( ioption(_name="dim.x",_prefix=this->prefix() ) ),
    M_dimY( ioption(_name="dim.y",_prefix=this->prefix() ) ),
    M_dimZ( ioption(_name="dim.z",_prefix=this->prefix() ) ),
    M_dimSpacing( doption(_name="dim.spacing",_prefix=this->prefix() ) ),
    M_radiusArrayName( soption(_name="radius-array-name",_prefix=this->prefix() ) )
{
    if ( !this->inputCenterlinesPath().empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}

void
ImageFromCenterlines::updateOutputPathFromInputFileName()
{
    CHECK( !this->inputCenterlinesPath().empty() ) << "input centerlines path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( M_outputDirectory.empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(M_outputDirectory).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(M_outputDirectory);
    else
        meshesdirectories = fs::path(M_outputDirectory);

    // get filename without extension
    fs::path gp = this->inputCenterlinesPath();
    std::string nameMeshFile = gp.stem().string();

    std::string dimXTag = (M_dimX > 0 )? (boost::format("%1%")%M_dimX).str() : "x";
    std::string dimYTag = (M_dimY > 0 )? (boost::format("%1%")%M_dimY).str() : "y";
    std::string dimZTag = (M_dimZ > 0 )? (boost::format("%1%")%M_dimZ).str() : "z";

    std::string dimComponentTag;
    if ( M_dimX > 0 || M_dimY > 0 || M_dimZ > 0 )
        dimComponentTag = "_" + dimXTag + "_" + dimYTag  + "_" + dimZTag;

    std::string dimSpacingTag;
    if ( std::abs(M_dimSpacing) > 1e-12 )
        dimSpacingTag = (boost::format("_spacing%1%")%M_dimSpacing).str();

    //std::string newFileName = (boost::format("%1%_%2%-%3%-%4%.mha")%nameMeshFile %M_dimX %M_dimY %M_dimZ ).str();
    std::string newFileName = (boost::format("%1%%2%%3%.mha")%nameMeshFile %dimComponentTag %dimSpacingTag ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath = outputPath.string();
}

void
ImageFromCenterlines::run()
{
    if ( !fs::exists( this->inputCenterlinesPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : image building not done because this input path not exist :" << this->inputCenterlinesPath() << "\n";
        return;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run ImageFromCenterlines \n"
            << "---------------------------------------\n";
    coutStr << "input centerlines   : " << this->inputCenterlinesPath() << "\n";
    if ( std::abs(M_dimSpacing) > 1e-12 )
        coutStr << "dimensions spacing  : " << M_dimSpacing << "\n";
    if ( M_dimX > 0 || M_dimY > 0 || M_dimZ > 0 )
        coutStr << "dimensions          : [" << M_dimX << "," << M_dimY << "," << M_dimZ << "]\n";
    coutStr << "radius array name   : " << M_radiusArrayName << "\n";
    coutStr << "output path         : " << this->outputPath() << "\n"
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
#if 0
        std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
        std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );

        std::ostringstream __str;
        __str << pythonExecutable << " ";
        __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtkcenterlinemodeller ";
        //__str << "-ifile " << this->inputPath() << " -radiusarray MaximumInscribedSphereRadius "
        __str << "-ifile " << this->inputCenterlinesPath() << " -radiusarray " << M_radiusArrayName << " "
              << "-negate 1 -dimensions " << M_dimX << " " << M_dimY << " " << M_dimZ << " ";
        __str << "-ofile " << this->outputPath();

        std::cout << "---------------------------------------\n"
                  << "run in system : \n" << __str.str() << "\n"
                  << "---------------------------------------\n";
        auto err = ::system( __str.str().c_str() );
#else

        vtkSmartPointer<vtkPolyDataReader> readerVTK = vtkSmartPointer<vtkPolyDataReader>::New();
        readerVTK->SetFileName(this->inputCenterlinesPath().c_str());
        readerVTK->Update();

        //vtkSmartPointer<vtkvmtkPolyBallModeller> modeller = vtkvmtkPolyBallModeller::New();
        vtkSmartPointer<angiotkPolyBallModeller> modeller = angiotkPolyBallModeller::New();
        modeller->SetInput( (vtkDataObject*)readerVTK->GetOutput() );
        modeller->SetRadiusArrayName(M_radiusArrayName.c_str());
        modeller->UsePolyBallLineOn();

        if ( std::abs(M_dimSpacing) > 1e-12 )
        {
            vtkPolyData* inputCenterlines = vtkPolyData::SafeDownCast(readerVTK->GetOutput());
            vtkDataArray* radiusArray = inputCenterlines->GetPointData()->GetArray(M_radiusArrayName.c_str());
            double maxRadius = radiusArray->GetRange()[1];
            double* bounds = inputCenterlines->GetBounds();
            double maxDist = 2.0 * maxRadius;
            if ( M_dimX == 0 )
                M_dimX = static_cast<int>( std::floor( (bounds[1] - bounds[0] + 2*maxDist )/M_dimSpacing) ) + 1;
            if ( M_dimY == 0 )
                M_dimY = static_cast<int>( std::floor( (bounds[3] - bounds[2] + 2*maxDist )/M_dimSpacing) ) + 1;
            if ( M_dimZ == 0 )
                M_dimZ = static_cast<int>( std::floor( (bounds[5] - bounds[4] + 2*maxDist )/M_dimSpacing) ) + 1;
        }
        CHECK( M_dimX > 0 && M_dimY > 0 && M_dimZ > 2 ) << "M_dimX,M_dimY,M_dimZ must be > 0";
        int sampleDimensions[3] = { M_dimX,M_dimY,M_dimZ };
        modeller->SetSampleDimensions( sampleDimensions );
        modeller->SetNegateFunction(1);
        modeller->Update();

        vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
        writer->SetInput(modeller->GetOutput());
        writer->SetFileName(this->outputPath().c_str());
        writer->Write();
#endif
    }

}

po::options_description
ImageFromCenterlines::options( std::string const& prefix )
{
    po::options_description myImageFromCenterlinesOptions( "Create Image from Centerlines options" );

    myImageFromCenterlinesOptions.add_options()
        (prefixvm(prefix,"input.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        (prefixvm(prefix,"dim.x").c_str(), Feel::po::value<int>()->default_value(0), "(int) dim.x")
        (prefixvm(prefix,"dim.y").c_str(), Feel::po::value<int>()->default_value(0), "(int) dim.y")
        (prefixvm(prefix,"dim.z").c_str(), Feel::po::value<int>()->default_value(0), "(int) dim.z")
        (prefixvm(prefix,"dim.spacing").c_str(), Feel::po::value<double>()->default_value(0.), "(double) dim.spacing")
        (prefixvm(prefix,"radius-array-name").c_str(), Feel::po::value<std::string>()->default_value("MaximumInscribedSphereRadius"), "(std::string) radius-array-name")
        ;
    return myImageFromCenterlinesOptions;
}

SurfaceFromImage::SurfaceFromImage( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputPath( AngioTkEnvironment::expand( soption(_name="input.filename",_prefix=this->prefix()) ) ),
    M_outputDirectory( AngioTkEnvironment::expand( soption(_name="output.directory",_prefix=this->prefix()) ) ),
    M_method( soption(_name="method",_prefix=this->prefix()) ),
    M_hasThresholdLower(false), M_hasThresholdUpper(false),
    M_thresholdLower(0.0),M_thresholdUpper(0.0),
    M_applyConnectivityLargestRegion( boption(_name="apply-connectivity.largest-region",_prefix=this->prefix()) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) )
{

    CHECK( M_method == "threshold" || M_method == "isosurface" ) << "invalid method : " << M_method << " -> must be threshold or isosurface";

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
            << "run SurfaceFromImage \n"
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
        __str << "-ifile " << this->inputPath() << " -interactive 0 ";
        if ( M_method == "threshold" )
        {
            __str << "-method threshold ";
            if ( M_hasThresholdLower )
                __str << " -lowerthreshold " << M_thresholdLower << " ";
            if ( M_hasThresholdUpper )
                __str << " -upperthreshold " << M_thresholdUpper << " ";
        }
        else if ( M_method == "isosurface" )
        {
            __str << "-method isosurface -isosurface 0 ";
        }

        __str << "-olevelsetsfile " << outputPathImageInit;

        std::cout << "---------------------------------------\n"
                  << "run in system : \n" << __str.str() << "\n"
                  << "---------------------------------------\n";
        auto err = ::system( __str.str().c_str() );
#if 0
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
#else
        std::string nameImageLevelSet = fs::path(this->outputPath()).stem().string()+"_levelset.vti";
        std::string outputPathImageLevelSetInit = (fs::path(this->outputPath()).parent_path()/fs::path(nameImageLevelSet)).string();

        std::ostringstream __str2;
        __str2 << pythonExecutable << " ";
        __str2 << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtklevelsetsegmentation ";
        __str2 << "-ifile " << this->inputPath() << " ";
        //__str2 << "-levelsetstype isosurface -isosurfacevalue 0 ";
        __str2 << "-initiallevelsetsfile " << outputPathImageInit << " ";
        //__str2 << "-initializationimagefile " << outputPathImageInit << " ";
        __str2 << "-iterations 1 -ofile " << outputPathImageLevelSetInit;
        std::cout << "---------------------------------------\n"
                  << "run in system : \n" << __str2.str() << "\n"
                  << "---------------------------------------\n";
        auto err2 = ::system( __str2.str().c_str() );


        // read image obtained by levelset segmentation
        vtkSmartPointer<vtkXMLImageDataReader> readerSegmentation = vtkSmartPointer<vtkXMLImageDataReader>::New();
        readerSegmentation->SetFileName(outputPathImageLevelSetInit.c_str());
        readerSegmentation->Update();
        //int* boundIMG = reader->GetWholeExtent();
        //vtkDataSet * hhh = reader->GetOutputAsDataSet();

        // read initial image
        vtkSmartPointer<vtkMetaImageReader> readerInitialImage = vtkSmartPointer<vtkMetaImageReader>::New();
        readerInitialImage->SetFileName(this->inputPath().c_str());
        readerInitialImage->Update();
#if 0
        double boundIMG[6];
        readerSegmentation->GetOutput()->GetBounds(boundIMG);
        std::cout << " boundIMG : \n"
                  << boundIMG[0] << " , " << boundIMG[1] << "\n"
                  << boundIMG[2] << " , " << boundIMG[3] << "\n"
                  << boundIMG[4] << " , " << boundIMG[5] << "\n";
        double boundIMG2[6];
        readerInitialImage->GetOutput()->GetBounds(boundIMG2);
        std::cout << " boundIMG2 : \n"
                  << boundIMG2[0] << " , " << boundIMG2[1] << "\n"
                  << boundIMG2[2] << " , " << boundIMG2[3] << "\n"
                  << boundIMG2[4] << " , " << boundIMG2[5] << "\n";
        std::cout << " diff boundIMG : \n"
                  << boundIMG[1]-boundIMG[0] << " vs " << boundIMG2[1]-boundIMG2[0] << "\n"
                  << boundIMG[3]-boundIMG[2] << " vs " << boundIMG2[3]-boundIMG2[2] << "\n"
                  << boundIMG[5]-boundIMG[4] << " vs " << boundIMG2[5]-boundIMG2[4] << "\n";
#endif
        double originSegmentationImage[3];
        readerSegmentation->GetOutput()->GetOrigin(originSegmentationImage);
        std::cout << "originSegmentationImage " << originSegmentationImage[0] << "," << originSegmentationImage[1] << "," << originSegmentationImage[2] << "\n";

        double originInitialImage[3];
        readerInitialImage->GetOutput()->GetOrigin(originInitialImage);
        std::cout << "originInitialImage " << originInitialImage[0] << "," << originInitialImage[1] << "," << originInitialImage[2] << "\n";
        readerSegmentation->GetOutput()->SetOrigin(originInitialImage);
#if 0
        double bary1[3] = { (boundIMG[1]+boundIMG[0])/2.,(boundIMG[3]+boundIMG[2])/2.,(boundIMG[5]+boundIMG[4])/2. };
        double bary2[3] = { (boundIMG2[1]+boundIMG2[0])/2.,(boundIMG2[3]+boundIMG2[2])/2.,(boundIMG2[5]+boundIMG2[4])/2. };
        double trans[3] = { bary2[0]-bary1[0],bary2[1]-bary1[1],bary2[2]-bary1[2] };
        std::cout << "trans " << trans[0] << ","<<trans[1] << "," << trans[2] << "\n";
        vtkSmartPointer<vtkImageTranslateExtent> imageTranslated = vtkSmartPointer<vtkImageTranslateExtent>::New();
        imageTranslated->SetInput(reader->GetOutput());
        imageTranslated->SetTranslation(trans[0],trans[1],trans[2]);
        imageTranslated->Update();

        double boundIMGFinal[6];
        imageTranslated->GetOutput()->GetBounds(boundIMGFinal);
        std::cout << " boundIMGFinal : \n"
                  << boundIMGFinal[0] << " , " << boundIMGFinal[1] << "\n"
                  << boundIMGFinal[2] << " , " << boundIMGFinal[3] << "\n"
                  << boundIMGFinal[4] << " , " << boundIMGFinal[5] << "\n";
#endif

        // marching cube
        vtkSmartPointer<vtkMarchingCubes> surface = vtkSmartPointer<vtkMarchingCubes>::New();
        //#if VTK_MAJOR_VERSION <= 5
        surface->SetInput(readerSegmentation->GetOutput());
        //surface->SetInput(imageTranslated->GetOutput());
        //#else
        //surface->SetInputData(volume);
        //#endif
        //surface->ComputeNormalsOn();
        double isoValue=0.0;
        surface->SetValue(0, isoValue);
        surface->Update();
        //std::cout<<"Number of points: " << surface->GetNumberOfPoints() << std::endl;

#if 0
        // create polydata from iso-surface
        vtkSmartPointer<vtkPolyData> marched = vtkSmartPointer<vtkPolyData>::New();
        surface->Update();
        marched->DeepCopy(surface->GetOutput());
        std::cout<<"Number of points: " << marched->GetNumberOfPoints() << std::endl;
#endif

         // clean surface : ensure that mesh has no duplicate entities
        vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
        //cleanPolyData->SetInputConnection(surface->GetOutputPort());
        cleanPolyData->SetInput(surface->GetOutput());
        cleanPolyData->Update();

        // keep only largest region
        vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
        if ( M_applyConnectivityLargestRegion )
        {
            connectivityFilter->SetInput(/*surface*/cleanPolyData->GetOutput());
            connectivityFilter->SetExtractionModeToLargestRegion();
            connectivityFilter->Update();
        }

        // save stl on disk
        vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
        stlWriter->SetFileName(this->outputPath().c_str());
        if ( M_applyConnectivityLargestRegion )
            stlWriter->SetInput(connectivityFilter->GetOutput());
        else
            stlWriter->SetInput(/*surface*/cleanPolyData->GetOutput());
        //stlWriter->SetInputConnection(surface->GetOutputPort());
        stlWriter->Write();

#endif
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
        (prefixvm(prefix,"method").c_str(), Feel::po::value<std::string>()->default_value("threshold"), "(string) method : threshold, isosurface")
        (prefixvm(prefix,"threshold.lower").c_str(), Feel::po::value<double>(), "(double) threshold lower")
        (prefixvm(prefix,"threshold.upper").c_str(), Feel::po::value<double>(), "(double) threshold upper")
        (prefixvm(prefix,"apply-connectivity.largest-region").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) apply-connectivity.largest-region")
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        ;
    return mySurfaceFromImageOptions;
}



SubdivideSurface::SubdivideSurface( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.filename",_prefix=this->prefix()) ) ),
    M_outputDirectory( AngioTkEnvironment::expand( soption(_name="output.directory",_prefix=this->prefix()) ) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
    M_method( soption(_name="method",_prefix=this->prefix()) ),
    M_nSubdivisions( ioption(_name="subdivisions",_prefix=this->prefix()) )
{
    CHECK( M_method == "linear" || M_method == "butterfly" || M_method == "loop" ) << "invalid method " << M_method << "\n";
    if ( !this->inputSurfacePath().empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}



void
SubdivideSurface::updateOutputPathFromInputFileName()
{
    CHECK( !this->inputSurfacePath().empty() ) << "input path is empty";

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
    std::string nameMeshFile = fs::path(this->inputSurfacePath()).stem().string();

    std::string newFileName = (boost::format("%1%_subdivide%2%%3%.stl")%nameMeshFile %M_nSubdivisions %M_method ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    M_outputPath = outputPath.string();
}

void
SubdivideSurface::run()
{
    if ( !fs::exists( this->inputSurfacePath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : smoothing surface not done because this input surface path not exist :" << this->inputSurfacePath() << "\n";
        return;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run SubdivideSurface \n"
            << "---------------------------------------\n";
    coutStr << "input surface : " << this->inputSurfacePath() << "\n";
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
        __str << "-ifile " << this->inputSurfacePath() << " "
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
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.filename",_prefix=this->prefix()) ) ),
    M_inputCenterlinesPath( AngioTkEnvironment::expand( soption(_name="input.centerlines.filename",_prefix=this->prefix()) ) ),
    M_outputDirectory( AngioTkEnvironment::expand( soption(_name="output.directory",_prefix=this->prefix()) ) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
    M_method( soption(_name="method",_prefix=this->prefix()) ),
    M_nIterations( ioption(_name="iterations",_prefix=this->prefix()) ),
    M_taubinPassBand( doption(_name="taubin.passband",_prefix=this->prefix()) ),
    M_laplaceRelaxationFactor( doption(_name="laplace.relaxation",_prefix=this->prefix()) )
{
    CHECK( M_method == "taubin" || M_method == "laplace" || M_method == "centerlines" ) << "invalid method " << M_method << "\n";
    if ( !this->inputSurfacePath().empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}

void
SmoothSurface::updateOutputPathFromInputFileName()
{
    CHECK( !this->inputSurfacePath().empty() ) << "input path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( M_outputDirectory.empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(M_outputDirectory).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(M_outputDirectory);
    else
        meshesdirectories = fs::path(M_outputDirectory);

    // get filename without extension
    fs::path gp = this->inputSurfacePath();
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
    if ( !fs::exists( this->inputSurfacePath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : smoothing surface not done because this input surface path not exist :" << this->inputSurfacePath() << "\n";
        return;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run SmoothSurface \n"
            << "---------------------------------------\n";
    coutStr << "input surface : " << this->inputSurfacePath() << "\n";
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
        __str << "-ifile " << this->inputSurfacePath() << " "
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
        (prefixvm(prefix,"input.centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
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
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.surface.filename",_prefix=this->prefix()) ) ),
    M_inputCenterlinesPath( AngioTkEnvironment::expand( soption(_name="input.centerlines.filename",_prefix=this->prefix()) ) ),
    M_outputDirectory( AngioTkEnvironment::expand( soption(_name="output.directory",_prefix=this->prefix()) ) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
    M_distanceClipScalingFactor( doption(_name="distance-clip.scaling-factor",_prefix=this->prefix() ) ),
    M_saveOutputSurfaceBinary( boption(_name="output.save-binary",_prefix=this->prefix() ) )
{
    if ( !M_inputSurfacePath.empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}


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
            << "output file type     : " << std::string((M_saveOutputSurfaceBinary)? "binary" : "ascii") << "\n"
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
    std::string method = "gmsh";
    if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )
    {
        if ( method == "gmsh" )
            this->runGMSH();
        else if ( method == "gmsh-executable" )
            this->runGMSHwithExecutable();
        else if ( method == "vmtk" )
            this->runVMTK();
    }

}


void
OpenSurface::runGMSH()
{
    CHECK( !this->inputSurfacePath().empty() && fs::exists(this->inputSurfacePath()) ) << "inputSurfacePath is empty or not exist";
    CHECK( !this->inputCenterlinesPath().empty() && fs::exists(this->inputCenterlinesPath()) ) << "inputCenterlinesPath is empty or not exist";

    if ( true )
        {
            GmshInitialize();
            CTX::instance()->terminal = 1;
            int verbosityLevel = 5;
            Msg::SetVerbosity( verbosityLevel );
            CTX::instance()->geom.tolerance=1e-8;
        }
    GModel * gmodel = new GModel();
    // add new model as current (important if this function is called more than 1 time)
    GModel::current(GModel::list.size() - 1);

    std::shared_ptr<AngioTkCenterline> centerlinesTool( new AngioTkCenterline );
    centerlinesTool->setModeClipMesh(true);
    centerlinesTool->setClipMeshScalingFactor( M_distanceClipScalingFactor );
    centerlinesTool->importSurfaceFromFile( this->inputSurfacePath() );
    centerlinesTool->importFile( this->inputCenterlinesPath() );
    centerlinesTool->runClipMesh();
    //CTX::instance()->mesh.binary = M_saveOutputSurfaceBinary;
    centerlinesTool->saveClipMeshSTL(this->outputPath(), M_saveOutputSurfaceBinary );

    delete gmodel;
}

void
OpenSurface::runGMSHwithExecutable()
{
    std::ostringstream geodesc;
    geodesc << "Mesh.Algorithm = 6; //(1=MeshAdapt, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) \n"
            << "Mesh.Algorithm3D = 1; //(1=tetgen, 4=netgen, 7=MMG3D, 9=R-tree) \n"
            << "Mesh.Binary = " << M_saveOutputSurfaceBinary <<";\n";

    geodesc << "Merge \""<< this->inputSurfacePath() <<"\";\n";
    geodesc << "Field[1] = AngioTkCenterline;\n";
    geodesc << "Field[1].FileName = \"" << this->inputCenterlinesPath() << "\";\n";
    geodesc << "Field[1].clipMesh =1;\n"
            << "Field[1].clipMeshScalingFactor = " << M_distanceClipScalingFactor << ";\n"
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

void
OpenSurface::runVMTK()
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



po::options_description
OpenSurface::options( std::string const& prefix )
{
    po::options_description myOpenSurfaceOptions( "Open the surface options" );

    myOpenSurfaceOptions.add_options()
        (prefixvm(prefix,"input.centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"input.surface.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        (prefixvm(prefix,"distance-clip.scaling-factor").c_str(), Feel::po::value<double>()->default_value(0.0), "(double) scaling-factor")
        (prefixvm(prefix,"output.save-binary").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) save-binary")
        ;
    return myOpenSurfaceOptions;
}




RemeshSTL::RemeshSTL( std::string prefix )
    :
    M_prefix( prefix ),
    M_packageType(soption(_name="package-type",_prefix=this->prefix())),
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.filename",_prefix=this->prefix())) ),
    M_inputCenterlinesPath( AngioTkEnvironment::expand(soption(_name="centerlines.filename",_prefix=this->prefix())) ),
    M_gmshRemeshNbPointsInCircle( ioption(_name="nb-points-in-circle",_prefix=this->prefix()) ),
    M_gmshRemeshPartitionForceRebuild( boption(_name="gmsh.remesh-partition.force-rebuild",_prefix=this->prefix()) ),
    M_vmtkArea( doption(_name="area",_prefix=this->prefix()) ),
    M_vmtkNumberOfIteration( ioption(_name="vmtk.n-iteration",_prefix=this->prefix()) ),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
    M_saveOutputSurfaceBinary( boption(_name="output.save-binary",_prefix=this->prefix() ) )
{
    CHECK( M_packageType == "gmsh" || M_packageType == "gmsh-executable" || M_packageType == "vmtk" ) << "error on packageType : " << M_packageType;
    if ( !this->inputSurfacePath().empty() && M_outputPathGMSH.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}


void
RemeshSTL::updateOutputPathFromInputFileName()
{
    CHECK( !this->inputSurfacePath().empty() ) << "input surface path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( M_outputDirectory.empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(M_outputDirectory).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(M_outputDirectory);
    else
        meshesdirectories = fs::path(M_outputDirectory);

    // get filename without extension
    fs::path gp = this->inputSurfacePath();
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
              << "inputSurfacePath     : " << this->inputSurfacePath() << "\n";
    if ( this->packageType() == "gmsh" || this->packageType() == "gmsh-executable" )
        std::cout << "inputCenterlinesPath : " << this->inputCenterlinesPath() << "\n";
    else if ( this->packageType() == "vmtk" )
        std::cout << "area                 : " << this->area() << "\n";
    std::cout << "output path          : " << this->outputPath() << "\n";
    if ( this->packageType() == "gmsh" || this->packageType() == "gmsh-executable" )
        std::cout << "output file type     : " << std::string((M_saveOutputSurfaceBinary)? "binary" : "ascii") << "\n";
    std::cout << "---------------------------------------\n"
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
    else if ( this->packageType() == "gmsh-executable" )
        this->runGMSHwithExecutable();
    else if ( this->packageType() == "vmtk" )
        this->runVMTK();

    std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n";
}

void
RemeshSTL::runVMTK()
{
    CHECK( !this->inputSurfacePath().empty() ) << "inputSurfacePath is empty";

    std::ostringstream __str;
    // source ~/packages/vmtk/vmtk.build2/Install/vmtk_env.sh
    std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
    __str << pythonExecutable << " ";
    //std::string dirBaseVmtk = "/Users/vincentchabannes/packages/vmtk/vmtk.build2/Install/bin/";
    std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );
    __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtksurfaceremeshing ";
    __str << "-ifile " << this->inputSurfacePath() << " ";
    __str << "-ofile " << this->outputPath() << " ";
    __str << "-area " << this->area() << " ";
    __str << "-iterations " << M_vmtkNumberOfIteration << " ";
    auto err = ::system( __str.str().c_str() );

    //std::cout << "hola\n"<< __str.str() <<"\n";
}


void
RemeshSTL::runGMSHwithExecutable()
{
    CHECK( !this->inputSurfacePath().empty() ) << "inputSurfacePath is empty";
    CHECK( !this->inputCenterlinesPath().empty() ) << "inputCenterlinesPath is empty";

    std::ostringstream geodesc;
    geodesc << "Mesh.Algorithm = 6; //(1=MeshAdapt, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) \n"
            << "Mesh.Algorithm3D = 1; //(1=tetgen, 4=netgen, 7=MMG3D, 9=R-tree) \n"
            << "Mesh.Binary = " << M_saveOutputSurfaceBinary <<";\n";

    //geodesc << "Merge \"stl_remesh_vmtk/fluidskin3.stl\"\n";
    geodesc << "Merge \""<< this->inputSurfacePath() <<"\";\n";

    //geodesc << "Field[1] = Centerline;\n";
    geodesc << "Field[1] = AngioTkCenterline;\n";


    //geodesc << "Field[1].FileName = \"../centerlines/fluidskin3.vtk\"\n";
    geodesc << "Field[1].FileName = \"" << this->inputCenterlinesPath() << "\";\n";
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

void
RemeshSTL::runGMSH()
{
    CHECK( !this->inputSurfacePath().empty() && fs::exists(this->inputSurfacePath()) ) << "inputSurfacePath is empty or not exist";
    CHECK( !this->inputCenterlinesPath().empty() && fs::exists(this->inputCenterlinesPath()) ) << "inputCenterlinesPath is empty or not exist";

    if ( true )
        {
            GmshInitialize();
            CTX::instance()->terminal = 1;
            int verbosityLevel = 5;
            Msg::SetVerbosity( verbosityLevel );
            CTX::instance()->geom.tolerance=1e-8;
        }
    //CTX::instance()->mesh.randFactor = 1e-10;//1e-8;//1e-14;//1e-10;
    CTX::instance()->mesh.algo2d = ALGO_2D_FRONTAL;//ALGO_2D_DELAUNAY;//ALGO_2D_FRONTAL //ALGO_2D_BAMG

    GModel * gmodel = new GModel();
    // add new model as current (important if this function is called more than 1 time)
    GModel::current(GModel::list.size() - 1);

    std::shared_ptr<AngioTkCenterline> centerlinesTool( new AngioTkCenterline );
    centerlinesTool->setIsCut(true);
    centerlinesTool->setRemeshNbPoints( this->remeshNbPointsInCircle() );
    centerlinesTool->importSurfaceFromFile( this->inputSurfacePath() );
    centerlinesTool->importFile( this->inputCenterlinesPath() );

    fs::path gp = this->inputSurfacePath();
    std::string remeshPartitionMeshFileName = gp.stem().string() + "_remeshpartition.msh";
    fs::path directory = fs::path(this->outputPath()).parent_path();
    std::string remeshPartitionMeshFilePath = (directory/fs::path(remeshPartitionMeshFileName)).string();
    centerlinesTool->runSurfaceRemesh(remeshPartitionMeshFilePath,M_gmshRemeshPartitionForceRebuild);
    //CTX::instance()->mesh.binary = M_saveOutputSurfaceBinary;
    centerlinesTool->saveSurfaceRemeshSTL(this->outputPath(), M_saveOutputSurfaceBinary );

    delete gmodel;
}


po::options_description
RemeshSTL::options( std::string const& prefix )
{
    po::options_description myMeshSurfaceOptions( "Mesh surface blood flow from STL and Centerlines options" );
    myMeshSurfaceOptions.add_options()

        ( prefixvm(prefix,"package-type").c_str(), po::value<std::string>()->default_value( "vmtk" ), "force-remesh" )
        ( prefixvm(prefix,"input.filename").c_str(), po::value<std::string>()->default_value( "" ), "stl.filename" )
        ( prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        ( prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        ( prefixvm(prefix,"output.save-binary").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) save-binary")
        // gmsh options
        ( prefixvm(prefix,"centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "centerlines.filename" )
        ( prefixvm(prefix,"nb-points-in-circle").c_str(), po::value<int>()->default_value( 15 ), "nb-points-in-circle" )
        ( prefixvm(prefix,"gmsh.remesh-partition.force-rebuild").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) force-rebuild")
        // vmtk options
        ( prefixvm(prefix,"area").c_str(), po::value<double>()->default_value( 0.5 ), "area" )
        ( prefixvm(prefix,"vmtk.n-iteration").c_str(), po::value<int>()->default_value( 10 ), "maxit" )
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

    //CTX::instance()->mesh.randFactor = 1e-10;//1e-8;//1e-14;//1e-10;
    //CTX::instance()->mesh.algo2d = ALGO_2D_DELAUNAY;//ALGO_2D_BAMG

#if 0
  if(GModel::current()->empty()){
      // if the current model is empty, make sure it's reaaally cleaned-up, and
      // reuse it
      GModel::current()->destroy();
      GModel::current()->getGEOInternals()->destroy();
  }
  else{
    // if the current model is not empty make it invisible and add a new model
    new GModel();
    GModel::current(GModel::list.size() - 1);
  }
#endif


    GModel * M_gmodel = new GModel();
    // add new model as current (important if this function is called more than 1 time)
    GModel::current(GModel::list.size() - 1);

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
        M_gmodel->writeSTL( outputMeshName, 1/*CTX::instance()->mesh.binary*/ );
    else
        CHECK( false ) << "error \n";


    delete M_gmodel;
}

} //  namespace detail



VolumeMeshing::VolumeMeshing( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputSTLPath( AngioTkEnvironment::expand( soption(_name="input.stl.filename",_prefix=this->prefix()) ) ),
    M_inputCenterlinesPath( AngioTkEnvironment::expand( soption(_name="input.centerlines.filename",_prefix=this->prefix())) ),
    M_inputInletOutletDescPath( AngioTkEnvironment::expand( soption(_name="input.desc.filename",_prefix=this->prefix())) ),
    M_remeshNbPointsInCircle( ioption(_name="nb-points-in-circle",_prefix=this->prefix() ) ),
    M_extrudeWall( boption(_name="extrude-wall",_prefix=this->prefix() ) ),
    M_extrudeWallNbElemLayer( ioption(_name="extrude-wall.nb-elt-layer",_prefix=this->prefix() ) ),
    M_extrudeWallhLayer( doption(_name="extrude-wall.h-layer",_prefix=this->prefix()) ),
    M_outputPath(),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) ),
    M_saveOutputVolumeBinary( boption(_name="output.save-binary",_prefix=this->prefix() ) )
{
    if ( !M_inputSTLPath.empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}


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
    std::cout << "input path              : " << this->inputSTLPath() << "\n"
              << "centerlines path        : " << this->inputCenterlinesPath() << "\n"
              << "inletoutlet desc path   : " << this->inputInletOutletDescPath() << "\n"
              << "extrude arterial wall   : " << std::boolalpha << M_extrudeWall << "\n";
    if ( M_extrudeWall )
        std::cout << "arterial wall thickness (radius percent) : " << M_extrudeWallhLayer << "\n"
                  << "arterial wall number of layer            : " << M_extrudeWallNbElemLayer <<"\n";
    std::cout << "output path             : " << this->outputPath() << "\n"
              << "output file type        : " << std::string((M_saveOutputVolumeBinary)? "binary" : "ascii") << "\n"
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
            << "Mesh.Algorithm3D = 1; //(1=tetgen, 4=netgen, 7=MMG3D, 9=R-tree) \n"
            << "Mesh.Binary = " << M_saveOutputVolumeBinary <<";\n";

    //CTX::instance()->mesh.binary = M_saveOutputSurfaceBinary;

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
        ( prefixvm(prefix,"output.save-binary").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) save-binary")
        ( prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        ( prefixvm(prefix,"nb-points-in-circle").c_str(), po::value<int>()->default_value( 15 ), "nb-points-in-circle" )
        ( prefixvm(prefix,"extrude-wall").c_str(),po::value<bool>()->default_value( false ), "extrude-wall" )
        ( prefixvm(prefix,"extrude-wall.nb-elt-layer").c_str(), po::value<int>()->default_value( 2 ), "nb-elt-layer" )
        ( prefixvm(prefix,"extrude-wall.h-layer").c_str(), po::value<double>()->default_value( 0.2 ), "h-layer" )
        ;
    return myMeshVolumeOptions;
}





} // namespace Feel
