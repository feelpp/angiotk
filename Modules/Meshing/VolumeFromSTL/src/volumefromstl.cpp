/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <volumefromstl.hpp>
#include <AngioTkCenterlineField.h>

#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkMetaImageReader.h> // .mha
#include <vtkImageChangeInformation.h>
#include <vtkXMLImageDataReader.h> // .vti
#include <vtkMarchingCubes.h>
#include <vtkSTLWriter.h>
#include <vtkImageTranslateExtent.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPointData.h>
//#include <vtkCellData.h>
#include <vtkCleanPolyData.h>

//#include <vtkvmtkPolyBallModeller.h>
#include <angiotkPolyBallModeller.h>

#include <vtkvmtkPolyDataCenterlines.h>
#include <vtkSTLReader.h>
#include <vtkPointLocator.h>
#include <vtkIdList.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>


Feel::fs::path AngioTkEnvironment::S_pathInitial;
boost::shared_ptr<Feel::Environment> AngioTkEnvironment::S_feelEnvironment;

namespace AngioTk
{
    using namespace Feel;

AngioTkFilterBase::AngioTkFilterBase( std::string const& prefix )
    :
    M_prefix( prefix ),
    M_outputDirectory( AngioTkEnvironment::expand( soption(_name="output.directory",_prefix=this->prefix()) ) ),
    M_outputPath( AngioTkEnvironment::expand( soption(_name="output.path",_prefix=this->prefix()) ) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) )
{
    if ( !this->outputPath().empty() )
        {
            //if ( !fs::exists( this->outputPath() ) )
            //     M_outputPath.clear();
            if ( fs::path(M_outputPath).is_relative() )
                M_outputPath = (fs::path(Feel::Environment::rootRepository())/fs::path(M_outputPath)).string();
            this->updateOutputDirFromOutputPath();
        }
}
po::options_description
AngioTkFilterBase::options( std::string const& prefix )
{
    po::options_description myFilterBaseOptions( "Modified an Image" );

    myFilterBaseOptions.add_options()
        (prefixvm(prefix,"output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"output.path").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output filename")
        (prefixvm(prefix,"force-rebuild").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) force-rebuild")
        ;
    return myFilterBaseOptions;
}


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

int
InletOutletDesc::loadFromSTL( std::string inputPath )
{

    std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
    std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );

    //if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )
    if ( !fs::exists(inputPath) )
    {
        std::cout << "WARNING!, the inputPath does not exist : " << inputPath <<"\n";
        return 1;
    }

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

    if ( !fs::exists(outputPath) )
    {
        std::cout << "WARNING!, outputPath does not exist : " << outputPath <<"\n";
        return 1;
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

    return 0;
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

void
InletOutletDesc::saveJSON( std::string outputPath )
{
    std::cout << "save json in : " << outputPath << "\n";
    fs::path dir = fs::path(outputPath).parent_path();
    if ( !fs::exists( dir ) )
        fs::create_directories( dir );
    
    boost::property_tree::ptree root, entries;
    for ( auto const& desc : *this )
    {
        boost::property_tree::ptree entry, marker, point, pdata;
        entry.put("marker", desc.markerLumen());
        for( int i = 0; i < desc.node().size(); i++)
        {
            pdata.put("", desc.node().at(i));
            point.push_back(std::make_pair("", pdata));
        }
        entry.add_child("point", point);
        root.add_child("fluid", entry);
    }
    write_json(outputPath, root);
}


CenterlinesFromSurface::CenterlinesFromSurface( std::string const& prefix )
    :
    super_type( prefix ),
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.surface.filename",_prefix=this->prefix()) ) ),
    M_inputCenterlinesPointSetPath( AngioTkEnvironment::expand( soption(_name="input.pointset.filename",_prefix=this->prefix()) ) ),
    M_inputCenterlinesPointPairPath( AngioTkEnvironment::expand( soption(_name="input.pointpair.filename",_prefix=this->prefix()) ) ),
    M_inputInletOutletDescPath( AngioTkEnvironment::expand( soption(_name="input.desc.filename",_prefix=this->prefix()) ) ),
    M_inputGeoCenterlinesPath( AngioTkEnvironment::expand( soption(_name="input.geo-centerlines.filename",_prefix=this->prefix()) ) ),
    M_costFunctionExpr( soption(_name="cost-function.expression",_prefix=this->prefix()) ),
    M_useInteractiveSelection( boption(_name="use-interactive-selection",_prefix=this->prefix()) ),
    M_viewResults( boption(_name="view-results",_prefix=this->prefix() ) ),
    M_viewResultsWithSurface( boption(_name="view-results.with-surface",_prefix=this->prefix() ) ),
    M_delaunayTessellationOutputDirectory( AngioTkEnvironment::expand( soption(_name="delaunay-tessellation.output.directory",_prefix=this->prefix()) ) ),
    M_delaunayTessellationForceRebuild( boption(_name="delaunay-tessellation.force-rebuild",_prefix=this->prefix() ) )
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

    if ( !M_inputSurfacePath.empty() && this->outputPath().empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}


void
CenterlinesFromSurface::updateOutputPathFromInputFileName()
{
    CHECK( !this->inputSurfacePath().empty() ) << "input path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( this->outputDirectory().empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(this->outputDirectory()).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(this->outputDirectory());
    else
        meshesdirectories = fs::path(this->outputDirectory());

    // get filename without extension
    fs::path gp = this->inputSurfacePath();
    std::string nameMeshFile = gp.stem().string();

    std::string newFileName = (boost::format("%1%_centerlines.vtk")%nameMeshFile ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    this->setOutputPath( outputPath.string() );
}

std::tuple< std::vector<std::vector<double> >, std::vector<std::vector<double> > >
CenterlinesFromSurface::loadFromCenterlinesPointSetFile()
{
    if ( !fs::exists( this->inputCenterlinesPointSetPath() ) )
        return std::tuple< std::vector<std::vector<double> >, std::vector<std::vector<double> > >();

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
          if ( typePt == 0 )
              sourcePts.push_back(pt);
          else if ( typePt == 1 )
              targetPts.push_back(pt);
      }
    fileLoaded.close();

    auto res = std::make_tuple( sourcePts, targetPts );
    return res;
}

std::vector< std::pair< std::vector<double>,std::vector<double> > >
CenterlinesFromSurface::loadFromCenterlinesPointPairFile()
{
    std::vector< std::pair< std::vector<double>,std::vector<double> > > res;
    if ( !fs::exists( this->inputCenterlinesPointPairPath() ) )
        return res;

    std::ifstream fileLoaded( this->inputCenterlinesPointPairPath(), std::ios::in);
    while ( !fileLoaded.eof() )
    {
        std::vector<double> pt1(3);
        std::vector<double> pt2(3);
        double radius1=0,radius2=0;
        int typePt1 = -1,typePt2 = -1;
        fileLoaded >> typePt1;
        if ( fileLoaded.eof() ) break;
        fileLoaded >> pt1[0] >> pt1[1] >> pt1[2] >> radius1;

        fileLoaded >> typePt2;
        if ( fileLoaded.eof() ) break;
        fileLoaded >> pt2[0] >> pt2[1] >> pt2[2] >> radius2;

        res.push_back( std::make_pair(pt1,pt2) );
    }
    fileLoaded.close();
    return res;
}


int
CenterlinesFromSurface::run()
{
    if ( !fs::exists( this->inputSurfacePath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : centerlines computation will not be done, because the input path does not exist :" << this->inputSurfacePath() << std::endl
                      << "Please set it with the \"input.surface.filename\" option." << std::endl;
        return 1;
    }
    if ( this->inputCenterlinesPointSetPath().empty() && this->inputCenterlinesPointPairPath().empty() &&
         this->sourceids().empty() && this->targetids().empty() && this->inputGeoCenterlinesPath().empty() )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : centerlines computation not done because this sourceids and targetids are empty\n";
        return 1;
    }


    fs::path directory = fs::path(this->outputPath()).parent_path();
    fs::path outputFileNamePath = fs::path(this->outputPath());
    std::string name = outputFileNamePath.stem().string();

    fs::path directoryDelaunayTessellation = M_delaunayTessellationOutputDirectory;
    if ( M_delaunayTessellationOutputDirectory.empty() )
        directoryDelaunayTessellation = directory;
    else if ( fs::path(M_delaunayTessellationOutputDirectory).is_relative() )
        directoryDelaunayTessellation = fs::path(Environment::rootRepository())/fs::path(M_delaunayTessellationOutputDirectory);
    std::string fileNameDelaunayTessellation = (boost::format("%1%_DelaunayTessellation" )%name ).str();
    std::string pathDelaunayTessellation = (directoryDelaunayTessellation/fs::path(fileNameDelaunayTessellation+".vtk")).string();
    bool rebuildDelaunayTessellation = M_delaunayTessellationForceRebuild || !fs::exists(pathDelaunayTessellation);


    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run CenterlinesFromSurface \n"
            << "---------------------------------------\n";
    coutStr << "input surface path   : " << this->inputSurfacePath() << "\n";
    if ( !M_useInteractiveSelection )
    {
        if ( !this->inputGeoCenterlinesPath().empty() )
            coutStr << "input GeoCenterlines : " << this->inputGeoCenterlinesPath() << "\n";
        else if ( !this->inputCenterlinesPointPairPath().empty() )
            coutStr << "input point pair file : " << this->inputCenterlinesPointPairPath() << "\n";
        else if ( !this->inputCenterlinesPointSetPath().empty() )
            coutStr << "input point set file : " << this->inputCenterlinesPointSetPath() << "\n";
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
        coutStr << "DelaunayTessellation path  : " << pathDelaunayTessellation << "\n"
                << "DelaunayTessellation build : " << std::string(( rebuildDelaunayTessellation )? "true":"false") << "\n";
    }

    coutStr << "output path          : " << this->outputPath() << "\n"
            << "---------------------------------------\n"
            << "---------------------------------------\n";
    std::cout << coutStr.str();


    // build directories if necessary
    if ( this->worldComm().isMasterRank() )
    {
        if ( !fs::exists( directory ) )
            fs::create_directories( directory );
        if ( !fs::exists( directoryDelaunayTessellation ) )
            fs::create_directories( directoryDelaunayTessellation );
    }
    // // wait all process
    this->worldComm().globalComm().barrier();



    //fs::path stlNamePath = fs::path(this->inputSurfacePath());
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
            centerlinesTool.createFromFile( this->inputGeoCenterlinesPath(),this->inputSurfacePath() );
            centerlinesTool.addFieldBranchIds();
            centerlinesTool.addFieldRadiusMin("MaximumInscribedSphereRadius");
            centerlinesTool.writeCenterlinesVTK( this->outputPath() );
        }
        else if ( !this->inputCenterlinesPointPairPath().empty() || !this->inputCenterlinesPointSetPath().empty() ||
                  !this->inputInletOutletDescPath().empty() ) // use c++ vmtkcenterlines
        {
            vtkSmartPointer<vtkSTLReader> readerSTL = vtkSmartPointer<vtkSTLReader>::New();
            readerSTL->SetFileName( this->inputSurfacePath().c_str());
            readerSTL->Update();

            vtkSmartPointer<vtkvmtkPolyDataCenterlines> centerlineFilter = vtkvmtkPolyDataCenterlines::New();
            centerlineFilter->SetInput( readerSTL->GetOutput() );

            vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
            pointLocator->SetDataSet( readerSTL->GetOutput() );
            pointLocator->BuildLocator();

            std::vector<std::pair<vtkSmartPointer<vtkIdList>,vtkSmartPointer<vtkIdList> > > vecSourceTargetIdList;
            if ( !this->inputCenterlinesPointPairPath().empty() )
            {
                auto vecpointpair = loadFromCenterlinesPointPairFile();
                CHECK( !vecpointpair.empty() ) << "PointPairFile is empty or invalid or not exist";
                for (auto const& pointpair : vecpointpair)
                {
                    vtkSmartPointer<vtkIdList> sourceIdList = vtkSmartPointer<vtkIdList>::New();
                    vtkSmartPointer<vtkIdList> targetIdList = vtkSmartPointer<vtkIdList>::New();
                    auto const& sourcePt = pointpair.first;
                    auto const& targetPt = pointpair.second;
                    double sourcePoint[3] = { sourcePt[0], sourcePt[1], sourcePt[2] };
                    double targetPoint[3] = { targetPt[0], targetPt[1], targetPt[2] };
                    vtkIdType vtkSourceId = pointLocator->FindClosestPoint(sourcePoint);
                    vtkIdType vtkTargetId = pointLocator->FindClosestPoint(targetPoint);
                    sourceIdList->InsertNextId(vtkSourceId);
                    targetIdList->InsertNextId(vtkTargetId);
                    vecSourceTargetIdList.push_back( std::make_pair(sourceIdList,targetIdList) );
                }
            }
            else if ( !this->inputCenterlinesPointSetPath().empty() )
            {
                auto pointset = loadFromCenterlinesPointSetFile();
                CHECK( !std::get<0>(pointset).empty() && !std::get<1>(pointset).empty() ) << "PointSetFile has no source or target points";
                vtkSmartPointer<vtkIdList> sourceIdList = vtkSmartPointer<vtkIdList>::New();
                vtkSmartPointer<vtkIdList> targetIdList = vtkSmartPointer<vtkIdList>::New();
                for ( std::vector<double> const& pt : std::get<0>(pointset) ) // source
                {
                    double point[3] = { pt[0], pt[1], pt[2] };
                    vtkIdType vtkid = pointLocator->FindClosestPoint(point);
                    sourceIdList->InsertNextId(vtkid);
                }
                for ( std::vector<double> const& pt : std::get<1>(pointset) ) // target
                {
                    double point[3] = { pt[0], pt[1], pt[2] };
                    vtkIdType vtkid = pointLocator->FindClosestPoint(point);
                    targetIdList->InsertNextId(vtkid);
                }
                vecSourceTargetIdList.push_back( std::make_pair(sourceIdList,targetIdList) );
            }
            else if ( !this->inputInletOutletDescPath().empty() )
            {
                InletOutletDesc iodesc( this->inputInletOutletDescPath() );
                vtkSmartPointer<vtkIdList> sourceIdList = vtkSmartPointer<vtkIdList>::New();
                vtkSmartPointer<vtkIdList> targetIdList = vtkSmartPointer<vtkIdList>::New();
                for ( int id : this->sourceids() )
                {
                    CHECK( id < iodesc.size() ) << "id : " << id << "not valid! must be < " << iodesc.size();
                    auto iodata = iodesc[id];
                    double point[3] = { iodata.nodeX(), iodata.nodeY(), iodata.nodeZ() };
                    vtkIdType vtkid = pointLocator->FindClosestPoint(point);
                    sourceIdList->InsertNextId(vtkid);
                }
                for ( int id : this->targetids() )
                {
                    CHECK( id < iodesc.size() ) << "id : " << id << "not valid! must be < " << iodesc.size();
                    auto iodata = iodesc[id];
                    double point[3] = { iodata.nodeX(), iodata.nodeY(), iodata.nodeZ() };
                    vtkIdType vtkid = pointLocator->FindClosestPoint(point);
                    targetIdList->InsertNextId(vtkid);
                }
                vecSourceTargetIdList.push_back( std::make_pair(sourceIdList,targetIdList) );
            }
            else
                CHECK(false) << "TODO get id from inlet/outlet cap";

            // load DelaunayTessellation if given
            if ( !rebuildDelaunayTessellation )
            {
                vtkSmartPointer<vtkUnstructuredGridReader> delaunayTessellationReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
                delaunayTessellationReader->SetFileName( pathDelaunayTessellation.c_str() );
                delaunayTessellationReader->Update();
                centerlineFilter->GenerateDelaunayTessellationOff();
                centerlineFilter->SetDelaunayTessellation( delaunayTessellationReader->GetOutput() );
            }

            int nCenterlinesComputed = vecSourceTargetIdList.size();
            std::vector<std::string> pathCenterlinesToFusion;
            for ( int k=0;k<nCenterlinesComputed;++k )
            {
                std::string outputPathUsed;// = this->outputPath();
                if ( nCenterlinesComputed > 1 )
                {
                    std::string myfilename = (boost::format("%1%_vmtkformat_part%2%.vtk")%name %k).str();
                    outputPathUsed = (directory / fs::path(myfilename) ).string();
                }
                else
                {
                    std::string myfilename = (boost::format("%1%_vmtkformat.vtk")%name).str();
                    outputPathUsed = (directory / fs::path(myfilename) ).string();
                }
                pathCenterlinesToFusion.push_back(outputPathUsed);
                if ( !this->forceRebuild() && fs::exists( outputPathUsed ) )
                    continue;

                if ( nCenterlinesComputed == 1 )
                {
                    if ( rebuildDelaunayTessellation )
                        std::cout << "Computing centerline with DelaunayTessellation ..." << std::endl
                                  << "Be patient because this step can take a long time to complete." << std::endl;
                    else
                        std::cout << "Computing centerline ..." << std::endl;
                }
                else
                {
                    if ( rebuildDelaunayTessellation )
                        std::cout << "Computing centerline " << k+1 << "/" << nCenterlinesComputed << " with DelaunayTessellation ..." << std::endl
                                  << "Be patient because this step can take a long time to complete." << std::endl;
                    else
                        std::cout << "Computing centerline " << k+1 << "/" << nCenterlinesComputed << " ..." << std::endl;
                }

                auto const& sourceTargetIdList = vecSourceTargetIdList[k];
                centerlineFilter->SetSourceSeedIds( sourceTargetIdList.first );
                centerlineFilter->SetTargetSeedIds( sourceTargetIdList.second );

                if ( rebuildDelaunayTessellation )
                    centerlineFilter->GenerateDelaunayTessellationOn();
                else
                    centerlineFilter->GenerateDelaunayTessellationOff();

                std::string radiusArrayName = "MaximumInscribedSphereRadius";
                centerlineFilter->SetRadiusArrayName( radiusArrayName.c_str() );
                centerlineFilter->SetCostFunction( M_costFunctionExpr.c_str() );
                centerlineFilter->SetFlipNormals( 0 );
                centerlineFilter->SetAppendEndPointsToCenterlines( 0 );
                centerlineFilter->SetSimplifyVoronoi( 0 );
                centerlineFilter->SetCenterlineResampling( 0 );
                centerlineFilter->SetResamplingStepLength( 1.0 );
                centerlineFilter->Update();

                vtkSmartPointer<vtkPolyDataWriter> centerlineWriter  = vtkSmartPointer<vtkPolyDataWriter>::New();
                centerlineWriter->SetInput( centerlineFilter->GetOutput() );
                centerlineWriter->SetFileName( outputPathUsed.c_str() );
                //centerlineWriter->SetFileTypeToBinary();
                centerlineWriter->SetFileTypeToASCII();
                centerlineWriter->Write();

                if ( rebuildDelaunayTessellation )
                {
                    vtkSmartPointer<vtkUnstructuredGridWriter> delaunayTessellationWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
                    delaunayTessellationWriter->SetInput( centerlineFilter->GetDelaunayTessellation() );
                    delaunayTessellationWriter->SetFileName( pathDelaunayTessellation.c_str() );
                    delaunayTessellationWriter->SetFileTypeToBinary();
                    delaunayTessellationWriter->Write();
                    rebuildDelaunayTessellation = false;
                }
            } // for (int k ... )


            if ( true )
            {
                std::shared_ptr<AngioTkCenterline> centerlinesTool;
                for ( std::string const& pathCenterline : pathCenterlinesToFusion )
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
                        centerlinesTool->importSurfaceFromFile( this->inputSurfacePath() );
                    }
                    centerlinesTool->importFile( pathCenterline );
                }
                centerlinesTool->addFieldBranchIds();
                centerlinesTool->writeCenterlinesVTK( this->outputPath() );
            }
        }
        else // use script python vmtkcenterlines
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
            //__str << " -delaunaytessellationfile " << "toto.vtk ";
#endif

            __str << " -ifile " << this->inputSurfacePath() << " ";
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
        std::cout << "The output file already exists in " << this->outputPath() << ", skipping step." << std::endl;
    }


    if ( M_viewResults )
    {
        std::ostringstream ostrView;
        ostrView << pythonExecutable << " ";// << dirBaseVmtk << "/vmtk "
        if ( M_viewResultsWithSurface )
        {
            ostrView << dirBaseVmtk << "/vmtksurfacereader -ifile " << this->inputSurfacePath() << " --pipe "
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

    return 0;
}


po::options_description
CenterlinesFromSurface::options( std::string const& prefix )
{
    po::options_description myCenterlinesOptions( "Centerlines from STL for blood flow mesh options" );
    myCenterlinesOptions.add_options()
        (prefixvm(prefix,"input.surface.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input filename" )
        (prefixvm(prefix,"input.pointset.filename").c_str(), po::value<std::string>()->default_value( "" ), "input.pointset.filename" )
        (prefixvm(prefix,"input.pointpair.filename").c_str(), po::value<std::string>()->default_value( "" ), "input.pointpair.filename" )
        (prefixvm(prefix,"input.desc.filename").c_str(), po::value<std::string>()->default_value( "" ), "inletoutlet-desc.filename" )
        (prefixvm(prefix,"input.geo-centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "geo-centerlines.filename" )

        (prefixvm(prefix,"use-interactive-selection").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) use-interactive-selection")
        (prefixvm(prefix,"cost-function.expression").c_str(), Feel::po::value<std::string>()->default_value("1/R"), "(string) cost-function")
        (prefixvm(prefix,"source-ids").c_str(), po::value<std::vector<int> >()->multitoken(), "(vector of int) source ids" )
        (prefixvm(prefix,"target-ids").c_str(), po::value<std::vector<int> >()->multitoken(), "(vector of int) target ids" )

        (prefixvm(prefix,"view-results").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) view-results")
        (prefixvm(prefix,"view-results.with-surface").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) view-results with surface")

        (prefixvm(prefix,"delaunay-tessellation.output.directory").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        (prefixvm(prefix,"delaunay-tessellation.force-rebuild").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) force-rebuild")
        ;
    return myCenterlinesOptions.add( super_type::options( prefix ) );
}


CenterlinesManager::CenterlinesManager( std::string const& prefix )
    :
    super_type( prefix ),
    M_inputCenterlinesPath(),// 1, soption(_name="input.centerlines.filename",_prefix=this->prefix()) ),
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.surface.filename",_prefix=this->prefix()) ) ),
    M_inputPointSetPath( AngioTkEnvironment::expand( soption(_name="input.point-set.filename",_prefix=this->prefix()) ) ),
    M_inputPointPairPath( AngioTkEnvironment::expand( soption(_name="input.point-pair.filename",_prefix=this->prefix()) ) ),
    M_applyThresholdMinRadius( doption(_name="threshold-radius.min",_prefix=this->prefix() ) ),
    M_applyThresholdMaxRadius( doption(_name="threshold-radius.max",_prefix=this->prefix() ) ),
    M_applyThresholdZonePointSetPath( AngioTkEnvironment::expand( soption(_name="thresholdzone-radius.point-set.filename",_prefix=this->prefix()) ) ),
    M_applyThresholdZoneMinRadius( doption(_name="thresholdzone-radius.min",_prefix=this->prefix() ) ),
    M_applyThresholdZoneMaxRadius( doption(_name="thresholdzone-radius.max",_prefix=this->prefix() ) ),
    M_avoidTubularColision( boption(_name="avoid-tubular-colision.apply",_prefix=this->prefix() ) ),
    M_avoidTubularColisionDistanceMin( doption(_name="avoid-tubular-colision.distance-min",_prefix=this->prefix() ) ),
    M_avoidTubularColisionRadiusMin( doption(_name="avoid-tubular-colision.radius-min",_prefix=this->prefix() ) ),
    M_avoidTubularColisionInputPointPairPath( AngioTkEnvironment::expand( soption(_name="avoid-tubular-colision.point-pair.filename",_prefix=this->prefix() ) ) ),
    M_smoothResample( boption(_name="smooth-resample.apply",_prefix=this->prefix() ) ),
    M_smoothResampleMeshSize( doption(_name="smooth-resample.mesh-size",_prefix=this->prefix() ) ),
    M_smoothResampleGeoPointSpacing( doption(_name="smooth-resample.geo-points-spacing",_prefix=this->prefix() ) )
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

    if ( !this->inputCenterlinesPath().empty() && !this->inputCenterlinesPath(0).empty() && this->outputPath().empty() )
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
    if ( this->outputDirectory().empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(this->outputDirectory()).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(this->outputDirectory());
    else
        meshesdirectories = fs::path(this->outputDirectory());

    // get filename without extension
    fs::path gp = this->inputCenterlinesPath(0);
    std::string nameMeshFile = gp.stem().string();

    std::string newFileName = (boost::format("%1%_up.vtk")%nameMeshFile ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    this->setOutputPath( outputPath.string() );
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


    fs::path directory = fs::path(this->outputPath()).parent_path();
    std::string filename = fs::path(this->outputPath()).stem().string();
    // build directories if necessary
    if ( !this->outputPath().empty() && this->worldComm().isMasterRank() )
    {
        //directory = fs::path(this->outputPath()).parent_path();
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
                    GmshInitialize();
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
        centerlinesTool->cleanBranch();

        centerlinesTool->addFieldBranchIds();
        if ( !centerlinesTool->hasField("RadiusMin") )
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

        if ( M_avoidTubularColision )
        {
            auto avoidTubularColisionParam = std::make_pair( M_avoidTubularColisionDistanceMin,M_avoidTubularColisionRadiusMin );
            if ( !M_avoidTubularColisionInputPointPairPath.empty() && fs::exists( M_avoidTubularColisionInputPointPairPath ) )
            {
                auto pointPairData = AngioTk::loadFromPointPairFile( M_avoidTubularColisionInputPointPairPath );
                centerlinesTool->applyTubularColisionFix( pointPairData, avoidTubularColisionParam );
            }
            else
                centerlinesTool->applyTubularColisionFix( avoidTubularColisionParam );
        }

        if ( M_smoothResample )
        {
            fs::path outputOriginalPath = directory / fs::path(filename+"_original.vtk");
            centerlinesTool->writeCenterlinesVTK( outputOriginalPath.string() );

            fs::path outputSmoothGeoPath = directory / fs::path(filename+".geo");
            AngioTkCenterline centerlinesSmooth;
            centerlinesSmooth.createFromCenterlines( *centerlinesTool, outputSmoothGeoPath.string(),
                                                     M_smoothResampleMeshSize, M_smoothResampleGeoPointSpacing );
            centerlinesSmooth.addFieldBranchIds();
            // redo TubularColision if asked
            if ( M_avoidTubularColision )
            {
                auto avoidTubularColisionParam = std::make_pair( M_avoidTubularColisionDistanceMin,M_avoidTubularColisionRadiusMin );
                if ( !M_avoidTubularColisionInputPointPairPath.empty() && fs::exists( M_avoidTubularColisionInputPointPairPath ) )
                {
                    auto pointPairData = AngioTk::loadFromPointPairFile( M_avoidTubularColisionInputPointPairPath );
                    centerlinesSmooth.applyTubularColisionFix( pointPairData, avoidTubularColisionParam );
                }
                else
                    centerlinesSmooth.applyTubularColisionFix( avoidTubularColisionParam );
            }

            centerlinesSmooth.writeCenterlinesVTK( this->outputPath() );
        }
        else
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
        (prefixvm(prefix,"input.point-pair.filename").c_str(), po::value<std::string>()->default_value( "" ), "input.pointpair.filename" )
        (prefixvm(prefix,"remove-branch-ids").c_str(), po::value<std::vector<int> >()->multitoken(), "(vector of int) remove branch ids" )
        (prefixvm(prefix,"use-window-interactor").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) use-window-interactor")
        (prefixvm(prefix,"window-width").c_str(), Feel::po::value<int>()->default_value(1024), "(int) window width")
        (prefixvm(prefix,"window-height").c_str(), Feel::po::value<int>()->default_value(768), "(int) window height")
        (prefixvm(prefix,"threshold-radius.min").c_str(), Feel::po::value<double>()->default_value(-1), "(double) threshold-radius.min")
        (prefixvm(prefix,"threshold-radius.max").c_str(), Feel::po::value<double>()->default_value(-1), "(double) threshold-radius.max")
        (prefixvm(prefix,"thresholdzone-radius.point-set.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input point-set filename" )
        (prefixvm(prefix,"thresholdzone-radius.max").c_str(), Feel::po::value<double>()->default_value(-1), "(double) threshold-radius.max")
        (prefixvm(prefix,"thresholdzone-radius.min").c_str(), Feel::po::value<double>()->default_value(-1), "(double) threshold-radius.max")
        (prefixvm(prefix,"avoid-tubular-colision.apply").c_str(), po::value<bool>()->default_value( false ), "(string) input point-set filename" )
        (prefixvm(prefix,"avoid-tubular-colision.distance-min").c_str(), po::value<double>()->default_value( 0.4 ), "(double) distance-min" )
        (prefixvm(prefix,"avoid-tubular-colision.radius-min").c_str(), po::value<double>()->default_value( 0.4 ), "(double) radius-min" )
        (prefixvm(prefix,"avoid-tubular-colision.point-pair.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input point-set filename" )
        (prefixvm(prefix,"smooth-resample.apply").c_str(), po::value<bool>()->default_value( false ), "(bool) smooth-resample.apply" )
        (prefixvm(prefix,"smooth-resample.mesh-size").c_str(), po::value<double>()->default_value( 1.0 ), "(double) mesh size" )
        (prefixvm(prefix,"smooth-resample.geo-points-spacing").c_str(), po::value<double>()->default_value( 4.0 ), "(double) geo-points-spacing" )
        ;
    return myCenterlinesManagerOptions.add( super_type::options( prefix ) );
}



ImageFromCenterlines::ImageFromCenterlines( std::string const& prefix )
    :
    super_type( prefix ),
    M_inputCenterlinesPath( AngioTkEnvironment::expand( soption(_name="input.centerlines.filename",_prefix=this->prefix()) ) ),
    M_dimX( ioption(_name="dim.x",_prefix=this->prefix() ) ),
    M_dimY( ioption(_name="dim.y",_prefix=this->prefix() ) ),
    M_dimZ( ioption(_name="dim.z",_prefix=this->prefix() ) ),
    M_dimSpacing( doption(_name="dim.spacing",_prefix=this->prefix() ) ),
    M_radiusArrayName( soption(_name="radius-array-name",_prefix=this->prefix() ) )
{
    if ( !this->inputCenterlinesPath().empty() && this->outputPath().empty() )
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
    if ( this->outputDirectory().empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(this->outputDirectory()).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(this->outputDirectory());
    else
        meshesdirectories = fs::path(this->outputDirectory());

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
    this->setOutputPath( outputPath.string() );
}

int
ImageFromCenterlines::run()
{
    if ( !fs::exists( this->inputCenterlinesPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : image building not done because this input path for centerlines does not exist :" << this->inputCenterlinesPath() << std::endl
                      << "Please set it with the \"input.centerlines.filename\" option." << std::endl;
        return 1;
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

    return 0;
}

po::options_description
ImageFromCenterlines::options( std::string const& prefix )
{
    po::options_description myImageFromCenterlinesOptions( "Create Image from Centerlines options" );

    myImageFromCenterlinesOptions.add_options()
        (prefixvm(prefix,"input.centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"dim.x").c_str(), Feel::po::value<int>()->default_value(0), "(int) dim.x")
        (prefixvm(prefix,"dim.y").c_str(), Feel::po::value<int>()->default_value(0), "(int) dim.y")
        (prefixvm(prefix,"dim.z").c_str(), Feel::po::value<int>()->default_value(0), "(int) dim.z")
        (prefixvm(prefix,"dim.spacing").c_str(), Feel::po::value<double>()->default_value(0.), "(double) dim.spacing")
        (prefixvm(prefix,"radius-array-name").c_str(), Feel::po::value<std::string>()->default_value("MaximumInscribedSphereRadius"), "(std::string) radius-array-name")
        ;
    return myImageFromCenterlinesOptions.add( super_type::options( prefix ) );
}

SurfaceFromImage::SurfaceFromImage( std::string const& prefix )
    :
    super_type( prefix ),
    M_inputImagesPath( /*AngioTkEnvironment::expand( soption(_name="input.filename",_prefix=this->prefix()) )*/ ),
    M_imageFusionOperator( soption(_name="image-fusion.operator",_prefix=this->prefix()) ),
    M_resizeFromRefImagePath( AngioTkEnvironment::expand( soption(_name="pre-process.resize-from-reference-image.path",_prefix=this->prefix()) ) ),
    M_method( soption(_name="method",_prefix=this->prefix()) ),
    M_hasThresholdLower(false), M_hasThresholdUpper(false),
    M_thresholdLower(0.0),M_thresholdUpper(0.0),
    M_applyConnectivityLargestRegion( boption(_name="apply-connectivity.largest-region",_prefix=this->prefix()) ),
    M_applyConnectivityNumberOfRegion( ioption(_name="apply-connectivity.number",_prefix=this->prefix()) )
{
    if ( Environment::vm().count(prefixvm(this->prefix(),"input.image.filename").c_str()) )
        M_inputImagesPath = Environment::vm()[prefixvm(this->prefix(),"input.image.filename").c_str()].as<std::vector<std::string> >();
    for ( int k=0;k<M_inputImagesPath.size();++k)
        M_inputImagesPath[k] = AngioTkEnvironment::expand( M_inputImagesPath[k] );

    CHECK( M_method == "threshold" || M_method == "isosurface" ) << "invalid method : " << M_method << " -> must be threshold or isosurface";

    CHECK( M_imageFusionOperator == "max" || M_imageFusionOperator == "min" || M_imageFusionOperator == "multiply" || M_imageFusionOperator == "subtract" )
        << "invalid fusion operator : " << M_imageFusionOperator << " -> must be max,min,multiply,subtract";

    for ( int k=0;k<M_inputImagesPath.size();++k)
        if ( !M_inputImagesPath[k].empty() && fs::path(M_inputImagesPath[k]).is_relative() )
            M_inputImagesPath[k] = (AngioTkEnvironment::pathInitial()/fs::path(M_inputImagesPath[k]) ).string();

    if ( !M_resizeFromRefImagePath.empty() && !fs::exists(M_resizeFromRefImagePath) )
    {
        //M_resizeFromRefImagePath.clear();
        std::cout << "SurfaceFromImage::Initialization" << std::endl
                  << "The pre-process.resize-from-reference-image.path option has been set, but the specified path does not exist ("
                  << M_resizeFromRefImagePath
                  << "). Please set it to a correct path or unset the option." << std::endl;
        exit(1);
    }

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

    if ( this->outputPath().empty() )
    {
        if ( !M_inputImagesPath.empty() && !M_inputImagesPath[0].empty() )
            this->updateOutputPathFromInputFileName();
    }
}

void
SurfaceFromImage::updateOutputPathFromInputFileName()
{
    CHECK( !M_inputImagesPath.empty() ) << "input path is empty";
    CHECK( !M_inputImagesPath[0].empty() ) << "input path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( this->outputDirectory().empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(this->outputDirectory()).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(this->outputDirectory());
    else
        meshesdirectories = fs::path(this->outputDirectory());

    // get filename without extension
    fs::path gp = M_inputImagesPath[0];
    std::string nameMeshFile = gp.stem().string();

    std::string newFileName = (boost::format("%1%.stl")%nameMeshFile ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    this->setOutputPath( outputPath.string() );
}

int
SurfaceFromImage::run()
{
    if ( this->inputImagesPath().empty() )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : surface segmentation will not be done, because you specified an empty input path." << std::endl
                      << "Please set it with the \"input.image.filename\" option." << std::endl;
        return 1;
    }

    if( !fs::exists( this->inputImagesPath(0) ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : surface segmentation will not be done because this input path not exist: " << this->inputImagesPath() << std::endl
                      << "Please set it with the \"input.image.filename\" option." << std::endl;
        return 1;
    }

    if ( !M_hasThresholdLower && !M_hasThresholdUpper )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : surface segmentation will not be done because no threshold has been given.\n" << std::endl
                      << "Please set the threshold values with the \"threshold.lower\" and \"threshold.upper\" options." << std::endl;
        return 1;
    }

    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run SurfaceFromImage \n"
            << "---------------------------------------\n";
    for ( int k=0;k<M_inputImagesPath.size();++k)
        coutStr << "inputImagesPath         : " << this->inputImagesPath(k) << "\n";
    if ( M_hasThresholdLower )
        coutStr << "threshold lower   : " << M_thresholdLower << "\n";
    if ( M_hasThresholdUpper )
        coutStr << "threshold upper   : " << M_thresholdUpper << "\n";
    coutStr << "output path       : " << this->outputPath() << "\n"
            << "---------------------------------------\n"
            << "---------------------------------------\n";
    std::cout << coutStr.str();


    fs::path directory = fs::path(this->outputPath()).parent_path();
    std::string outputNameWithoutExt = fs::path(this->outputPath()).stem().string();

    // build directories if necessary
    if ( !this->outputPath().empty() && this->worldComm().isMasterRank() )
    {
        if ( !fs::exists( directory ) )
            fs::create_directories( directory );
    }
    // // wait all process
    this->worldComm().globalComm().barrier();

    if ( !fs::exists( this->outputPath() ) || this->forceRebuild() )
    {
        std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
        std::string dirBaseVmtk = BOOST_PP_STRINGIZE( VMTK_BINARY_DIR );

        std::string nameImageInit = outputNameWithoutExt+"_levelsetInit.vti";
        std::string outputPathImageInit = (directory/fs::path(nameImageInit)).string();

        std::string inputImagePath = this->inputImagesPath(0);

        CHECK( M_inputImagesPath.size() <= 2 ) << "TODO image fusion with more than 2 images";
        if ( M_inputImagesPath.size() == 2 )
        {
            std::string nameImageFusion = outputNameWithoutExt+"_fusion.mha";
            std::string outputPathImageFusion = (directory/fs::path(nameImageFusion)).string();
            std::ostringstream ostrImageFusion;
            ostrImageFusion << pythonExecutable << " " << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtkimagecompose ";
            ostrImageFusion << "-ifile " << this->inputImagesPath(0) << " -i2file " << this->inputImagesPath(1) << " -operation " <<  M_imageFusionOperator << " ";
            ostrImageFusion << "-ofile " << outputPathImageFusion;
            // run vmtk script
            auto err = ::system( ostrImageFusion.str().c_str() );
            // use new image generated as input
            inputImagePath = outputPathImageFusion;
        }

        if ( !M_resizeFromRefImagePath.empty() )
        {
            std::string nameImageResized = outputNameWithoutExt+"_resized.mha";
            std::string outputPathImageResized = (directory/fs::path(nameImageResized)).string();
            ImagesManager myImagesManager;
            myImagesManager.setInputPath( inputImagePath );
            myImagesManager.setOutputPath( outputPathImageResized );
            myImagesManager.setResizeFromRefImageApply( true );
            myImagesManager.setForceRebuild( true );
            myImagesManager.setResizeFromRefImagePath( M_resizeFromRefImagePath );
            myImagesManager.run();
            inputImagePath = myImagesManager.outputPath();
        }

        std::ostringstream __str;
        __str << pythonExecutable << " ";
        __str << dirBaseVmtk << "/vmtk " << dirBaseVmtk << "/vmtkimageinitialization ";
        __str << "-ifile " << inputImagePath << " -interactive 0 ";
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
        __str2 << "-ifile " << inputImagePath << " -initiallevelsetsfile " << outputPathImageInit << " -iterations 1 ";
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
        __str2 << "-ifile " << inputImagePath << " ";
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
        readerInitialImage->SetFileName(inputImagePath/*this->inputImagesPath()*/.c_str());
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
            if ( M_applyConnectivityNumberOfRegion == 1 )
            {
                connectivityFilter->SetInput(/*surface*/cleanPolyData->GetOutput());
                connectivityFilter->SetExtractionModeToLargestRegion();
                connectivityFilter->Update();
            }
            else if ( M_applyConnectivityNumberOfRegion > 1 )
            {
                connectivityFilter->SetInput(/*surface*/cleanPolyData->GetOutput());
                connectivityFilter->SetExtractionModeToAllRegions();
                connectivityFilter->ColorRegionsOn();
                connectivityFilter->Update();
                int nRegion = connectivityFilter->GetNumberOfExtractedRegions();

                if ( connectivityFilter->GetOutput()->GetPointData()->HasArray("RegionId") == 0 )
                    CHECK( false ) << "RegionId must be present";

                std::vector<int> counterCellInRegion(nRegion,0);
                vtkIdType dataRegionIdSize = connectivityFilter->GetOutput()->GetPointData()->GetArray("RegionId")->GetSize();
                auto arrayRegionId = connectivityFilter->GetOutput()->GetPointData()->GetArray("RegionId");
                for (vtkIdType k=0;k<dataRegionIdSize;++k)
                {
                    double* val = arrayRegionId->GetTuple(k);
                    vtkIdType curRegionId = static_cast<vtkIdType>( val[0] );
                    CHECK( curRegionId >= 0 && curRegionId <= nRegion ) << "invalid curRegionId" << curRegionId;
                    ++counterCellInRegion[curRegionId];
                }

                // thanks to std::set with tuple, we can have automaticaly an ordering with respect to connectivity size
                std::set<std::tuple<int,int> > cellInRegionIdSet;
                for (int k=0;k<nRegion;++k )
                {
                    cellInRegionIdSet.insert( std::make_pair( counterCellInRegion[k],k ) );
                }

                std::set<int> regionIdExtracted;
                for ( auto rit=cellInRegionIdSet.rbegin(); rit != cellInRegionIdSet.rend(); ++rit)
                {
                    regionIdExtracted.insert(std::get<1>(*rit));
                    if ( regionIdExtracted.size() == M_applyConnectivityNumberOfRegion ) break;
                }

                connectivityFilter->SetExtractionModeToSpecifiedRegions();
                connectivityFilter->InitializeSpecifiedRegionList();
                for ( vtkIdType _regionId : regionIdExtracted )
                    connectivityFilter->AddSpecifiedRegion(_regionId);
                connectivityFilter->Update();

            }
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
    
    return 0;
}

po::options_description
SurfaceFromImage::options( std::string const& prefix )
{
    po::options_description mySurfaceFromImageOptions( "Create Surface from Image options" );

    mySurfaceFromImageOptions.add_options()
        (prefixvm(prefix,"input.image.filename").c_str(), po::value<std::vector<std::string> >()->multitoken(), "(vector of string) input image filename" )
        (prefixvm(prefix,"method").c_str(), Feel::po::value<std::string>()->default_value("threshold"), "(string) method : threshold, isosurface")
        (prefixvm(prefix,"image-fusion.operator").c_str(), Feel::po::value<std::string>()->default_value("max"), "(string) operator : max, min, multiply, subtract")
        (prefixvm(prefix,"pre-process.resize-from-reference-image.path").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) resize-from-reference-image")

        (prefixvm(prefix,"threshold.lower").c_str(), Feel::po::value<double>(), "(double) threshold lower")
        (prefixvm(prefix,"threshold.upper").c_str(), Feel::po::value<double>(), "(double) threshold upper")
        (prefixvm(prefix,"apply-connectivity.largest-region").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) apply-connectivity.largest-region")
        (prefixvm(prefix,"apply-connectivity.number").c_str(), Feel::po::value<int>()->default_value(1), "(bool) number of largest region")
        ;
    return mySurfaceFromImageOptions.add( super_type::options( prefix ) );
}


ImagesManager::ImagesManager()
    :
    super_type(),
    M_resizeFromRefImageApply( false )
{}

ImagesManager::ImagesManager( std::string const& prefix )
    :
    super_type( prefix ),
    M_inputPath(),
    M_resizeFromRefImageApply( boption(_name="resize-from-reference-image.apply",_prefix=this->prefix()) ),
    M_resizeFromRefImagePath( AngioTkEnvironment::expand( soption(_name="resize-from-reference-image.path",_prefix=this->prefix()) ) )
{
    if ( Environment::vm().count(prefixvm(this->prefix(),"input.image.path").c_str()) )
        M_inputPath = Environment::vm()[prefixvm(this->prefix(),"input.image.path").c_str()].as<std::vector<std::string> >();
    for ( int k=0;k<M_inputPath.size();++k)
        M_inputPath[k] = AngioTkEnvironment::expand( M_inputPath[k] );

    if ( this->outputPath().empty() && !this->inputPath().empty() && !this->inputPath(0).empty() )
        this->updateOutputPathFromInputFileName();
}

void
ImagesManager::updateOutputPathFromInputFileName()
{
    CHECK( !this->inputPath().empty() ) << "input path is empty";
    CHECK( !this->inputPath(0).empty() ) << "input path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( this->outputDirectory().empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(this->outputDirectory()).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(this->outputDirectory());
    else
        meshesdirectories = fs::path(this->outputDirectory());

    // get filename without extension
    fs::path gp = M_inputPath[0];
    std::string nameMeshFile = gp.stem().string();

    std::string newFileName = (boost::format("%1%_cvrt.mha")%nameMeshFile ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    this->setOutputPath( outputPath.string() );
}

void
ImagesManager::printInfo() const
{
    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run ImagesManager \n"
            << "---------------------------------------\n";
    for ( int k=0;k<this->inputPath().size();++k)
        coutStr << "inputPath         : " << this->inputPath(k) << "\n";
    coutStr << "output path       : " << this->outputPath() << "\n";
    if ( this->resizeFromRefImageApply() )
        coutStr << "input reference image path : " << this->resizeFromRefImagePath() << "\n";
    coutStr << "---------------------------------------\n"
            << "---------------------------------------\n";
    std::cout << coutStr.str();
}

void
ImagesManager::run()
{
    this->printInfo();

    if ( this->inputPath().empty() )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : ImagesManager do nothing because input path is empty\n";
        return;
    }
    if ( !fs::exists( this->inputPath(0) ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : ImagesManager do nothing done because this input path does not exist :" << this->inputPath(0) << "\n";
        return;
    }
    if ( this->resizeFromRefImageApply() && !fs::exists( this->resizeFromRefImagePath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : ImagesManager do nothing because resizeFromRefImage path does not exist :" << this->resizeFromRefImagePath() << "\n";
        return;
    }

    this->createOutputDirectory();

    if ( this->resizeFromRefImageApply() )
        this->updateResizeFromRefImage();
}

void
ImagesManager::updateResizeFromRefImage()
{
    std::string extInitialImage = fs::path(this->inputPath(0)).extension().string();
    std::string extRefImage = fs::path(this->resizeFromRefImagePath()).extension().string();
    CHECK( extInitialImage == ".mha" /*|| extInitialImage == ".nii"*/ ) << "image formart not support " << extInitialImage;
    CHECK( extRefImage == ".mha" /*|| extInitialImage == ".nii"*/ ) << "image formart not support " << extRefImage;

    vtkSmartPointer<vtkMetaImageReader> readerRefImage = vtkSmartPointer<vtkMetaImageReader>::New();
    readerRefImage->SetFileName(this->resizeFromRefImagePath().c_str());
    readerRefImage->Update();
    double originRefImage[3];
    readerRefImage->GetOutput()->GetOrigin(originRefImage);
    int dimRefImage[3];
    readerRefImage->GetOutput()->GetDimensions( dimRefImage );
    double spacingRefImage[3];
    readerRefImage->GetOutput()->GetSpacing( spacingRefImage );
    double lengthRefImage[3] = { spacingRefImage[0]*(dimRefImage[0]-1), spacingRefImage[1]*(dimRefImage[1]-1),spacingRefImage[2]*(dimRefImage[2]-1) };

    vtkSmartPointer<vtkMetaImageReader> readerInitialImage = vtkSmartPointer<vtkMetaImageReader>::New();
    readerInitialImage->SetFileName(this->inputPath(0).c_str());
    readerInitialImage->Update();
    double originInitialImage[3];
    readerInitialImage->GetOutput()->GetOrigin(originInitialImage);
    int dimInitialImage[3];
    readerInitialImage->GetOutput()->GetDimensions( dimInitialImage );
    double spacingInitialImage[3];
    readerInitialImage->GetOutput()->GetSpacing( spacingInitialImage );
    double lengthInitialImage[3] = { spacingInitialImage[0]*(dimInitialImage[0]-1), spacingInitialImage[1]*(dimInitialImage[1]-1),spacingInitialImage[2]*(dimInitialImage[2]-1) };

    if ( true )
    {
        std::cout << "originRefImage " << originRefImage[0] << "," << originRefImage[1] << "," << originRefImage[2] << "\n";
        std::cout << "dimRefImage " << dimRefImage[0] << "," << dimRefImage[1] << "," << dimRefImage[2] << "\n";
        std::cout << "lengthRefImage " << lengthRefImage[0] << "," << lengthRefImage[1] << "," << lengthRefImage[2] << "\n";
        std::cout << "spacingRefImage " << spacingRefImage[0] << "," << spacingRefImage[1] << "," << spacingRefImage[2] << "\n";
        std::cout << "originInitialImage " << originInitialImage[0] << "," << originInitialImage[1] << "," << originInitialImage[2] << "\n";
        std::cout << "dimInitialImage " << dimInitialImage[0] << "," << dimInitialImage[1] << "," << dimInitialImage[2] << "\n";
        std::cout << "spacingInitialImage " << spacingInitialImage[0] << "," << spacingInitialImage[1] << "," << spacingInitialImage[2] << "\n";
        std::cout << "lengthInitialImage " << lengthInitialImage[0] << "," << lengthInitialImage[1] << "," << lengthInitialImage[2] << "\n";
    }

    double spacingNewImage[3] = { lengthRefImage[0]/(dimInitialImage[0]-1), lengthRefImage[1]/(dimInitialImage[1]-1), lengthRefImage[2]/(dimInitialImage[2]-1) };

    vtkSmartPointer<vtkImageChangeInformation> imageChangeInfo = vtkSmartPointer<vtkImageChangeInformation>::New();
    imageChangeInfo->SetInput(readerInitialImage->GetOutput());
    imageChangeInfo->SetOutputOrigin(originRefImage);
    imageChangeInfo->SetOutputSpacing(spacingNewImage);
    imageChangeInfo->Update();

    vtkSmartPointer<vtkMetaImageWriter> metaImageWriter = vtkSmartPointer<vtkMetaImageWriter>::New();
    metaImageWriter->SetFileName(this->outputPath().c_str());
    metaImageWriter->SetInput(imageChangeInfo->GetOutput());

    metaImageWriter->Write();

}

po::options_description
ImagesManager::options( std::string const& prefix )
{
    po::options_description myImagesManagerOptions( "Modified an Image" );

    myImagesManagerOptions.add_options()
        (prefixvm(prefix,"input.image.path").c_str(), po::value<std::vector<std::string> >()->multitoken(), "(vector of string) input image filename" )
        (prefixvm(prefix,"resize-from-reference-image.path").c_str(), Feel::po::value<std::string>()->default_value(""), "(string) reference-image path")
        (prefixvm(prefix,"resize-from-reference-image.apply").c_str(), Feel::po::value<bool>()->default_value(false), "(bool) resize-from-reference-image")
        ;
    return myImagesManagerOptions.add( super_type::options( prefix ) );
}



SubdivideSurface::SubdivideSurface( std::string const& prefix )
    :
    super_type( prefix ),
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.surface.filename",_prefix=this->prefix()) ) ),
    M_method( soption(_name="method",_prefix=this->prefix()) ),
    M_nSubdivisions( ioption(_name="subdivisions",_prefix=this->prefix()) )
{
    CHECK( M_method == "linear" || M_method == "butterfly" || M_method == "loop" ) << "invalid method " << M_method << "\n";
    if ( !this->inputSurfacePath().empty() && this->outputPath().empty() )
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
    if ( this->outputDirectory().empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(this->outputDirectory()).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(this->outputDirectory());
    else
        meshesdirectories = fs::path(this->outputDirectory());

    // get filename without extension
    //fs::path gp = M_inputPath;
    std::string nameMeshFile = fs::path(this->inputSurfacePath()).stem().string();

    std::string newFileName = (boost::format("%1%_subdivide%2%%3%.stl")%nameMeshFile %M_nSubdivisions %M_method ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    this->setOutputPath( outputPath.string() );
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
        (prefixvm(prefix,"input.surface.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"method").c_str(), Feel::po::value<std::string>()->default_value("butterfly"), "(string) linear, butterfly, loop")
        (prefixvm(prefix,"subdivisions").c_str(), Feel::po::value<int>()->default_value(1), "(int) number of subdivisions")
        ;
    return mySubdivideSurfaceOptions.add( super_type::options( prefix ) );
}

SmoothSurface::SmoothSurface( std::string const& prefix )
    :
    super_type( prefix ),
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.surface.filename",_prefix=this->prefix()) ) ),
    M_inputCenterlinesPath( AngioTkEnvironment::expand( soption(_name="input.centerlines.filename",_prefix=this->prefix()) ) ),
    M_method( soption(_name="method",_prefix=this->prefix()) ),
    M_nIterations( ioption(_name="iterations",_prefix=this->prefix()) ),
    M_taubinPassBand( doption(_name="taubin.passband",_prefix=this->prefix()) ),
    M_laplaceRelaxationFactor( doption(_name="laplace.relaxation",_prefix=this->prefix()) )
{
    CHECK( M_method == "taubin" || M_method == "laplace" || M_method == "centerlines" ) << "invalid method " << M_method << "\n";
    if ( !this->inputSurfacePath().empty() && this->outputPath().empty() )
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
    if ( this->outputDirectory().empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(this->outputDirectory()).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(this->outputDirectory());
    else
        meshesdirectories = fs::path(this->outputDirectory());

    // get filename without extension
    fs::path gp = this->inputSurfacePath();
    std::string nameMeshFile = gp.stem().string();

    std::string newFileName;
    if ( M_method == "taubin")
        newFileName = (boost::format("%1%_smooth%2%Taubin%3%.stl")%nameMeshFile %M_nIterations %M_taubinPassBand ).str();
    else
        newFileName = (boost::format("%1%_smooth%2%Laplace%3%.stl")%nameMeshFile %M_nIterations %M_laplaceRelaxationFactor ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    this->setOutputPath( outputPath.string() );
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
        (prefixvm(prefix,"input.surface.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"input.centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        (prefixvm(prefix,"method").c_str(), Feel::po::value<std::string>()->default_value("taubin"), "(string) taubin or laplace")
        (prefixvm(prefix,"iterations").c_str(), Feel::po::value<int>()->default_value(30), "(int) number of iterations")
        (prefixvm(prefix,"taubin.passband").c_str(), Feel::po::value<double>()->default_value(0.1), "(double) taubin passband")
        (prefixvm(prefix,"laplace.relaxation").c_str(), Feel::po::value<double>()->default_value(0.01), "(double) laplace.relaxation")
        ;
    return mySmoothSurfaceOptions.add( super_type::options( prefix ) );
}



OpenSurface::OpenSurface( std::string const& prefix )
    :
    super_type( prefix ),
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.surface.filename",_prefix=this->prefix()) ) ),
    M_inputCenterlinesPath( AngioTkEnvironment::expand( soption(_name="input.centerlines.filename",_prefix=this->prefix()) ) ),
    M_distanceClipScalingFactor( doption(_name="distance-clip.scaling-factor",_prefix=this->prefix() ) ),
    M_radiusUncertainty( doption(_name="radius-uncertainty",_prefix=this->prefix()) ),
    M_saveOutputSurfaceBinary( boption(_name="output.save-binary",_prefix=this->prefix() ) )
{
    if ( !M_inputSurfacePath.empty() && this->outputPath().empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}


void
OpenSurface::updateOutputPathFromInputFileName()
{
    CHECK( !M_inputSurfacePath.empty() ) << "input path is empty";
    if(M_inputSurfacePath.empty())
    {
        std::cout << "WARNING : opening surface step will not be done, because this input surface path for centerlines does not exist: " << this->inputSurfacePath() << std::endl
                  << "Please set it with the input.surface.filename option" << std::endl; 
        exit(1);
    }

    // define output directory
    fs::path meshesdirectories;
    if ( this->outputDirectory().empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(this->outputDirectory()).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(this->outputDirectory());
    else
        meshesdirectories = fs::path(this->outputDirectory());

    // get filename without extension
    fs::path gp = M_inputSurfacePath;
    std::string nameMeshFile = gp.stem().string();

    std::string newFileName = (boost::format("%1%_open.stl")%nameMeshFile ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    this->setOutputPath( outputPath.string() );
}

int
OpenSurface::run()
{
    if ( !fs::exists( this->inputSurfacePath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : opening surface step will not be done, because this input surface path for centerlines does not exist: " << this->inputSurfacePath() << std::endl
                      << "Please set it with the \"input.surface.filename\" option" << std::endl; 
        return 1;
    }
    if ( !fs::exists( this->inputCenterlinesPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : opening surface not done because this input centerlines path for centerlines not exist :" << this->inputCenterlinesPath() << std::endl
                      << "Please set it with the \"input.centerlines.filename\" option" << std::endl; 
        return 1;
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

    return 0;
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
    centerlinesTool->setCutMeshRadiusUncertainty( M_radiusUncertainty );
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
        (prefixvm(prefix,"distance-clip.scaling-factor").c_str(), Feel::po::value<double>()->default_value(0.0), "(double) scaling-factor")
        (prefixvm(prefix,"radius-uncertainty").c_str(), Feel::po::value<double>()->default_value(0.), "(double) radius-uncertainty")
        (prefixvm(prefix,"output.save-binary").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) save-binary")
        ;
    return myOpenSurfaceOptions.add( super_type::options( prefix ) );
}




RemeshSurface::RemeshSurface( std::string const& prefix )
    :
    super_type( prefix ),
    M_packageType(soption(_name="package-type",_prefix=this->prefix())),
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.surface.filename",_prefix=this->prefix())) ),
    M_inputCenterlinesPath( AngioTkEnvironment::expand(soption(_name="gmsh.centerlines.filename",_prefix=this->prefix())) ),
    M_gmshRemeshNbPointsInCircle( ioption(_name="gmsh.nb-points-in-circle",_prefix=this->prefix()) ),
    M_gmshRemeshRadiusUncertainty( doption(_name="gmsh.radius-uncertainty",_prefix=this->prefix()) ),
    M_gmshRemeshPartitionForceRebuild( boption(_name="gmsh.remesh-partition.force-rebuild",_prefix=this->prefix()) ),
    M_vmtkArea( doption(_name="vmtk.area",_prefix=this->prefix()) ),
    M_vmtkNumberOfIteration( ioption(_name="vmtk.n-iteration",_prefix=this->prefix()) ),
    M_saveOutputSurfaceBinary( boption(_name="output.save-binary",_prefix=this->prefix() ) )
{
    CHECK( M_packageType == "gmsh" || M_packageType == "gmsh-executable" || M_packageType == "vmtk" ) << "error on packageType : " << M_packageType;
    if ( !this->inputSurfacePath().empty() && this->outputPath().empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}


void
RemeshSurface::updateOutputPathFromInputFileName()
{
    CHECK( !this->inputSurfacePath().empty() ) << "input surface path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( this->outputDirectory().empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(this->outputDirectory()).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(this->outputDirectory());
    else
        meshesdirectories = fs::path(this->outputDirectory());

    // get filename without extension
    fs::path gp = this->inputSurfacePath();
    std::string nameMeshFile = gp.stem().string();

    std::string specStr;
    if ( this->packageType() =="gmsh" || this->packageType() == "gmsh-executable" )
        specStr = (boost::format( "remeshGMSHpt%1%") %this->remeshNbPointsInCircle() ).str();
    else if ( this->packageType() == "vmtk" )
        specStr = (boost::format( "remeshVMTKarea%1%") %this->area() ).str();

    std::string newFileName = (boost::format("%1%_%2%.stl")%nameMeshFile %specStr ).str();
    this->setOutputPath( ( meshesdirectories / fs::path(newFileName) ).string() );
}

int
RemeshSurface::run()
{
    if ( !fs::exists( this->inputSurfacePath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : remesh surface not done because this input surface path for centerlines not exist :" << this->inputSurfacePath() << std::endl
                      << "Please set it with the \"input.surface.filename\" option." << std::endl;
        return 1;
    }
    if ( ( this->packageType() == "gmsh" || this->packageType() == "gmsh-executable" ) && !fs::exists( this->inputCenterlinesPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : remesh surface not done because this input centerlines path for centerlines not exist :" << this->inputCenterlinesPath() << std::endl
                      << "Please set it with the \"gmsh.centerlines.filename\" option." << std::endl;
        return 1;
    }

    std::cout << "\n"
              << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
              << "---------------------------------------\n"
              << "run RemeshSurface \n"
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
        std::cout << "The output file already exists in " << this->outputPath() << ", skipping remeshSTL." << std::endl
                  << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
        return 0;
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

    return 0;
}

void
RemeshSurface::runVMTK()
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
RemeshSurface::runGMSHwithExecutable()
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
RemeshSurface::runGMSH()
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
    centerlinesTool->setCutMeshRadiusUncertainty( this->gmshRemeshRadiusUncertainty() );
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
RemeshSurface::options( std::string const& prefix )
{
    po::options_description myMeshSurfaceOptions( "Mesh surface blood flow from STL and Centerlines options" );
    myMeshSurfaceOptions.add_options()

        ( prefixvm(prefix,"package-type").c_str(), po::value<std::string>()->default_value( "vmtk" ), "force-remesh" )
        ( prefixvm(prefix,"input.surface.filename").c_str(), po::value<std::string>()->default_value( "" ), "surface filename" )
        ( prefixvm(prefix,"output.save-binary").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) save-binary")
        // gmsh options
        ( prefixvm(prefix,"gmsh.centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "centerlines.filename" )
        ( prefixvm(prefix,"gmsh.nb-points-in-circle").c_str(), po::value<int>()->default_value( 15 ), "nb-points-in-circle" )
        ( prefixvm(prefix,"gmsh.remesh-partition.force-rebuild").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) force-rebuild")
        ( prefixvm(prefix,"gmsh.radius-uncertainty").c_str(), Feel::po::value<double>()->default_value(0.), "(double) radius-uncertainty")
        // vmtk options
        ( prefixvm(prefix,"vmtk.area").c_str(), po::value<double>()->default_value( 0.5 ), "area" )
        ( prefixvm(prefix,"vmtk.n-iteration").c_str(), po::value<int>()->default_value( 10 ), "maxit" )
        ;
    return myMeshSurfaceOptions.add( super_type::options( prefix ) );
}


TubularExtension::TubularExtension( std::string const& prefix )
    :
    super_type( prefix ),
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.surface.filename",_prefix=this->prefix()) ) ),
    M_inputCenterlinesPath( AngioTkEnvironment::expand( soption(_name="input.centerlines.filename",_prefix=this->prefix()) ) ),
    M_saveOutputSurfaceBinary( boption(_name="output.save-binary",_prefix=this->prefix() ) )
{
    if ( !M_inputSurfacePath.empty() && this->outputPath().empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}

void
TubularExtension::updateOutputPathFromInputFileName()
{
    CHECK( !M_inputSurfacePath.empty() ) << "input path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( this->outputDirectory().empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(this->outputDirectory()).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(this->outputDirectory());
    else
        meshesdirectories = fs::path(this->outputDirectory());

    // get filename without extension
    fs::path gp = M_inputSurfacePath;
    std::string nameMeshFile = gp.stem().string();

    std::string newFileName = (boost::format("%1%_tubeext.stl")%nameMeshFile ).str();
    fs::path outputPath = meshesdirectories / fs::path(newFileName);
    this->setOutputPath( outputPath.string() );
}

void
TubularExtension::run()
{

    if ( !fs::exists( this->inputSurfacePath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : opening surface not done because this input surface path for centerlines not exist :" << this->inputSurfacePath() << "\n";
        return;
    }
#if 0
    if ( !fs::exists( this->inputCenterlinesPath() ) )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "WARNING : opening surface not done because this input centerlines path for centerlines not exist :" << this->inputCenterlinesPath() << "\n";
        return;
    }
#endif
    std::ostringstream coutStr;
    coutStr << "\n"
            << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
            << "---------------------------------------\n"
            << "run TubularExtension \n"
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
    // // wait all process
    this->worldComm().globalComm().barrier();


    if ( true )
    {
        GmshInitialize();
        CTX::instance()->terminal = 1;
        int verbosityLevel = 5;
        Msg::SetVerbosity( verbosityLevel );
        CTX::instance()->geom.tolerance=1e-8;
        CTX::instance()->mesh.algo2d = ALGO_2D_FRONTAL;//ALGO_2D_DELAUNAY;//ALGO_2D_FRONTAL //ALGO_2D_BAMG
    }
    GModel * gmodel = new GModel();
    // add new model as current (important if this function is called more than 1 time)
    GModel::current(GModel::list.size() - 1);

    std::shared_ptr<AngioTkCenterline> centerlinesTool( new AngioTkCenterline );
    centerlinesTool->importSurfaceFromFile( this->inputSurfacePath() );
    centerlinesTool->importFile( this->inputCenterlinesPath() );
    centerlinesTool->runTubularExtension();
    centerlinesTool->saveTubularExtensionSTL( this->outputPath(), M_saveOutputSurfaceBinary );

    delete gmodel;
}

po::options_description
TubularExtension::options( std::string const& prefix )
{
    po::options_description myTubularExtensionOptions( "Extension of inlet/outlet options" );
    myTubularExtensionOptions.add_options()
        ( prefixvm(prefix,"input.centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        ( prefixvm(prefix,"input.surface.filename").c_str(), po::value<std::string>()->default_value( "" ), "(string) input centerline filename" )
        ( prefixvm(prefix,"output.save-binary").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) save-binary")
        ;
    return myTubularExtensionOptions.add( super_type::options( prefix ) );
}

namespace detail
{
void
generateMeshFromGeo( std::string const& inputGeoName,std::string const& outputMeshName,int dim )
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



VolumeMeshing::VolumeMeshing( std::string const& prefix )
    :
    super_type( prefix ),
    M_inputSurfacePath( AngioTkEnvironment::expand( soption(_name="input.surface.filename",_prefix=this->prefix()) ) ),
    M_inputCenterlinesPath( AngioTkEnvironment::expand( soption(_name="input.centerlines.filename",_prefix=this->prefix())) ),
    M_inputInletOutletDescPath( AngioTkEnvironment::expand( soption(_name="input.desc.filename",_prefix=this->prefix())) ),
    M_remeshNbPointsInCircle( ioption(_name="nb-points-in-circle",_prefix=this->prefix() ) ),
    M_extrudeWall( boption(_name="extrude-wall",_prefix=this->prefix() ) ),
    M_extrudeWallNbElemLayer( ioption(_name="extrude-wall.nb-elt-layer",_prefix=this->prefix() ) ),
    M_extrudeWallhLayer( doption(_name="extrude-wall.h-layer",_prefix=this->prefix()) ),
    M_saveOutputVolumeBinary( boption(_name="output.save-binary",_prefix=this->prefix() ) )
{
    if ( !this->inputSurfacePath().empty() && this->outputPath().empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}


void
VolumeMeshing::updateOutputPathFromInputFileName()
{
    CHECK( !this->inputSurfacePath().empty() ) << "input path is empty";

    // define output directory
    fs::path meshesdirectories;
    if ( this->outputDirectory().empty() )
        meshesdirectories = fs::current_path();
    else if ( fs::path(this->outputDirectory()).is_relative() )
        meshesdirectories = fs::path(Environment::rootRepository())/fs::path(this->outputDirectory());
    else
        meshesdirectories = fs::path(this->outputDirectory());

    // get filename without extension
    fs::path gp = this->inputSurfacePath();
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
    this->setOutputPath( outputPath.string() );
}

void
VolumeMeshing::run()
{
    std::cout << "\n"
              << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
              << "---------------------------------------\n"
              << "run VolumeMeshing \n"
              << "---------------------------------------\n";
    std::cout << "input surface path      : " << this->inputSurfacePath() << "\n"
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

    CHECK( !this->inputSurfacePath().empty() ) << "inputSurfacePath is empty";
    CHECK( !this->inputCenterlinesPath().empty() ) << "centerlinesFileName is empty";
    CHECK( !this->inputInletOutletDescPath().empty() ) << "inletOutletDescFileName is empty";

    if ( fs::exists( this->outputPath() ) && !this->forceRebuild() )
    {
        std::cout << "The output file already exists in " << this->outputPath() << ", skipping meshVolume" << std::endl
                  << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
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
VolumeMeshing::generateGeoFor3dVolumeFromSTLAndCenterlines(std::string const& geoname)
{
    std::ostringstream geodesc;
    geodesc << "General.ExpertMode=1;\n";
    geodesc << "Mesh.Algorithm = 6; //(1=MeshAdapt, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) \n"
            << "Mesh.Algorithm3D = 1; //(1=tetgen, 4=netgen, 7=MMG3D, 9=R-tree) \n"
            << "Mesh.Binary = " << M_saveOutputVolumeBinary <<";\n";

    //CTX::instance()->mesh.binary = M_saveOutputSurfaceBinary;

    //Mesh.CharacteristicLengthFactor=0.015;

    //Merge "stl_remesh_gmsh/fluidskin_P15.stl";
    geodesc << "Merge \""<< this->inputSurfacePath() <<"\";\n";

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
        ( prefixvm(prefix,"input.surface.filename").c_str(), po::value<std::string>()->default_value( "" ), "stl.filename" )
        ( prefixvm(prefix,"input.centerlines.filename").c_str(), po::value<std::string>()->default_value( "" ), "centerlines.filename" )
        ( prefixvm(prefix,"input.desc.filename").c_str(), po::value<std::string>()->default_value( "" ), "inletoutlet-desc.filename" )
        ( prefixvm(prefix,"output.save-binary").c_str(), Feel::po::value<bool>()->default_value(true), "(bool) save-binary")
        ( prefixvm(prefix,"nb-points-in-circle").c_str(), po::value<int>()->default_value( 15 ), "nb-points-in-circle" )
        ( prefixvm(prefix,"extrude-wall").c_str(),po::value<bool>()->default_value( false ), "extrude-wall" )
        ( prefixvm(prefix,"extrude-wall.nb-elt-layer").c_str(), po::value<int>()->default_value( 2 ), "nb-elt-layer" )
        ( prefixvm(prefix,"extrude-wall.h-layer").c_str(), po::value<double>()->default_value( 0.2 ), "h-layer" )
        ;
    return myMeshVolumeOptions.add( super_type::options( prefix ) );
}





} // namespace Feel
