/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#ifndef __VOLUMEFROMSTL_H
#define __VOLUMEFROMSTL_H 1

#include <feel/feelcore/environment.hpp>

#include <angiotkMeshingConfig.h>

#include <feel/feelfilters/gmsh.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#if defined( FEELPP_HAS_GMSH_H )

#include <GmshConfig.h>
#include <Gmsh.h>
#include <GModel.h>
#include <OpenFile.h>
#include <GmshDefines.h>
#include <Context.h>
#endif

/*#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"*/

class AngioTkEnvironment : boost::noncopyable
{
public :
    template <class ArgumentPack>
    AngioTkEnvironment( ArgumentPack const& args )
    {
        S_pathInitial = Feel::fs::current_path();
        S_feelEnvironment.reset(new Feel::Environment( Feel::_argc=args[Feel::_argc], Feel::_argv=args[Feel::_argv],
                                                       Feel::_desc=args[Feel::_desc|Feel::feel_nooptions()],
                                                       Feel::_about=args[Feel::_about| Feel::makeAboutDefault/*detail::makeAbout*/( args[Feel::_argv][0] )] ) );
    }
    ~AngioTkEnvironment()
    {
        //std::cout << "use_count() " << S_feelEnvironment.use_count()  <<"\n";
        S_feelEnvironment.reset();
    }
    static Feel::fs::path pathInitial() { return S_pathInitial; }
    static Feel::Environment const& feelEnvironment() { return *S_feelEnvironment; }

    BOOST_PARAMETER_CONSTRUCTOR(
    AngioTkEnvironment, ( AngioTkEnvironment ), Feel::tag,
        ( required
          ( argc,* )
          ( argv,* ) )
        ( optional
          ( desc,* )
          //( desc_lib,* )
          ( about,* )
          //( directory,( std::string ) )
        ) ) // no semicolon
    //{}
    static std::string expand( std::string const& expr ) { return Feel::Environment::expand( expr ); }

    static Feel::po::variables_map const& vm() { return Feel::Environment::vm(); }

private :
    static Feel::fs::path S_pathInitial;
    static boost::shared_ptr<Feel::Environment> S_feelEnvironment;
};



namespace AngioTk
{
//namespace fs = boost::filesystem;
//namespace fs = Feel::fs;
//namespace po = Feel::po;
using namespace Feel;

class AngioTkFilterBase
{
public :
    AngioTkFilterBase()
        :
        M_forceRebuild( true )
    {}
    AngioTkFilterBase( std::string const& prefix );
    AngioTkFilterBase( AngioTkFilterBase const& e ) = default;
    Feel::WorldComm const& worldComm() const { return Feel::Environment::worldComm(); }
    std::string const& prefix() const { return M_prefix; }
    std::string const& outputDirectory() const { return M_outputDirectory; }
    std::string const& outputPath() const { return M_outputPath; }
    bool forceRebuild() const { return M_forceRebuild; }
    void setOutputPath( std::string const& path ) { M_outputPath=path; this->updateOutputDirFromOutputPath(); }
    void setOutputDirectory( std::string const& path ) { M_outputDirectory=path; }
    void setForceRebuild( bool b ) { M_forceRebuild = b; }

    void createOutputDirectory()
    {
        // build directories if necessary
        if ( !this->outputDirectory().empty() && this->worldComm().isMasterRank() )
            {
                if ( !fs::exists( this->outputDirectory() ) )
                    fs::create_directories( this->outputDirectory() );
            }
        // // wait all process
        this->worldComm().globalComm().barrier();
    }
    static po::options_description options( std::string const& prefix );

private :
    void updateOutputDirFromOutputPath()
    {
        if ( this->outputPath().empty() ) return;
        if ( fs::path(this->outputPath()).is_absolute() )
            M_outputDirectory = fs::path(this->outputPath()).parent_path().string();
    }
    std::string M_prefix;
    std::string M_outputDirectory, M_outputPath;
    bool M_forceRebuild;
};

enum BCType { BC_INLET=0,BC_OUTLET=1 };

class InletOutletData : public boost::tuple<BCType,std::string,std::string,std::vector<double> >
{
    typedef boost::tuple<BCType,std::string,std::string,std::vector<double> > super_type;
public :
    InletOutletData( BCType bctype,std::string markerLumen,std::string markerArterialWall,double x, double y, double z)
        :
        super_type(bctype,markerLumen,markerArterialWall,{x,y,z} )
    {}
    InletOutletData( InletOutletData const& e )
        :
        super_type( e )
    {}
public :
    BCType bcType() const { return this->get<0>(); }
    std::string bcTypeStr() const { return ( this->bcType() == BCType::BC_INLET ) ? "INLET" : "OUTLET"; }
    std::string markerLumen() const { return this->get<1>(); }
    std::string markerArterialWall() const { return this->get<2>(); }
    std::vector<double> const& node() const { return this->get<3>(); }
    double nodeX() const { return this->node()[0]; }
    double nodeY() const { return this->node()[1]; }
    double nodeZ() const { return this->node()[2]; }

};

class InletOutletDesc : public std::vector<InletOutletData>
{
    typedef std::vector<InletOutletData> super_type;
public :
    InletOutletDesc()
        :
        super_type()
    {}
    InletOutletDesc( std::string const& path );

    InletOutletDesc( InletOutletDesc const& e )
        :
        super_type( e )
        //M_path( e.M_path )
    {}

    void add( InletOutletData const& data );

    void loadFromSTL( std::string inputPath );
    void save( std::string outputPath );
    void saveJSON( std::string outputPath );

 private :
//std::string M_path;

};

class CenterlinesFromSurface : public AngioTkFilterBase
{
    typedef AngioTkFilterBase super_type;
public :

    CenterlinesFromSurface( std::string const& prefix );
    CenterlinesFromSurface( CenterlinesFromSurface const& e ) = default;

    std::string inputSurfacePath() const { return M_inputSurfacePath; }
    std::string inputCenterlinesPointSetPath() const { return M_inputCenterlinesPointSetPath; }
    std::string inputCenterlinesPointPairPath() const { return M_inputCenterlinesPointPairPath; }
    std::string inputInletOutletDescPath() const { return M_inputInletOutletDescPath; }
    std::string inputGeoCenterlinesPath() const { return M_inputGeoCenterlinesPath; }

    std::set<int> const& targetids() const { return M_targetids; }
    std::set<int> const& sourceids() const { return M_sourceids; }
    bool viewResults() const { return M_viewResults; }

    void setInputSurfacePath(std::string s) { M_inputSurfacePath=s; }
    void setTargetids( std::set<int> const& s ) { M_targetids=s; }
    void setSourceids( std::set<int> const& s ) { M_sourceids=s; }
    void setViewResults(bool b) { M_viewResults=b; }

    void updateOutputPathFromInputFileName();

    void run();

    static po::options_description options( std::string const& prefix );
private :
    //std::tuple< std::vector<double>, std::vector<double> >
    std::tuple< std::vector<std::vector<double> >, std::vector<std::vector<double> > >
    loadFromCenterlinesPointSetFile();

    std::vector< std::pair< std::vector<double>,std::vector<double> > >
    loadFromCenterlinesPointPairFile();

private :
    std::string M_inputSurfacePath, M_inputCenterlinesPointSetPath, M_inputCenterlinesPointPairPath, M_inputInletOutletDescPath;
    std::string M_inputGeoCenterlinesPath;
    std::set<int> M_targetids, M_sourceids;
    std::string M_costFunctionExpr;

    bool M_useInteractiveSelection;
    bool M_viewResults,M_viewResultsWithSurface;

    std::string M_delaunayTessellationOutputDirectory;
    bool M_delaunayTessellationForceRebuild;
};

class CenterlinesManager : public AngioTkFilterBase
{
    typedef AngioTkFilterBase super_type;
public :

    CenterlinesManager( std::string const& prefix );
    CenterlinesManager( CenterlinesManager const& e ) = default;

    void updateOutputPathFromInputFileName();

    std::map<int,std::vector<std::tuple<double,double,double> > >
    loadPointSetFile( std::string const& filepath );


    void run();

    static po::options_description options( std::string const& prefix );

    std::vector<std::string> const& inputCenterlinesPath() const { return M_inputCenterlinesPath; }
    std::string const& inputCenterlinesPath(int k) const { return M_inputCenterlinesPath[k]; }
    std::string const& inputSurfacePath() const { return M_inputSurfacePath; }
    std::string const& inputPointSetPath() const { return M_inputPointSetPath; }
    std::string const& inputPointPairPath() const { return M_inputPointPairPath; }

    void setInputCenterlinesPath(std::string const& path) { M_inputCenterlinesPath = { path }; }
    void setInputSurfacePath(std::string const& path) { M_inputSurfacePath=path; }
    void setInputPointSetPath(std::string const& path) { M_inputPointSetPath=path; }
    void setInputPointPairPath(std::string const& path) { M_inputPointPairPath=path; }

private :
    std::vector<std::string> M_inputCenterlinesPath;
    std::string M_inputSurfacePath, M_inputPointSetPath, M_inputPointPairPath;
    std::set<int> M_removeBranchIds;
    double M_applyThresholdMinRadius,M_applyThresholdMaxRadius;
    std::string M_applyThresholdZonePointSetPath;
    double M_applyThresholdZoneMinRadius,M_applyThresholdZoneMaxRadius;
    bool M_avoidTubularColision;
    double M_avoidTubularColisionDistanceMin, M_avoidTubularColisionRadiusMin;
    std::string M_avoidTubularColisionInputPointPairPath;
    bool M_smoothResample;
    double M_smoothResampleMeshSize, M_smoothResampleGeoPointSpacing;
};

class ImageFromCenterlines : public AngioTkFilterBase
{
    typedef AngioTkFilterBase super_type;
public :

    ImageFromCenterlines( std::string const& prefix );
    ImageFromCenterlines( ImageFromCenterlines const& e ) = default;

    void updateOutputPathFromInputFileName();

    void run();

    static po::options_description options( std::string const& prefix );

    std::string inputCenterlinesPath() const { return M_inputCenterlinesPath; }

private :
    std::string M_inputCenterlinesPath;
    int M_dimX,M_dimY,M_dimZ;
    double M_dimSpacing;
    std::string M_radiusArrayName;
};

class SurfaceFromImage : public AngioTkFilterBase
{
    typedef AngioTkFilterBase super_type;
public :

    SurfaceFromImage( std::string const& prefix );
    SurfaceFromImage( SurfaceFromImage const& e ) = default;

    void updateOutputPathFromInputFileName();

    int run();

    static po::options_description options( std::string const& prefix );

    std::vector<std::string> const& inputImagesPath() const { return M_inputImagesPath; }
    std::string const& inputImagesPath( int k ) const { return M_inputImagesPath[k]; }
    std::string const& method() const { return M_method; }

private :
    std::vector<std::string> M_inputImagesPath;
    std::string M_imageFusionOperator;
    std::string M_resizeFromRefImagePath;
    std::string M_method;
    double M_thresholdLower,M_thresholdUpper;
    bool M_hasThresholdLower,M_hasThresholdUpper;
    bool M_applyConnectivityLargestRegion;
    int M_applyConnectivityNumberOfRegion;
};


class ImagesManager : public AngioTkFilterBase
{
    typedef AngioTkFilterBase super_type;
public :

    ImagesManager();
    ImagesManager( std::string const& prefix );
    ImagesManager( ImagesManager const& e ) = default;
    void updateOutputPathFromInputFileName();
    void run();
    void printInfo() const;
    static po::options_description options( std::string const& prefix );

    std::vector<std::string> const& inputPath() const { return M_inputPath; }
    std::string const& inputPath(int k) const { return M_inputPath[k]; }
    bool resizeFromRefImageApply() const { return M_resizeFromRefImageApply; }
    std::string const& resizeFromRefImagePath() const { return M_resizeFromRefImagePath; }

    void setInputPath( std::string const& path ) { M_inputPath.clear(); M_inputPath.push_back( path ); }
    void setResizeFromRefImageApply( bool b ) { M_resizeFromRefImageApply = b; }
    void setResizeFromRefImagePath( std::string const& path ) { M_resizeFromRefImagePath = path; }

private :
    void updateResizeFromRefImage();
private :
    std::vector<std::string> M_inputPath;
    bool M_resizeFromRefImageApply;
    std::string M_resizeFromRefImagePath;
};


class SubdivideSurface : public AngioTkFilterBase
{
    typedef AngioTkFilterBase super_type;
public :

    SubdivideSurface( std::string const& prefix );
    SubdivideSurface( SubdivideSurface const& e ) = default;

    void updateOutputPathFromInputFileName();

    void run();

    std::string inputSurfacePath() const { return M_inputSurfacePath; }

    void setInputSurfacePath(std::string const& path) { M_inputSurfacePath=path; }

    static po::options_description options( std::string const& prefix );

private :
    std::string M_inputSurfacePath;
    std::string M_method;
    int M_nSubdivisions;
};

class SmoothSurface : public AngioTkFilterBase
{
    typedef AngioTkFilterBase super_type;
public :

    SmoothSurface( std::string const& prefix );
    SmoothSurface( SmoothSurface const& e ) = default;

    void updateOutputPathFromInputFileName();

    void run();

    std::string inputSurfacePath() const { return M_inputSurfacePath; }

    void setInputSurfacePath(std::string const& path) { M_inputSurfacePath=path; }

    static po::options_description options( std::string const& prefix );

private :
    std::string M_inputSurfacePath, M_inputCenterlinesPath;
    std::string M_method;
    int M_nIterations;
    double M_taubinPassBand;
    double M_laplaceRelaxationFactor;
};

class OpenSurface : public AngioTkFilterBase
{
    typedef AngioTkFilterBase super_type;
public :

    OpenSurface( std::string const& prefix );
    OpenSurface( OpenSurface const& e ) = default;

    void updateOutputPathFromInputFileName();

    void run();
    void runGMSH();
    void runGMSHwithExecutable();
    void runVMTK();

    static po::options_description options( std::string const& prefix );

    std::string inputSurfacePath() const { return M_inputSurfacePath; }
    std::string inputCenterlinesPath() const { return M_inputCenterlinesPath; }

    void setInputSurfacePath(std::string const& s) { M_inputSurfacePath=s; }
    void setInputCenterlinesPath(std::string const& s) { M_inputCenterlinesPath=s; }

private :
    std::string M_inputSurfacePath, M_inputCenterlinesPath;
    double M_distanceClipScalingFactor, M_radiusUncertainty;
    bool M_saveOutputSurfaceBinary;
};

namespace detail
{
void generateMeshFromGeo( std::string const& inputGeoName,std::string const& outputMeshName,int dim);
}

class RemeshSurface : public AngioTkFilterBase
{
    typedef AngioTkFilterBase super_type;
public :

    RemeshSurface( std::string const& prefix );
    RemeshSurface( RemeshSurface const& e ) = default;

    std::string packageType() const { return M_packageType; }
    std::string inputSurfacePath() const { return M_inputSurfacePath; }
    std::string inputCenterlinesPath() const { return M_inputCenterlinesPath; }

    // gmsh
    int remeshNbPointsInCircle() const { return M_gmshRemeshNbPointsInCircle; }
    double gmshRemeshRadiusUncertainty() const { return M_gmshRemeshRadiusUncertainty; }
    // vmtk
    double area() const { return M_vmtkArea; }

    void setPackageType( std::string type)
    {
        CHECK( type == "gmsh" || type == "gmsh-executable" || type == "vmtk" ) << "error on packageType : " << type;
        M_packageType=type;
    }
    void setInputSurfacePath(std::string const& s) { M_inputSurfacePath=s; }
    void setInputCenterlinesPath(std::string const& s) { M_inputCenterlinesPath=s; }

    void updateOutputPathFromInputFileName();

    void run();
    void runVMTK();
    void runGMSH();
    void runGMSHwithExecutable();

    static po::options_description options( std::string const& prefix );

private :
    std::string M_packageType;
    std::string M_inputSurfacePath;
    // with gmsh
    std::string M_inputCenterlinesPath;
    int M_gmshRemeshNbPointsInCircle;
    bool M_gmshRemeshPartitionForceRebuild;
    double M_gmshRemeshRadiusUncertainty;
    // with vmtk
    double M_vmtkArea;
    int M_vmtkNumberOfIteration;

    bool M_saveOutputSurfaceBinary;
};

class TubularExtension : public AngioTkFilterBase
{
    typedef AngioTkFilterBase super_type;
public :

    TubularExtension( std::string const& prefix );
    TubularExtension( TubularExtension const& e ) = default;

    void updateOutputPathFromInputFileName();
    void run();

    std::string inputSurfacePath() const { return M_inputSurfacePath; }
    std::string inputCenterlinesPath() const { return M_inputCenterlinesPath; }

    void setInputSurfacePath(std::string const& path) { M_inputSurfacePath=path; }
    void setInputCenterlinesPath(std::string const& s) { M_inputCenterlinesPath=s; }

    static po::options_description options( std::string const& prefix );

private :
    std::string M_inputSurfacePath, M_inputCenterlinesPath;
    bool M_saveOutputSurfaceBinary;
};

class VolumeMeshing : public AngioTkFilterBase
{
    typedef AngioTkFilterBase super_type;
public :

    VolumeMeshing( std::string const& prefix );
    VolumeMeshing( VolumeMeshing const& e ) = default;

    std::string const& inputSurfacePath() const { return M_inputSurfacePath; }
    std::string const& inputCenterlinesPath() const { return M_inputCenterlinesPath; }
    std::string const& inputInletOutletDescPath() const { return M_inputInletOutletDescPath; }
    int remeshNbPointsInCircle() const { return M_remeshNbPointsInCircle; }
    bool extrudeWall() const { return M_extrudeWall; }
    int extrudeWall_nbElemLayer() const { return M_extrudeWallNbElemLayer; }
    double extrudeWall_hLayer() const { return M_extrudeWallhLayer; }

    void setInputSurfacePath(std::string const& s) { M_inputSurfacePath=s; }
    void setInputCenterlinesPath(std::string const& s) { M_inputCenterlinesPath=s; }
    void setInputInletOutletDescPath(std::string const& s) { M_inputInletOutletDescPath=s; }

    void updateOutputPathFromInputFileName();

    void run();
    void generateGeoFor3dVolumeFromSTLAndCenterlines(std::string const& geoname);

    static po::options_description options( std::string const& prefix );


private :

    std::string M_inputSurfacePath;
    std::string M_inputCenterlinesPath;

    std::string M_inputInletOutletDescPath;

    int M_remeshNbPointsInCircle;

    bool M_extrudeWall;
    int M_extrudeWallNbElemLayer;
    double M_extrudeWallhLayer;

    bool M_saveOutputVolumeBinary;
};

} // namespace AngioTk

#endif // __VOLUMEFROMSTL_H
