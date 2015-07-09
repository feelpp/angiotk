/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#ifndef __VOLUMEFROMSTL_H
#define __VOLUMEFROMSTL_H 1

#include <feel/feelcore/environment.hpp>

#include <angiotkMeshingConfig.h>

#include <feel/feelfilters/gmsh.hpp>
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


namespace detail
{
class AngioTkEnvironment
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
    static Feel::fs::path pathInitial() { return S_pathInitial; }
    static Feel::Environment const& feelEnvironment() { return *S_feelEnvironment; }

private :
    static Feel::fs::path S_pathInitial;
    static boost::shared_ptr<Feel::Environment> S_feelEnvironment;
};
} // namespace detail

class AngioTkEnvironment : public detail::AngioTkEnvironment
{
public:
    BOOST_PARAMETER_CONSTRUCTOR(
    AngioTkEnvironment, ( detail::AngioTkEnvironment ), Feel::tag,
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
};
namespace Feel
{

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

 private :
//std::string M_path;

};

class CenterlinesFromSTL
{
public :

    CenterlinesFromSTL( std::string prefix );
    CenterlinesFromSTL( CenterlinesFromSTL const& e ) = default;

    std::string prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }

    std::string inputPath() const { return M_inputPath; }
    std::string inputCenterlinesPointSetPath() const { return M_inputCenterlinesPointSetPath; }
    std::string inputInletOutletDescPath() const { return M_inputInletOutletDescPath; }
    std::string outputPath() const { return M_outputPath; }
    void setOutputPath(std::string const& path) { M_outputPath=path; }
    std::set<int> const& targetids() const { return M_targetids; }
    std::set<int> const& sourceids() const { return M_sourceids; }
    bool forceRebuild() const { return M_forceRebuild; }
    bool viewResults() const { return M_viewResults; }

    void setStlFileName(std::string s) { M_inputPath=s; }
    void setTargetids( std::set<int> const& s ) { M_targetids=s; }
    void setSourceids( std::set<int> const& s ) { M_sourceids=s; }
    void setForceRebuild( bool b ) { M_forceRebuild=b; }
    void setViewResults(bool b) { M_viewResults=b; }

    void updateOutputPathFromInputFileName();

    void run();

    static po::options_description options( std::string const& prefix );
private :
    //std::tuple< std::vector<double>, std::vector<double> >
    std::tuple< std::vector<std::vector<double> >, std::vector<std::vector<double> > >
    loadFromCenterlinesPointSetFile();

private :
    std::string M_prefix;
    std::string M_inputPath, M_inputCenterlinesPointSetPath, M_inputInletOutletDescPath, M_outputPath;
    std::string M_outputDirectory;
    std::set<int> M_targetids, M_sourceids;
    bool M_forceRebuild;
    bool M_useInteractiveSelection;
    bool M_viewResults,M_viewResultsWithSurface;
};

class CenterlinesManager
{
public :

    CenterlinesManager( std::string prefix );
    CenterlinesManager( CenterlinesManager const& e ) = default;

    void updateOutputPathFromInputFileName();

    void run();

    static po::options_description options( std::string const& prefix );

    std::string const& prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }
    std::vector<std::string> const& inputCenterlinesPath() const { return M_inputCenterlinesPath; }
    std::string const& inputCenterlinesPath(int k) const { return M_inputCenterlinesPath[k]; }
    std::string const& inputSurfacePath() const { return M_inputSurfacePath; }
    std::string const& inputPointSetPath() const { return M_inputPointSetPath; }
    std::string const& outputPath() const { return M_outputPath; }
    bool forceRebuild() const { return M_forceRebuild; }

    void setInputCenterlinesPath(std::string const& path) { M_inputCenterlinesPath = { path }; }
    void setInputSurfacePath(std::string const& path) { M_inputSurfacePath=path; }
    void setInputPointSetPath(std::string const& path) { M_inputPointSetPath=path; }
    void setOutputPath(std::string const& path) { M_outputPath=path; }
    void setOutputDirectory(std::string const& path) { M_outputDirectory=path; }

private :
    std::string M_prefix;
    std::vector<std::string> M_inputCenterlinesPath;
    std::string M_inputSurfacePath, M_inputPointSetPath, M_outputDirectory, M_outputPath;
    bool M_forceRebuild;
    bool M_useWindowInteractor;
    std::set<int> M_removeBranchIds;
};

class ImageFromCenterlines
{
public :

    ImageFromCenterlines( std::string prefix );
    ImageFromCenterlines( ImageFromCenterlines const& e );

    void updateOutputPathFromInputFileName();

    void run();

    static po::options_description options( std::string const& prefix );

    std::string prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }
    std::string inputPath() const { return M_inputPath; }
    std::string outputPath() const { return M_outputPath; }
    bool forceRebuild() const { return M_forceRebuild; }

private :
    std::string M_prefix;
    std::string M_inputPath, M_outputDirectory,M_outputPath;
    double M_dimX,M_dimY,M_dimZ;
    bool M_forceRebuild;
};

class SurfaceFromImage
{
public :

    SurfaceFromImage( std::string prefix );
    SurfaceFromImage( SurfaceFromImage const& e );

    void updateOutputPathFromInputFileName();

    void run();

    static po::options_description options( std::string const& prefix );

    std::string prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }
    std::string inputPath() const { return M_inputPath; }
    std::string outputPath() const { return M_outputPath; }
    bool forceRebuild() const { return M_forceRebuild; }

private :
    std::string M_prefix;
    std::string M_inputPath, M_outputDirectory,M_outputPath;
    double M_thresholdLower,M_thresholdUpper;
    bool M_hasThresholdLower,M_hasThresholdUpper;
    bool M_forceRebuild;
};

class SubdivideSurface
{
public :

    SubdivideSurface( std::string prefix );
    SubdivideSurface( SubdivideSurface const& e );

    void updateOutputPathFromInputFileName();

    void run();

    std::string prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }
    std::string inputPath() const { return M_inputPath; }
    void setInputPath(std::string const& path) {M_inputPath=path; }
    std::string outputPath() const { return M_outputPath; }
    void setOutputPath(std::string const& path) { M_outputPath=path; }
    void setOutputDirectory(std::string const& path) { M_outputDirectory=path; }

    bool forceRebuild() const { return M_forceRebuild; }
    static po::options_description options( std::string const& prefix );

private :
    std::string M_prefix;
    std::string M_inputPath, M_outputDirectory, M_outputPath;
    std::string M_method;
    int M_nSubdivisions;
    bool M_forceRebuild;
};

class SmoothSurface
{
public :

    SmoothSurface( std::string prefix );
    SmoothSurface( SmoothSurface const& e );

    void updateOutputPathFromInputFileName();

    void run();

    std::string prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }
    std::string inputPath() const { return M_inputPath; }
    std::string outputPath() const { return M_outputPath; }
    bool forceRebuild() const { return M_forceRebuild; }

    static po::options_description options( std::string const& prefix );

private :
    std::string M_prefix;
    std::string M_inputPath, M_outputDirectory, M_outputPath;
    std::string M_method;
    int M_nIterations;
    double M_taubinPassBand;
    double M_laplaceRelaxationFactor;
    bool M_forceRebuild;
};

class OpenSurface
{
public :

    OpenSurface( std::string prefix );
    OpenSurface( OpenSurface const& e );

    void updateOutputPathFromInputFileName();

    void run();
    void runGMSH();

    static po::options_description options( std::string const& prefix );

    std::string prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }
    std::string inputSurfacePath() const { return M_inputSurfacePath; }
    std::string inputCenterlinesPath() const { return M_inputCenterlinesPath; }
    std::string outputPath() const { return M_outputPath; }
    bool forceRebuild() const { return M_forceRebuild; }

private :
    std::string M_prefix;
    std::string M_inputSurfacePath, M_inputCenterlinesPath, M_outputDirectory, M_outputPath;
    bool M_forceRebuild;
};

namespace detail
{
void generateMeshFromGeo( std::string inputGeoName,std::string outputMeshName,int dim);
}

class RemeshSTL
{
public :

    RemeshSTL( std::string prefix );
    RemeshSTL( RemeshSTL const& e );

    std::string prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }
    std::string packageType() const { return M_packageType; }
    std::string inputPath() const { return M_inputPath; }
    std::string centerlinesFileName() const { return M_centerlinesFileName; }
    int remeshNbPointsInCircle() const { return M_remeshNbPointsInCircle; }
    double area() const { return M_area; }
    std::string outputPath() const { if ( this->packageType() =="gmsh" ) return M_outputPathGMSH; else return M_outputPathVMTK; }

    void setPackageType( std::string type)
    {
        CHECK( type == "gmsh" || type == "vmtk" ) << "error on packageType : " << type;
        M_packageType=type;
    }
    void setInputPath(std::string s) { M_inputPath=s; }
    void setCenterlinesFileName(std::string s) { M_centerlinesFileName=s; }

    bool forceRebuild() const { return M_forceRebuild; }
    void setForceRebuild( bool b ) { M_forceRebuild=b; }

    void updateOutputPathFromInputFileName();

    void run();
    void runVMTK();
    void runGMSH();

    static po::options_description options( std::string const& prefix );

private :
    std::string M_prefix;
    std::string M_packageType;
    std::string M_inputPath;
    // with gmsh
    std::string M_centerlinesFileName;
    int M_remeshNbPointsInCircle;
    // with vmtk
    double M_area;
    int M_nIterationVMTK;

    std::string M_outputPathGMSH, M_outputPathVMTK;
    std::string M_outputDirectory;
    bool M_forceRebuild;

}; // class RemeshSTL

class VolumeMeshing
{
public :

    VolumeMeshing( std::string prefix );
    VolumeMeshing( VolumeMeshing const& e );

    std::string prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }

    std::string inputSTLPath() const { return M_inputSTLPath; }
    std::string inputCenterlinesPath() const { return M_inputCenterlinesPath; }
    std::string inputInletOutletDescPath() const { return M_inputInletOutletDescPath; }
    std::string outputPath() const { return M_outputPath; }
    bool forceRebuild() const { return M_forceRebuild; }
    int remeshNbPointsInCircle() const { return M_remeshNbPointsInCircle; }
    bool extrudeWall() const { return M_extrudeWall; }
    int extrudeWall_nbElemLayer() const { return M_extrudeWallNbElemLayer; }
    double extrudeWall_hLayer() const { return M_extrudeWallhLayer; }

    void setInputSTLPath(std::string s) { M_inputSTLPath=s; }
    void setInputCenterlinesPath(std::string s) { M_inputCenterlinesPath=s; }
    void setInputInletOutletDescPath(std::string s) { M_inputInletOutletDescPath=s; }
    void setForceRebuild( bool b ) { M_forceRebuild=b; }

    void updateOutputPathFromInputFileName();

    void run();
    void generateGeoFor3dVolumeFromSTLAndCenterlines(std::string geoname);

    static po::options_description options( std::string const& prefix );


private :

    std::string M_prefix;
    std::string M_inputSTLPath;
    std::string M_inputCenterlinesPath;

    std::string M_inputInletOutletDescPath;

    int M_remeshNbPointsInCircle;

    bool M_extrudeWall;
    int M_extrudeWallNbElemLayer;
    double M_extrudeWallhLayer;

    std::string M_outputPath;
    std::string M_outputDirectory;
    bool M_forceRebuild;

};

} // namespace Feel

#endif // __VOLUMEFROMSTL_H
