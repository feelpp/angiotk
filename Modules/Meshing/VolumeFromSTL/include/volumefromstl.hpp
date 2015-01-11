/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#ifndef __VOLUMEFROMSTL_H
#define __VOLUMEFROMSTL_H 1

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
namespace Feel
{

class InletOutletData : public boost::tuple<std::string,std::string,std::vector<double> >
{
    typedef boost::tuple<std::string,std::string,std::vector<double> > super_type;
public :
    InletOutletData( std::string markerLumen,std::string markerArterialWall,double x, double y, double z)
        :
        super_type(markerLumen,markerArterialWall,{x,y,z} )
    {}
    InletOutletData( InletOutletData const& e )
        :
        super_type( e )
    {}
public :
    std::string markerLumen() const { return this->get<0>(); }
    std::string markerArterialWall() const { return this->get<1>(); }
    std::vector<double> const& node() const { return this->get<2>(); }
    double nodeX() const { return this->node()[0]; }
    double nodeY() const { return this->node()[1]; }
    double nodeZ() const { return this->node()[2]; }

};

class InletOutletDesc : public std::vector<InletOutletData>
{
    typedef std::vector<InletOutletData> super_type;
public :

    InletOutletDesc( std::string const& path );

    InletOutletDesc( InletOutletDesc const& e )
        :
        super_type( e ),
        M_path( e.M_path )
    {}

    void add( InletOutletData const& data );

 private :
    std::string M_path;

};

class CenterlinesFromSTL
{
public :

    CenterlinesFromSTL( std::string prefix );
    CenterlinesFromSTL( CenterlinesFromSTL const& e );

    std::string prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }

    std::string inputPath() const { return M_inputPath; }
    std::string inputInletOutletDescPath() const { return M_inputInletOutletDescPath; }
    std::string outputPath() const { return M_outputPath; }
    std::set<int> const& targetids() const { return M_targetids; }
    std::set<int> const& sourceids() const { return M_sourceids; }
    bool forceRebuild() const { return M_forceRebuild; }

    void setStlFileName(std::string s) { M_inputPath=s; }
    void setTargetids( std::set<int> const& s ) { M_targetids=s; }
    void setSourceids( std::set<int> const& s ) { M_sourceids=s; }
    void setForceRebuild( bool b ) { M_forceRebuild=b; }

    void updateOutputPathFromInputFileName();

    void run();

    static po::options_description options( std::string const& prefix );

private :
    std::string M_prefix;
    std::string M_inputPath, M_inputInletOutletDescPath, M_outputPath;
    std::string M_outputDirectory;
    std::set<int> M_targetids, M_sourceids;
    bool M_forceRebuild;
    bool M_viewResults,M_viewResultsWithSurface;
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
    std::string M_centerlinesFileName;
    int M_remeshNbPointsInCircle;
    double M_area;
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
