/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#ifndef __toolboxbloodflowmesh_H
#define __toolboxbloodflowmesh_H 1

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

//#define MY_PYTHON_EXECUTABLE /opt/local/bin/python2.7

namespace Feel
{

class ToolBoxCenterlines
{
public :
    ToolBoxCenterlines()
        :
        M_stlFileName( soption(_name="centerlines.stl.filename") )
        {
            this->updateOutputFileName();
        }

    ToolBoxCenterlines( ToolBoxCenterlines const& e )
        :
        M_stlFileName( e.M_stlFileName ),
        M_outputFileName( e.M_outputFileName ),
        M_targetids( e.M_targetids ),
        M_sourceids( e.M_sourceids )
        {}

    std::string stlFileName() const { return M_stlFileName; }
    //std::string setFileNameSTL(std::string s) { return M_stlFileName=s;this->updateOutputFileName(); }
    void setStlFileName(std::string s) { M_stlFileName=s; this->updateOutputFileName(); }

    std::string outputFileName() const { return M_outputFileName; }

    std::set<int> const& targetids() const { return M_targetids; }
    std::set<int> const& sourceids() const { return M_sourceids; }
    void setTargetids( std::set<int> const& s ) { M_targetids=s; }
    void setSourceids( std::set<int> const& s ) { M_sourceids=s; }

    void updateOutputFileName()
    {
        std::string outputName = fs::path(this->stlFileName()).stem().string()+"_centerlines.vtk";
        //std::string outputName = fs::path(this->stlFileName()).string()+"_centerlines.vtk";
        M_outputFileName = (fs::current_path()/outputName).string();
    }

    void run()
    {

        std::ostringstream coutStr;

        coutStr << "\n"
                  << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
                  << "---------------------------------------\n"
                  << "run ToolBoxCenterlines \n"
                  << "---------------------------------------\n";
        coutStr << "stlFileName          : " << this->stlFileName() << "\n";
        coutStr << "targetids : ";
        for ( int id : this->targetids() )
            coutStr << id << " ";
        coutStr << "\n";
        coutStr << "sourceids : ";
        for ( int id : this->sourceids() )
            coutStr << id << " ";
        coutStr << "\n"
                << "outputFileName       : " << this->outputFileName() << "\n"
                << "---------------------------------------\n"
                << "---------------------------------------\n";
        std::cout << coutStr.str();

        if ( fs::exists( this->outputFileName() ) )
        {
            std::cout << "already file exist, ignore centerline\n"
                      << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
            return;
        }

        //fs::path stlNamePath = fs::path(this->stlFileName());
        fs::path outputFileNamePath = fs::path(this->outputFileName());
        std::string name = outputFileNamePath.stem().string();

        std::ostringstream __str;

        // source ~/packages/vmtk/vmtk.build2/Install/vmtk_env.sh
        //std::string pythonExecutable="/opt/local/bin/python2.7";
        std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );

        __str << pythonExecutable << " ";
        //std::string dirBaseVmtk = "/Users/vincentchabannes/packages/vmtk/new/vmtk.build/Install/bin/";
        std::string dirBaseVmtk = "/Users/vincentchabannes/packages/vmtk/vmtk.build2/Install/bin/";
        __str << dirBaseVmtk << "vmtk " << dirBaseVmtk << "vmtkcenterlines "
              <<"-seedselector profileidlist ";
        __str << "-targetids ";
        for ( int id : this->targetids() )
            __str << id << " ";
        __str << "-sourceids ";
        for ( int id : this->sourceids() )
            __str << id << " ";

        __str << "-ifile " << this->stlFileName() << " ";
        __str << "-ofile " << name << ".vtp "
              << "--pipe " << dirBaseVmtk << "vmtksurfacewriter "
              << "-ifile " << name << ".vtp "
              << "-ofile " << name << ".vtk "
              << " -mode ascii ";

        std::cout << "---------------------------------------\n"
                  << "run in system : \n" << __str.str() << "\n"
                  << "---------------------------------------\n";
        auto err = ::system( __str.str().c_str() );

        std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n";

    }

    static
    po::options_description
    options()
    {
        po::options_description myCenterlinesOptions( "Centerlines from STL for blood flow mesh options" );
        myCenterlinesOptions.add_options()
            ( "centerlines.stl.filename", po::value<std::string>()->default_value( "" ), "stl.filename" )
            ;
        return myCenterlinesOptions;
    }

private :
    std::string M_stlFileName, M_outputFileName;
    std::set<int> M_targetids, M_sourceids;

}; // class ToolBoxCenterlines

namespace detail
{
    void generateMeshFromGeo( std::string inputGeoName,std::string outputMeshName,int dim/*,std::string formatExportMesh = "msh"*/)
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


}

class ToolBoxBloodFlowReMeshSTL
{
public :

    ToolBoxBloodFlowReMeshSTL( std::string type = "gmsh" )
        :
        M_packageType( type ),
        M_stlFileName(soption(_name="mesh-surface.stl.filename")),
        M_forceRemesh( boption(_name="mesh-surface.force-remesh") ),
        M_centerlinesFileName(soption(_name="mesh-surface.centerlines.filename")),
        M_remeshNbPointsInCircle( ioption(_name="mesh-surface.nb-points-in-circle") ),
        M_area( doption(_name="mesh-surface.area") ),
        M_outputFileName()
        {
            CHECK( M_packageType == "gmsh" || M_packageType == "vmtk" ) << "error on packageType : " << M_packageType;
            this->updateOutputFileName();
        }

    ToolBoxBloodFlowReMeshSTL( ToolBoxBloodFlowReMeshSTL const& e )
        :
        M_packageType( e.M_packageType ),
        M_stlFileName( e.M_stlFileName ),
        M_forceRemesh( e.M_forceRemesh ),
        M_centerlinesFileName( e.M_centerlinesFileName ),
        M_remeshNbPointsInCircle( e.M_remeshNbPointsInCircle ),
        M_area( e.M_area ),
        M_outputFileName( e.M_outputFileName )
        {}


    std::string packageType() const { return M_packageType; }
    std::string stlFileName() const { return M_stlFileName; }
    bool forceRemesh() const { return M_forceRemesh; }
    std::string centerlinesFileName() const { return M_centerlinesFileName; }
    int remeshNbPointsInCircle() const { return M_remeshNbPointsInCircle; }
    double area() const { return M_area; }
    std::string outputFileName() const { return M_outputFileName; }

    void setStlFileName(std::string s) { M_stlFileName=s; this->updateOutputFileName(); }
    void setCenterlinesFileName(std::string s) { M_centerlinesFileName=s; }


    void updateOutputFileName()
    {
        std::string specStr;
        if ( this->packageType() == "gmsh" )
            specStr += (boost::format( "remeshGMSHpt%1%") %this->remeshNbPointsInCircle() ).str();
        else if ( this->packageType() == "vmtk" )
            specStr += (boost::format( "remeshVMTKarea%1%") %this->area() ).str();

        std::string outputName = (boost::format( "%1%_%2%.stl")
                                  %fs::path(this->stlFileName()).stem().string()
                                  % specStr
                                  ).str();


        //std::string outputName = fs::path(this->stlFileName()).stem().string()+"_remeshGMSH.stl";
        M_outputFileName = (fs::current_path()/outputName).string();
    }

    void run()
    {
        std::cout << "\n"
                  << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
                  << "---------------------------------------\n"
                  << "run ToolBoxBloodFlowReMeshSTL \n"
                  << "---------------------------------------\n";
        std::cout << "type                 : " << M_packageType << "\n"
                  << "stlFileName          : " << this->stlFileName() << "\n";
        if ( this->packageType() == "gmsh" )
            std::cout << "centerlinesFileName  : " << this->centerlinesFileName() << "\n";
        else if ( this->packageType() == "vmtk" )
            std::cout << "area  : " << this->area() << "\n";
        std::cout << "outputFileName       : " << this->outputFileName() << "\n"
                  << "---------------------------------------\n"
                  << "---------------------------------------\n";

        if ( !this->forceRemesh() && fs::exists( this->outputFileName() ) )
        {
            std::cout << "already file exist, ignore remeshSTL\n"
                      << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
            return;
        }

        if ( this->packageType() == "gmsh" )
            this->runGMSH();
        else if ( this->packageType() == "vmtk" )
            this->runVMTK();

        std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n";

    }
    void runVMTK()
    {
        CHECK( !this->stlFileName().empty() ) << "stlFileName is empty";



        std::ostringstream __str;
        // source ~/packages/vmtk/vmtk.build2/Install/vmtk_env.sh
        std::string pythonExecutable = BOOST_PP_STRINGIZE( PYTHON_EXECUTABLE );
        __str << pythonExecutable << " ";
        std::string dirBaseVmtk = "/Users/vincentchabannes/packages/vmtk/vmtk.build2/Install/bin/";
        __str << dirBaseVmtk << "vmtk " << dirBaseVmtk << "vmtksurfaceremeshing ";
        __str << "-ifile " << this->stlFileName() << " ";
        __str << "-ofile " << this->outputFileName() << " ";
        __str << "-area " << this->area() << " ";
        auto err = ::system( __str.str().c_str() );

        //std::cout << "hola\n"<< __str.str() <<"\n";
    }
    void runGMSH()
    {
        CHECK( !this->stlFileName().empty() ) << "stlFileName is empty";
        CHECK( !this->centerlinesFileName().empty() ) << "centerlinesFileName is empty";

        std::ostringstream geodesc;
        geodesc << "Mesh.Algorithm = 6; //(1=MeshAdapt, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) \n"
                << "Mesh.Algorithm3D = 1; //(1=tetgen, 4=netgen, 7=MMG3D, 9=R-tree) \n";

        //geodesc << "Merge \"stl_remesh_vmtk/fluidskin3.stl\"\n";
        geodesc << "Merge \""<< this->stlFileName() <<"\";\n";

        geodesc << "Field[1] = Centerline;\n";
        //geodesc << "Field[1].FileName = \"../centerlines/fluidskin3.vtk\"\n";
        geodesc << "Field[1].FileName = \"" << this->centerlinesFileName() << "\";\n";
        geodesc << "Field[1].nbPoints = "<< this->remeshNbPointsInCircle() << ";//15//25 //number of mesh elements in a circle\n";
        //geodesc << "Field[1].nbPoints = 15;//25 //number of mesh elements in a circle\n";

        //Remesh the initial stl
        geodesc << "Field[1].reMesh =1;\n"
                << "Field[1].run;\n"
                << "Background Field = 1;\n";

        fs::path outputMeshNamePath = fs::path(this->outputFileName());
        std::string _name = outputMeshNamePath.stem().string();
        std::string geoname=_name+".geo";

        std::ofstream geofile( geoname.c_str() );
        geofile << geodesc.str();
        geofile.close();

        //std::cout << " this->outputFileName() " << this->outputFileName() << "\n";
        detail::generateMeshFromGeo(geoname/*remeshGeoFileName*/,this->outputFileName()/*remeshMshFileName*/,2);
    }
private :
    std::string M_packageType;
    std::string M_stlFileName;
    bool M_forceRemesh;
    std::string M_centerlinesFileName;
    int M_remeshNbPointsInCircle;
    double M_area;
    std::string M_outputFileName;

}; // class ToolBoxBloodFlowReMeshSTL

//template< class FluidMeshType, class SolidMeshType >
class ToolBoxBloodFlowMesh
{
public :

    ToolBoxBloodFlowMesh()
        :
        M_stlFileName(soption(_name="mesh-volume.stl.filename")),
        M_centerlinesFileName(soption(_name="mesh-volume.centerlines.filename")),
        M_inletOutletDescFileName(soption(_name="mesh-volume.inletoutlet-desc.filename")),
        M_remeshNbPointsInCircle( ioption(_name="mesh-volume.nb-points-in-circle" ) ),
        M_extrudeWall( boption(_name="mesh-volume.extrude-wall") ),
        M_extrudeWallNbElemLayer( ioption(_name="mesh-volume.extrude-wall.nb-elt-layer") ),
        M_extrudeWallhLayer( doption(_name="mesh-volume.extrude-wall.h-layer") ),
        M_outputFileName()
        {
            this->updateOutputFileName();
        }
    ToolBoxBloodFlowMesh( ToolBoxBloodFlowMesh const& e )
        :
        M_stlFileName( e.M_stlFileName ),
        M_centerlinesFileName( e.M_centerlinesFileName ),
        M_inletOutletDescFileName( e.M_inletOutletDescFileName ),
        M_remeshNbPointsInCircle( e.M_remeshNbPointsInCircle ),
        M_extrudeWall( e.M_extrudeWall ),
        M_extrudeWallNbElemLayer( e.M_extrudeWallNbElemLayer ),
        M_extrudeWallhLayer( e.M_extrudeWallhLayer ),
        M_outputFileName( e.M_outputFileName )
        {}

    std::string stlFileName() const { return M_stlFileName; }
    std::string centerlinesFileName() const { return M_centerlinesFileName; }

    void setStlFileName(std::string s) { M_stlFileName=s;this->updateOutputFileName(); }
    void setCenterlinesFileName(std::string s) { M_centerlinesFileName=s; }

    std::string inletOutletDescFileName() const { return M_inletOutletDescFileName; }
    void setInletOutletDescFileName(std::string s) { M_inletOutletDescFileName=s; }

    int remeshNbPointsInCircle() const { return M_remeshNbPointsInCircle; }

    bool extrudeWall() const { return M_extrudeWall; }
    int extrudeWall_nbElemLayer() const { return M_extrudeWallNbElemLayer; }
    double extrudeWall_hLayer() const { return M_extrudeWallhLayer; }

    std::string outputFileName() const { return M_outputFileName; }

    void updateOutputFileName()
    {
        std::string arterialWallSpecStr = (this->extrudeWall())? "WithArterialWall" : "WithoutArterialWall";
        if ( this->extrudeWall() )
            arterialWallSpecStr += (boost::format("hLayer%1%nEltLayer%2%")%this->extrudeWall_hLayer() %this->extrudeWall_nbElemLayer() ).str();

        //std::string outputName = fs::path(this->stlFileName()).stem().string()+"_volumeMesh.vtk";
        std::string outputName = (boost::format( "%1%_%2%pt%3%%4%.msh")
                                  %fs::path(this->stlFileName()).stem().string()
                                  %"volumeMesh"
                                  %this->remeshNbPointsInCircle()
                                  %arterialWallSpecStr
                                  ).str();
        //std::string outputName = fs::path(this->stlFileName()).string()+"_centerlines.vtk";
        M_outputFileName = (fs::current_path()/outputName).string();
    }

    void run()
    {
        std::cout << "\n"
                  << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
                  << "---------------------------------------\n"
                  << "run ToolBoxBloodFlowMesh \n"
                  << "---------------------------------------\n";
        std::cout << "stlFileName             : " << this->stlFileName() << "\n"
                  << "centerlinesFileName     : " << this->centerlinesFileName() << "\n"
                  << "inletOutletDescFileName : " << this->inletOutletDescFileName() << "\n"
                  << "outputFileName          : " << this->outputFileName() << "\n"
                  << "---------------------------------------\n"
                  << "---------------------------------------\n";

        CHECK( !this->stlFileName().empty() ) << "stlFileName is empty";
        CHECK( !this->centerlinesFileName().empty() ) << "centerlinesFileName is empty";
        CHECK( !this->inletOutletDescFileName().empty() ) << "inletOutletDescFileName is empty";

        if ( fs::exists( this->outputFileName() ) )
        {
            std::cout << "already file exist, ignore meshVolume\n"
                      << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
            return;
        }


#if 0
        if ( true )
        {
            std::string remeshGeoFileName = "firstStep.geo";
            std::string remeshMshFileName = "firstStep.stl";
            this->generateGeoForRemeshSurfaceFromSTLAndCenterlines(remeshGeoFileName);
            /*this->*/detail::generateMeshFromGeo(remeshGeoFileName,remeshMshFileName,2);

            this->setStlFileName( remeshMshFileName );
        }
        else
        {
            std::string remeshMshFileName = "firstStep.stl";
            this->setStlFileName( remeshMshFileName );
        }
#endif
        std::string volumeGeoFileName = "secondStep.geo";
        std::string volumeMshFileName = "secondStep.msh";
        this->generateGeoFor3dVolumeFromSTLAndCenterlines(volumeGeoFileName);
        detail::generateMeshFromGeo(volumeGeoFileName,this->outputFileName()/*volumeMshFileName*/,3);

        std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";

    }
#if 0
    void generateMeshFromGeo( std::string inputGeoName,std::string outputMeshName,int dim/*,std::string formatExportMesh = "msh"*/)
    {
        //fs::path outputMeshNamePath=fs::path(Environment::findFile(outputMeshName));
        fs::path outputMeshNamePath=fs::path(outputMeshName);
        std::string _name = outputMeshNamePath.stem().string();
        std::cout << "\n _name " << _name << "\n";
        std::cout << "\n _ext " << outputMeshNamePath.extension() << "\n";

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

#endif
#if 0
    void generateGeoForRemeshSurfaceFromSTLAndCenterlines(std::string geoname)
    {
        std::ostringstream geodesc;
        geodesc << "Mesh.Algorithm = 6; //(1=MeshAdapt, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) \n"
                << "Mesh.Algorithm3D = 1; //(1=tetgen, 4=netgen, 7=MMG3D, 9=R-tree) \n";

        //geodesc << "Merge \"stl_remesh_vmtk/fluidskin3.stl\"\n";
        geodesc << "Merge \""<< this->stlFileName() <<"\";\n";

        geodesc << "Field[1] = Centerline;\n";
        //geodesc << "Field[1].FileName = \"../centerlines/fluidskin3.vtk\"\n";
        geodesc << "Field[1].FileName = \"" << this->centerlinesFileName() << "\";\n";
        geodesc << "Field[1].nbPoints = "<< this->remeshNbPointsInCircle() << ";//15//25 //number of mesh elements in a circle\n";
        //geodesc << "Field[1].nbPoints = 15;//25 //number of mesh elements in a circle\n";

        //Remesh the initial stl
        geodesc << "Field[1].reMesh =1;\n"
                << "Field[1].run;\n"
                << "Background Field = 1;\n";

        std::ofstream geofile( geoname.c_str() );
        geofile << geodesc.str();
        geofile.close();
    }
#endif
    void generateGeoFor3dVolumeFromSTLAndCenterlines(std::string geoname)
    {
        std::ostringstream geodesc;
        geodesc << "General.ExpertMode=1;\n";
        geodesc << "Mesh.Algorithm = 6; //(1=MeshAdapt, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) \n"
                << "Mesh.Algorithm3D = 1; //(1=tetgen, 4=netgen, 7=MMG3D, 9=R-tree) \n";

        //Mesh.CharacteristicLengthFactor=0.015;

        //Merge "stl_remesh_gmsh/fluidskin_P15.stl";
        geodesc << "Merge \""<< this->stlFileName() <<"\";\n";

        geodesc << "Field[1] = Centerline;\n";
        // centerline file in vtk format
        geodesc << "Field[1].FileName = \"" << this->centerlinesFileName() << "\";\n";

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
        geodesc << "Field[1].descInletOutlet = \"" << this->inletOutletDescFileName() << "\";\n";

        // apply field
        geodesc << "Field[1].run;\n"
                << "Background Field = 1;\n";


        std::ofstream geofile( geoname.c_str() );
        geofile << geodesc.str();
        geofile.close();
    }

private :
    std::string M_stlFileName;
    std::string M_centerlinesFileName;

    std::string M_inletOutletDescFileName;

    int M_remeshNbPointsInCircle;

    bool M_extrudeWall;
    int M_extrudeWallNbElemLayer;
    double M_extrudeWallhLayer;

    std::string M_outputFileName;

};

} // namespace Feel

#endif // __toolboxbloodflowmesh_H
