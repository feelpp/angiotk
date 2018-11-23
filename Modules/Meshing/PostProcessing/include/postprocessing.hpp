/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#ifndef __toolboxfsimesh_H
#define __toolboxfsimesh_H 1


#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/mesh.hpp>


namespace Feel
{

class MeshPartitioner
{
public :
    MeshPartitioner( std::string prefix );
    MeshPartitioner( MeshPartitioner const& m );

    std::string prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }
    worldcomm_ptr_t const& worldCommPtr() const { return Environment::worldCommPtr(); }
    std::string inputPath() const { return M_inputPath; }
    std::vector<std::string> const& outputPath() const { return M_outputPath; }
    std::string const& outputPath(int k) const { return M_outputPath[k]; }
    bool forceRebuild() const { return M_forceRebuild; }
    std::vector<int> const& nPartitions() const { return M_nPartitions; }
    int nPartitions(int k) const { return M_nPartitions[k]; }
    int partitioner() const { return M_partitioner; }

    void setInputPath( std::string s ) { M_inputPath=s; }
    void setOutputPath( std::string s ) { M_outputPath[0]=s; }
    void setNumberOfPartitions(int p) { M_nPartitions.resize(1);M_nPartitions[0]=p; }
    void setPartitioner(int p) { M_partitioner=p; }

    void updateOutputPathFromInputPath();
    void updateOutputPathFromInputFileName();
    void updateOutputPathFromPath( std::string path );

    void run();

    static po::options_description options( std::string const& prefix );
private :
    std::string M_prefix;
    std::string M_inputPath;//,M_outputPath;
    std::vector<std::string> M_outputPath;
    std::string M_outputDirectory;
    bool M_forceRebuild;
    std::vector<int> M_nPartitions;
    int M_partitioner;
};

template< class MeshType=Mesh<Simplex<3> > >
class ExtractSubMeshFromFSIMesh
{
public :

    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef typename mesh_type::shape_type shape_type;

    ExtractSubMeshFromFSIMesh( std::string prefix );
    ExtractSubMeshFromFSIMesh( ExtractSubMeshFromFSIMesh const& m );

    std::string prefix() const { return M_prefix; }
    WorldComm const& worldComm() const { return Environment::worldComm(); }
    std::string inputPath() const { return M_inputBloodFlowMeshFilename; }
    std::string outputPathLumen() const { return M_outputPathLumen; }
    std::string outputPathArterialWall() const { return M_outputPathArterialWall; }
    bool forceRebuild() const { return M_forceRebuild; }
    std::string markerNameLumenVolume() const { return M_markerNameLumenVolume; }
    std::string markerNameArterialWallVolume() const { return M_markerNameArterialWallVolume; }

    void updateOutputPathFromInputFileName();

    void run();

    static po::options_description options( std::string const& prefix );

private :
    std::string M_prefix;
    std::string M_inputBloodFlowMeshFilename;
    std::string M_outputDirectory;
    std::string M_outputPathLumen,M_outputPathArterialWall;
    bool M_forceRebuild;

    std::string M_markerNameLumenVolume,M_markerNameArterialWallVolume;

    mesh_ptrtype M_meshLumen, M_meshArterialWall;
};


} //Feel


#endif // __toolboxfsimesh_H
