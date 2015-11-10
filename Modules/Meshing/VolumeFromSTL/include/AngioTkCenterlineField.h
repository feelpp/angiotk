// Gmsh - Copyright (C) 1997-2014 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@geuz.org>.
//
// Contributor(s):
//   Emilie Marchandise

#ifndef _CENTERLINEFIELD_H_
#define _CENTERLINEFIELD_H_

#include "feel/feelconfig.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <memory> // std::shared_ptr

#include <gmshHeadersMissing/Field.h>
#include <MEdge.h>

#include <meshGFaceDelaunayInsertion.h>

#include <centerlinesmanageriodata.hpp>

class GModel;
class GFace;
class GPoint;
class MLine;
class MVertex;
class GEntity;
class MTriangle;   
class discreteEdge;
class discreteFace;
class MElement;
class SPoint3;


// A branch of a 1D tree
struct BranchDesc {
  int tag;
  std::vector<MLine*> lines;
  double length;
  MVertex *vB;
  MVertex *vE;
  //std::vector<Branch> children;
  double minRad;
  double maxRad;
  double boundMinX,boundMaxX,boundMinY,boundMaxY,boundMinZ,boundMaxZ;

  BranchDesc( std::vector<MLine*> _lines, int _tag, double _length,  MVertex *_vB, MVertex *_vE,double _minRad,double _maxRad )
  :
  lines(  _lines ), tag(_tag), length( _length ), vB(_vB),vE(_vE),minRad(_minRad), maxRad(_maxRad),
    boundMinX(0.),boundMaxX(0.),boundMinY(0.),boundMaxY(0.),boundMinZ(0.),boundMaxZ(0.)
  {}
  BranchDesc( BranchDesc const& _branch ) = default;
  //~BranchDesc(){}
  bool
  isInsideBox( MVertex * vertex, double dist )
  {
    double lengthAdded = maxRad+dist;
    return ( ( vertex->x() > ( boundMinX - lengthAdded ) ) && ( vertex->x() < ( boundMaxX + lengthAdded ) ) &&
	     ( vertex->y() > ( boundMinY - lengthAdded ) ) && ( vertex->y() < ( boundMaxY + lengthAdded ) ) &&
	     ( vertex->z() > ( boundMinZ - lengthAdded ) ) && ( vertex->z() < ( boundMaxZ + lengthAdded ) ) );
  }

};

#if defined(FEELPP_HAS_ANN_H)
class ANNkd_tree;


int readVTKPolyDataFields( const std::string &name, std::map<std::string,std::vector< std::vector<double> > > & fieldsPointData, bool bigEndian=false );

// This class takes as input A 1D mesh which is the centerline
// of a tubular 2D surface mesh
// It computes a mesh size field function of the distance to the centerlines
// and a cross field that is the direction of the centerline 
// It splits the tubular structure in many mesh partitions
// along planes perpendicuar to the centerlines 

class AngioTkCenterline : public Field{

 protected: 
  GModel *current; //current GModel
  std::shared_ptr<GModel> mod;
  //GModel *mod; //centerline GModel
  GModel *split; //split GModel
  ANNkd_tree *kdtree, *kdtreeR; 
  std::string fileName;
  int nbPoints;
  double recombine;
  int NF, NV, NE, NR;
  bool useGmshExecutable;
  int is_cut, is_closed, is_extruded;
  int is_clip_mesh;
  double M_clipMeshScalingFactor;
  double hLayer;
  double hSecondLayer;
  int nbElemLayer;
  int nbElemSecondLayer;
  std::string descInletOutlet;
  std::map<int,std::string> mapBoundEdgeIdToPhysicalMarkerLumen;
  std::map<int,std::string> mapBoundEdgeIdToPhysicalMarkerArerialWall;
  //inlet point
  SPoint3 ptin;
  //all (unique) lines of centerlines
  std::vector<MLine*> lines;
  //the stuctured tree of the centerlines
  std::vector<BranchDesc> edges;
  //the radius of the surface mesh at a given line
  std::map<MLine*,double> radiusl;
  //the junctions of the tree
  std::set<MVertex*> junctions;
  //some colors (int) for all points and lines
  std::map<MVertex*,int> colorp;
  std::map<MLine*,int> colorl;

  std::vector<GEdge*> modEdges;

  // save map from vertex to set of ( branchId,lineId )
  std::map<MVertex*,std::set<std::pair<int,int> > > M_vertexToLinesId;
  // save junction points and extremity points 
  //std::set<MVertex*> M_junctionsVertex;
  std::map<MVertex*,std::set<int> > M_junctionsVertex;
  //std::set<MVertex*> M_extremityVertex;
  // vertex -> ( branchId, lineIdInBranch )
  std::map<MVertex*, std::pair<int,int> > M_extremityVertex;


  //the tubular surface mesh
  std::vector<MTriangle*> triangles;
  std::vector<MVertex*> vertices;
  
  //the lines cut of the tubular mesh by planes
  std::set<MEdge,Less_Edge> theCut;
  std::map<int, std::set<MEdge,Less_Edge> > theCutMap;
  std::set<MVertex*> theCutV;

  //discrete edes and faces created by the cut
  std::vector<discreteEdge*> discEdges;
  std::vector<discreteFace*> discFaces;

 public:
  std::map<std::string,std::vector<std::vector<double> > >  M_centerlinesFieldsPointData;

  std::vector<std::vector<double> > const& centerlinesFieldsPointData(std::string const& key ) const
    {
      return M_centerlinesFieldsPointData.find(key)->second;
    }
  std::vector<double> const& centerlinesFieldsPointData(std::string const& key, MVertex* vertex ) const
    {
      return this->centerlinesFieldsPointData(key, this->mapVertexGmshIdToVtkId(vertex->getIndex()));
    }
  std::vector<double> const& centerlinesFieldsPointData(std::string const& key,int vtkId ) const
    {
      return M_centerlinesFieldsPointData.find(key)->second[vtkId];
    }

 public:
  AngioTkCenterline(std::string fileName);
  AngioTkCenterline();
  ~AngioTkCenterline();

  virtual bool isotropic () const {return false;}
  virtual const char *getName()
  {
    return "centerline Field";
  }
  virtual std::string getDescription()
  {
    return "The value of this field is the distance to the centerline.\n\n"
" You should specify a fileName that contains the centerline."
" The centerline of a surface can be obtained with the open source software vmtk (http://www.vmtk.org/)"
" using the following script:\n\n"
"vmtk vmtkcenterlines -seedselector openprofiles -ifile mysurface.stl -ofile centerlines.vtp --pipe vmtksurfacewriter -ifile centerlines.vtp -ofile centerlines.vtk\n";
  }

  //isotropic operator for mesh size field function of distance to centerline
  double operator() (double x, double y, double z, GEntity *ge=0);
  //anisotropic operator
  void operator() (double x, double y, double z, SMetric3 &metr, GEntity *ge=0);

  void computeCrossField(double x,double y,double z,
			 SVector3 &d1, SVector3 &d2,  SVector3 &d3);



  void createFromFile( std::string const& fileName, std::string const& inputSurfacePath );
  void createFromGeo( std::tuple< std::vector< std::pair<GPoint,double> >, std::vector<std::vector<int> > > const& geoDesc, std::string const& outputSmoothGeoPath );
  void createFromCenterlines( AngioTkCenterline const& inputCenterlines, std::string const& outputSmoothGeoPath,
			      double meshSizeUniform=1., double resampleGeoPointSpacing = 4. );
 private :
  void createFromGeoCenterlinesFile( std::string const& fileName, std::string const& inputSurfacePath );
  void updateCenterlinesForUse(std::vector<GEdge*> const& _modEdges);
  void updateCenterlinesForUse(std::map<int,std::vector<MLine*> > const& _modEdges);
  //void fixBranchConnectivity();
  void checkCenterlinesConnectivity();
  void updateFieldsDataAfterReduction( std::map<int,int> const& _previousMapVertexGmshIdToVtkId );

 public:
  void updateCenterlinesFieldsFromFile(std::string fileName);

  void removeBranchIds( std::set<int> const& _removeBranchIds );
  void cleanBranch();
  void removeDuplicateBranch();
  void addFieldBranchIds( std::string const& fieldName = "BranchIds" );
  void addFieldRadiusMin( std::string const& fieldName = "RadiusMin" );
  bool hasField( std::string const& fieldName ) const { return M_centerlinesFieldsPointData.find( fieldName ) != M_centerlinesFieldsPointData.end(); }

  void applyFieldThresholdMin( std::string const& fieldName,double value );
  void applyFieldThresholdMax( std::string const& fieldName,double value );
  void applyFieldThresholdMin( std::vector<std::string> const& fieldName,double value );
  void applyFieldThresholdMax( std::vector<std::string> const& fieldName,double value );
  void applyFieldThresholdZoneMin(std::string const& fieldName, double value, std::map<int,std::vector<std::tuple<double,double,double> > > const& mapPointPair );
  void applyFieldThresholdZoneMax(std::string const& fieldName, double value, std::map<int,std::vector<std::tuple<double,double,double> > > const& mapPointPair );
  void applyFieldThresholdZoneMin(std::vector<std::string> const& fieldName, double value, std::map<int,std::vector<std::tuple<double,double,double> > > const& mapPointPair );
  void applyFieldThresholdZoneMax(std::vector<std::string> const& fieldName, double value, std::map<int,std::vector<std::tuple<double,double,double> > > const& mapPointPair );
  void applyFieldThresholdZoneImpl(std::vector<std::string> const& fieldName, double value, std::map<int,std::vector<std::tuple<double,double,double> > > const& mapPointPair, int type );

  void applyTubularColisionFix( std::vector<MVertex*> const& vTestedSet, double distMin );
  void applyTubularColisionFix( std::map< MVertex*,std::set<MVertex*> > const& mapVertexTested, double distMin, int method /*= 0*/, int maxrecurrence/*=-1*/,int nrecurrence = 0 );
  void applyTubularColisionFix( AngioTk::pointpair_data_type const& pointPair, double distMin );
  void applyTubularColisionFix( double distMin );


  void writeCenterlinesVTK( std::string fileName );

  std::map<MVertex*, std::pair<int,int> > const&
    centerlinesExtremities() const { return M_extremityVertex; }
  std::vector<BranchDesc> const& centerlinesBranch() const { return edges; }
  BranchDesc const& centerlinesBranch(int k) const { return edges[k]; }

  std::tuple<MVertex*,double> foundClosestPointInCenterlines( double ptToLocalize[3]);
  //std::tuple< std::vector<MLine*> , std::map<MVertex*,int>, double > pathBetweenVertex( MVertex* vertexA, MVertex* vertexB );
  std::tuple< std::vector<MLine*> , double > pathBetweenVertex( MVertex* vertexA, MVertex* vertexB );
  int vertexOnSameBranch( MVertex* vertexA, MVertex* vertexB );
  //std::set<int>
  std::tuple<MVertex*, std::vector<int> > vertexOnNeighboringBranch( MVertex* vertexA, MVertex* vertexB );
  bool canFindPathBetweenVertex( MVertex* vertexA, MVertex* vertexB );

  double maxScalarValueInPath( std::vector<MLine*> const& path, std::string const& fieldName ) const;

  std::map<MLine*,double> const& centerlinesRadiusl() const { return radiusl; }
  double minRadiusAtVertex( MVertex* myvertex ) const;

  //double centerlinesRadiusFromLine() const { 
    //auto itr = M_angioTkCenterlines->centerlinesRadiusl().find( mylines[lineIdInBranch]/*mylines.front()*/);
  //}
  //std::pair<std::map<int,int>,std::map<int,int> > computeRelationVertex() const;
  void updateRelationMapVertex();
  void updateRelationMapVertex(std::map<int,int> & _mapVertexGmshIdToVtkId,
			       std::map<int,int> & _mapVertexVtkIdToGmshId );
  std::map<int,int> M_mapVertexGmshIdToVtkId, M_mapVertexVtkIdToGmshId;
  int mapVertexGmshIdToVtkId(int k) const { return M_mapVertexGmshIdToVtkId.find(k)->second; }
  int mapVertexVtkIdToGmshId(int k) const { return M_mapVertexVtkIdToGmshId.find(k)->second; }


  std::set<std::pair<int,int> > M_registerLinesToRemoveFromPointIdPairInModelEdge;
  std::set<MLine*> M_registerLinesDuplicatedToIgnore;

  std::vector<std::shared_ptr<AngioTkCenterline> > M_attachAngioTkCenterline;
  void attachAngioTkCenterline( std::shared_ptr<AngioTkCenterline> const& obj ) { if ( obj ) M_attachAngioTkCenterline.push_back(obj); }

  void importSurfaceFromFile(std::string const& fileName );
  //import the 1D mesh of the centerlines (in vtk format)
  //and fill the vector of lines
  void importFile(std::string fileName);

  //refine the 1D mesh to have many points on the 1D centerlines 
  //for annKDTree distance computations
  void buildKdTree();

  //Creates the branch structure (topology, connectivity) from the 
  //vector of lines
  void createBranches(int maxN);

  //Computes for the Branches the min and maxRadius of the tubular structure
  //this function needs the current GModel
  void computeRadii();

  //Computes for each MLine the minRadius
  void distanceToSurface();

  //actions
  void run();

  //load desc file and define marker which are close to a point
  void initPhysicalMarkerFromDescFile( std::vector<GEdge*> boundEdges );

  // Cut the mesh in different parts of small aspect ratio
  void cutMesh( std::string const& remeshPartitionMeshFile="" );
  //Create In and Outlet Planar Faces
  void createClosedVolume(GEdge *gin, std::vector<GEdge*> boundEdges);
  //extrude outer wall
  void extrudeBoundaryLayerWall(GEdge *gin, std::vector<GEdge*> boundEdges);

  // Cut the tubular structure with a disk
  // perpendicular to the tubular structure
  bool cutByDisk(SVector3 &pt, SVector3 &dir, double &maxRad, int tag = -1);

  //create discrete faces
  void createFaces();
  void createFacesFromClip();
  void createSplitCompounds();

  //Print for debugging
  void printSplit() const;
 
  SMetric3 metricBasedOnSurfaceCurvature(SVector3 dMin, SVector3 dMax, double cMin, double cMax,
					  double lc_n, double lc_t1, double lc_t2);

  // mode remesh
  void setIsCut(bool b) { is_cut = b; }
  void setRemeshNbPoints( int i ) { nbPoints=i; }
  void setSurfaceRemeshRadiusUncertainty(double d) { M_surfaceRemeshRadiusUncertainty = d; }
  void runSurfaceRemesh( std::string const& remeshPartitionMeshFile="", bool forceRebuildPartition=true );
  void saveSurfaceRemeshSTL(std::string const outputPath, bool binary );

  // mode clip mesh
  void setModeClipMesh(bool b) { is_clip_mesh = b; }
  void setClipMeshScalingFactor( double d ) { M_clipMeshScalingFactor = d; }
  void runClipMesh();
  void saveClipMeshSTL(std::string const outputPath, bool binary );

  // tubular extension
  void runTubularExtension();
  void saveTubularExtensionSTL(std::string const outputPath, bool binary );


 private :
  void saveCurrentGModelSTL(std::string const outputPath, bool binary );
  void updateMergeFromExtremities( AngioTkCenterline const& centerlinesMerged,
				   std::map<MVertex*,std::pair< std::vector<std::pair<MLine*,int> >, MVertex*> > & indexVertexReplaced );

  // represent a length where the radius can vary ( e.g. can be equal to the image accuracy )
  double M_surfaceRemeshRadiusUncertainty;
};
#else
class AngioTkCenterline : public Field{

 public:
  AngioTkCenterline(std::string fileName){ Msg::Error("Gmsh has to be compiled with ANN support to use CenterlineFields");}
  AngioTkCenterline(){ Msg::Error("Gmsh has to be compiled with ANN support to use CenterlineFields");}
  ~AngioTkCenterline();

  virtual bool isotropic () const {return false;}
  virtual const char *getName()
  {
    return "centerline Field";
  }
  virtual std::string getDescription()
  {
    return "The value of this field is the distance to the centerline.\n\n"
" You should specify a fileName that contains the centerline."
" The centerline of a surface can be obtained with the open source software vmtk (http://www.vmtk.org/)"
" using the following script:\n\n"
"vmtk vmtkcenterlines -seedselector openprofiles -ifile mysurface.stl -ofile centerlines.vtp --pipe vmtksurfacewriter -ifile centerlines.vtp -ofile centerlines.vtk\n";
  }

  //isotropic operator for mesh size field function of distance to centerline
  double operator() (double x, double y, double z, GEntity *ge=0);
  //anisotropic operator
  void operator() (double x, double y, double z, SMetric3 &metr, GEntity *ge=0);

  //temporary operator where v1, v2 and v3 are three orthonormal directions
  void operator()(double x,double y,double z,SVector3& v1,SVector3& v2,SVector3& v3,GEntity* ge=0);

};

#endif

#endif
