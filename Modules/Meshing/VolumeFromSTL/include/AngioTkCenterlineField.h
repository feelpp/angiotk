// Gmsh - Copyright (C) 1997-2014 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@geuz.org>.
//
// Contributor(s):
//   Emilie Marchandise

#ifndef _CENTERLINEFIELD_H_
#define _CENTERLINEFIELD_H_

#include <vector>
#include <map>
#include <set>
#include <string>
#include <memory> // std::shared_ptr

#include <gmshHeadersMissing/Field.h>
#include <MEdge.h>

#include <meshGFaceDelaunayInsertion.h>
class GModel;
class GFace;
class MLine;
class MVertex;
class GEntity;
class MTriangle;   
class discreteEdge;
class discreteFace;
class MElement;
class SPoint3;

// A branch of a 1D tree
struct Branch{
  int tag;
  std::vector<MLine*> lines;
  double length;
  MVertex *vB;
  MVertex *vE;
  std::vector<Branch> children;
  double minRad;
  double maxRad;
};

#if defined(HAVE_ANN)
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
  int is_cut, is_closed, is_extruded;
  int is_clip_mesh;
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
  std::vector<Branch> edges;
  //the radius of the surface mesh at a given line
  std::map<MLine*,double> radiusl;
  //the junctions of the tree
  std::set<MVertex*> junctions;
  //some colors (int) for all points and lines
  std::map<MVertex*,int> colorp;
  std::map<MLine*,int> colorl;

  std::vector<GEdge*> modEdges;

  // save junction points and extremity points 
  std::set<MVertex*> M_junctionsVertex;
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
  std::map<std::string,std::vector<std::vector<double> > >  centerlinesFieldsPointData;

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


  void updateCenterlinesFromFile( std::string fileName );
  void updateCenterlinesFieldsFromFile(std::string fileName);
  void removeBranchIds( std::set<int> const& _removeBranchIds );
  void addFieldBranchIds();
  void writeCenterlinesVTK( std::string fileName );

  std::map<MVertex*, std::pair<int,int> > const&
    centerlinesExtremities() const { return M_extremityVertex; }
  std::vector<Branch> const& centerlinesBranch() const { return edges; }
  Branch const& centerlinesBranch(int k) const { return edges[k]; }

  std::map<MLine*,double> /*const&*/ centerlinesRadiusl() const { return radiusl; }
  //double centerlinesRadiusFromLine() const { 
    //auto itr = M_angioTkCenterlines->centerlinesRadiusl().find( mylines[lineIdInBranch]/*mylines.front()*/);
  //}
  //std::pair<std::map<int,int>,std::map<int,int> > computeRelationVertex() const;
  void updateRelationMapVertex();
  void updateRelationMapVertex(std::map<int,int> & _mapVertexGmshIdToVtkId,
			       std::map<int,int> & _mapVertexVtkIdToGmshId );
  std::map<int,int> M_mapVertexGmshIdToVtkId, M_mapVertexVtkIdToGmshId;
  std::set<std::pair<int,int> > M_registerLinesToRemoveFromPointIdPairInModelEdge;


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

  void runClipMesh();
  void createFacesFromClip();

  //load desc file and define marker which are close to a point
  void initPhysicalMarkerFromDescFile( std::vector<GEdge*> boundEdges );

  // Cut the mesh in different parts of small aspect ratio
  void cutMesh();
  //Create In and Outlet Planar Faces
  void createClosedVolume(GEdge *gin, std::vector<GEdge*> boundEdges);
  //extrude outer wall
  void extrudeBoundaryLayerWall(GEdge *gin, std::vector<GEdge*> boundEdges);

  // Cut the tubular structure with a disk
  // perpendicular to the tubular structure
  bool cutByDisk(SVector3 &pt, SVector3 &dir, double &maxRad, int tag = -1);

  //create discrete faces
  void createFaces();
  void createSplitCompounds();

  //Print for debugging
  void printSplit() const;
 
  SMetric3 metricBasedOnSurfaceCurvature(SVector3 dMin, SVector3 dMax, double cMin, double cMax,
					  double lc_n, double lc_t1, double lc_t2);

};
#else
class AngioTkCenterline : public Field{

 public:
  Centerline(std::string fileName){ Msg::Error("Gmsh has to be compiled with ANN support to use CenterlineFields");}
  Centerline(){ Msg::Error("Gmsh has to be compiled with ANN support to use CenterlineFields");}
  ~Centerline();

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
