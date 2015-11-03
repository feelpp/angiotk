// Gmsh - Copyright (C) 1997-2014 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@geuz.org>.
//
// Contributor(s):
//   Emilie Marchandise
//

#include <AngioTkCenterlineField.h>

#include <vector>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <OS.h>
#include <GModel.h>
#include <MElement.h>
#include <MPoint.h>
#include <MTriangle.h>
#include <MVertex.h>
#include <MLine.h>
#include <StringUtils.h>
#include <GEntity.h>
//#include <Field.h>
#include <GFace.h>
#include <discreteEdge.h>
#include <discreteFace.h>
#include <GEdgeCompound.h>
#include <GFaceCompound.h>
#include <gmshHeadersMissing/BackgroundMesh.h>
#include <meshGFace.h>
#include <meshGEdge.h>
#include <MQuadrangle.h>
#include <gmshHeadersMissing/Curvature.h>
#include <MElement.h>
#include <Context.h>
#include <directions3D.h>
#include <gmshHeadersMissing/meshGRegion.h>
#include <gmshHeadersMissing/MVertexPositionSet.h>


#if defined(FEELPP_HAS_ANN_H)
#include <ANN/ANN.h>


static void erase(std::vector<MLine*>& lines, MLine* l)
{
  std::vector<MLine*>::iterator it = std::remove(lines.begin(), lines.end(), l);
  lines.erase(it, lines.end());
}

static double computeLength(std::vector<MLine*> lines)
{
  double length= 0.0;
  for (unsigned int i = 0; i< lines.size(); i++){
    length += lines[i]->getLength();
  }
  return length;
}

static bool isClosed(std::set<MEdge, Less_Edge> &theCut)
{
  std::multiset<MVertex*> boundV;
  std::set<MEdge,Less_Edge>::iterator it = theCut.begin();
  for (; it!= theCut.end(); it++){
    boundV.insert(it->getVertex(0));
    boundV.insert(it->getVertex(1));
  }
  std::multiset<MVertex*>::iterator itb = boundV.begin();
  for ( ; itb != boundV.end(); ++itb){
    if (boundV.count(*itb) != 2) {
      return false;
    }
  }
  return true;
}

static void orderMLines(std::vector<MLine*> &lines, MVertex *vB, MVertex *vE)
{
  std::vector<MLine*> _m;
  std::list<MLine*> segments;
  std::map<MVertex*, MLine*> boundv;
  std::vector<int> _orientation;

  // store all lines in a list : segments
  for (unsigned int i = 0; i < lines.size(); i++){
    segments.push_back(lines[i]);
  }

  // find a lonely MLine
  for (std::list<MLine*>::iterator it = segments.begin();
       it != segments.end(); ++it){
    MVertex *vL = (*it)->getVertex(0);
    MVertex *vR = (*it)->getVertex(1);
    std::map<MVertex*,MLine*>::iterator it1 = boundv.find(vL);
    if (it1 == boundv.end()) boundv.insert(std::make_pair(vL,*it));
    else  boundv.erase(it1);
    std::map<MVertex*,MLine*>::iterator it2 = boundv.find(vR);
    if (it2 == boundv.end()) boundv.insert(std::make_pair(vR,*it));
    else boundv.erase(it2);
  }

  // find the first MLine and erase it from the list segments
  MLine *firstLine;
  if (boundv.size() == 2){   // non periodic
    MVertex *v = (boundv.begin())->first;
    if ( v == vB) firstLine = (boundv.begin())->second;
    else{
      MVertex *v = (boundv.rbegin())->first;
      if (v == vB) firstLine = (boundv.rbegin())->second;
      else{
        Msg::Error("begin vertex not found for branch");
        return;
      }
    }
    for (std::list<MLine*>::iterator it = segments.begin();
         it != segments.end(); ++it){
      if (*it == firstLine){
        segments.erase(it);
        break;
      }
    }
  }
  else if (boundv.size() == 0){   // periodic
    if ( vB != vE )
      Msg::Error("vB and vE must be identical in periodic centerline");

    std::list<MLine*>::iterator it = segments.begin();
    for ( ; it != segments.end() ; ++it) {
      if ( (*it)->getVertex(0) == vB )
	break;
      }
    if ( it != segments.end() )
      {
	firstLine = *it;
	segments.erase(it);
      }
    else
      Msg::Error("segment not find");
  }
  else{
    Msg::Error("line is wrong (it has %d end points)",  boundv.size());
  }

  // loop over all segments to order segments and store it in the list _m
  _m.push_back(firstLine);
  _orientation.push_back(1);
  MVertex *first = _m[0]->getVertex(0);
  MVertex *last = _m[0]->getVertex(1);
  while (first != last){
    if (segments.empty())break;
    bool found = false;
    for (std::list<MLine*>::iterator it = segments.begin();
         it != segments.end(); ++it){
      MLine *e = *it;
      if (e->getVertex(0) == last){
        _m.push_back(e);
        segments.erase(it);
        _orientation.push_back(1);
        last = e->getVertex(1);
        found = true;
        break;
      }
      else if (e->getVertex(1) == last){
        _m.push_back(e);
        segments.erase(it);
        _orientation.push_back(0);
        last = e->getVertex(0);
        found = true;
        break;
      }
    }
    if (!found  && _orientation[0]==1){ //reverse orientation of first Line
      Msg::Warning("reverse orientation of first Line (%d,%d)",boundv.size(),segments.size());
      if (_m.size() == 1 ){
        MVertex *temp = first;
        first = last;
        last = temp;
        _orientation[0] = 0;
      }
      else {
        Msg::Error("lines is wrong");
        return;
      }
    }
  }

  //lines is now a list of ordered MLines
  lines = _m;

  //special case reverse orientation
  if (lines.size() < 2) return;
  if (_orientation[0] && lines[0]->getVertex(1) != lines[1]->getVertex(1)
      && lines[0]->getVertex(1) != lines[1]->getVertex(0)){
    for (unsigned int i = 0; i < lines.size(); i++)
      _orientation[i] = !_orientation[i];
  }
}

static void recurConnectByMEdge(const MEdge &e,
				std::multimap<MEdge, MTriangle*, Less_Edge> const&e2e,
				std::set<MTriangle*> &group,
				std::set<MEdge, Less_Edge> &touched,
				std::set<MEdge, Less_Edge> const&theCut)
{
  if (touched.find(e) != touched.end()) return;
  touched.insert(e);
  for (std::multimap <MEdge, MTriangle*, Less_Edge>::const_iterator it = e2e.lower_bound(e);
       it != e2e.upper_bound(e); ++it){
    if ( !it->second )
      {
	Msg::Error("triangle pointer is null");
	continue;
      }
    // if triangle already register, go to the next
    if ( group.find( it->second ) != group.end() )
      continue;
    group.insert(it->second);

    for (int i = 0; i < it->second->getNumEdges(); ++i){
      MEdge me = it->second->getEdge(i);
      if (theCut.find(me) != theCut.end()) {
	touched.insert(me); //break;
      }
      else recurConnectByMEdge(me, e2e, group, touched, theCut);
    }
  }
}

static MVertex* getCommonVertexInEdge( MEdge const& edge1, MEdge const& edge2 )
{
  if ( edge1.getVertex(0) == edge2.getVertex(0) || edge1.getVertex(0) == edge2.getVertex(1) )
    return edge1.getVertex(0);
  else if ( edge1.getVertex(1) == edge2.getVertex(0) || edge1.getVertex(1) == edge2.getVertex(1) )
    return edge1.getVertex(1);

  Msg::Error("not common edge");
  return edge1.getVertex(0); 
}
static
std::tuple<MVertex*,std::pair<MVertex*,MVertex*> >
getVertexDistributionFromEdgePair( MEdge const& edge1, MEdge const& edge2)
{
    MVertex* commonVertex = getCommonVertexInEdge(edge1,edge2);
    std::vector<MVertex*> otherV;
    if ( commonVertex == edge1.getVertex(0) ) otherV.push_back( edge1.getVertex(1) );
    else otherV.push_back( edge1.getVertex(0) );
    if ( commonVertex == edge2.getVertex(0) ) otherV.push_back( edge2.getVertex(1) );
    else otherV.push_back( edge2.getVertex(0) );
    return std::make_tuple( commonVertex ,std::make_pair(otherV[0],otherV[1]) );
}


static void cutTriangle(MTriangle *tri,
                        std::map<MEdge,MVertex*,Less_Edge> &cutEdges,
                        std::vector<MVertex*> &cutVertices,
                        std::vector<MTriangle*> &newTris,
                        std::set<MEdge,Less_Edge> &newCut)
{
  MVertex *c[3] = {0,0,0};
  for (int j=0;j<3;j++){
    MEdge ed = tri->getEdge(j);
    std::map<MEdge,MVertex*,Less_Edge> :: iterator it = cutEdges.find(ed);
    if (it != cutEdges.end()){
      c[j] = it->second;

    }
  }
  MVertex *old_v0  = tri->getVertex(0);
  MVertex *old_v1  = tri->getVertex(1);
  MVertex *old_v2  = tri->getVertex(2);

#if 0 // vtk
  // tri.eddge0 -> (tri.pt0,tri.pt1)  ; tri.eddge1 -> (tri.pt1,tri.pt2)  ;   tri.eddge2 -> (tri.pt0,tri.pt2) 
  if ( tri->getEdge(0).getVertex(0) == old_v2 || tri->getEdge(0).getVertex(1) == old_v2 ) Msg::Error("Error cut triangle0");
  if ( tri->getEdge(1).getVertex(0) == old_v0 || tri->getEdge(1).getVertex(1) == old_v0 ) Msg::Error("Error cut triangle1");
  if ( tri->getEdge(2).getVertex(0) == old_v1 || tri->getEdge(2).getVertex(1) == old_v1 ) Msg::Error("Error cut triangle2");
#endif



#if 1
  std::set<int> edgeMustBeCut;
  for ( int eId : std::vector<int>({0,1,2}) )
    if ( c[eId] ) edgeMustBeCut.insert(eId);

  if (edgeMustBeCut.size()>2){
    Msg::Error( "invalid cut : 3 edges has an intersection");
  }
  else if (edgeMustBeCut.size()==2){
    int eId0 = *edgeMustBeCut.begin(), eId1 = *edgeMustBeCut.rbegin();
    auto resDist = getVertexDistributionFromEdgePair( tri->getEdge(eId0),tri->getEdge(eId1) );
    MVertex* commonVertex = std::get<0>(resDist);
    MVertex* otherVertexInEdge0 = std::get<1>(resDist).first;
    MVertex* otherVertexInEdge1 = std::get<1>(resDist).second;
    if ( tri->getEdge(eId0).getVertex(0) == commonVertex )
      {
	newTris.push_back(new MTriangle (commonVertex,c[eId0],c[eId1]));
	newTris.push_back(new MTriangle (c[eId0],otherVertexInEdge0,otherVertexInEdge1));
	newTris.push_back(new MTriangle (otherVertexInEdge1,c[eId1],c[eId0] ));
      }
    else
      {
	newTris.push_back(new MTriangle (c[eId0],commonVertex,c[eId1]));
	newTris.push_back(new MTriangle (otherVertexInEdge0, c[eId0],otherVertexInEdge1));
	newTris.push_back(new MTriangle (otherVertexInEdge1,c[eId0],c[eId1] ));
      }
    newCut.insert(MEdge(c[eId0],c[eId1]));
  }
  else if (edgeMustBeCut.size()==1) {

    int eIdCut = *edgeMustBeCut.begin();
    int eId0 = (eIdCut==0)?1:0, eId1 = (eIdCut==0)? 2 : (eIdCut==1)? 2 : 1;
    auto resDist = getVertexDistributionFromEdgePair( tri->getEdge(eId0),tri->getEdge(eId1) );
    MVertex* commonVertex = std::get<0>(resDist);
    MVertex* otherVertexInEdge0 = std::get<1>(resDist).first;
    MVertex* otherVertexInEdge1 = std::get<1>(resDist).second;

    if (std::find(cutVertices.begin(), cutVertices.end(), commonVertex) != cutVertices.end())
      {
	if ( tri->getEdge(eIdCut).getVertex(0) == otherVertexInEdge0 )
	  {
	    newTris.push_back(new MTriangle (otherVertexInEdge0,c[eIdCut],commonVertex));
	    newTris.push_back(new MTriangle (commonVertex,c[eIdCut],otherVertexInEdge1));
	  }
	else
	  {
	    newTris.push_back(new MTriangle (otherVertexInEdge1,c[eIdCut],commonVertex));
	    newTris.push_back(new MTriangle (commonVertex,c[eIdCut],otherVertexInEdge0));
	  }
	    newCut.insert(MEdge(c[eIdCut],commonVertex));
      }
    else
      {
	//Msg::Error( "invalid cut0");
      }
  }
  else if ( edgeMustBeCut.empty() ) {
    newTris.push_back(tri);

    if (std::find(cutVertices.begin(), cutVertices.end(), old_v0) != cutVertices.end() &&
	std::find(cutVertices.begin(), cutVertices.end(), old_v1) != cutVertices.end())
      newCut.insert(MEdge(old_v0,old_v1));
    else if (std::find(cutVertices.begin(), cutVertices.end(), old_v1) != cutVertices.end() &&
	std::find(cutVertices.begin(), cutVertices.end(), old_v2) != cutVertices.end())
      newCut.insert(MEdge(old_v1,old_v2));
    else if (std::find(cutVertices.begin(), cutVertices.end(), old_v2) != cutVertices.end() &&
	std::find(cutVertices.begin(), cutVertices.end(), old_v0) != cutVertices.end())
      newCut.insert(MEdge(old_v2,old_v0));
  }







#else

  if (c[0] && c[1]){
    newTris.push_back(new MTriangle (c[0],old_v1,c[1]));
    newTris.push_back(new MTriangle (old_v0,c[0],old_v2));
    newTris.push_back(new MTriangle (old_v2,c[0],c[1]));
    newCut.insert(MEdge(c[0],c[1]));
  }
  else if (c[0] && c[2]){
    newTris.push_back(new MTriangle (old_v0,c[0],c[2]));
    newTris.push_back(new MTriangle (c[0],old_v1,old_v2));
    newTris.push_back(new MTriangle (old_v2,c[2],c[0]));
    newCut.insert(MEdge(c[0],c[2]));
  }
  else if (c[1] && c[2]){
    newTris.push_back(new MTriangle (old_v2,c[2],c[1]));
    newTris.push_back(new MTriangle (old_v0,old_v1,c[2]));
    newTris.push_back(new MTriangle (c[2],old_v1,c[1]));
    newCut.insert(MEdge(c[1],c[2]));
  }
  else if (c[0]){
    newTris.push_back(new MTriangle (old_v0,c[0],old_v2));
    newTris.push_back(new MTriangle (old_v2,c[0],old_v1));
    if (std::find(cutVertices.begin(), cutVertices.end(), old_v0) != cutVertices.end()){
      newCut.insert(MEdge(c[0],old_v0));
    }
    else if (std::find(cutVertices.begin(), cutVertices.end(), old_v1) != cutVertices.end()) {
      newCut.insert(MEdge(c[0],old_v1));
    }
    else if (std::find(cutVertices.begin(), cutVertices.end(), old_v2) != cutVertices.end()){
      newCut.insert(MEdge(c[0],old_v2));
    }
  }
  else if (c[1]){
    newTris.push_back(new MTriangle (old_v1,c[1],old_v0));
    newTris.push_back(new MTriangle (old_v0,c[1],old_v2));
    if (std::find(cutVertices.begin(), cutVertices.end(), old_v0) != cutVertices.end()){
      newCut.insert(MEdge(c[1],old_v0));
    }
    else if (std::find(cutVertices.begin(), cutVertices.end(), old_v1) != cutVertices.end()) {
      newCut.insert(MEdge(old_v1, c[1]));
    }
    else if (std::find(cutVertices.begin(), cutVertices.end(), old_v2) != cutVertices.end()){
      newCut.insert(MEdge(c[1],old_v2));
    }
  }
  else if (c[2]){
      newTris.push_back(new MTriangle (old_v0,old_v1, c[2]));
      newTris.push_back(new MTriangle (old_v1,old_v2, c[2]));
    if (std::find(cutVertices.begin(), cutVertices.end(), old_v0) != cutVertices.end()){
      newCut.insert(MEdge(c[2],old_v0));
    }
    else if (std::find(cutVertices.begin(), cutVertices.end(), old_v1) != cutVertices.end()) {
      newCut.insert(MEdge(c[2], old_v1));
    }
    else if (std::find(cutVertices.begin(), cutVertices.end(), old_v2) != cutVertices.end()){
      newCut.insert(MEdge(c[2], old_v2));
    }
  }
  else {
    newTris.push_back(tri);
    if (std::find(cutVertices.begin(), cutVertices.end(), old_v0) != cutVertices.end() &&
	std::find(cutVertices.begin(), cutVertices.end(), old_v1) != cutVertices.end())
      newCut.insert(MEdge(old_v0,old_v1));
    else if (std::find(cutVertices.begin(), cutVertices.end(), old_v1) != cutVertices.end() &&
	std::find(cutVertices.begin(), cutVertices.end(), old_v2) != cutVertices.end())
      newCut.insert(MEdge(old_v1,old_v2));
    else if (std::find(cutVertices.begin(), cutVertices.end(), old_v2) != cutVertices.end() &&
	std::find(cutVertices.begin(), cutVertices.end(), old_v0) != cutVertices.end())
      newCut.insert(MEdge(old_v2,old_v0));
  }
#endif

}





int numberOfVertices( std::vector<Branch> const& edges )
{
  // get the number of vertices in centerlines
  std::set<int> thePtIds;
  unsigned int nBranch = edges.size();
  for(unsigned int i = 0; i < nBranch; ++i)
  {
    std::vector<MLine*> mylines = edges[i].lines;
    for(unsigned int k = 0; k < mylines.size(); ++k)
    {
      MLine *l = mylines[k];
      MVertex *v0 = l->getVertex(0);
      MVertex *v1 = l->getVertex(1);
      thePtIds.insert(v0->getIndex());
      thePtIds.insert(v1->getIndex());
    }
  }
  return thePtIds.size();
}

int maxVerticesIndex( std::vector<Branch> const& edges )
{
  int nBranch = edges.size();
  int maxId = 0;
  for( int i = 0; i < nBranch; ++i )
  {
    std::vector<MLine*> mylines = edges[i].lines;
    for(unsigned int k = 0; k < mylines.size(); ++k)
    {
      MLine *l = mylines[k];
      MVertex *v0 = l->getVertex(0);
      MVertex *v1 = l->getVertex(1);
      maxId = std::max( maxId,std::max( v0->getIndex(),v1->getIndex() ) );
    }
  }
  return maxId;
}





AngioTkCenterline::AngioTkCenterline(std::string fileName): kdtree(0), kdtreeR(0)
{
  recombine = (CTX::instance()->mesh.recombineAll) || (CTX::instance()->mesh.recombine3DAll);
  nbPoints = 25;
  hLayer = 0.3;
  hSecondLayer = 0.3;
  nbElemLayer = 3;
  nbElemSecondLayer = 0;
  useGmshExecutable = false;
  is_cut = 0;
  is_closed = 0;
  is_extruded = 0;

  importFile(fileName);
  buildKdTree();
  update_needed = false;
}

AngioTkCenterline::AngioTkCenterline(): kdtree(0), kdtreeR(0)
{

  recombine = (CTX::instance()->mesh.recombineAll) || (CTX::instance()->mesh.recombine3DAll);
  fileName = "centerlines.vtk";//default
  nbPoints = 25;
  hLayer = 0.3;
  hSecondLayer = 0.3;
  nbElemLayer = 3;
  nbElemSecondLayer = 0;
  useGmshExecutable = false;
  is_cut = 0;
  is_closed = 0;
  is_extruded = 0;
  is_clip_mesh = 0;
  M_clipMeshScalingFactor=2.0;
  descInletOutlet = "";

  options["closeVolume"] = new FieldOptionInt
    (is_closed, "Action: Create In/Outlet planar faces");
  options["extrudeWall"] = new FieldOptionInt
    (is_extruded, "Action: Extrude wall");
  options["reMesh"] = new FieldOptionInt
    (is_cut, "Action: Cut the initial mesh in different mesh partitions using the "
     "centerlines");
  options["clipMesh"] = new FieldOptionInt
    (is_clip_mesh, "Action: Cut the initial mesh in different mesh partitions using thecenterlines");
  options["clipMeshScalingFactor"] = new FieldOptionDouble
    (M_clipMeshScalingFactor, "scaling factor (distance proportional to radius)");

  callbacks["run"] = new FieldCallbackGeneric<AngioTkCenterline>
    (this, &AngioTkCenterline::run, "Run actions (closeVolume, extrudeWall, cutMesh) \n");

  options["FileName"] = new FieldOptionString
    (fileName, "File name for the centerlines", &update_needed);
  options["nbPoints"] = new FieldOptionInt
    (nbPoints, "Number of mesh elements in a circle");
  options["nbElemLayer"] = new FieldOptionInt
    (nbElemLayer, "Number of mesh elements the extruded layer");
  options["hLayer"] = new FieldOptionDouble
    (hLayer, "Thickness (% of radius) of the extruded layer");
  options["nbElemSecondLayer"] = new FieldOptionInt
    (nbElemSecondLayer, "Number of mesh elements the second extruded layer");
  options["hSecondLayer"] = new FieldOptionDouble
    (hSecondLayer, "Thickness (% of radius) of the second extruded layer");

  //double * outlet0= new double(0);//outlet0CenterX
  //double SavehLayer = hLayer;
  options["descInletOutlet"] = new FieldOptionString
    (descInletOutlet, "descInletOutlet");


}

AngioTkCenterline::~AngioTkCenterline()
{

  // WARNING  memory leaks here (mod must be destroy but can be problematic in particular case)
  // -> TODO : need to introduce mod as std::shared_ptr
  //if (mod) delete mod;
  
  if(kdtree){
    ANNpointArray nodes = kdtree->thePoints();
    if(nodes) annDeallocPts(nodes);
    delete kdtree;
  }
  if(kdtreeR){
    ANNpointArray nodesR = kdtreeR->thePoints();
    if(nodesR) annDeallocPts(nodesR);
    delete kdtreeR;
  }
}

//#include <OpenFile.h>
void AngioTkCenterline::importSurfaceFromFile(std::string const& fileName )
{
  if ( fileName.empty() || !Feel::fs::exists( fileName ) )
       return;

  current = GModel::current();
  current->load(fileName);
  current->removeDuplicateMeshVertices(1.e-8);

  std::vector<GFace*> currentFaces(current->firstFace(), current->lastFace());
  for (unsigned int i = 0; i < currentFaces.size(); i++){
    GFace *gf = currentFaces[i];
     if (gf->geomType() == GEntity::DiscreteSurface){
     	for(unsigned int j = 0; j < gf->triangles.size(); j++)
     	  triangles.push_back(gf->triangles[j]);
	if (is_cut || is_clip_mesh ){
	  gf->triangles.clear();
	  gf->deleteVertexArrays();
	  current->remove(gf);
	}
     }
  }

  //std::cout << "number of triangle : " << triangles.size() << "\n";
#if 0
  if(triangles.empty()){
    Msg::Error("Current GModel has no triangles ...");
    return;
  }
#endif

}


void
AngioTkCenterline::updateMergeFromExtremities( AngioTkCenterline const& centerlinesMerged,
					       std::map<MVertex*,std::pair< std::vector<std::pair<MLine*,int> >, MVertex*> > & indexVertexReplaced )
{
  for ( auto const& extremityPair : centerlinesMerged.centerlinesExtremities() )
    {
      double ptToLocalize[3] = { extremityPair.first->x(),extremityPair.first->y(),extremityPair.first->z() };

      int branchId = extremityPair.second.first;
      int lineIdInBranch = extremityPair.second.second;
      std::vector<MLine*> mylines = centerlinesMerged.centerlinesBranch(branchId).lines;
      MLine* myline = mylines[lineIdInBranch];
      if ( centerlinesMerged.centerlinesRadiusl().find(myline) == centerlinesMerged.centerlinesRadiusl().end() )
	Msg::Error("radius not find for this line \n");
      double radius = centerlinesMerged.centerlinesRadiusl().find(myline)->second;

      auto ptFoundData = this->foundClosestPointInCenterlines( ptToLocalize );
      MVertex* pointFound = std::get<0>(ptFoundData);
      double dist = std::get<1>(ptFoundData);
      if ( dist > 2*(radius+this->minRadiusAtVertex(pointFound) ) )
	pointFound = NULL;

      if ( pointFound == NULL )
	continue;

      // if pointFound is already replaced, use replace point
      if ( indexVertexReplaced.find( pointFound ) != indexVertexReplaced.end() )
	pointFound = indexVertexReplaced.find( pointFound )->second.second;

      // define point to replace as the extremity point
      MVertex* _ptToReplaced = extremityPair.first;
      MLine* _lineWhichHasToReplaced = myline;
      int _idVertexInLine = ( _ptToReplaced == myline->getVertex(0) )? 0 : 1;
      if ( myline->getVertex(_idVertexInLine) != _ptToReplaced )
	Msg::Error("Error localToGlobalVertex");

      MVertex* _mylink = ( _ptToReplaced == myline->getVertex(0) )? myline->getVertex(1) : myline->getVertex(0);

      // search other point closer (close to the extremity point) 
      double curMinDist = dist/*[0]*/;
      bool isForwardSearch = ( lineIdInBranch == 0 );
      double lcTotal = 0;
      int indexSearch = 0;
      for ( int q=1; ( lcTotal < 4*radius && q < (mylines.size()-1) ) ; ++q )
	{
	  int newLineId = (isForwardSearch)?lineIdInBranch+q : lineIdInBranch-q;
	  MLine* mylineSearch = mylines[newLineId];
	  int vIdInLine = (isForwardSearch)?0 : 1;
	  MVertex *newV = mylineSearch->getVertex( vIdInLine );
	  if ( _mylink != newV ) Msg::Error("link not valid");
	  _mylink = mylineSearch->getVertex( (int)(vIdInLine+1)%2 );
	  double newdist = pointFound->distance(newV);
	  if ( newdist < curMinDist )
	    {
	      curMinDist = newdist;
	      _ptToReplaced = newV;
	      _lineWhichHasToReplaced = mylineSearch;
	      _idVertexInLine = vIdInLine;
	      indexSearch=q;
	    }
	  lcTotal += mylineSearch->getLength();
	}
      //std::cout << "curMinDist " << curMinDist << " VS dist[0] " << dist[0] << " indexSearch "<< indexSearch << "\n";
      // store lines which use this point 
      std::vector<std::pair<MLine*,int> > myvecLine;
      myvecLine.push_back( std::make_pair(_lineWhichHasToReplaced,_idVertexInLine) );
      if ( indexSearch > 0 )
	{
	  for ( int q=0; q<indexSearch ; ++q )
	    {
	      int newLineId = (isForwardSearch)?lineIdInBranch+q : lineIdInBranch-q;
	      MLine* mylineSearch = mylines[newLineId];
	      int v0Id = mylineSearch->getVertex(0)->getIndex();
	      int v1Id = mylineSearch->getVertex(1)->getIndex();
	      if ( mylineSearch->getVertex(0) == _ptToReplaced )
		myvecLine.push_back( std::make_pair(mylineSearch,0) );
	      if ( mylineSearch->getVertex(1) == _ptToReplaced )
		myvecLine.push_back( std::make_pair(mylineSearch,1) );
	    }
	}

      // insert new vertex to replace
      auto itFindVertexReplaced = indexVertexReplaced.find(_ptToReplaced);
      if ( itFindVertexReplaced == indexVertexReplaced.end() )
	indexVertexReplaced[_ptToReplaced] = std::make_pair( myvecLine ,pointFound );
      else
	{
	  MVertex* ptUsed = itFindVertexReplaced->second.second;
	  if ( pointFound != ptUsed )
	    Msg::Error("points must identical");
	  indexVertexReplaced[_ptToReplaced].first.insert( itFindVertexReplaced->second.first.end(), myvecLine.begin(), myvecLine.end() );
	}

    } // for ( auto const& extremityPair ... )
}




void AngioTkCenterline::importFile(std::string fileName)
{
  current = GModel::current();
  //current->removeDuplicateMeshVertices(1.e-8);
  
  //if ( !current ) current = new GModel();
  if ( triangles.empty() )
    {
      std::vector<GFace*> currentFaces(current->firstFace(), current->lastFace());
      for (unsigned int i = 0; i < currentFaces.size(); i++){
	GFace *gf = currentFaces[i];
	if (gf->geomType() == GEntity::DiscreteSurface){
	  for(unsigned int j = 0; j < gf->triangles.size(); j++)
	    triangles.push_back(gf->triangles[j]);
	  if (is_cut || is_clip_mesh ){
	    gf->triangles.clear();
	    gf->deleteVertexArrays();
	    current->remove(gf);
	  }
	}
      }
    }
  //std::cout << "number of triangle : " << triangles.size() << "\n";
#if 0
  if(triangles.empty()){
    Msg::Error("Current GModel has no triangles ...");
    return;
  }
#endif

  // if current centerlines if empty, create new one from file
  if ( modEdges.empty() )
    {
      this->createFromFile( fileName,"" );
      return;
    }

  // merge new centerlines file into current centerlines
  Msg::Info("AngioTkCenterline: importFile : merge centerlines start");

  std::shared_ptr<AngioTkCenterline> newCenterlines;
  //int previousMaxVertexIndex=0;

  newCenterlines.reset( new AngioTkCenterline );
  newCenterlines->importFile( fileName );

  this->attachAngioTkCenterline(newCenterlines);

  // copy previous mapping
  std::map<int,int> _previousMapVertexGmshIdToVtkId;
  _previousMapVertexGmshIdToVtkId.insert( this->M_mapVertexGmshIdToVtkId.begin(),this->M_mapVertexGmshIdToVtkId.end() );

  // search points which can be connected between centerlines
  std::map<MVertex*,std::pair< std::vector<std::pair<MLine*,int> >, MVertex*> > indexVertexReplaced;
  this->updateMergeFromExtremities(*newCenterlines,indexVertexReplaced);
  newCenterlines->updateMergeFromExtremities(*this,indexVertexReplaced);

  // search max tag already register
  int previousMaxTag = 0;
  for (unsigned int i = 0; i < modEdges.size(); i++){
    GEdge *ge = modEdges[i];
    previousMaxTag = std::max(previousMaxTag, ge->tag());
    //std::cout << "ge->tag() " << ge->tag() << "\n";
  }
  // max vertex index register in current centerlines
  int previousMaxVertexIndex = maxVerticesIndex( this->edges );

  // replace duplicate points between current and merged centerlines
  std::set<int> ptDone;
  int nDuplicatePt=0;
  auto itEdge=newCenterlines->mod->firstEdge();
  auto enEdge=newCenterlines->mod->lastEdge();
  for ( ; itEdge!=enEdge; ++itEdge )
    {
      (*itEdge)->setTag( previousMaxTag+(*itEdge)->tag() );
      GEdge *ge = *itEdge;
      for(unsigned int j = 0; j < ge->lines.size(); j++)
	{
	  MLine *l = ge->lines[j];
	  MVertex *v0 = l->getVertex(0);
	  MVertex *v1 = l->getVertex(1);
	  std::vector<MVertex*> myVerticesInLine = { v0,v1 };
	  for ( MVertex* myvertex : myVerticesInLine )
	    {
	      if ( ptDone.find(myvertex->getNum()) == ptDone.end() )
		{
		  ptDone.insert(myvertex->getNum());
		  // change vertex index (by shifting with previoux max index )
		  myvertex->setIndex( previousMaxVertexIndex + myvertex->getIndex() );

		  double ptToLocalize[3] = { myvertex->x(),myvertex->y(),myvertex->z() };
		  auto ptFoundData = this->foundClosestPointInCenterlines( ptToLocalize );
		  MVertex* pointFound = std::get<0>(ptFoundData);
		  double dist = std::get<1>(ptFoundData);
		  if ( dist < 1e-8 )
		    {
		      ++nDuplicatePt;
		      for ( MLine* mylineSearch : newCenterlines->lines )
			{
			  MVertex *v0S = mylineSearch->getVertex(0);
			  MVertex *v1S = mylineSearch->getVertex(1);
			  if ( v0S == myvertex || v1S == myvertex )
			    {
			      int locVertexId = (v0S == myvertex)? 0 : 1;
			      std::vector<std::pair<MLine*,int> > myvecLine(1,std::make_pair(mylineSearch,locVertexId));
			      if ( indexVertexReplaced.find(myvertex) == indexVertexReplaced.end() )
				{
				  indexVertexReplaced[myvertex] = std::make_pair( myvecLine ,pointFound );
				}
			      else indexVertexReplaced[myvertex].first.push_back( std::make_pair(mylineSearch,locVertexId) );
			    }
			}
		    } // if ( dist ... )
		}
	    } // for ( MVertex* myvertex ... )

	}
    }
  if ( nDuplicatePt>0 )
    Msg::Info("AngioTkCenterline: find and remove %d duplicate points in merging",nDuplicatePt);

  // detect if some lines in newCenterlines are duplicated (due to the fact that the 2 pts in line are replaced) 
  itEdge=newCenterlines->mod->firstEdge();
  for ( ; itEdge!=enEdge; ++itEdge )
    {
      (*itEdge)->setTag( previousMaxTag+(*itEdge)->tag() );
      GEdge *ge = *itEdge;
      for(unsigned int j = 0; j < ge->lines.size(); j++)
	{
	  MLine *l = ge->lines[j];
	  MVertex *v0 = l->getVertex(0);
	  MVertex *v1 = l->getVertex(1);
	  if ( indexVertexReplaced.find(v0) != indexVertexReplaced.end() &&
	       indexVertexReplaced.find(v1) != indexVertexReplaced.end() )
	    {
	      M_registerLinesDuplicatedToIgnore.insert( l );
	    }
	}
    }

  // replace points in line
  for ( auto const& idReplacePair : indexVertexReplaced )
    {
      MVertex* newPt = idReplacePair.second.second;
      for ( auto const& vecLine : idReplacePair.second.first)
	{
	  MLine* lineModified = vecLine.first;
	  int ptIdInLine = vecLine.second;
	  lineModified->setVertex(ptIdInLine,newPt);
#if 0
	  std::cout << "lineModified->getLength() "<< lineModified->getLength() << " with ptIdInLine " << ptIdInLine <<"\n";
	  std::cout << " initPt : " << idReplacePair.first->x() << "," << idReplacePair.first->y() << "," << idReplacePair.first->z()
		    << "with index " << idReplacePair.first->getIndex()<<"\n";
	  std::cout << " newPt; " << newPt->x() << "," << newPt->y() << "," << newPt->z() << "with index " << newPt->getIndex()<<"\n";
	  std::cout << " newline vertex0 " << lineModified->getVertex(0)->x() << "," << lineModified->getVertex(0)->y() << "," << lineModified->getVertex(0)->z() <<"\n"
		    << " newline vertex1 " << lineModified->getVertex(1)->x() << "," << lineModified->getVertex(1)->y() << "," << lineModified->getVertex(1)->z() <<"\n";
#endif
	}
    }

  //---------------------------------------------//

  //modEdges.insert(modEdges.end(),mod->firstEdge(), mod->lastEdge());
  modEdges.insert(modEdges.end(),newCenterlines->mod->firstEdge(), newCenterlines->mod->lastEdge());
  //std::cout << "new distance" << std::distance( newCenterlines->mod->firstEdge(), newCenterlines->mod->lastEdge())<<"\n";
  // update centerlines : data branch, extremities, kdtree,...
  this->updateCenterlinesForUse(modEdges);

  //---------------------------------------------------------------------------------------------//

  // in case of merging centerline : update M_centerlinesFieldsPointData
  if ( newCenterlines )
    {
      std::map<std::string,std::vector<std::vector<double> > >  newCenterlinesFieldsPointData;

      // delete datas when fieldname is not presence in all centerlines
      std::set<std::string> fieldsToErase;
      for (auto const& fieldsPointDataPair : this->M_centerlinesFieldsPointData )
	{
	  if ( newCenterlines->M_centerlinesFieldsPointData.find(fieldsPointDataPair.first) == newCenterlines->M_centerlinesFieldsPointData.end() )
	    fieldsToErase.insert( fieldsPointDataPair.first );
	}
      for ( std::string const fieldName : fieldsToErase )
	this->M_centerlinesFieldsPointData.erase(fieldName);

      for (auto const& fieldsPointDataPair : M_centerlinesFieldsPointData )
	{
	  newCenterlinesFieldsPointData[fieldsPointDataPair.first].clear();
	  newCenterlinesFieldsPointData[fieldsPointDataPair.first].resize(this->M_mapVertexVtkIdToGmshId.size());
	  //std::cout << "fieldsPointDataPair.first " << fieldsPointDataPair.first << " : " << newCenterlinesFieldsPointData[fieldsPointDataPair.first].size() << "\n";
	  for (int k=0;k<this->M_mapVertexVtkIdToGmshId.size();++k)
	    {
	      int gmshId = this->M_mapVertexVtkIdToGmshId[k];
	      if ( gmshId > previousMaxVertexIndex )
		{
		  auto itFindGmshIdInNew = newCenterlines->M_mapVertexGmshIdToVtkId.find( gmshId-previousMaxVertexIndex );
		  if ( itFindGmshIdInNew == newCenterlines->M_mapVertexGmshIdToVtkId.end() )
		    Msg::Error("Error in mapping 1\n");
		  int vtkId = itFindGmshIdInNew->second;
		  newCenterlinesFieldsPointData[fieldsPointDataPair.first][k] = newCenterlines->M_centerlinesFieldsPointData[fieldsPointDataPair.first][vtkId];
		}
	      else
		{
		  auto itFindGmshIdInPrevious = _previousMapVertexGmshIdToVtkId.find( gmshId );
		  if ( itFindGmshIdInPrevious == _previousMapVertexGmshIdToVtkId.end() )
		    Msg::Error("Error in mapping 2\n");
		  int vtkId = itFindGmshIdInPrevious->second;
		  newCenterlinesFieldsPointData[fieldsPointDataPair.first][k] = this->M_centerlinesFieldsPointData[fieldsPointDataPair.first][vtkId];
		}
	    }
	}
      M_centerlinesFieldsPointData.clear();
      M_centerlinesFieldsPointData = newCenterlinesFieldsPointData;
    }

  Msg::Info("AngioTkCenterline: finish importFile");

}

#if 0
#include <cstdio>
#if defined __cplusplus
extern "C" {
#endif
FILE *fmemopen(void *buf, size_t size, const char *mode);
#ifdef __cplusplus
}
#endif

//#include <feel/feelcore/fmemopen.h>

int gmsh_yyparse();
int gmsh_yylex();
void gmsh_yyflush();
int gmsh_yy_scan_string ( const char *str );
extern std::string gmsh_yyname;
extern int gmsh_yylineno;
extern FILE* gmsh_yyin;
extern int gmsh_yyerrorstate;
extern int gmsh_yyviewindex;
extern char *gmsh_yytext;
int
ParseGeoFromMemory( GModel* model, std::string const& name, std::string const& geo  )
{
    gmsh_yyname = name;
    gmsh_yylineno = 1;
    gmsh_yyerrorstate = 0;
    gmsh_yyviewindex = 0;
    gmsh_yyin = fmemopen( (void*)geo.c_str(), geo.size()*sizeof(char), "r");

    while(!feof(gmsh_yyin))
    {
        gmsh_yyparse();
        if(gmsh_yyerrorstate > 20)
        {
          if(gmsh_yyerrorstate != 999) // 999 is a volontary exit
            std::cerr << "Too many errors: aborting parser...\n";
          gmsh_yyflush();
          break;
        }
    }
    fclose(gmsh_yyin);

    int imported = model->importGEOInternals();
    return imported;
}

#endif

void AngioTkCenterline::createFromGeo( std::tuple< std::vector< std::pair<GPoint,double> >, std::vector<std::vector<int> > > const& geoDesc,
				       std::string const& outputSmoothGeoPath )
{
  auto const& ptDesc = std::get<0>( geoDesc );
  auto const& lineDesc = std::get<1>( geoDesc );
  if ( lineDesc.empty() ) return;

  std::ostringstream ostr;
  int curPtId = 1;
  int curLineId = 1;

  for ( int p=0 ; p<ptDesc.size() ; ++p )
    {
      GPoint const& gp = std::get<0>(ptDesc[p]);
      double pointMeshSize = std::get<1>(ptDesc[p]);
      ostr << "Point(" << p+1 << ")={"
	   << gp.x() << "," << gp.y() << "," << gp.z() << ","
	   << pointMeshSize <<"};\n"; 
    }

  for ( int branchId=0 ; branchId<lineDesc.size() ; ++branchId )
    {
      auto const& ptIdInBranch = lineDesc[branchId];
      int nPtInBranch = ptIdInBranch.size();
      if ( nPtInBranch < 2 )
	continue;
      ostr << "Spline(" << curLineId++ << ")={";
      for ( int k = 0 ; k < nPtInBranch ; ++k )
	{
	  ostr << ptIdInBranch[k]+1;
	  if ( k < (nPtInBranch-1) )
	    ostr << ",";
	}
      ostr << "};\n";
    }

  //std::cout << ostr.str();
#if 0
  mod.reset( new GModel() );
  mod->setFactory("Gmsh");

  ParseGeoFromMemory(mod.get(),"hola",ostr.str());
  mod->mesh(1);

  mod->removeDuplicateMeshVertices(1.e-8);
  // call indexMeshVertices allow to assign id on vertices from 1
  bool saveAll = true;//false;
  int numVerticesInCenterlinesGModel = mod->indexMeshVertices(saveAll);

  modEdges.insert(modEdges.end(),mod->firstEdge(), mod->lastEdge());
  //this->importSurfaceFromFile( inputSurfacePath );
  // update centerlines : data branch, extremities, kdtree,...
  this->updateCenterlinesForUse(modEdges);
#else
  std::ofstream fileWrited( outputSmoothGeoPath, std::ios::out | std::ios::trunc);
  fileWrited << ostr.str();
  fileWrited.close();
  this->createFromFile( outputSmoothGeoPath, "" );
#endif
  Msg::Info("AngioTkCenterline: createFromGeoCenterlinesFile : finish");
}

void AngioTkCenterline::createFromCenterlines( AngioTkCenterline const& inputCenterlines, std::string const& outputSmoothGeoPath,
					       double meshSizeUniform, double resampleGeoPointSpacing )
{
  if ( inputCenterlines.centerlinesBranch().empty() )
    return;

  auto const& inputBranchCollection = inputCenterlines.centerlinesBranch();
  int inputNBranch = inputBranchCollection.size();

  std::vector< std::pair<GPoint,double> > ptDesc;
  std::vector<std::vector<int> > lineDesc(inputNBranch);
  std::map<MVertex*,int> mapInputVertexToPointDescId;
  for (int i=0;i<inputNBranch;++i )
    {
      auto const& inputBranch = inputBranchCollection[i];
      double minRad = inputBranch.minRad;//+inputBranch[i].maxRad/2.
      double curLength = 0;
      MVertex* vB = inputBranch.vB;
      if ( mapInputVertexToPointDescId.find( vB ) == mapInputVertexToPointDescId.end() )
	{
	  ptDesc.push_back( std::make_pair( GPoint(vB->x(),vB->y(),vB->z()), meshSizeUniform ) );
	  mapInputVertexToPointDescId[ vB ] = ptDesc.size()-1;
	}
      lineDesc[i].push_back( mapInputVertexToPointDescId.find( vB )->second );
      for ( int k=0;k<(inputBranch.lines.size()-1);++k )
	{
	  MLine* curInputLine = inputBranch.lines[k];
	  curLength += curInputLine->getLength();
	  if ( curLength < resampleGeoPointSpacing )//4 )//resampleLengthScalingFactor*minRad )
	    continue;

	  MVertex* v1 = curInputLine->getVertex(1);
	  if ( mapInputVertexToPointDescId.find( v1 ) == mapInputVertexToPointDescId.end() )
	    {
	      ptDesc.push_back( std::make_pair( GPoint(v1->x(),v1->y(),v1->z()), meshSizeUniform ) );
	      mapInputVertexToPointDescId[ v1 ] = ptDesc.size()-1;
	    }
	  lineDesc[i].push_back( mapInputVertexToPointDescId.find( v1 )->second );
	  curLength = 0;
	}
      MVertex* vE = inputBranch.vE;
      if ( mapInputVertexToPointDescId.find( vE ) == mapInputVertexToPointDescId.end() )
	{
	  ptDesc.push_back( std::make_pair( GPoint(vE->x(),vE->y(),vE->z()), meshSizeUniform ) );
	  mapInputVertexToPointDescId[ vE ] = ptDesc.size()-1;
	}
      lineDesc[i].push_back( mapInputVertexToPointDescId.find( vE )->second );
    }
  this->createFromGeo( std::make_tuple(ptDesc,lineDesc), outputSmoothGeoPath );

  std::string fieldName = "RadiusMin";
  if ( inputCenterlines.hasField(fieldName) )
    {
      // init fields
      int numVertices = numberOfVertices(edges);
      M_centerlinesFieldsPointData[fieldName].resize( numVertices );
      std::vector<bool> ptDone(numVertices,false);
      for ( auto& ptPair : colorp )
	{
	  MVertex* myvertex = ptPair.first;
	  int vId = myvertex->getIndex();
	  int vVtkId = this->M_mapVertexGmshIdToVtkId[vId];
	  M_centerlinesFieldsPointData[fieldName][vVtkId] = { 1. };
	}

      // init localisation tool
      std::vector<MVertex*> verticesPosition;
      for ( auto const& vertexPair : mapInputVertexToPointDescId )
	{
	  MVertex* myvertex = vertexPair.first;
	  verticesPosition.push_back(new MVertex(myvertex->x(), myvertex->y(), myvertex->z(), NULL, myvertex->getIndex() ));
	}

      MVertexPositionSet pos(verticesPosition);
#if 0
      double tolerance = 1e-1;//1e-8;//1.0e-8;
      SBoundingBox3d bbox = mod->bounds();
      double lc = bbox.empty() ? 1. : norm(SVector3(bbox.max(), bbox.min()));
      //double eps = lc * tolerance;
#endif
      double eps = (2./3.)*meshSizeUniform;

      for (int i=0;i<edges.size();++i )
	{
	  auto const& outputBranch = edges[i];
	  if ( outputBranch.lines.empty() ) continue;
	  double curLength = 0;
	  MVertex* vStart = outputBranch.lines.front()->getVertex(0);
	  MVertex *vertexFindStart = pos.find(vStart->x(), vStart->y(), vStart->z(), eps,true);
	  if ( !vertexFindStart )
	    Msg::Error("not find begin vertex");
	  std::set<MVertex*> vertexDone;
	  vertexDone.insert( vStart );
	  std::vector<std::tuple<MVertex*,double> > vertexOrdered;
	  vertexOrdered.push_back( std::make_tuple(vStart,curLength) );

	  int inputVertexIdStart = vertexFindStart->getNum();
	  int inputVtkIdStart = inputCenterlines.mapVertexGmshIdToVtkId(inputVertexIdStart);
	  double inputRadiusStart = inputCenterlines.centerlinesFieldsPointData(fieldName,inputVtkIdStart)[0];

	  //MVertex* vEnd = theBranch.lines.front()->getVertex(1);
	  for ( int k=0;k<outputBranch.lines.size();++k )
	    {
	      MLine* myline = outputBranch.lines[k];
	      curLength += myline->getLength();

	      MVertex* v0 = myline->getVertex(0);
	      MVertex* v1 = myline->getVertex(1);
	      MVertex *vertexFindCur0 = NULL;
	      if ( vertexDone.find(v0) == vertexDone.end() )
		{
		  Msg::Warning("maybe wrong");
		  vertexFindCur0 = pos.find(v0->x(), v0->y(), v0->z(), eps,false);
		  vertexDone.insert(v0);
		  vertexOrdered.push_back( std::make_tuple(v0,curLength ) );
		}
	      MVertex *vertexFindCur1 = NULL;
	      if ( vertexDone.find(v1) == vertexDone.end() )
		{
		  vertexFindCur1 = pos.find(v1->x(), v1->y(), v1->z(), eps,false);
		  vertexDone.insert(v1);
		  vertexOrdered.push_back( std::make_tuple(v1,curLength) );
		}

	      if ( vertexFindCur0 && vertexFindCur1 )
		Msg::Error("impossible : something is wrong");

	      if ( vertexFindCur0 )
		Msg::Error("aieiei");

	      if ( vertexFindCur1 )
		{
		  //std::cout << "find in branch " << i << "for k=" << k << "/"<< outputBranch.lines.size() << " and " <<vertexOrdered.size() << "\n";
		  int inputVertexIdEnd = vertexFindCur1->getNum();
		  int inputVtkIdEnd = inputCenterlines.mapVertexGmshIdToVtkId(inputVertexIdEnd);
		  double inputRadiusEnd = inputCenterlines.centerlinesFieldsPointData(fieldName,inputVtkIdEnd)[0];

		  if ( vertexOrdered.size() < 2 ) 
		    Msg::Error("a branch must have at least 2 vertex");

		  for ( auto const&/*MVertex**/ outputVertex : vertexOrdered )
		    {
		      int outputVertexId = std::get<0>(outputVertex)->getIndex();
		      double outputLengthFromStart = std::get<1>(outputVertex);
		      double outputRadius = ( (inputRadiusEnd-inputRadiusStart)/curLength )*outputLengthFromStart + inputRadiusStart;
		      int outputVtkIdEnd = this->mapVertexGmshIdToVtkId(outputVertexId);
		      this->M_centerlinesFieldsPointData[fieldName][outputVtkIdEnd] = { outputRadius };
		    }

		  curLength=0.;
		  vertexOrdered.clear();
		  vStart = v1;
		  vertexOrdered.push_back( std::make_tuple(vStart,curLength) );
		  inputVertexIdStart = inputVertexIdEnd;
		  inputVtkIdStart = inputVtkIdEnd;
		  inputRadiusStart = inputRadiusEnd;
		}

	    }
	}
      }

}



void AngioTkCenterline::createFromGeoCenterlinesFile( std::string const& fileName, std::string const& inputSurfacePath )
{
  Msg::Info("AngioTkCenterline: createFromGeoCenterlinesFile : %s",fileName.c_str());

  std::vector<std::tuple<std::pair<GPoint,double>,std::pair<GPoint,double> > > geoLines;

  std::ifstream fileLoaded( fileName, std::ios::in);
  while ( !fileLoaded.eof() )
    {
      std::vector<double> pt(3);
      std::string typeEntity;
      double pointMeshSize = 0.;
      fileLoaded >> typeEntity;
      if ( fileLoaded.eof() ) break;

      if ( typeEntity == "LINE" )
	{
	  fileLoaded >> pt[0] >> pt[1] >> pt[2] >> pointMeshSize;
	  auto p1 = std::make_pair( GPoint(pt[0],pt[1],pt[2]),pointMeshSize );
	  fileLoaded >> pt[0] >> pt[1] >> pt[2] >> pointMeshSize;
	  auto p2 = std::make_pair( GPoint(pt[0],pt[1],pt[2]),pointMeshSize );
	  geoLines.push_back( std::make_tuple(p1,p2) );
	}
      else
	Msg::Error( "invalid typeEntity : %s",typeEntity.c_str());
    }
  fileLoaded.close();

  mod.reset( new GModel() );
  mod->setFactory("Gmsh");
      
  int curPhysicalTag = 1;
  for ( auto const& linePts : geoLines )
    {
      auto const& p1 =  std::get<0>(linePts);
      auto const& p2 =  std::get<1>(linePts);
      GPoint const& gp1 = p1.first;
      double pointMeshSize1 = p1.second;
      GPoint const& gp2 = p2.first;
      double pointMeshSize2 = p2.second;

      GVertex * v1 = mod->addVertex(gp1.x(), gp1.y(), gp1.z(), pointMeshSize1 );
      GVertex * v2 = mod->addVertex(gp2.x(), gp2.y(), gp2.z(), pointMeshSize2 );
      GEdge *myLine = mod->addLine(v1,v2);
      myLine->addPhysicalEntity(curPhysicalTag);
      ++curPhysicalTag;
    }
  mod->mesh(1);

  mod->removeDuplicateMeshVertices(1.e-8);
  // call indexMeshVertices allow to assign id on vertices from 1
  bool saveAll = false;//true;//false;
  int numVerticesInCenterlinesGModel = mod->indexMeshVertices(saveAll);

  modEdges.insert(modEdges.end(),mod->firstEdge(), mod->lastEdge());
  this->importSurfaceFromFile( inputSurfacePath );

  // update centerlines : data branch, extremities, kdtree,...
  this->updateCenterlinesForUse(modEdges);

  Msg::Info("AngioTkCenterline: createFromGeoCenterlinesFile : finish");
}

void AngioTkCenterline::createFromFile( std::string const& fileName, std::string const& inputSurfacePath )
{
  std::string inputExt = Feel::fs::path( fileName ).extension().string();
  if ( inputExt == ".data" )
    {
      this->createFromGeoCenterlinesFile( fileName, inputSurfacePath );
      return;
    }
  Msg::Info("AngioTkCenterline: createFromFile : %s",fileName.c_str());

  mod.reset( new GModel() );
  mod->load( fileName );
  if ( inputExt == ".geo" )
    mod->mesh(1);

  mod->removeDuplicateMeshVertices(1.e-8);
  // call indexMeshVertices allow to assign id on vertices from 1
  bool saveAll = true;//false;
  int numVerticesInCenterlinesGModel = mod->indexMeshVertices(saveAll);

  //Msg::Info("AngioTkCenterline: createFromGmshFile : nModEdge %d",int(std::distance( mod->firstEdge(), mod->lastEdge())) );

  modEdges.insert(modEdges.end(),mod->firstEdge(), mod->lastEdge());
  if ( !inputSurfacePath.empty() )
    this->importSurfaceFromFile( inputSurfacePath );

  // update centerlines : data branch, extremities, kdtree,...
  this->updateCenterlinesForUse(modEdges);

  if ( inputExt == ".vtk" )
    this->updateCenterlinesFieldsFromFile(fileName);

  Msg::Info("AngioTkCenterline: createFromFile : finish");
}

void AngioTkCenterline::updateCenterlinesForUse(std::vector<GEdge*> const& _modEdges)
{
  std::map<int,std::vector<MLine*> > _modEdgesConvert;
  for (int i = 0; i < modEdges.size(); i++)
    {
      GEdge *ge = modEdges[i];
      int geTag = ge->tag();
      for(unsigned int j = 0; j < ge->lines.size(); j++)
	{
	  MLine *myline = ge->lines[j];
	  _modEdgesConvert[geTag].push_back(myline);
	}
    }
  this->updateCenterlinesForUse(_modEdgesConvert);
}
void AngioTkCenterline::updateCenterlinesForUse(std::map<int,std::vector<MLine*> > const& _modEdges)
{
  // create coloring
  lines.clear();
  colorp.clear();
  colorl.clear();
  int maxN = 0.0;

  // import first build previous centerlines (with maybe some mod) and after add new centerlines
  std::set<std::pair<int,int> > registerLinesBuildFromPointPair;
  std::map<MVertex*,int> vertexToNumberUse;
  std::map<MVertex*,std::vector<std::pair<MLine*,int> > > vertexToLines;
  std::map<MVertex*,std::set<int> > checkJunctionVertex;

  std::set<int> allTag;
  for ( auto const& _modEdgesPair : _modEdges )
    {
      int geTag = _modEdgesPair.first;
      allTag.insert( geTag );
      auto const& linesVec = _modEdgesPair.second;
      //Msg::Info("AngioTkCenterline:updateForUse in tag %d (%d) ",geTag,int(linesVec.size()));
      for(unsigned int j = 0; j < linesVec.size(); j++)
	{
	  MLine *l = linesVec[j];
	  MVertex *v0 = l->getVertex(0); MVertex *v1 = l->getVertex(1);
	  int v0Id = v0->getIndex(), v1Id = v1->getIndex();

	  // after remove duplacte vertices, it is possible to have a wrong line with 2 same pts
	  if ( v0 == v1 )
	    {
	      //std::cout << "ignore a wrong line\n";
	      continue;
	    }
	  if ( M_registerLinesDuplicatedToIgnore.find( l ) != M_registerLinesDuplicatedToIgnore.end() )
	    continue;
	  if ( M_registerLinesToRemoveFromPointIdPairInModelEdge.find( std::make_pair(v0Id,v1Id ) ) != M_registerLinesToRemoveFromPointIdPairInModelEdge.end() ||
	       M_registerLinesToRemoveFromPointIdPairInModelEdge.find( std::make_pair(v1Id,v0Id ) ) != M_registerLinesToRemoveFromPointIdPairInModelEdge.end() )
	    {
	      //std::cout << "ignore a line\n";
	      continue;
	    }

	  if ( registerLinesBuildFromPointPair.find( std::make_pair(v0Id,v1Id) ) == registerLinesBuildFromPointPair.end() &&
	       registerLinesBuildFromPointPair.find( std::make_pair(v1Id,v0Id) ) == registerLinesBuildFromPointPair.end() )
	    {
	      lines.push_back(l);
	      colorl.insert(std::make_pair(l, geTag));
	      maxN = std::max(maxN, geTag);
	      registerLinesBuildFromPointPair.insert( std::make_pair(v1Id,v0Id) );
	      registerLinesBuildFromPointPair.insert( std::make_pair(v0Id,v1Id) );
	    }
	  std::map<MVertex*, int>::iterator it0 = colorp.find(v0);
	  if (it0 == colorp.end())
	    {
	      colorp.insert(std::make_pair(v0, geTag));
	      vertexToNumberUse[v0]=0;
	    }
	  else if ( it0->second != geTag )
	    {
	      checkJunctionVertex[v0].insert( geTag );
	      checkJunctionVertex[v0].insert( it0->second );
	    }
	  ++vertexToNumberUse[v0];
	  vertexToLines[v0].push_back(std::make_pair(l,0));

	  std::map<MVertex*, int>::iterator it1 = colorp.find(v1);
	  if (it1 == colorp.end())
	    {
	      colorp.insert(std::make_pair(v1, geTag));
	      vertexToNumberUse[v1]=0;
	    }
	  else if ( it1->second != geTag )
	    {
	      checkJunctionVertex[v1].insert( geTag );
	      checkJunctionVertex[v1].insert( it1->second );
	    }
	  ++vertexToNumberUse[v1];
	  vertexToLines[v1].push_back(std::make_pair(l,1));
	}
    }
#if 0
  std::cout << "lines.size() " << lines.size() << " colorl.size() " << colorl.size() << "\n" << "colorp.size() " << colorp.size() << "\n";
  std::cout << "TAG LIST :";
  for ( int mytag : allTag )
    std::cout << " " << mytag;
  std::cout << "\n";
#endif

  std::vector<std::set<int> > tagGroups;
  std::map<int,int> tagToChange;
  for ( auto const& checkJunctionVertexPair : checkJunctionVertex )
    {
      MVertex* vertex = checkJunctionVertexPair.first;
      auto const& tagSet = checkJunctionVertexPair.second;
      if ( tagSet.size() == 2 && vertexToNumberUse[vertex] == 2 )
	{
	  tagGroups.push_back( tagSet );
	}
    }

  int cptIt=0;
  while ( true )
    {
      std::map<int,int> fusionGroup;
      std::set<int> groupDone;
      for (int k1=0;k1<tagGroups.size();++k1)
	{
	  if ( groupDone.find(k1) != groupDone.end() ) continue;
	  for ( int j1 : tagGroups[k1] )
	    {
	      for (int k2=0;k2<tagGroups.size();++k2)
		{
		  if (k1==k2) continue;
		  if ( groupDone.find(k2) != groupDone.end() ) continue;
		  for ( int j2 : tagGroups[k2] )
		    {
		      if ( j1 == j2 )
			{
			  groupDone.insert(k1);
			  groupDone.insert(k2);
			  fusionGroup[k1]=k2;
			}
		      
		    }
		}

	    }
	}
      // finish
      if ( fusionGroup.empty() ) break;
      if ( ++cptIt > 10000 ) { Msg::Error("PROBLEM");break; }

      for ( auto const& fusionPair : fusionGroup )
	{
	  int k1 = fusionPair.first, k2 = fusionPair.second;
	  tagGroups[k1].insert( tagGroups[k2].begin(),tagGroups[k2].end() );
	  tagGroups[k2].clear();
	}
    }

      for (int k1=0;k1<tagGroups.size();++k1)
	{
	  if ( tagGroups[k1].empty() ) continue;
	  int newGrougTag = *tagGroups[k1].begin();//smallest id
	  for ( int j1 : tagGroups[k1] )
	    {
	      tagToChange[j1] = newGrougTag;
	    }
	}


  // update new tag
  std::set<int> tagPrintedDone;
  for ( auto & myl : colorl)
    {
      int tag = myl.second;
      if ( tagToChange.find(tag) != tagToChange.end() )
	{
	  if ( false )
	    if ( tagPrintedDone.find( myl.second ) == tagPrintedDone.end() ) {
	      std::cout << "CHANGE TAG " << myl.second << " -> " << tagToChange.find(tag)->second << "\n";
	      tagPrintedDone.insert(myl.second);
	    }
	  myl.second = tagToChange.find(tag)->second;

	  MVertex *v0 = myl.first->getVertex(0);
	  MVertex *v1 = myl.first->getVertex(1);
	  colorp[v0] = myl.second;
	  colorp[v1] = myl.second;

	}
    }
  // update also GEdge in modEdges
  for (int i = 0; i < modEdges.size(); i++)
    {
      GEdge *ge = modEdges[i];
      int geTag = ge->tag();
      auto itNewTag = tagToChange.find(geTag);
      if ( itNewTag != tagToChange.end() )
	{
	  ge->setTag( itNewTag->second );
	}
    }
      

  std::map<MLine *,bool> lineDone;

  std::set<int> pointReverseDone;
  std::set<int> lineReverseDone;

  std::set<int> tagReverseDone;
  bool finishAllBranch = false;
  while( !finishAllBranch )
    {
      // search a first line
      MVertex * currentPoint = NULL;
      int branchTag = -1;
      for ( auto const& mylinePair : colorl)
	{
	  MLine* myline = mylinePair.first;
	  if ( lineDone.find(myline) != lineDone.end() ) continue;
	  //if (  lineReverseDone.find(myline) != lineReverseDone.end() ) continue;
	  if ( tagReverseDone.find(mylinePair.second) != tagReverseDone.end() ) continue;
	  MVertex *v0 = myline->getVertex(0);
	  MVertex *v1 = myline->getVertex(1);
	  if ( vertexToNumberUse.find(v0)->second != 2 )
	      currentPoint = v1;
	  else if ( vertexToNumberUse.find(v1)->second != 2 )
	      currentPoint = v0;

	  if ( currentPoint != NULL )
	    {
 	      lineDone[myline]=true;
	      branchTag = mylinePair.second;
	      break;
	    }
	}

      if ( currentPoint == NULL )
	{
	  finishAllBranch = true;
	  break;
	}

      tagReverseDone.insert(branchTag);
      //std::cout << "Tag "<<branchTag<<" with "<< currentPoint->getIndex() << "\n" << std::flush;
      bool finishBranch = finishAllBranch;//false;
      while ( !finishBranch )
	{
	  if (  vertexToNumberUse.find(currentPoint) == vertexToNumberUse.end() ) Msg::Error("bad");
	  int nVertexUse = vertexToNumberUse.find(currentPoint)->second;
	  if ( nVertexUse != 2 )
	    {
	      finishBranch=true;
	      break;
	    }

	  //pointReverseDone.insert( currentPoint->getIndex() );

	  MLine * line1 = vertexToLines[currentPoint][0].first;

	  MLine * line2 = (nVertexUse>1)?vertexToLines[currentPoint][1].first : line1;
	  int colorl1 = colorl.find(line1)->second;
	  int colorl2 = colorl.find(line2)->second;
	  //if ( branchTag == 11 )
	  //std::cout<< "colorl1 " << colorl1 << " colorl2 " << colorl2 <<"\n";
	  if ( colorl1 != colorl2 ) { finishBranch=true;/*break;*/}
	  int idInLine1 = vertexToLines[currentPoint][0].second;
	  int idInLine2 = vertexToLines[currentPoint][1].second;

	  if ( false)//branchTag == 11 )
	    std::cout << "LINE1 " << line1->getVertex(0)->getIndex() << "," <<line1->getVertex(1)->getIndex() << "\n"
		      << "LINE2 " << line2->getVertex(0)->getIndex() << "," <<line2->getVertex(1)->getIndex() << "\n";
	  if ( colorl1 == branchTag && lineDone.find( line1 ) == lineDone.end() )
	    {
	      //std::cout << "OK1 with " << line1->getVertex( idInLine1 )->getIndex() << "\n"<<std::flush;
	      if ( lineDone.find( line2 ) == lineDone.end() ) Msg::Error("not allow1");
	      if ( idInLine1 == idInLine2 )
		{
		  line1->reverse();
		  for (auto & v2l : vertexToLines[line1->getVertex(0)] )
		    {
		      if ( v2l.first == line1 )
			v2l.second = 0;
		    }
		  for (auto & v2l : vertexToLines[line1->getVertex(1)] )
		    {
		      if ( v2l.first == line1 )
			v2l.second = 1;
		    }
		  currentPoint = line1->getVertex( idInLine1 );
		}
	      else
		currentPoint = line1->getVertex( (int)((idInLine1+1)%2) );
	      //if ( branchTag == 11 )
	      //std::cout << "OK1 with " << currentPoint->getIndex() << "\n"<<std::flush;

	      lineDone[line1]=true;
	    }
	  else if ( colorl2 == branchTag && lineDone.find( line2 ) == lineDone.end() )
	    {
	      //std::cout << "OK2 with " << line2->getVertex( idInLine2 )->getIndex() << "\n"<<std::flush;
	      if ( lineDone.find( line1 ) == lineDone.end() ) Msg::Error("not allow2");
	      if ( idInLine1 == idInLine2 )
		{
		  line2->reverse();
		  for (auto & v2l : vertexToLines[line2->getVertex(0)] )
		    {
		      if ( v2l.first == line2 )
			v2l.second = 0;
		    }
		  for (auto & v2l : vertexToLines[line2->getVertex(1)] )
		    {
		      if ( v2l.first == line2 )
			v2l.second = 1;
		    }
		  currentPoint = line2->getVertex( idInLine2 );
		}
	      else
		currentPoint = line2->getVertex( (int)((idInLine2+1)%2) );

	      //if ( branchTag == 11 )
	      //	std::cout << "OK2 with " << currentPoint->getIndex() << "\n"<<std::flush;

	      lineDone[line2]=true;
	    }
	  else {/*std::cout <<"BAD\n"; */Msg::Error("error in connectivity fix"); break;}//
	}
      //tagReverseDone.insert(branchTag);
    }

  //---------------------------------------------------------------------------------------------//
  // create centerlines branch
  this->createBranches(maxN);

  this->updateRelationMapVertex();

  this->checkCenterlinesConnectivity();

  this->buildKdTree();

  //this->cleanBranch();

}

void AngioTkCenterline::cleanBranch()
{
  Msg::Info("AngioTkCenterline: cleanBranch start");

  // check if a branch connected to an extremity has a length enough large
  // or not a small branch created after an approximate connection
  // else delete this branch
  if(triangles.empty() )
    return;

  std::set<int> branchIdsRemove;
  for(int i = 0; i < edges.size(); ++i)
    {
      if ( branchIdsRemove.find( i ) != branchIdsRemove.end() )
	continue;
      /*if ( edges[i].vB == edges[i].vE )
	continue;*/

      if ( this->centerlinesExtremities().find( edges[i].vB ) != this->centerlinesExtremities().end() ||
	   this->centerlinesExtremities().find( edges[i].vE ) != this->centerlinesExtremities().end() )
	{
	  //if ( 3*edges[i].length < (edges[i].minRad+edges[i].maxRad)/2. )
	  if ( edges[i].length < 2*edges[i].maxRad )
	    {
	      branchIdsRemove.insert(i);
	    }

	}
      else if ( M_junctionsVertex.find( edges[i].vB ) != M_junctionsVertex.end() &&
		M_junctionsVertex.find( edges[i].vE ) != M_junctionsVertex.end() )
	{
	  for(int i2 = 0; i2 < edges.size(); ++i2)
	    {
	      if ( i == i2 ) continue;
	      /*if ( edges[i2].vB == edges[i2].vE )
		continue;*/

	      if ( branchIdsRemove.find( i2 ) != branchIdsRemove.end() )
		continue;

	      if ( ( edges[i].vB == edges[i2].vB && edges[i].vE == edges[i2].vE ) ||
		   ( edges[i].vB == edges[i2].vE && edges[i].vE == edges[i2].vB ) )
		{
		  if ( edges[i2].length < 2*(edges[i2].minRad+edges[i2].maxRad)/2. )
		    {
		      branchIdsRemove.insert(i2);
		    }
		}
	    }
	}
    }

  // re-run updateCenterlinesForUse if there are branchs to remove
  if ( !branchIdsRemove.empty() )
    {
      for ( int branchId : branchIdsRemove )
	{
	  //std::cout << "find NEW remove : edges[i2].lines.size() " << edges[branchId].lines.size() << "\n";
	  std::vector<MLine*> mylinesToRemove = edges[branchId].lines;
	  for(int j = 0; j < mylinesToRemove.size(); ++j)
	    {
	      int v0Id = mylinesToRemove[j]->getVertex(0)->getIndex();
	      int v1Id = mylinesToRemove[j]->getVertex(1)->getIndex();
	      M_registerLinesToRemoveFromPointIdPairInModelEdge.insert( std::make_pair(v0Id,v1Id) );
	      M_registerLinesToRemoveFromPointIdPairInModelEdge.insert( std::make_pair(v1Id,v0Id) );
	    }
	}
      // copy previous mapping
      std::map<int,int> _previousMapVertexGmshIdToVtkId;
      _previousMapVertexGmshIdToVtkId.insert( this->M_mapVertexGmshIdToVtkId.begin(),this->M_mapVertexGmshIdToVtkId.end() );
      // update for use
      this->updateCenterlinesForUse(modEdges);
      // update fields data
      this->updateFieldsDataAfterReduction(_previousMapVertexGmshIdToVtkId);
    }

}

void AngioTkCenterline::createBranches(int maxN)
{
  Msg::Info("AngioTkCenterline: start createBranches (maxNin =%d)",maxN);

  //sort colored lines and create edges
  std::vector<std::vector<MLine*> > color_edges;
  color_edges.resize(maxN);
  std::multiset<MVertex*> allV;
  std::map<MLine*, int>::iterator itl = colorl.begin();
  while (itl != colorl.end()){
    int color = itl->second;
    if ( color > color_edges.size() ) Msg::Error( "invalid coloring %g but max is %g",color, color_edges.size());
    MLine* l = itl->first;
    allV.insert(l->getVertex(0));
    allV.insert(l->getVertex(1));
    color_edges[color-1].push_back(l);
    itl++;
  }
  Msg::Info("AngioTkCenterline: color_edges size = %d",(int)color_edges.size());

  //detect junctions
  junctions.clear();
  std::multiset<MVertex*>::iterator it = allV.begin();
  for ( ; it != allV.end(); ++it){
    if (allV.count(*it) != 2) {
      junctions.insert(*it);
    }
  }
  Msg::Info("AngioTkCenterline: junction+extrermity size = %d", junctions.size() );

  //split edges
  edges.clear();
  int tag = 0;
  for(unsigned int i = 0; i < color_edges.size(); ++i){
    std::vector<MLine*> lines = color_edges[i];
    while (!lines.empty()) {
      std::vector<MLine*> myLines;
      std::vector<MLine*>::iterator itl = lines.begin();
      MVertex *vB = (*itl)->getVertex(0);
      MVertex *vE = (*itl)->getVertex(1);
      myLines.push_back(*itl);
      erase(lines, *itl);
      itl = lines.begin();
      //std::cout << "i="<<i << "("<< color_edges.size() <<")" <<" lines.size() "<<lines.size()<<"\n";

      while ( lines.size() > 0 && itl != lines.end() &&  !( junctions.find(vE) != junctions.end() &&
							    junctions.find(vB) != junctions.end()) ) {
  	MVertex *v1 = (*itl)->getVertex(0);
  	MVertex *v2 = (*itl)->getVertex(1);
	bool goVE = (junctions.find(vE) == junctions.end()) ? true : false ;
	bool goVB = (junctions.find(vB) == junctions.end()) ? true : false;
      	if (v1 == vE && goVE){
    	  myLines.push_back(*itl);
  	  erase(lines, *itl);
	  itl = lines.begin();
  	  vE = v2;
  	}
  	else if ( v2 == vE && goVE){
    	  myLines.push_back(*itl);
  	  erase(lines, *itl);
	  itl = lines.begin();
  	  vE = v1;
  	}
  	else if ( v1 == vB && goVB){
    	  myLines.push_back(*itl);
  	  erase(lines, *itl);
	  itl = lines.begin();
  	  vB = v2;
  	}
  	else if ( v2 == vB && goVB){
   	  myLines.push_back(*itl);
	  erase(lines, *itl);
	  itl = lines.begin();
  	  vB = v1;
  	}
	else {itl++;}
      }
      if (vB == vE && false ) {
        Msg::Error("Begin and end points branch are the same \n");
	if ( myLines.back()->getVertex(0) == vE )
	  vE = myLines.back()->getVertex(1);
	else if ( myLines.back()->getVertex(1) == vE )
	  vE = myLines.back()->getVertex(0);
	else
	  Msg::Error("Not good");
	myLines.pop_back();
        //break;
      }

      orderMLines(myLines, vB, vE);
      std::vector<Branch> children;

      Branch newBranch ={ tag++, myLines, computeLength(myLines), vB, vE,
                          children, 1.e6, 0.0};
      edges.push_back(newBranch) ;
    }
  }

  Msg::Info("AngioTkCenterline: build %d branches in centerlines",(int)edges.size());



  //create children
  for(unsigned int i = 0; i < edges.size(); ++i) {
    //std::cout << "edges["<<i<<"].lines.size()" << edges[i].lines.size() << "\n";
    MVertex *vE = edges[i].vE;
    std::vector<Branch> myChildren;
    for (std::vector<Branch>::iterator it = edges.begin(); it != edges.end(); ++it){
      Branch myBranch = *it;
      if (myBranch.vB == vE) myChildren.push_back(myBranch);
    }
    edges[i].children = myChildren;
  }

  // compute junction and extremity points and branch bounds
  M_vertexToLinesId.clear();
  M_junctionsVertex.clear();
  M_extremityVertex.clear();
  double _boundMinX=0,_boundMaxX=0,_boundMinY=0,_boundMaxY=0,_boundMinZ=0,_boundMaxZ=0;
  for(unsigned int i = 0; i < edges.size(); i++){
    std::vector<MLine*> mylines = edges[i].lines;
    bool initBounds = false;
    for (unsigned int j= 0; j < mylines.size(); j++){
      MVertex *v1 = mylines[j]->getVertex(0);
      MVertex *v2 = mylines[j]->getVertex(1);
      M_vertexToLinesId[v1].insert( std::make_pair(i,j) );
      M_vertexToLinesId[v2].insert( std::make_pair(i,j) );

      if ( !initBounds ) {
	_boundMinX=v1->x();_boundMaxX=v1->x(),_boundMinY=v1->y(),_boundMaxY=v1->y(),_boundMinZ=v1->z(),_boundMaxZ=v1->z();
	initBounds=true;
      }
      else {
	_boundMinX=std::min(v1->x(),_boundMinX);_boundMaxX=std::max(v1->x(),_boundMaxX);
	_boundMinY=std::min(v1->y(),_boundMinY);_boundMaxY=std::max(v1->y(),_boundMaxY);
	_boundMinZ=std::min(v1->z(),_boundMinZ);_boundMaxZ=std::max(v1->z(),_boundMaxZ);
      }
      _boundMinX=std::min(v2->x(),_boundMinX);_boundMaxX=std::max(v2->x(),_boundMaxX);
      _boundMinY=std::min(v2->y(),_boundMinY);_boundMaxY=std::max(v2->y(),_boundMaxY);
      _boundMinZ=std::min(v2->z(),_boundMinZ);_boundMaxZ=std::max(v2->z(),_boundMaxZ);
    }
    edges[i].boundMinX = _boundMinX; edges[i].boundMaxX = _boundMaxX;
    edges[i].boundMinY = _boundMinY; edges[i].boundMaxY = _boundMaxY;
    edges[i].boundMinZ = _boundMinZ; edges[i].boundMaxZ = _boundMaxZ;
  }
  for ( auto const& nPtsPair : M_vertexToLinesId )
    {
      auto const& vertexSetInfo = nPtsPair.second;
      if ( vertexSetInfo.size() == 1 )
	M_extremityVertex[nPtsPair.first] = *(vertexSetInfo.begin());// std::make_pair(0,0);
      else if ( vertexSetInfo.size() > 2 )
	{
	  for ( std::pair<int,int> const& vertexInfo : vertexSetInfo )
	    M_junctionsVertex[nPtsPair.first].insert( vertexInfo.first );
	}
    }

  Msg::Info("AngioTkCenterline: number of branch = %d, extremity = %d, junction = %d",
            (int) edges.size(), (int)M_extremityVertex.size(), (int)M_junctionsVertex.size());
  if ( false )
    for(int i = 0; i < edges.size(); ++i) {
      Msg::Info("AngioTkCenterline: number of discrete lines in edge %d = %d",i,(int)edges[i].lines.size());
    }

  // compute measures between surface and centerlines
  if(!triangles.empty())
  {
    //compute radius
    distanceToSurface();

    computeRadii();

    //print for debug
    printSplit();
  }
  Msg::Info("AngioTkCenterline: createBranches done");

}



#if 0
void AngioTkCenterline::fixBranchConnectivity()
{

  Msg::Info("AngioTkCenterline: fixBranchConnectivity");

  std::map<int,int> tagMap;
  std::map<int,std::vector<int> > mapNewTagToBranchIds;
  for(int i = 0; i < edges.size(); i++)
    {
      //int geTag = i;
      if ( tagMap.find( i ) == tagMap.end() )
	{
	  tagMap[i]=i+1;
	  auto itInsert = mapNewTagToBranchIds[ tagMap[i] ].begin();
	  mapNewTagToBranchIds[ tagMap[i] ].insert( itInsert,i );

	  bool branchIsClosedAtJunction = false;
	  MVertex* currentVertexB = edges[i].vB;
	  MVertex* currentVertexE = edges[i].vE;

	  while ( !branchIsClosedAtJunction)
	    {
	      bool vBisJE = ( M_junctionsVertex.find( currentVertexB ) != M_junctionsVertex.end() ||
			      M_extremityVertex.find( currentVertexB ) != M_extremityVertex.end() );
	      bool vEisJE = ( M_junctionsVertex.find( currentVertexE ) != M_junctionsVertex.end() ||
			      M_extremityVertex.find( currentVertexE ) != M_extremityVertex.end() );
	      if ( vBisJE && vEisJE )
		branchIsClosedAtJunction = true;
	      else
		{
		  for(int ii = 0; ii < edges.size(); ii++)
		    {
		      if ( ii == i ) continue;
		      if ( !vBisJE )
			{
			  if ( currentVertexB == edges[ii].vB )
			    {
			      tagMap[ii] = tagMap[i];
			      auto itInsert = mapNewTagToBranchIds[ tagMap[ii] ].begin();
			      mapNewTagToBranchIds[ tagMap[ii] ].insert( itInsert,ii );
			      currentVertexB = edges[ii].vE;
			    }
			  else if ( currentVertexB == edges[ii].vE )
			    {
			      tagMap[ii] = tagMap[i];
			      auto itInsert = mapNewTagToBranchIds[ tagMap[ii] ].begin();
			      mapNewTagToBranchIds[ tagMap[ii] ].insert( itInsert,ii );
			      currentVertexB = edges[ii].vB;
			    }
			}
		      if ( !vEisJE )
			{
			  if ( currentVertexE == edges[ii].vB )
			    {
			      tagMap[ii] = tagMap[i];
			      auto itInsert = mapNewTagToBranchIds[ tagMap[ii] ].end();
			      mapNewTagToBranchIds[ tagMap[ii] ].insert( itInsert,ii );
			      currentVertexE = edges[ii].vE;
			    }
			  else if ( currentVertexE == edges[ii].vE )
			    {
			      tagMap[ii] = tagMap[i];
			      auto itInsert = mapNewTagToBranchIds[ tagMap[ii] ].end();
			      mapNewTagToBranchIds[ tagMap[ii] ].insert( itInsert,ii );
			      currentVertexE = edges[ii].vB;
			    }
			}

		    }
		}
	    }
	}
    }

  // add some information about vertex extremities
  std::vector<std::pair<int,bool> > vertexExtrimityInLines;
  for(int i = 0; i < edges.size(); i++)
    {
      MVertex* currentVertexB = edges[i].vB;
      MVertex* currentVertexE = edges[i].vE;

      int idPtInListB = -1;
      bool isOnFrontLinesB = false;
      if ( currentVertexB->getIndex() == edges[i].lines.front()->getVertex(0)->getIndex() )
	{
	  idPtInListB=0;
	  isOnFrontLinesB=true;
	}
      else if ( currentVertexB->getIndex() == edges[i].lines.front()->getVertex(1)->getIndex() )
	{
	  idPtInListB=1;
	  isOnFrontLinesB=true;
	}
      if ( idPtInListB == -1 ) Msg::Error("TODO");

      vertexExtrimityInLines.push_back(std::make_pair(idPtInListB,isOnFrontLinesB));
    }

  // define if the line ordering must be reverse (in a composite branch )
  std::set<int> branchIdToReverse;
  for ( auto const& mapNewTagToBranchIdsPair : mapNewTagToBranchIds )
    {
#if 0
      std::cout << "newTag " << mapNewTagToBranchIdsPair.first << " : ";
      for ( int k=0;k<mapNewTagToBranchIdsPair.second.size();++k )
	std::cout << " " << mapNewTagToBranchIdsPair.second[k];
      std::cout << "\n";
#endif
      int lastVertexIdE = -1,lastVertexIdB = -1;
      int lastVertexIdInLineE = -1, lastVertexIdInLineB = -1;
      bool lastBranchIsReverted = false;
      for ( int k=0;k<mapNewTagToBranchIdsPair.second.size();++k )
	{
	  int idBranch = mapNewTagToBranchIdsPair.second[k];
	  int vertexIdInLineB = vertexExtrimityInLines[idBranch].first;
	  int vertexIdInLineE = (int)((vertexIdInLineB+1)%2);
	  bool isOnFrontLine = vertexExtrimityInLines[idBranch].second;
	  MLine* myLineB = (isOnFrontLine)? edges[idBranch].lines.front() : edges[idBranch].lines.back();
	  int vertexIdB = myLineB->getVertex(vertexIdInLineB)->getIndex();
	  MLine* myLineE = (!isOnFrontLine)? edges[idBranch].lines.front() : edges[idBranch].lines.back();
	  int vertexIdE = myLineE->getVertex( vertexIdInLineE )->getIndex();
	  if ( lastVertexIdE > 0 )
	    {
	      if ( lastVertexIdE == vertexIdB ) 
		{
		  if ( (!lastBranchIsReverted && lastVertexIdInLineE == vertexIdInLineB ) ||
		       ( lastBranchIsReverted && lastVertexIdInLineE != vertexIdInLineB ) )
		    { branchIdToReverse.insert( idBranch ); lastBranchIsReverted=true; }
		  lastVertexIdInLineE = vertexIdInLineE;
		  lastVertexIdE = vertexIdE;
		}
	      else if ( lastVertexIdE == vertexIdE )
		{
		  if ( (!lastBranchIsReverted && lastVertexIdInLineE == vertexIdInLineE ) ||
		       ( lastBranchIsReverted && lastVertexIdInLineE != vertexIdInLineE ) )
		    { branchIdToReverse.insert( idBranch );lastBranchIsReverted=true; }
		  lastVertexIdInLineE = vertexIdInLineB;
		  lastVertexIdE = vertexIdB;
		}
	      else if ( lastVertexIdB ==  vertexIdE )
		{
		  if ( (!lastBranchIsReverted && lastVertexIdInLineB == vertexIdInLineE ) ||
		       ( lastBranchIsReverted && lastVertexIdInLineB != vertexIdInLineE ) )
		    { branchIdToReverse.insert( idBranch );lastBranchIsReverted=true; }
		  lastVertexIdInLineB = vertexIdInLineB;
		  lastVertexIdB = vertexIdB;
		}
	      else if ( lastVertexIdB ==  vertexIdB )
		{
		  if ( (!lastBranchIsReverted && lastVertexIdInLineB == vertexIdInLineB ) ||
		       ( lastBranchIsReverted && lastVertexIdInLineB != vertexIdInLineB ) )
		    { branchIdToReverse.insert( idBranch );lastBranchIsReverted=true; }
		  lastVertexIdInLineB = vertexIdInLineE;
		  lastVertexIdB = vertexIdE;
		}
	      else Msg::Error("invalid case");
	    }
	  else
	    {
	      lastVertexIdB = vertexIdB;
	      lastVertexIdE = vertexIdE;
	      lastVertexIdInLineB = vertexIdInLineB;
	      lastVertexIdInLineE = vertexIdInLineE;
	    }
	}

    }

  // reverse lines necessary to have each new branch with a full line ordered
  for ( int k : branchIdToReverse )
    {
      std::cout << "branchIdToReverse : " << k << "\n";
      for ( MLine* myline : edges[k].lines )
	{
	  myline->reverse();
	}
    }

#if 0
  if ( edges.size() == 40 )
  for(int i = 0; i < edges.size(); i++)
    {
      if ( i!=20 ) continue;
      std::cout << "REVERSE myline\n";
      for ( MLine* myline : edges[i].lines )
	{
	  myline->reverse();
	}
    }
#endif
}
#endif



void AngioTkCenterline::distanceToSurface()
{
  Msg::Info("AngioTkCenterline: computing distance to surface mesh ");

  //COMPUTE WITH REVERSE ANN TREE (SURFACE POINTS IN TREE)
  std::set<MVertex*> allVS;
  for(unsigned int j = 0; j < triangles.size(); j++)
    for(int k = 0; k<3; k++) allVS.insert(triangles[j]->getVertex(k));
  int nbSNodes = allVS.size();
  ANNpointArray nodesR = annAllocPts(nbSNodes, 3);
  vertices.resize(nbSNodes);
  int ind = 0;
  std::set<MVertex*>::iterator itp = allVS.begin();
  while (itp != allVS.end()){
    MVertex *v = *itp;
    nodesR[ind][0] = v->x();
    nodesR[ind][1] = v->y();
    nodesR[ind][2] = v->z();
    vertices[ind] = v;
    itp++; ind++;
  }
  kdtreeR = new ANNkd_tree(nodesR, nbSNodes, 3);

  for(unsigned int i = 0; i < lines.size(); i++){
    MLine *l = lines[i];
    MVertex *v1 = l->getVertex(0);
    MVertex *v2 = l->getVertex(1);
    double midp[3] = {0.5*(v1->x()+v2->x()), 0.5*(v1->y()+v1->y()),0.5*(v1->z()+v2->z())};
    ANNidx index[1];
    ANNdist dist[1];
    kdtreeR->annkSearch(midp, 1, index, dist);
    double minRad = sqrt(dist[0]);
    radiusl.insert(std::make_pair(lines[i], minRad));
    //if ( minRad <1e-3 ) std::cout << minRad << "\n";
  }

}

void AngioTkCenterline::computeRadii()
{
  Msg::Info("AngioTkCenterline: computeRadii");

  for(unsigned int i = 0; i < edges.size(); ++i) {
    std::vector<MLine*> lines = edges[i].lines;
    for(unsigned int j = 0; j < lines.size(); ++j) {
      MLine *l = lines[j];
      std::map<MLine*,double>::iterator itr = radiusl.find(l);
      if (itr != radiusl.end()){
	edges[i].minRad = std::min(itr->second, edges[i].minRad);
	edges[i].maxRad = std::max(itr->second, edges[i].maxRad);
      }
      else printf("ARGG line not found \n");
    }
  }

}

void AngioTkCenterline::buildKdTree()
{
  FILE * f = Fopen("myPOINTS.pos","w");
  fprintf(f, "View \"\"{\n");

  int nbPL = 3;  //10 points per line
  //int nbNodes  = (lines.size()+1) + (nbPL*lines.size());
  int nbNodes  = (colorp.size()) + (nbPL*lines.size());

  ANNpointArray nodes = annAllocPts(nbNodes, 3);
  int ind = 0;
  std::map<MVertex*, int>::iterator itp = colorp.begin();
  while (itp != colorp.end()){
     MVertex *v = itp->first;
     nodes[ind][0] = v->x();
     nodes[ind][1] = v->y();
     nodes[ind][2] = v->z();
     itp++; ind++;
  }
  for(unsigned int k = 0; k < lines.size(); ++k){
   MVertex *v0 = lines[k]->getVertex(0);
   MVertex *v1 = lines[k]->getVertex(1);
   SVector3 P0(v0->x(),v0->y(), v0->z());
   SVector3 P1(v1->x(),v1->y(), v1->z());
   for (int j= 1; j < nbPL+1; j++){
     double inc = (double)j/(double)(nbPL+1);
     SVector3 Pj = P0+inc*(P1-P0);
     nodes[ind][0] = Pj.x();
     nodes[ind][1] = Pj.y();
     nodes[ind][2] = Pj.z();
     ind++;
   }
 }

  if ( kdtree != NULL )
    delete kdtree;
 kdtree = new ANNkd_tree(nodes, nbNodes, 3);

 for(int i = 0; i < nbNodes; ++i){
   fprintf(f, "SP(%g,%g,%g){%g};\n",
	   nodes[i][0], nodes[i][1],nodes[i][2],1.0);
 }
 fprintf(f,"};\n");
 fclose(f);
}

#include <meshGEdge.h>
#include <meshGFace.h>

void AngioTkCenterline::createSplitCompounds()
{
  //number of discrete vertices, edges, faces and regions for the mesh
  NV = current->getMaxElementaryNumber(0);
  NE = current->getMaxElementaryNumber(1);
  NF = current->getMaxElementaryNumber(2);
  NR = current->getMaxElementaryNumber(3);

  // Remesh new faces (Compound Lines and Compound Surfaces)
  Msg::Info("AngioTkCenterline: creating split compounds ...");

  //Parametrize Compound Lines
  for (int i=0; i < NE; i++){
    std::vector<GEdge*>e_compound;
    GEdge *pe = current->getEdgeByTag(i+1);//current edge

    e_compound.push_back(pe);
    int num_gec = NE+i+1;
    Msg::Info("Create Compound Line (%d) with tag = %d discrete edge",num_gec, pe->tag());
    GEdge *gec =  current->addCompoundEdge(e_compound,num_gec);
    //Msg::Info("Create Compound Line (%d) tag = %d discrete edge",num_gec, gec->tag());

    if (CTX::instance()->mesh.algo2d != ALGO_2D_BAMG){
      gec->meshAttributes.method = MESH_TRANSFINITE;
      gec->meshAttributes.nbPointsTransfinite = nbPoints+1;
      gec->meshAttributes.typeTransfinite = 0;
      gec->meshAttributes.coeffTransfinite = 1.0;
    }

    if ( !useGmshExecutable )
    {
      // mesh the compound line
      meshGEdge mesher;
      mesher(gec);
    }
  }
  // Parametrize Compound surfaces
  std::list<GEdge*> U0;
#if 0
  for (int i=0; i < NF; i++){
    std::vector<GFace*> f_compound;
    GFace *pf =  current->getFaceByTag(i+1);//current face
    f_compound.push_back(pf);
    //NF = current->getMaxElementaryNumber(2);
    int num_gfc = NF+i+1;
    Msg::Info("Create Compound Surface (%d)     = %d discrete face", num_gfc, pf->tag());

    //1=conf_spectral 4=convex_circle, 7=conf_fe
    GFace *gfc = current->addCompoundFace(f_compound, 7, 0, num_gfc);
    //GFace *gfc = current->addCompoundFace(f_compound, 1, 0, num_gfc);
    // 4 marche pas mal
    //GFace *gfc = current->addCompoundFace(f_compound, 4, 0, num_gfc);
    //GFace *gfc = current->addCompoundFace(f_compound, 5, 0, num_gfc);

    Msg::Info("Create Compound Surface (%d) tag = %d discrete face",num_gfc, gfc->tag());

    gfc->meshAttributes.recombine = recombine;
    gfc->addPhysicalEntity(1);
    current->setPhysicalName("wall", 2, 1);//tag 1
  }
#else
  std::map<int,std::vector<std::tuple<GFace *,int> > > mapGFaceCompound;
  std::set<int> gfaceByTagDone;
  std::vector<int> algoListToUse = { 7, 4 };
  if ( useGmshExecutable )
    algoListToUse = { 7 };
  for ( int typeCompoundFace : algoListToUse )
    {
      int startFaceNum = current->getMaxElementaryNumber(2) + 1;
      for (int i=0; i < NF; i++)
	{
	  std::vector<GFace*> f_compound;
	  GFace *pf =  current->getFaceByTag(i+1);//current face
	  f_compound.push_back(pf);
	  int num_gfc = startFaceNum + i;
	  if ( typeCompoundFace == algoListToUse.front() )
	    Msg::Info("Create Compound Surface (%d) with tag = %d ", num_gfc, pf->tag());

	  //1=conf_spectral 4=convex_circle, 7=conf_fe
	  GFace *gfc = current->addCompoundFace(f_compound, typeCompoundFace/*7*/, 0, num_gfc);
	  //GFace *gfc = current->addCompoundFace(f_compound, 7, 0, num_gfc);
	  // 4 marche pas mal
	  //GFace *gfc = current->addCompoundFace(f_compound, 4, 0, num_gfc);

	  int genus = ((GFaceCompound*)gfc)->genusGeom();
	  if ( genus != 0 )
	    {
	      if ( typeCompoundFace == algoListToUse.front() )
		Msg::Error("genus is not equal to 0 (genus=%d) : ignore compound surface %d with tag = %d",genus,num_gfc,pf->tag());
	      current->remove( gfc );
	      continue;
	    }

	  //gfc->setMeshingAlgo(ALGO_2D_DELAUNAY);
	  gfc->meshAttributes.recombine = recombine;
	  gfc->addPhysicalEntity(1);


	  mapGFaceCompound[pf->tag()].push_back( std::make_tuple(gfc,typeCompoundFace) );
	}
    }
  current->setPhysicalName("wall", 2, 1);//tag 1

  if ( !useGmshExecutable )
    {
      for ( auto itgf : mapGFaceCompound )
	{
	  bool meshingWasSuccessful = false;
	  int idAlgoSuccess = -1;
	  for ( int k=0;k< itgf.second.size() && !meshingWasSuccessful; ++k )
	    {
	      auto itgfc = itgf.second[k];
	      GFace * gfc = std::get<0>(itgfc);
	      int typeCompoundFace = std::get<1>(itgfc);
	      /*//gfc->setMeshingAlgo(ALGO_2D_DELAUNAY);
	      gfc->meshAttributes.recombine = recombine;
	      gfc->addPhysicalEntity(1);*/

	      Msg::Info("run meshing of face with tag=%d with algo %d",gfc->tag(),typeCompoundFace);

	      //backgroundMesh::current()->unset();
#if 0
	      meshGFace mesher(true);
	      mesher(gfc);
#else
	      gfc->model()->setCurrentMeshEntity(gfc);
	      bool repairSelfIntersecting1dMesh=false/*true*/,onlyInitialMesh=false,debug=false;
	      meshingWasSuccessful = meshGenerator(gfc, 0, repairSelfIntersecting1dMesh, onlyInitialMesh,debug);

	      if ( meshingWasSuccessful )
		idAlgoSuccess = k;
#endif
	    } // for ( int k ...

	  if ( !meshingWasSuccessful )
	    Msg::Error("meshing fail : we ignore this partition");

	  for ( int k=0;k< itgf.second.size(); ++k )
	    {
	      if ( k==idAlgoSuccess ) continue;
	      GFace * gfc = std::get<0>( itgf.second[k] );
	      current->remove( gfc );
	    }
	} //for ( auto itgf ...
    }
#endif
    

}

void AngioTkCenterline::createFaces()
{
  std::vector<std::vector<MTriangle*> > faces;

  std::multimap<MEdge, MTriangle*, Less_Edge> e2e;
  for(unsigned int i = 0; i < triangles.size(); ++i)
    for(int j = 0; j < 3; j++)
      e2e.insert(std::make_pair(triangles[i]->getEdge(j), triangles[i]));
  while(!e2e.empty()){
    std::set<MTriangle*> group;
    std::set<MEdge, Less_Edge> touched;
    group.clear();
    touched.clear();
    std::multimap<MEdge, MTriangle*, Less_Edge>::iterator ite = e2e.begin();
    MEdge me = ite->first;
    while (theCut.find(me) != theCut.end()){
      ite++;
      me = ite->first;
    }
    recurConnectByMEdge(me,e2e, group, touched, theCut);
    std::vector<MTriangle*> temp;
    temp.insert(temp.begin(), group.begin(), group.end());
    faces.push_back(temp);
    Msg::Info("AngioTkCenterline: add remesh partition with %d triangle",(int)temp.size());

    for(std::set<MEdge, Less_Edge>::iterator it = touched.begin();
        it != touched.end(); ++it)
      e2e.erase(*it);
  }
  Msg::Info("AngioTkCenterline: action (cutMesh) has cut surface mesh in %d faces ",
            (int)faces.size());

  //create discFaces
  for(unsigned int i = 0; i < faces.size(); ++i){
    int numF = current->getMaxElementaryNumber(2) + 1;
    discreteFace *f = new discreteFace(current, numF);
    current->add(f);
    discFaces.push_back(f);
    std::set<MVertex*> myVertices;
    std::vector<MTriangle*> myFace = faces[i];
    for(unsigned int j= 0; j< myFace.size(); j++){
      MTriangle *t = myFace[j];
      f->triangles.push_back(t);
      for (int k= 0; k< 3; k++){
	MVertex *v = t->getVertex(k);
	myVertices.insert(v);
	v->setEntity(f);
      }
    }
    f->mesh_vertices.insert(f->mesh_vertices.begin(),
  			    myVertices.begin(), myVertices.end());
  }
}



void AngioTkCenterline::createFacesFromClip()
{
  std::vector<std::vector<MTriangle*> > faces;

  std::multimap<MEdge, MTriangle*, Less_Edge> e2e;
  for(unsigned int i = 0; i < triangles.size(); ++i)
    for(int j = 0; j < triangles[i]->getNumEdges() /*3*/; j++)
      e2e.insert(std::make_pair(triangles[i]->getEdge(j), triangles[i]));

  //std::cout << "initial e2e.size() " << e2e.size() << "with triangles.size() "<<triangles.size()<<"\n";

  while(!e2e.empty()){
    std::set<MTriangle*> group;
    std::set<MEdge, Less_Edge> touched;
    group.clear();
    touched.clear();
    //std::cout << "e2e.size() " << e2e.size() << "\n";
    std::multimap<MEdge, MTriangle*, Less_Edge>::iterator ite = e2e.begin();
    MEdge me = ite->first;
    while (theCut.find(me) != theCut.end()){
      ite++;
      //if ( ite == e2e.end() ) break;
      me = ite->first;
    }
    //if ( ite == e2e.end() ) continue;
    //std::cout << "before recurConnectByMEdge\n"<<std::flush;
    recurConnectByMEdge(me,e2e, group, touched, theCut);
    //std::cout << "touched.size() " << touched.size() << "\n"<<std::flush;
    std::set<int> groupTouchClipTag;
    for ( MEdge const& edgeTouch : touched )
      {
	auto itCut = theCutMap.begin();
	auto enCut = theCutMap.end();
	bool find = false;
	for ( ; (itCut != enCut) && (!find) ; ++itCut )
	  {
	    int tagCut = itCut->first;
	    if ( itCut->second.find(edgeTouch) != itCut->second.end() )
	      {
		find = true;
		groupTouchClipTag.insert( tagCut );
	      }
	  }
      }
    //std::cout << "groupTouchClipTag.size() " << groupTouchClipTag.size() << "\n";

    if ( groupTouchClipTag.size() != 1 )
      {
	std::vector<MTriangle*> temp;
	temp.insert(temp.begin(), group.begin(), group.end());
	faces.push_back(temp);
      }
    else if ( groupTouchClipTag.size() == 1 )
      {
#if 0
	// compute barycenter of group
	SPoint3 baryGroup(0., 0., 0.);
	for ( MTriangle* theTriangle : group )
	  {
	    SPoint3 theBary = theTriangle->barycenter();
	    baryGroup[0] += theBary.x();
	    baryGroup[1] += theBary.y();
	    baryGroup[2] += theBary.z();
	  }
	baryGroup[0] /= (double)group.size();
	baryGroup[1] /= (double)group.size();
	baryGroup[2] /= (double)group.size();
	std::cout << "baryGroup : " << baryGroup[0] << " " << baryGroup[1] << " " << baryGroup[2] << "\n";
#endif
      }
    for(std::set<MEdge, Less_Edge>::iterator it = touched.begin();
        it != touched.end(); ++it)
      e2e.erase(*it);
  }
  Msg::Info("AngioTkCenterline: action (cutMesh) has cut surface mesh in %d faces ",(int)faces.size());

  //create discFaces
  for(unsigned int i = 0; i < faces.size(); ++i){
    int numF = current->getMaxElementaryNumber(2) + 1;
    discreteFace *f = new discreteFace(current, numF);
    current->add(f);
    discFaces.push_back(f);
    std::set<MVertex*> myVertices;
    std::vector<MTriangle*> myFace = faces[i];
    for(unsigned int j= 0; j< myFace.size(); j++){
      MTriangle *t = myFace[j];
      f->triangles.push_back(t);
      for (int k= 0; k< 3; k++){
	MVertex *v = t->getVertex(k);
	myVertices.insert(v);
	v->setEntity(f);
      }
    }
    f->mesh_vertices.insert(f->mesh_vertices.begin(),
  			    myVertices.begin(), myVertices.end());
  }

}



void AngioTkCenterline::initPhysicalMarkerFromDescFile( std::vector<GEdge*> boundEdges )
{
  if ( !descInletOutlet.empty() )
   {
      Msg::Info("AngioTkCenterline: action (initPhysicalMarkerFromDescFile) use descInletOutlet %s ",descInletOutlet.c_str());
      //std::cout << "descInletOutlet : " << descInletOutlet << "\n";
      std::map<std::pair<std::string,std::string>,SPoint3> mapPhysicalMarkerToPointLoc;

      std::ifstream fichier(descInletOutlet.c_str(), std::ios::in);  // on ouvre le fichier en lecture
      if(fichier)  // si l'ouverture a réussi
        {       
	  // instructions
	  std::string bcType,physicalMarkerLumen,physicalMarkerArterialWall;double ptx,pty,ptz;
	  while ( !fichier.eof() )
	    {
	      fichier >> bcType;
	      if ( fichier.eof() ) break;
	      fichier >> physicalMarkerLumen >> physicalMarkerArterialWall >> ptx >> pty >> ptz;
	      //std::cout << "load physicalMarker " << physicalMarkerLumen << "," << physicalMarkerArterialWall << " :" << ptx << " " << pty << " " << ptz << "\n";
	      //mapPhysicalMarkerToPointLoc[std::make_pair(physicalMarkerLumen,physicalMarkerLumen+"ring")] = SPoint3( ptx, pty, ptz );
	      mapPhysicalMarkerToPointLoc[std::make_pair(physicalMarkerLumen,physicalMarkerArterialWall)] = SPoint3( ptx, pty, ptz );

	    }
	  fichier.close();  // on ferme le fichier
        }
        else  // sinon
	  Msg::Fatal("Impossible d'ouvrir le fichier %s ", descInletOutlet.c_str());

      //std::map<std::string,SPoint3>::iterator itPM = mapPhysicalMarkerToPointLoc.begin();
      //std::map<std::string,SPoint3>::iterator enPM = mapPhysicalMarkerToPointLoc.end();
      std::map<std::pair<std::string,std::string>,SPoint3>::iterator itPM = mapPhysicalMarkerToPointLoc.begin();
      std::map<std::pair<std::string,std::string>,SPoint3>::iterator enPM = mapPhysicalMarkerToPointLoc.end();
      for ( ; itPM != enPM ; ++itPM )
	{
	  double minDist = 1.e6;
	  unsigned int idxE = boundEdges.size()+10; 
	  SPoint3 pt = itPM->second;
	  //for each boundary edges search 
	  for (unsigned int i = 0; i<  boundEdges.size(); i++)
	    {
	      //TODO trouver le plus proche
	      GEdge * gec;
	      if(is_cut) gec = current->getEdgeByTag(NE+boundEdges[i]->tag());
	      else gec = current->getEdgeByTag(boundEdges[i]->tag());

	      GVertex *gv = gec->getBeginVertex();
	      SPoint3 ptTest(gv->x(), gv->y(), gv->z());
	      double dist = pt.distance(ptTest);
	      if ( dist < minDist )
		{
		  minDist=dist;
		  idxE = i;
		}
	    }
	  //mapEdgeIdToPhysicalMarker[idxE] = itPM->first;
	  //mapBoundEdgeIdToPhysicalMarker[boundEdges[idxE]->tag()] = itPM->first;
	  mapBoundEdgeIdToPhysicalMarkerLumen[boundEdges[idxE]->tag()] = itPM->first.first;
	  mapBoundEdgeIdToPhysicalMarkerArerialWall[boundEdges[idxE]->tag()] = itPM->first.second;
	}
   } // if ( !descInletOutlet.empty() )
}


void AngioTkCenterline::createClosedVolume(GEdge *gin, std::vector<GEdge*> boundEdges)
{

  current->setFactory("Gmsh");
  std::vector<std::vector<GFace *> > myFaceLoops;
  std::vector<GFace *> myFaces;

  std::map<std::string,int> mapMarkerToTag;

  for (unsigned int i = 0; i<  boundEdges.size(); i++){
    std::vector<std::vector<GEdge *> > myEdgeLoops;
    std::vector<GEdge *> myEdges;
    GEdge * gec;
    if(is_cut) gec = current->getEdgeByTag(NE+boundEdges[i]->tag());
    else gec = current->getEdgeByTag(boundEdges[i]->tag());
    myEdges.push_back(gec);
    myEdgeLoops.push_back(myEdges);
    GFace *newFace = current->addPlanarFace(myEdgeLoops);

    if ( mapBoundEdgeIdToPhysicalMarkerLumen.find( boundEdges[i]->tag() ) != mapBoundEdgeIdToPhysicalMarkerLumen.end() )
      {
	std::string markerFace = mapBoundEdgeIdToPhysicalMarkerLumen.find( boundEdges[i]->tag() )->second;
	std::string markerEdge = "edge_"+markerFace;

	int thetagFace = current->getMaxPhysicalNumber(-1)+1;
	if ( mapMarkerToTag.find(markerFace) != mapMarkerToTag.end() )
	  thetagFace = mapMarkerToTag.find(markerFace)->second;
	else
	    mapMarkerToTag[markerFace] = thetagFace;
	newFace->addPhysicalEntity( thetagFace );
	current->setPhysicalName(markerFace, 2, thetagFace);

	int thetagEdge = current->getMaxPhysicalNumber(-1)+1;
	gec->addPhysicalEntity( thetagEdge );
	current->setPhysicalName(markerEdge, 1, thetagEdge);
      }
    else
      {

	if (gin==boundEdges[i]) {
	  newFace->addPhysicalEntity(2);
	  current->setPhysicalName("inlet", 2, 2);//tag 2
	}
	else{
	  newFace->addPhysicalEntity(3);
	  current->setPhysicalName("outlets", 2, 3);//tag 3
	}
      }
    myFaces.push_back(newFace);
  }

  Msg::Info("AngioTkCenterline: action (closeVolume) has created %d in/out planar faces ",
            (int)boundEdges.size());

  for (int i = 0; i < NF; i++){
    GFace * gf;
    if(is_cut) gf = current->getFaceByTag(NF+i+1);
    else gf = current->getFaceByTag(i+1);
    myFaces.push_back(gf);
  }
  myFaceLoops.push_back(myFaces);
  GRegion *reg = current->addVolume(myFaceLoops);
  int thePhysicalTagVolume = current->getMaxPhysicalNumber(-1)+1;
  reg->addPhysicalEntity(thePhysicalTagVolume);
  current->setPhysicalName("lumenVolume", 3, thePhysicalTagVolume);

  Msg::Info("AngioTkCenterline: action (closeVolume) has created volume %d ", reg->tag());

}

void AngioTkCenterline::extrudeBoundaryLayerWall(GEdge* gin, std::vector<GEdge*> boundEdges)
{
  Msg::Info("AngioTkCenterline: extrude boundary layer wall (%d, %g%%R) ", nbElemLayer,  hLayer);

  //orient extrude direction outward
#if 1
  int dir = 0;
  MElement *e = current->getFaceByTag(1)->getMeshElement(0);
  SVector3 ne = e->getFace(0).normal();
  SVector3 ps(e->getVertex(0)->x(), e->getVertex(0)->y(), e->getVertex(0)->z());
  double xyz[3] = {ps.x(), ps.y(), ps.z()};
  ANNidx index[1];
  ANNdist dist[1];
  kdtree->annkSearch(xyz, 1, index, dist);
  ANNpointArray nodes = kdtree->thePoints();
  SVector3 pc(nodes[index[0]][0], nodes[index[0]][1], nodes[index[0]][2]);
  SVector3 nc = ps-pc;
  if (dot(ne,nc) < 0) dir = 1;
  if (dir == 1 && hLayer > 0 ) hLayer *= -1.0;
#endif
  //int shift = 0;
  //if(is_cut) shift = NE;
  for (int i= 0; i< NF; i++){
    GFace *gfc ;
    if (is_cut) gfc = current->getFaceByTag(NF+i+1);
    else gfc = current->getFaceByTag(i+1);
    current->setFactory("Gmsh");
#if 0
    int dir = 0;
    MElement *e = current->getFaceByTag(i+1)->getMeshElement(0);
    SVector3 ne = e->getFace(0).normal();
    SVector3 ps(e->getVertex(0)->x(), e->getVertex(0)->y(), e->getVertex(0)->z());
    double xyz[3] = {ps.x(), ps.y(), ps.z()};
    ANNidx index[1];
    ANNdist dist[1];
    kdtree->annkSearch(xyz, 1, index, dist);
    ANNpointArray nodes = kdtree->thePoints();
    SVector3 pc(nodes[index[0]][0], nodes[index[0]][1], nodes[index[0]][2]);
    SVector3 nc = ps-pc;
    if (dot(ne,nc) < 0) dir = 1;
    if (dir == 1 && hLayer > 0 ) hLayer *= -1.0;
#endif

    //view -5 to scale hLayer by radius in BoundaryLayers.cpp
    std::vector<GEntity*> extrudedE = current->extrudeBoundaryLayer
      (gfc, nbElemLayer,  hLayer, dir, -5);
    GFace *eFace = (GFace*) extrudedE[0];
    int currentPhysicalTag = current->getMaxPhysicalNumber(-1)+1;
    eFace->addPhysicalEntity(currentPhysicalTag);
    current->setPhysicalName("outerWall", 2, currentPhysicalTag);//dim 2
    GRegion *eRegion = (GRegion*) extrudedE[1];
    currentPhysicalTag = current->getMaxPhysicalNumber(-1)+1;
    eRegion->addPhysicalEntity(currentPhysicalTag);
    current->setPhysicalName("wallVolume", 3, currentPhysicalTag);//dim 3

    //if double extruded layer
    if (nbElemSecondLayer > 0){
      std::vector<GEntity*> extrudedESec = current->extrudeBoundaryLayer
      	(eFace, nbElemSecondLayer,  hSecondLayer, dir, -5);
      GFace *eFaceSec = (GFace*) extrudedESec[0];
      eFaceSec->addPhysicalEntity(9);                    //tag 9
      current->setPhysicalName("outerSecondWall", 2, 9);//dim 2 tag 9
      GRegion *eRegionSec = (GRegion*) extrudedESec[1];
      eRegionSec->addPhysicalEntity(10);             //tag 10
      current->setPhysicalName("wallVolume", 3, 10);//dim 3 tag 10
    }
    //end double extrusion


    std::map<std::string,int> mapMarkerToTag;
    for (unsigned int j = 2; j < extrudedE.size(); j++){
      GFace *elFace = (GFace*) extrudedE[j];
      std::list<GEdge*> l_edges = elFace->edges();
      for(std::list<GEdge*>::iterator it = l_edges.begin(); it != l_edges.end(); it++){
	GEdge *myEdge = *it;
	if (is_cut) myEdge = current->getEdgeByTag((*it)->tag()-NE);
	if( std::find(boundEdges.begin(), boundEdges.end(), myEdge) != boundEdges.end() ){

	  if ( mapBoundEdgeIdToPhysicalMarkerArerialWall.find( myEdge->tag() ) != mapBoundEdgeIdToPhysicalMarkerArerialWall.end() )
	    {
	      std::string markerFace = mapBoundEdgeIdToPhysicalMarkerArerialWall.find( myEdge->tag() )->second;
	      int thetagFace = current->getMaxPhysicalNumber(-1)+1;

	      if ( mapMarkerToTag.find(markerFace) != mapMarkerToTag.end() )
		thetagFace = mapMarkerToTag.find(markerFace)->second;
	      else
		  mapMarkerToTag[markerFace] = thetagFace;

	      elFace->addPhysicalEntity( thetagFace);
	      current->setPhysicalName(markerFace, 2, thetagFace);
	    }
	  else
	    {
	      if (myEdge==gin){
		elFace->addPhysicalEntity(7);
		current->setPhysicalName("inletRing", 2, 7);//tag 7
	      }
	      else{
		elFace->addPhysicalEntity(8);
		current->setPhysicalName("outletRings", 2, 8);//tag 8
	      }
	    }
	}
      }
    }
  }

}

void AngioTkCenterline::run()
{
  //std::cout << "hLayer " << hLayer << "\n";
  useGmshExecutable = true;
  double t1 = Cpu();
  if (update_needed){
    std::ifstream input;
    //std::string pattern = FixRelativePath(fileName, "./");
    //Msg::StatusBar(true, "Reading TEST '%s'...", pattern.c_str());
    //input.open(pattern.c_str());
    input.open(fileName.c_str());
    if(StatFile(fileName))
      Msg::Fatal("AngioTkCenterline file '%s' does not exist ", fileName.c_str());
    importFile(fileName);
    buildKdTree();
    update_needed = false;
  }

  if ( is_clip_mesh )
  {
    this->runClipMesh();
    return;
  }


  if (is_cut) this->runSurfaceRemesh("");/*cutMesh();*/
  else{
    GFace *gf = current->getFaceByTag(1);
    int currentPhysicalTag = current->getMaxPhysicalNumber(-1)+1;
    gf->addPhysicalEntity(currentPhysicalTag);
    current->setPhysicalName("wall", 2, currentPhysicalTag);
    current->createTopologyFromMesh();
    ((GFaceCompound*)gf)->coherenceNormals();
    NV = current->getMaxElementaryNumber(0);
    NE = current->getMaxElementaryNumber(1);
    NF = current->getMaxElementaryNumber(2);
    NR = current->getMaxElementaryNumber(3);
  }

  //identify the boundary edges by looping over all discreteFaces
  std::vector<GEdge*> boundEdges;
  double dist_inlet = 1.e6;
  GEdge *gin = NULL;
  //int cptInOut=0;
  for (int i= 0; i< NF; i++){
    GFace *gf = current->getFaceByTag(i+1);
    std::list<GEdge*> l_edges = gf->edges();
    for(std::list<GEdge*>::iterator it = l_edges.begin(); it != l_edges.end(); it++){
      std::vector<GEdge*>::iterator ite = std::find(boundEdges.begin(),
                                                    boundEdges.end(), *it);
      if (ite != boundEdges.end()) boundEdges.erase(ite);
      else boundEdges.push_back(*it);
      GVertex *gv = (*it)->getBeginVertex();
      //std::cout << "inletoutlet " << cptInOut++ << " = " << gv->x() << " " << gv->y() << " " << gv->z() << "\n";
      SPoint3 pt(gv->x(), gv->y(), gv->z());
      double dist = pt.distance(ptin);
      if(dist < dist_inlet){
	dist_inlet = dist;
	gin = *it;
      }
    }
  }

  if (is_closed || is_extruded)
    this->initPhysicalMarkerFromDescFile( boundEdges );
  if (is_closed)   createClosedVolume(gin, boundEdges);
  if (is_extruded) extrudeBoundaryLayerWall(gin, boundEdges);

  double t2 = Cpu();
  Msg::Info("AngioTkCenterline operators computed in %g (s) ",t2-t1);
}

void AngioTkCenterline::runClipMesh()
{
  Msg::Info("runClipMesh");
  //std::cout << "runClipMesh \n";
  // point, direction, radius
  double scalingClip = M_clipMeshScalingFactor;//2;//1;
  //std::vector<std::tuple<SVector3,SVector3,double> > cutDesc; 
  std::map< MVertex *, std::vector< std::tuple<SVector3,SVector3,double> > > cutDesc2; 

  for(unsigned int i = 0; i < edges.size(); ++i)
    {
      MVertex *vB = edges[i].vB;
      MVertex *vE = edges[i].vE;
      //std::cout << "edge " << i << " vB : " << vB->x() << " " << vB->y() << " " << vB->z() << "\n";
      //std::cout << "edge " << i << " vE : " << vE->x() << " " << vE->y() << " " << vE->z() << "\n";

      auto itExtremityB = this->centerlinesExtremities().find( vB );
      auto itExtremityE = this->centerlinesExtremities().find( vE );
      bool vBeginConnected = ( itExtremityB == this->centerlinesExtremities().end() );
      bool vEndConnected = ( itExtremityE == this->centerlinesExtremities().end() );

      if ( !vBeginConnected )
	{
	  int branchIdB = itExtremityB->second.first;
	  if ( branchIdB != i )
	    Msg::Error("error branchIdB !!");

	  auto lines = edges[i].lines;
	  if (lines.size() > 0)
	    {
	      MVertex *v1;
	      MVertex *v2;
#if 0
	      std::map<MLine*,double>::iterator itrLineB = radiusl.find(lines.front());
	      double radiusLineB = itrLineB->second;//0.8;
	      //double radiusLineB = edges[i].maxRad;
#else
	      double radiusLineB = 0;
	      int nAvergage = std::min((int)lines.size(),(int)50);
	      for ( int ls=0;ls<nAvergage ;++ls)
		{
		  radiusLineB += radiusl.find(lines[ls])->second;
		}
	      radiusLineB /= nAvergage;
	      radiusLineB = std::max( radiusLineB, radiusl.find(lines.front())->second );
#endif
	      MLine* curLine = lines.front();

	      double curLength = 0.;

	      bool extremityIsPoint0 = ( lines.front()->getVertex(0) == vB ) || (lines.front()->getVertex(0) == vE);
	      bool extremityIsPoint1 = ( lines.front()->getVertex(1) == vB ) || (lines.front()->getVertex(1) == vE);
	      if ( !extremityIsPoint0 && !extremityIsPoint1 )
		Msg::Error("error !!");

	      SVector3 normalPlanDir(0.,0.,0.);
	      double minLengthBeforeStart=scalingClip*radiusLineB;
	      double maxLengthForTest=std::max(minLengthBeforeStart,4*radiusLineB);
	      for ( int k=0 ; k<lines.size() && curLength <= maxLengthForTest/*scalingClip*radiusLineB*/ ; ++k )
		{
		  curLine = lines[k];
		  curLength += curLine->getLength();

		  if ( curLength < minLengthBeforeStart )
		    continue;

		  if ( extremityIsPoint0 )
		    {
		      v1 = curLine->getVertex(1);
		      v2 = curLine->getVertex(0);
		      normalPlanDir = v2->point() - v1->point();
		      if ( (k+1) <lines.size() )
			{
			  MLine* nextLine = lines[k+1];
			  MVertex* v1Next = nextLine->getVertex(1);
			  MVertex* v2Next = nextLine->getVertex(0);
			  SVector3 dir2 = v2Next->point() - v1Next->point();
			  normalPlanDir = (1./2.)*(normalPlanDir+dir2);
			}
		    }
		  else
		    {
		      Msg::Error("NOT ALLOW");
		      v1 = curLine->getVertex(0);
		      v2 = curLine->getVertex(1);
		    }
		}
	      MVertex *thept = v1;//v2;
	      std::map<MLine*,double>::iterator itr = radiusl.find(curLine);
	      double radius = 1.1*itr->second;

	      cutDesc2[vB].push_back( std::make_tuple( SVector3(thept->x(),thept->y(),thept->z()),
						       normalPlanDir,
						       radius ) );

	    }
	}
      if ( !vEndConnected )
	{
	  int branchIdE = itExtremityE->second.first;
	  if ( branchIdE != i )
	    Msg::Error("error branchIdE !!");
	  auto lines = edges[i].lines;
	  if (lines.size() > 0)
	    {
	      MVertex *v1;
	      MVertex *v2;
#if 0
	      std::map<MLine*,double>::iterator itrLineB = radiusl.find(lines.back());
	      double radiusLineB = itrLineB->second;//0.8;
	      //double radiusLineB = edges[i].maxRad;
#else
	      double radiusLineE = 0;
	      int nAvergage = std::min((int)lines.size(),(int)50);
		for ( int ls=0;ls<nAvergage ;++ls)
		{
		  radiusLineE += radiusl.find(lines[lines.size()-1-ls])->second;
		}
		radiusLineE /= nAvergage;
		radiusLineE = std::max( radiusLineE, radiusl.find(lines.back())->second );
		//radiusLineE = edges[branchIdE].maxRad; //NEW
#endif


	      MLine* curLine = lines.back();
	      double curLength = 0.;

	      bool extremityIsPoint0 = ( lines.back()->getVertex(0) == vB ) || (lines.back()->getVertex(0) == vE);
	      bool extremityIsPoint1 = ( lines.back()->getVertex(1) == vB ) || (lines.back()->getVertex(1) == vE);
	      if ( !extremityIsPoint0 && !extremityIsPoint1 )
		Msg::Error("error !!");

	      SVector3 normalPlanDir(0.,0.,0.);
	      double minLengthBeforeStart=scalingClip*radiusLineE;
	      double maxLengthForTest=std::max(minLengthBeforeStart,4*radiusLineE);
	      for ( int k=0 ; k<lines.size() && curLength <= maxLengthForTest/*scalingClip*radiusLineE*/ ; ++k )
		{
		  curLine = lines[lines.size()-1-k];
		  curLength += curLine->getLength();

		  if ( curLength < minLengthBeforeStart )
		    continue;

		  if ( extremityIsPoint1 ) 
		    {
		      v1 = curLine->getVertex(0);
		      v2 = curLine->getVertex(1);
		      normalPlanDir = v2->point() - v1->point();
		      if ( int( lines.size()-1-(k+1) ) >= 0 )
			{
			  MLine* nextLine = lines[( lines.size()-1-(k+1))];
			  MVertex* v1Next = nextLine->getVertex(0);
			  MVertex* v2Next = nextLine->getVertex(1);
			  SVector3 dir2 = v2Next->point() - v1Next->point();
			  normalPlanDir = (1./2.)*(normalPlanDir+dir2);
			}
		    }
		  else
		    {
		      Msg::Error("NOT ALLOW");
		      v1 = curLine->getVertex(1);
		      v2 = curLine->getVertex(0);
		    }

		  MVertex *thept = v1;//v2;
		  std::map<MLine*,double>::iterator itr = radiusl.find(curLine);
		  double radius = 1.1*itr->second;
		  cutDesc2[vE].push_back( std::make_tuple( SVector3(thept->x(),thept->y(),thept->z()),
							   normalPlanDir, radius ) );

		}
	    }

	}

    }

  Msg::Info("AngioTkCenterline: prepare data for cutting done");


  std::ofstream fileWrited( "sphereremovebranch.data", std::ios::out | std::ios::trunc);
  int cptRealCut=0;
#if 0
  for ( int k=0;k<cutDesc.size(); ++k)
    {
      //if ( k==4 ) continue;
      auto & pt =  std::get<0>(cutDesc[k]);
      auto & dir =  std::get<1>(cutDesc[k]);
      double radius =  std::get<2>(cutDesc[k]);
      fileWrited << 0 << " " << pt[0] << " " << pt[1] << " " << pt[2] << " " << radius << "\n";
#if 0
      std::cout << "pt : " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
      std::cout << "dir : " << dir[0] << " " << dir[1] << " " << dir[2] << "\n";
      std::cout << "radius " << radius << "\n";
#endif
      bool cutSucess = cutByDisk(pt, dir, radius, k+1);
      if (!cutSucess)
	    Msg::Error("cut by disk fail");
      else
	++cptRealCut;
    }
#else
  //for ( int k=0;k<cutDesc.size(); ++k)
  for ( auto const& cutDescPair : cutDesc2 )
    {
      MVertex* vExtremity = cutDescPair.first;
      bool cutSucess = false;
      for ( auto const& cutDesc : cutDescPair.second )
	{
	  //if ( k==4 ) continue;
	  SVector3/*auto const&*/ pt =  std::get<0>(cutDesc);
	  SVector3/*auto const&*/ dir =  std::get<1>(cutDesc);
	  double radius =  std::get<2>(cutDesc);
	  fileWrited << 0 << " " << pt[0] << " " << pt[1] << " " << pt[2] << " " << radius << "\n";
	  cutSucess = cutByDisk(pt, dir, radius, cptRealCut+1);
	  if (cutSucess)
	    {
	      ++cptRealCut;
	      break;
	    }
	}
      if (!cutSucess)
	Msg::Error("open surface fail for the extermity pt : %g,%g,%g",vExtremity->x(),vExtremity->y(),vExtremity->z());
    }
#endif
  fileWrited.close();
  Msg::Info("AngioTkCenterline: all cuts done for open the surface (nb cut applied %d, planned %d)", cptRealCut,cutDesc2.size());

  //create discreteFaces
  //createFaces();
  this->createFacesFromClip();
  current->createTopologyFromFaces(discFaces);
  current->exportDiscreteGEOInternals();

#if 0
  //write
  Msg::Info("AngioTkCenterline: writing splitted mesh 'myCLIPPARTS.msh'");
  current->writeMSH("myCLIPPARTS.msh", 2.2, false, false);
#endif

}



void AngioTkCenterline::cutMesh(std::string const& remeshPartitionMeshFile)
{
  Msg::Info("AngioTkCenterline: action (cutMesh) splits surface mesh (%d tris) using %s ",
            triangles.size(), fileName.c_str());

  // i->j->(pt,radius)
  //std::vector< std::map<int, std::pair<SVector3,double> > > cutDiskToPerform(edges.size());
  std::vector< std::map<int, std::tuple<SVector3,double,double> > > cutDiskToPerform(edges.size());//pt,radius,lengthFromBranchBegin
  // first pass
  for(unsigned int i = 0; i < edges.size(); i++){
    std::vector<MLine*> lines = edges[i].lines;
    if ( lines.empty() ) continue;
    double L = edges[i].length;
    //double D = 2.*edges[i].minRad;  //(edges[i].minRad+edges[i].maxRad);
    double D = (edges[i].minRad+edges[i].maxRad);
    double AR = L/D;

    bool vBisJunc = M_junctionsVertex.find(edges[i].vB) != M_junctionsVertex.end();
    bool vEisJunc = M_junctionsVertex.find(edges[i].vE) != M_junctionsVertex.end();

    bool vBisInOut = M_extremityVertex.find(edges[i].vB) != M_extremityVertex.end();
    bool vEisInOut = M_extremityVertex.find(edges[i].vE) != M_extremityVertex.end();

    double radiusB = radiusl.find(lines.front())->second;
    double radiusE = radiusl.find(lines.back())->second;

#define USE_NEW_ANGIOTK_DEV 1


    // printf("*** AngioTkCenterline branch %d (AR=%.1f) \n", edges[i].tag, AR);
#if 0 //VINCENT
    int nbSplit = (int)ceil(AR/2 + 1.1); //AR/2 + 0.9
#else
    //int nbSplit = 4*(int)ceil(AR + 1.1); //AR/2 + 0.9
    int nbSplit = (int)ceil(AR + 1.1); //AR/2 + 0.9
    nbSplit = (int)ceil(1.5*AR + 1.1); // celui-ci a marcher
    nbSplit = (int)ceil(AR + 1.1);
    nbSplit = std::max( nbSplit, 3 );//3;//std::min( nbSplit, 3 );
#endif
    //int nbSplit = (int)ceil(AR/4 + 1.1); //AR/2 + 0.9
    //nbSplit = std::min( nbSplit, 8 );
    if( nbSplit > 1 ){
      double li  = L/nbSplit;
      if ( false ) Msg::Info("->> prepare cut branch %i in %d parts (L=%f,D=%f,AR=%f,li=%f,radiusB=%f,radiusE=%f)", i, nbSplit, L,D,AR,li,radiusB,radiusE);

#if !USE_NEW_ANGIOTK_DEV
      double lc = 0.0, lcTotal = 0.0;
#else
      double lc = lines[0]->getLength()/2., lcTotal = lines[0]->getLength()/2.;
#endif
      for (unsigned int j= 0; j < lines.size(); j++){
#if !USE_NEW_ANGIOTK_DEV
	lc += lines[j]->getLength();
#endif
	MVertex *v1 = lines[j]->getVertex(0);
	MVertex *v2 = lines[j]->getVertex(1);

	bool applyCut=false;

	if ( lc > li ) {
	  //SVector3 pt(v1->x(), v1->y(), v1->z());
	  SVector3 pt((v1->x()+v2->x())/2., (v1->y()+v2->y())/2., (v1->z()+v2->z())/2.);//new
	  SVector3 dir(v2->x()-v1->x(),v2->y()-v1->y(),v2->z()-v1->z());
	  std::map<MLine*,double>::iterator itr = radiusl.find(lines[j]);
	  double radius= itr->second;
	  //cutByDisk(pt, dir, itr->second);

	  /*bool*/ applyCut=true;

#if 0 //!USE_NEW_ANGIOTK_DEV
	  if ( vBisInOut && ( lcTotal < 2*radiusB ) )
	    applyCut=false;
	  if ( vEisInOut && ( lcTotal > (L-2*radiusE ) ) )
	    applyCut=false;
#else
	  double lenghVBjunc = (vBisInOut)? std::max(2*radiusB,li) : std::min(2*radiusB,li);
	  double lenghVEjunc = (vEisInOut)? std::max(2*radiusE,li) : std::min(2*radiusE,li);
	  if ( applyCut && lcTotal < lenghVBjunc )
	    applyCut=false;
	  if ( applyCut && lcTotal > (L-lenghVEjunc) )
	    applyCut=false;
#endif

	  for(unsigned int ii = 0; ii < edges.size(); ii++){
	    auto itj = cutDiskToPerform[ii].begin();
	    auto enj = cutDiskToPerform[ii].end();
	    for (; itj!=enj && applyCut;++itj)
	      {
		int jTest= itj->first;
		SVector3 ptTest = std::get<0>(itj->second);//itj->second.first;
		double radiusTest= std::get<1>(itj->second);//itj->second.second;
		double lengthFromBranchBeginTest= std::get<2>(itj->second);
		double distBetweenCenter = std::sqrt( std::pow( pt.x()-ptTest.x(),2)+std::pow( pt.y()-ptTest.y(),2)+std::pow( pt.z()-ptTest.z(),2) );

		// ignore cut done quite far of branch extremities 
#if !USE_NEW_ANGIOTK_DEV
		if ( ii != i &&
		    lengthFromBranchBeginTest > 3*edges[ii].maxRad &&
		    lengthFromBranchBeginTest < (edges[ii].length -3*edges[ii].maxRad) )
		  {
		    continue;
		  }
#elif 0
		if ( ii != i &&
	            lengthFromBranchBeginTest > 3*edges[ii].maxRad &&
		    lengthFromBranchBeginTest < (edges[ii].length -3*edges[ii].maxRad) )
		  {
		    continue;
		  }
#endif

		MVertex *v1Test = lines[j]->getVertex(0);
		MVertex *v2Test = lines[j]->getVertex(1);
		SVector3 dirTest(v2Test->x()-v1Test->x(),v2Test->y()-v1Test->y(),v2Test->z()-v1Test->z());


		double theScaling = 1;

#if USE_NEW_ANGIOTK_DEV
		if ( ii == i )
#endif
		  {
		if (norm(dir-dirTest) < 1e-2 )
		  theScaling = 1./2;
		if ( lc > 4*li )
		  theScaling = 1./2;
	}

#if 1
		if ( ii == i )
		  {
		    double theScalingCoarse = 2;
		    double lengthBeforeJunctionB = 4*radiusB;//4*radiusB;//2;//4;//edges[i].maxRad
		    double lengthBeforeJunctionE = 4*radiusE;
		    if ( !vBisJunc && vEisJunc && lcTotal < (L-lengthBeforeJunctionE) )
		      theScaling = theScalingCoarse;
		    if ( vBisJunc && !vEisJunc && lcTotal > lengthBeforeJunctionB )
		      theScaling = theScalingCoarse;
		    if ( vBisJunc && vEisJunc && lcTotal > lengthBeforeJunctionB && lcTotal < (L-lengthBeforeJunctionE) )
		      theScaling = theScalingCoarse;
		    if ( !vBisJunc && !vEisJunc )
		      theScaling = theScalingCoarse;
		  }
#endif
		double epsDist = (radius+radiusTest)/100.;
		if ( distBetweenCenter < (radius+radiusTest+epsDist)*theScaling ) //if ( distBetweenCenter < (radius+radiusTest) )
		  applyCut=false;
	      }
	} // for ( ii
	  if ( applyCut )
	    {
	      //cutDiskToPerform[i][j] = std::make_pair(pt,itr->second );
	      cutDiskToPerform[i][j] = std::make_tuple(pt,radius/*itr->second*/,lcTotal );
	      lc = 0.0;
	    }
	  //nbSplit--;
	  //lc = 0.0;
      } // if ( lc > li && nbSplit > 1)
#if USE_NEW_ANGIOTK_DEV
	  if ( !applyCut )
	    {
	  lc += lines[j]->getLength()/2.;
	  if ( (j+1) < lines.size() )
	    lc += lines[j+1]->getLength()/2.;
	}
	  lcTotal += lines[j]->getLength()/2.;
	  if ( (j+1) < lines.size() )
	    lcTotal += lines[j+1]->getLength()/2.;

#else
	lcTotal += lines[j]->getLength();
#endif
         } // for (unsigned int j= 0; j < lines.size(); j++)
     }
   } // end first for ( i

  for(unsigned int i = 0; i < edges.size(); i++){
    //if (i==2) continue;
    std::vector<MLine*> lines = edges[i].lines;
    auto/*std::map<int, std::pair<SVector3,double> >::iterator*/ itj = cutDiskToPerform[i].begin();
    auto/*std::map<int, std::pair<SVector3,double> >::iterator*/ enj = cutDiskToPerform[i].end();
    int cptRealCut=0;
    for (; itj!=enj;++itj)
      {
	int j = itj->first;
	SVector3 pt = std::get<0>(itj->second);
	double radius= std::get<1>(itj->second);
	MVertex *v1 = lines[j]->getVertex(0);
	MVertex *v2 = lines[j]->getVertex(1);
	SVector3 dir(v2->x()-v1->x(),v2->y()-v1->y(),v2->z()-v1->z());
	bool cutHasSucess = cutByDisk(pt, dir, radius);
	if ( cutHasSucess )
	  ++cptRealCut;
#if 0
	if ( cutHasSucess )
	  std::cout << "cut at pt " << pt.x() << "," << pt.y() << "," << pt.z() << " radius " << radius << "\n";
	else
	  std::cout << "fail cut at pt " << pt.x() << "," << pt.y() << "," << pt.z() << " radius " << radius << "\n";
#endif
      }
    Msg::Info("AngioTkCenterline: ->> cut branch %d/%d done in %d (planned %d) parts",i,(edges.size()-1), cptRealCut,cutDiskToPerform[i].size());

  }

#if 0
  //splitMesh
  for(unsigned int i = 0; i < edges.size(); i++){
    std::vector<MLine*> lines = edges[i].lines;
    double L = edges[i].length;
    double D = 2.*edges[i].minRad;  //(edges[i].minRad+edges[i].maxRad);
    double AR = L/D;
    // printf("*** AngioTkCenterline branch %d (AR=%.1f) \n", edges[i].tag, AR);

    int nbSplit = (int)ceil(AR/2 + 1.1); //AR/2 + 0.9
    //int nbSplit = (int)ceil(AR/4 + 1.1); //AR/2 + 0.9
    //nbSplit = std::min( nbSplit, 8 );
    if( nbSplit > 1 ){
      printf("->> cut branch in %d parts \n",  nbSplit);
      std::cout << "->> cut branch "<< i << "(L="<< L << " ,D="<< D << " ,AR=" << AR << ")"
		<< " in " << nbSplit << " parts\n";
      double li  = L/nbSplit;
      double lc = 0.0;
      for (unsigned int j= 0; j < lines.size(); j++){
      //for (unsigned int j= 2/*0*/; j < lines.size()-2; j++){
	lc += lines[j]->getLength();
	//if ( j < 3 || j >  lines.size()-3 ) continue;
	//if ( j != (int)ceil( lines.size()/2) ) continue;
	if (lc > li && nbSplit > 1) {
	  MVertex *v1 = lines[j]->getVertex(0);
	  MVertex *v2 = lines[j]->getVertex(1);
	  SVector3 pt(v1->x(), v1->y(), v1->z());
	  SVector3 dir(v2->x()-v1->x(),v2->y()-v1->y(),v2->z()-v1->z());
	  std::map<MLine*,double>::iterator itr = radiusl.find(lines[j]);
	  cutByDisk(pt, dir, itr->second);
	  nbSplit--;
	  lc = 0.0;
	}
      }
    }
#if 0
    if(edges[i].children.size() > 0.0 && AR > 1.0){
      MVertex *v1 = lines[lines.size()-1]->getVertex(1);//end vertex
      MVertex *v2;
      if(AR < 1.5) v2 = lines[0]->getVertex(0);
      else if (lines.size() > 4) v2 = lines[lines.size()-4]->getVertex(0);
      else v2 = lines[lines.size()-1]->getVertex(0);
      SVector3 pt(v1->x(), v1->y(), v1->z());
      SVector3 dir(v2->x()-v1->x(),v2->y()-v1->y(),v2->z()-v1->z());
      //printf("-->> cut branch at bifurcation \n");
      std::map<MLine*,double>::iterator itr = radiusl.find(lines[lines.size()-1]);
      //bool cutted =
      cutByDisk(pt, dir, itr->second);
      // if(!cutted){
      //   int l = lines.size()-1-lines.size()/(4*nbSplit); //chech this!
      //   v1 = lines[l]->getVertex(1);
      //   v2 = lines[l]->getVertex(0);
      //   pt = SVector3(v1->x(), v1->y(), v1->z());
      //   dir = SVector3(v2->x()-v1->x(),v2->y()-v1->y(),v2->z()-v1->z());
      //   printf("-->> cut bifurcation NEW \n");
      //   itr = radiusl.find(lines[l]);
      //   cutted = cutByDisk(pt, dir, itr->second);
      // }
    }
#endif
 }
#endif
  //create discreteFaces
  createFaces();
  current->createTopologyFromFaces(discFaces);
  current->exportDiscreteGEOInternals();

  //write
  std::string _remeshPartitionMeshFile = ( remeshPartitionMeshFile.empty() )? "myPARTS.msh" : remeshPartitionMeshFile;
  Msg::Info("AngioTkCenterline: writing splitted mesh '%s'",_remeshPartitionMeshFile.c_str());
  current->writeMSH(remeshPartitionMeshFile, 2.2, false, false);
  //exit(0);
#if 0
  //create compounds
  createSplitCompounds();
#endif
#if 0
  if ( !useGmshExecutable )
    current->writeSTL("myPARTS.stl",  false, false);
#endif
  Msg::Info("Done splitting mesh by centerlines");
}



void AngioTkCenterline::runSurfaceRemesh( std::string const& remeshPartitionMeshFile, bool forceRebuildPartition )
{
  if ( remeshPartitionMeshFile.empty() || StatFile(remeshPartitionMeshFile) || forceRebuildPartition )
    {
      this->cutMesh( remeshPartitionMeshFile );
    }
  else
    {
      GModel * gmodel = new GModel();
      // add new model as current
      GModel::current(GModel::list.size() - 1);

      current = GModel::current();
      current->load(remeshPartitionMeshFile);
      //current->removeDuplicateMeshVertices(1.e-8);

      // tolerance may be important in createTopologyFromMesh (the cutting in triangulation can generated point very close)
      // if we remove some points the mesh is not valid
      CTX::instance()->geom.tolerance=1e-8;
      current->createTopologyFromMesh();
    }

  // create compounds, apply surface parametrisation by each partition and get final mesh remeshed
  this->createSplitCompounds();
}

void AngioTkCenterline::saveCurrentGModelSTL(std::string const outputPath, bool binary )
{
  if (!current) { Msg::Error("not model available");return; }
  bool saveAll=false; double scalingFactor=1.0;
  current->writeSTL( outputPath,binary,saveAll,scalingFactor );
}

void AngioTkCenterline::saveSurfaceRemeshSTL(std::string const outputPath, bool binary )
{
  this->saveCurrentGModelSTL( outputPath, binary );
}

void AngioTkCenterline::saveClipMeshSTL(std::string const outputPath, bool binary )
{
  this->saveCurrentGModelSTL( outputPath, binary );
}

void AngioTkCenterline::saveTubularExtensionSTL(std::string const outputPath, bool binary )
{
  this->saveCurrentGModelSTL( outputPath, binary );
}



bool AngioTkCenterline::cutByDisk(SVector3 &PT, SVector3 &NORM, double &maxRad, int tag)
{
  double a = NORM.x();
  double b = NORM.y();
  double c = NORM.z();
  double d = -a * PT.x() - b * PT.y() - c * PT.z();

  int maxStep = 10;//20;
  //const double EPS = 0.007;
  // old value
  const double EPS = 1e-5;//1e-8;//1e-5;//0.007;
  //const double EPS = 1e-8;// new

  //std::set<MEdge,Less_Edge> allEdges;
  std::vector<MEdge> allEdges;
  for(unsigned int i = 0; i < triangles.size(); i++){
    for ( unsigned int j= 0; j <  3; j++){
      allEdges.push_back(triangles[i]->getEdge(j));
      //allEdges.insert(triangles[i]->getEdge(j));
    }
  }
  std::unique(allEdges.begin(), allEdges.end());

  bool closedCut = false;
  int step = 0;
  while (!closedCut && step < maxStep){
    double rad = 1.1*maxRad+0.05*step*maxRad;
    std::map<MEdge,MVertex*,Less_Edge> cutEdges;
    std::vector<MVertex*> cutVertices;
    std::vector<MTriangle*> newTris;
    std::set<MEdge,Less_Edge> newCut;
    cutEdges.clear();
    cutVertices.clear();
    newTris.clear();
    newCut.clear();
    // for (std::set<MEdge,Less_Edge>::iterator it = allEdges.begin();
    // 	 it != allEdges.end() ; ++it){
    // MEdge me = *it;
    for (unsigned int j = 0; j < allEdges.size(); j++){
      MEdge me = allEdges[j];
      SVector3 P1(me.getVertex(0)->x(),me.getVertex(0)->y(), me.getVertex(0)->z());
      SVector3 P2(me.getVertex(1)->x(),me.getVertex(1)->y(), me.getVertex(1)->z());
      double V1 = a * P1.x() + b * P1.y() + c * P1.z() + d;
      double V2 = a * P2.x() + b * P2.y() + c * P2.z() + d;
      //bool inters = (V1*V2<=0.0) ? true: false;
      bool inters = (V1*V2<(EPS*EPS)) ? true: false;
      //inters=false;
      //if ( (V1 < -EPS && V2 > -EPS) || (V1 > EPS && V2 < EPS) ) inters=true;

      bool inDisk = ((norm(P1-PT) < rad ) || (norm(P2-PT) < rad)) ? true : false;
      //bool inDisk2 = ((norm(P1-PT) < rad ) && (norm(P2-PT) < rad)) ? true : false;
      double rdist = -V1/(V2-V1);

      if ( std::abs(V1)<=EPS/*1e-4*/ && std::abs(V2)<=EPS/*1e-4*/ && (norm(P1-PT) < rad ) && (norm(P2-PT) < rad) ){
	inters = false;
	//Msg::Warning("les 2 pts sont dans le plan de coupe");
	cutVertices.push_back(me.getVertex(0));
	cutVertices.push_back(me.getVertex(1));
      }
      else if ( std::abs(V1)<=EPS && (norm(P1-PT) < rad ) )
	cutVertices.push_back(me.getVertex(0));
      else if ( std::abs(V2)<=EPS && (norm(P2-PT) < rad ) )
	cutVertices.push_back(me.getVertex(1));
      else if (inters && rdist > EPS && rdist < (1.-EPS) ){
      //else if (inters && rdist > 1e-5 && rdist < (1.-1e-5) ){
      //if (inters && rdist > EPS ){ // && rdist < (1.-EPS) ){
	SVector3 PZ = P1+rdist*(P2-P1);
	bool inDiskPZ = (norm(PZ-PT) < rad );
	if (inDiskPZ){
	  // cut edge me at point newv
          MVertex *newv = new MVertex (PZ.x(), PZ.y(), PZ.z());
          cutEdges.insert(std::make_pair(me,newv));
        }
      }
#if 0
      else if (inters && rdist <= EPS && inDisk && (norm(P1-PT) < rad ) )
	cutVertices.push_back(me.getVertex(0));
      else if (inters && rdist >= (1.-EPS) && inDisk && (norm(P2-PT) < rad) )
	cutVertices.push_back(me.getVertex(1));
#endif
    }
    //std::cout << "cutEdges.size() " << cutEdges.size() << "\n";
    //std::cout << "cutVertices.size() " << cutVertices.size() << "\n";
    for(unsigned int i = 0; i < triangles.size(); i++){
      cutTriangle(triangles[i], cutEdges,cutVertices, newTris, newCut);
    }
    if (!newCut.empty() && isClosed(newCut)) {
      triangles.clear();
      triangles = newTris;
      theCut.insert(newCut.begin(),newCut.end());
      if ( tag > 0 )
	theCutMap[tag].insert(newCut.begin(),newCut.end());
      break;
    }
    else {
      step++;
      //if (step == maxStep) {printf("no closed cut %d \n", (int)newCut.size()); };
      // // triangles = newTris;
      // // theCut.insert(newCut.begin(),newCut.end());
      // char name[256];
      // sprintf(name, "myCUT-%d.pos", step);
      // FILE * f2 = Fopen(name,"w");
      // fprintf(f2, "View \"\"{\n");
      // std::set<MEdge,Less_Edge>::iterator itp =  newCut.begin();
      // while (itp != newCut.end()){
      // 	MEdge l = *itp;
      // 	fprintf(f2, "SL(%g,%g,%g,%g,%g,%g){%g,%g};\n",
      // 		  l.getVertex(0)->x(), l.getVertex(0)->y(), l.getVertex(0)->z(),
      // 		  l.getVertex(1)->x(), l.getVertex(1)->y(), l.getVertex(1)->z(),
      // 		  1.0,1.0);
      // 	  itp++;
      // 	}
      // 	fprintf(f2,"};\n");
      // 	fclose(f2);
    }
  }


  if (step < maxStep){
    //printf("cutByDisk OK step =%d  \n", step);
    return true;
  }
  else {
    //printf("cutByDisk not succeeded \n");
    return false;
  }

}

double AngioTkCenterline::operator() (double x, double y, double z, GEntity *ge)
{

  if (update_needed){
     std::ifstream input;
     input.open(fileName.c_str());
     if(StatFile(fileName))
       Msg::Fatal("Centerline file '%s' does not exist", fileName.c_str());
     importFile(fileName);
     buildKdTree();
     update_needed = false;
   }

   double xyz[3] = {x,y,z};
   //take xyz = closest point on boundary in case we are on the planar in/out faces
   //or in the volume
   bool isCompound = false;
   if(ge){
     if (ge->dim() == 2 && ge->geomType() == GEntity::CompoundSurface) isCompound = true;
     std::list<GFace*> cFaces;
     if (isCompound) cFaces = ((GFaceCompound*)ge)->getCompounds();
     if ( ge->dim() == 3 || (ge->dim() == 2 && ge->geomType() == GEntity::Plane) ||
	  (isCompound && (*cFaces.begin())->geomType() == GEntity::Plane) ){
       const int num_neighbours = 1;
       ANNidx index[num_neighbours];
       ANNdist dist[num_neighbours];
       kdtreeR->annkSearch(xyz, num_neighbours, index, dist);
       ANNpointArray nodesR = kdtreeR->thePoints();
       xyz[0] = nodesR[index[0]][0];
       xyz[1] = nodesR[index[0]][1];
       xyz[2] = nodesR[index[0]][2];
     }
   }

   const int num_neighbours = 1;
   ANNidx index[num_neighbours];
   ANNdist dist[num_neighbours];
   kdtree->annkSearch(xyz, num_neighbours, index, dist);
   double rad = sqrt(dist[0]);

   //double cmax, cmin;
   //SVector3 dirMax,dirMin;
   //cmax = ge->curvatures(SPoint2(u, v),&dirMax, &dirMin, &cmax,&cmin);
   //cmax = ge->curvatureMax(SPoint2(u,v));
   //double radC = 1./cmax;

   double lc = 2*M_PI*rad/nbPoints;

   if(!ge) { return rad;}
   else  return lc;

}

void  AngioTkCenterline::operator() (double x, double y, double z, SMetric3 &metr, GEntity *ge)
{

   if (update_needed){
     std::ifstream input;
     input.open(fileName.c_str());
     if(StatFile(fileName))
       Msg::Fatal("Centerline file '%s' does not exist", fileName.c_str());
     importFile(fileName);
     buildKdTree();
     update_needed = false;
   }


   //take xyz = closest point on boundary in case we are on
   //the planar IN/OUT FACES or in VOLUME
   double xyz[3] = {x,y,z};
   bool onTubularSurface = true;
   double ds = 0.0;
   bool isCompound = (ge->dim() == 2 && ge->geomType() == GEntity::CompoundSurface) ?
     true : false;
   bool onInOutlets = (ge->geomType() == GEntity::Plane) ? true: false;
   std::list<GFace*> cFaces;
   if (isCompound) cFaces = ((GFaceCompound*)ge)->getCompounds();
   if ( ge->dim() == 3 || (ge->dim() == 2 && ge->geomType() == GEntity::Plane) ||
	(isCompound && (*cFaces.begin())->geomType() == GEntity::Plane) ){
     onTubularSurface = false;
   }

   ANNidx index[1];
   ANNdist dist[1];
   kdtreeR->annkSearch(xyz, 1, index, dist);
   if (! onTubularSurface){
     ANNpointArray nodesR = kdtreeR->thePoints();
     ds = sqrt(dist[0]);
     xyz[0] = nodesR[index[0]][0];
     xyz[1] = nodesR[index[0]][1];
     xyz[2] = nodesR[index[0]][2];
   }

   ANNidx index2[2];
   ANNdist dist2[2];
   kdtree->annkSearch(xyz, 2, index2, dist2);
   double radMax = sqrt(dist2[0]);
   ANNpointArray nodes = kdtree->thePoints();
   SVector3  p0(nodes[index2[0]][0], nodes[index2[0]][1], nodes[index2[0]][2]);
   SVector3  p1(nodes[index2[1]][0], nodes[index2[1]][1], nodes[index2[1]][2]);

   //dir_a = direction along the centerline
   //dir_n = normal direction of the disk
   //dir_t = tangential direction of the disk
   SVector3 dir_a = p1-p0; dir_a.normalize();
   SVector3 dir_n(xyz[0]-p0.x(), xyz[1]-p0.y(), xyz[2]-p0.z()); dir_n.normalize();
   SVector3 dir_cross  = crossprod(dir_a,dir_n); dir_cross.normalize();
   // if (ge->dim()==3){
   //   printf("coucou dim ==3 !!!!!!!!!!!!!!! \n");
   //   SVector3 d1,d2,d3;
   //   computeCrossField(x,y,z,d1,d2,d3);
   //   exit(1);
   // }

   //find discrete curvature at closest vertex
   Curvature& curvature = Curvature::getInstance();
   if( !Curvature::valueAlreadyComputed() ) {
     Msg::Info("Need to compute discrete curvature");
     Curvature::typeOfCurvature type = Curvature::RUSIN;
     curvature.computeCurvature(current, type);
   }
   double curv, cMin, cMax;
   SVector3 dMin, dMax;
   int isAbs = 1.0;
   curvature.vertexNodalValuesAndDirections(vertices[index[0]],&dMax, &dMin, &cMax, &cMin, isAbs);
   curvature.vertexNodalValues(vertices[index[0]], curv, 1);
   if (cMin == 0) cMin =1.e-12;
   if (cMax == 0) cMax =1.e-12;
   double rhoMin = 1./cMin;
   double rhoMax = 1./cMax;
   //double signMin = (rhoMin > 0.0) ? -1.0: 1.0;
   //double signMax = (rhoMax > 0.0) ? -1.0: 1.0;

   //user-defined parameters
   //define h_n, h_t1, and h_t2
   double thickness = radMax/3.;
   double h_far = radMax/5.;
   double beta = (ds <= thickness) ? 1.2 : 2.1; //CTX::instance()->mesh.smoothRatio;
   double ddist = (ds <= thickness) ? ds: thickness;

   double h_n_0 = thickness/20.;
   double h_n   = std::min( (h_n_0+ds*log(beta)), h_far);

   double betaMin = 10.0;
   double betaMax = 3.1;
   double oneOverD2_min = 1./(2.*rhoMin*rhoMin*(betaMin*betaMin-1)) *
     (sqrt(1+ (4.*rhoMin*rhoMin*(betaMin*betaMin-1))/(h_n*h_n))-1.);
   double oneOverD2_max = 1./(2.*rhoMax*rhoMax*(betaMax*betaMax-1)) *
    (sqrt(1+ (4.*rhoMax*rhoMax*(betaMax*betaMax-1))/(h_n*h_n))-1.);
   double h_t1_0 = sqrt(1./oneOverD2_min);
   double h_t2_0 = sqrt(1./oneOverD2_max);
   //double h_t1 =  h_t1_0*(rhoMin+signMin*ddist)/rhoMin ;
   //double h_t2 =  h_t2_0*(rhoMax+signMax*ddist)/rhoMax ;
   double h_t1  = std::min( (h_t1_0+(ddist*log(beta))), radMax);
   double h_t2  = std::min( (h_t2_0+(ddist*log(beta))), h_far);

   double dCenter = radMax-ds;
   double h_a_0 = 0.5*radMax;
   double h_a = h_a_0 - (h_a_0-h_t1_0)/(radMax)*dCenter;

   //length between min and max
   double lcMin = ((2 * M_PI *radMax) /( 50*nbPoints )); //CTX::instance()->mesh.lcMin;
   double lcMax =  lcMin*2000.; //CTX::instance()->mesh.lcMax;
   h_n = std::max(h_n, lcMin);    h_n = std::min(h_n, lcMax);
   h_t1 = std::max(h_t1, lcMin);  h_t1 = std::min(h_t1, lcMax);
   h_t2 = std::max(h_t2, lcMin);  h_t2 = std::min(h_t2, lcMax);

   //curvature metric
   SMetric3 curvMetric, curvMetric1, curvMetric2;
   SMetric3 centMetric1, centMetric2, centMetric;
   if (onInOutlets){
     metr = buildMetricTangentToCurve(dir_n,h_n,h_t2);
   }
   else {
     //on surface and in volume boundary layer
     if ( ds < thickness ){
       metr = metricBasedOnSurfaceCurvature(dMin, dMax, cMin, cMax, h_n, h_t1, h_t2);
     }
     //in volume
     else {
       //curvMetric = metricBasedOnSurfaceCurvature(dMin, dMax, cMin, cMax, h_n, h_t1, h_t2);
       metr = SMetric3( 1./(h_a*h_a), 1./(h_n*h_n), 1./(h_n*h_n), dir_a, dir_n, dir_cross);

       //metr = intersection_conserveM1_bis(metr, curvMetric);
       //metr = intersection_conserveM1(metr,curvMetric);
       //metr = intersection_conserve_mostaniso(metr, curvMetric);
       //metr = intersection(metr,curvMetric);
     }
   }

   return;

}

SMetric3 AngioTkCenterline::metricBasedOnSurfaceCurvature(SVector3 dirMin, SVector3 dirMax,
                                                   double cmin, double cmax,
						   double h_n, double h_t1, double h_t2)
{

  SVector3 dirNorm = crossprod(dirMax,dirMin);
  SMetric3 curvMetric (1./(h_t1*h_t1),1./(h_t2*h_t2),1./(h_n*h_n), dirMin, dirMax, dirNorm);

  return curvMetric;
}

void AngioTkCenterline::computeCrossField(double x,double y,double z,
				   SVector3 &d1, SVector3 &d2,  SVector3 &d3)
{

  //Le code suivant permet d'obtenir les trois vecteurs orthogonaux en un point.
  // int NumSmooth = 10;//CTX::instance()->mesh.smoothCrossField
  // Matrix m2;
  // if(NumSmooth)
  //   m2 = Frame_field::findCross(x,y,z);
  // else
  //   m2 = Frame_field::search(x,y,z);

}

//int readVTKPolyDataFields( const std::string &name, std::map<std::string,std::vector< std::vector<double> > > & fieldsPointData, bool bigEndian=false )
int readVTKPolyDataFields( const std::string &name, std::map<std::string,std::vector< std::vector<double> > > & fieldsPointData, bool bigEndian )
{
  FILE *fp = Fopen(name.c_str(), "rb");
  if(!fp){
    Msg::Error("Unable to open file '%s'", name.c_str());
    return 0;
  }

  char buffer[256], buffer2[256];
  std::map<int, std::map<int, std::string> > physicals[4];

  if(!fgets(buffer, sizeof(buffer), fp)){ fclose(fp); return 0; } // version line
  if(!fgets(buffer, sizeof(buffer), fp)){ fclose(fp); return 0; } // title

  if(fscanf(fp, "%s", buffer) != 1) // ASCII or BINARY
    Msg::Error("Failed reading buffer");
  bool binary = false;
  if(!strcmp(buffer, "BINARY")) binary = true;

  if(fscanf(fp, "%s %s", buffer, buffer2) != 2){ fclose(fp); return 0; }

  bool unstructured = false;
  if(!strcmp(buffer, "DATASET") && !strcmp(buffer2, "UNSTRUCTURED_GRID"))
    unstructured = true;

  if((strcmp(buffer, "DATASET") &&  strcmp(buffer2, "UNSTRUCTURED_GRID")) ||
     (strcmp(buffer, "DATASET") &&  strcmp(buffer2, "POLYDATA"))){
    Msg::Error("VTK reader can only read unstructured or polydata datasets");
    fclose(fp);
    return 0;
  }

  // read mesh vertices
  int numVertices;
  if(fscanf(fp, "%s %d %s\n", buffer, &numVertices, buffer2) != 3) return 0;
  if(strcmp(buffer, "POINTS") || !numVertices){
    Msg::Warning("No points in dataset");
    fclose(fp);
    return 0;
  }
  int datasize;
  if(!strcmp(buffer2, "double"))
    datasize = sizeof(double);
  else if(!strcmp(buffer2, "float"))
    datasize = sizeof(float);
  else{
    Msg::Warning("VTK reader only accepts float or double datasets");
    fclose(fp);
    return 0;
  }
  Msg::Info("Reading %d points", numVertices);
  std::vector<MVertex*> vertices(numVertices);
  for(int i = 0 ; i < numVertices; i++){
    double xyz[3];
    if(binary){
      if(datasize == sizeof(float)){
        float f[3];
        if(fread(f, sizeof(float), 3, fp) != 3){ fclose(fp); return 0; }
        if(!bigEndian) SwapBytes((char*)f, sizeof(float), 3);
        for(int j = 0; j < 3; j++) xyz[j] = f[j];
      }
      else{
        if(fread(xyz, sizeof(double), 3, fp) != 3){ fclose(fp); return 0; }
        if(!bigEndian) SwapBytes((char*)xyz, sizeof(double), 3);
      }
    }
    else{
      if(fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]) != 3){
        fclose(fp);
        return 0;
      }
    }
    vertices[i] = new MVertex(xyz[0], xyz[1], xyz[2]);
  }

  // read mesh elements
  int numElements, totalNumInt;
  if(fscanf(fp, "%s %d %d\n", buffer, &numElements, &totalNumInt) != 3){
    fclose(fp);
    return 0;
  }

  bool haveCells = true;
  bool haveLines = false;
  if( !strcmp(buffer, "CELLS") && numElements > 0)
    Msg::Info("Reading %d cells", numElements);
  else if (!strcmp(buffer, "POLYGONS") && numElements > 0)
    Msg::Info("Reading %d polygons", numElements);
  else if (!strcmp(buffer, "LINES") && numElements > 0){
    haveCells = false;
    haveLines = true;
    Msg::Info("Reading %d lines", numElements);
  }
  else{
    Msg::Warning("No cells or polygons in dataset");
    fclose(fp);
    return 0;
  }

  if (haveCells){
	Msg::Error("Cell not implement ");
    }
  else if ( haveLines )
    {
      if(!binary){
	int v0, v1;
	char line[100000], *p, *pEnd, *pEnd2;
	int iLine = 1;
	for (int k= 0; k < numElements; k++){
	  physicals[1][iLine][1] = "centerline";
	  if(!fgets(line, sizeof(line), fp)){ fclose(fp); return 0; }
	  v0 = (int)strtol(line, &pEnd, 10); //ignore first line
	  v0 = (int)strtol(pEnd, &pEnd2, 10);
	  p=pEnd2;
	  while(1){
	    v1 = strtol(p, &pEnd, 10);
	    if (p == pEnd )  break;
	    //elements[1][iLine].push_back(new MLine(vertices[v0],vertices[v1]));
	    p = pEnd;
	    v0 = v1;
	  }
	  iLine++;
	}
      }
      else{
	Msg::Error("TODO: implement reading lines for binary files \n");
      }
    }

  int numVertices2;
  if(fscanf(fp, "%s %d\n", buffer, &numVertices2) != 2){
    std::cout << "AIE\n";
    fclose(fp);
    return 0;
  }
  //else std::cout << "buffer "<< buffer <<"\n"
  int numFields;
  if(fscanf(fp, "%s %s %d\n", buffer, buffer2, &numFields) != 3){
    std::cout << "AIE2\n";
    fclose(fp);
    return 0;
  }

  for ( int ff =0;ff < numFields;++ff )
    {
  int nComp,numVertices3;
  if(fscanf(fp, "%s %d %d %s\n", buffer, &nComp,&numVertices3, buffer2) != 4){
    std::cout << "AIE3\n";
    fclose(fp);
    return 0;
  }
  std::string nameField( buffer );
  //std::cout <<  "nameField " << nameField << "\n";

  int datasizeField;
  if(!strcmp(buffer2, "double"))
    datasizeField = sizeof(double);
  else if(!strcmp(buffer2, "float"))
    datasizeField = sizeof(float);
  /*else{
    Msg::Warning("VTK reader only accepts float or double datasets");
    fclose(fp);
    return 0;
    }*/


  fieldsPointData[nameField].resize( numVertices3 );
  //std::vector<double> valuesInField( numVertices3 );
  for(int i = 0 ; i < numVertices3; i++)
    {
      fieldsPointData[nameField]/*valuesInFiel*/[i].resize(nComp);
      for(int comp = 0 ; comp < nComp; comp++)
      fscanf(fp, "%lf", &fieldsPointData[nameField]/*valuesInFiel*/[i][comp] );
      //fscanf(fp, "%lf", &valuesInField[i] );
    }

    }

  fclose(fp);


  return 1;



}


int writeVTKPolyData( std::shared_ptr<GModel>/*GModel **/ gmodel,std::vector<Branch> const& edges,
		      const std::string &name,
		      std::map<std::string,std::vector<std::vector<double> > > const& fieldsPointData,
		      bool binary=false,
		      bool saveAll=false, double scalingFactor=1.0,
		      bool bigEndian=false)
{
  unsigned int nBranch = edges.size();
  if ( nBranch == 0 )
    {
      Msg::Warning("centerline file is not exported because centerline is empty (path=%s)", name.c_str());
      return 1;
    }

  FILE *fp = Fopen(name.c_str(), binary ? "wb" : "w");
  if(!fp){
    Msg::Error("Unable to open file '%s'", name.c_str());
    return 0;
  }

  if(gmodel->noPhysicalGroups()) saveAll = true;

  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "%s, Created by Gmsh\n", gmodel->getName().c_str());
  if(binary)
    fprintf(fp, "BINARY\n");
  else
    fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET POLYDATA\n");


  // get the number of vertices
  int numVertices = numberOfVertices(edges);
  fprintf(fp, "POINTS %d double\n", numVertices);
  //std::cout << "WRITEVTK nBranch " << nBranch << "\n"; 
  //std::cout << "WRITEVTK numVertices " << numVertices << "\n"; 

  // write points
  int cptPtId = 0;
  std::map<int,int> mapInputIdToVtkId;
  std::map<int,int> mapVtkIdToInputId;
  for(unsigned int i = 0; i < nBranch; ++i)
  {
    std::vector<MLine*> mylines = edges[i].lines;
    for(unsigned int k = 0; k < mylines.size(); ++k)
    {
      MLine *l = mylines[k];
      MVertex *v0 = l->getVertex(0);
      MVertex *v1 = l->getVertex(1);
      int v0Id = v0->getIndex();
      int v1Id = v1->getIndex();
      auto itFind0 = mapInputIdToVtkId.find( v0Id );
      if ( itFind0 == mapInputIdToVtkId.end() )
      {
	v0->writeVTK(fp, binary, scalingFactor, bigEndian);
	mapInputIdToVtkId[v0Id] = cptPtId;
	mapVtkIdToInputId[cptPtId] = v0Id;
	++cptPtId;
      }
      auto itFind1 = mapInputIdToVtkId.find( v1Id );
      if ( itFind1 == mapInputIdToVtkId.end() )
      {
	v1->writeVTK(fp, binary, scalingFactor, bigEndian);
	mapInputIdToVtkId[v1Id] = cptPtId;
	mapVtkIdToInputId[cptPtId] = v1Id;
	++cptPtId;
      }
    }
  }
  fprintf(fp, "\n");

  // write lines
  unsigned int cptLines=0;
  for(unsigned int i = 0; i < nBranch; ++i){
    std::vector<MLine*> lines = edges[i].lines;
    cptLines +=lines.size()+2;
  }
  //std::cout << "WRITEVTK cptLines " << cptLines << "\n"; 

  fprintf(fp, "LINES %d %d\n", nBranch ,cptLines);
  for(unsigned int i = 0; i < nBranch; ++i){
    std::vector<MLine*> lines = edges[i].lines;
    bool firstPtDone = false;

    int idPtInListStart = -1, idPtInListFinish=-1;
    if ( edges[i].vB->getIndex() == lines.front()->getVertex(0)->getIndex() )
      {
       idPtInListStart=0;
       idPtInListFinish=1;
      }
    else if ( edges[i].vB->getIndex() == lines.front()->getVertex(1)->getIndex() )
      {
       idPtInListStart=1;
       idPtInListFinish=0;
      }
    if ( idPtInListStart==-1 ) Msg::Error("not found a good point to start");

    for(unsigned int k = 0; k < lines.size(); ++k){
      MLine *l = lines[k];
      if ( !firstPtDone )
	{
	  int nPointInBranch = lines.size() +1;
	  fprintf(fp, "%d", nPointInBranch );
	  if ( mapInputIdToVtkId.find(l->getVertex(idPtInListStart)->getIndex()) == mapInputIdToVtkId.end() ) Msg::Error("error mapInputIdToVtkId did not contains this entry");
	  fprintf(fp, " %d", mapInputIdToVtkId[l->getVertex(idPtInListStart)->getIndex()] );
	  firstPtDone = true;
	}

      if ( mapInputIdToVtkId.find(l->getVertex(idPtInListFinish)->getIndex()) == mapInputIdToVtkId.end() ) Msg::Error("error mapInputIdToVtkId did not contains this entry");
      fprintf(fp, " %d", mapInputIdToVtkId[l->getVertex(idPtInListFinish)->getIndex()] );
      //std::cout << "l->getVertex(1)->getIndex() " << l->getVertex(1)->getIndex() << "\n";
    }
    //std::cout << "edges[i].vB->getIndex() " << edges[i].vB->getIndex() << "\n";
    //std::cout << "edges[i].vE->getIndex() " << edges[i].vE->getIndex() << "\n";

    fprintf(fp, "\n");
  }

  // write fields
  if ( fieldsPointData.size() > 0 )
    {
    fprintf(fp, "\nPOINT_DATA %d\n", numVertices);
    fprintf(fp, "FIELD FieldData %d\n", (int)fieldsPointData.size());
    Msg::Info("AngioTkCenterline: writeVTKPolyData : fieldsPointData.size() : %d ", (int)fieldsPointData.size());
    Msg::Info("AngioTkCenterline: writeVTKPolyData : mapVtkIdToInputId.size() : %d ", (int)mapVtkIdToInputId.size());
    }
  Msg::Info("AngioTkCenterline: writeVTKPolyData : fieldsPointData start ");

  for ( auto const& thefield : fieldsPointData )
    {
      //std::cout << "thefield.second.size() " << thefield.second.size()  << "\n";
      //int nComp = thefield.second[mapVtkIdToInputId[0]-1].size();
      if ( thefield.second.empty() ) continue;
      int nComp = thefield.second.begin()->size();
      /*std::cout << "thefield.first" << thefield.first << " nComp=" << nComp
	<< " and mapVtkIdToInputId.size() " << mapVtkIdToInputId.size() << "\n";*/
      fprintf(fp, "%s %d %d double\n", thefield.first.c_str(), nComp, numVertices);

      int cptVal=0;
      for ( int k=0 ; k< mapVtkIdToInputId.size() ; ++k )
	{
	  //auto const& valAllComp = thefield.second[mapVtkIdToInputId[k]-1];
	  auto const& valAllComp = thefield.second[k];// NEW use vtk numeroation
	  //std::vector<double> valAllComp(nComp,0.);
	  if ( valAllComp.empty() ) Msg::Error("valAllComp is empty");
	  for ( double val : valAllComp )
	    {
	      fprintf(fp, "%.16g ", val );
	      if ( cptVal%6==5 ) 
		fprintf(fp, "\n" );
	      ++cptVal;
	    }	  
	}
      fprintf(fp, "\n" );
    }

  fclose(fp);
  return 1;

}

void AngioTkCenterline::writeCenterlinesVTK(std::string fileName)
{
  writeVTKPolyData(mod,edges,fileName,M_centerlinesFieldsPointData);
}

void AngioTkCenterline::updateCenterlinesFieldsFromFile(std::string fileName)
{

  std::vector<MVertex*> verticesPosition;
  int numVertices = numberOfVertices(edges);
  //std::cout << " numVertices " << numVertices << " and colorp " << colorp.size() << "\n";
  for ( auto& ptPair : colorp )
    {
      MVertex* myvertex = ptPair.first;
      //std::cout << myvertex->x() << " " << myvertex->y() << " "<< myvertex->z() << "\n";
      verticesPosition.push_back(new MVertex(myvertex->x(), myvertex->y(), myvertex->z(), NULL, myvertex->getIndex() ));
    }
  // init localisation tool
  MVertexPositionSet pos(verticesPosition);
  double tolerance = 1.0e-8;
  SBoundingBox3d bbox = mod->bounds();
  double lc = bbox.empty() ? 1. : norm(SVector3(bbox.max(), bbox.min()));
  double eps = lc * tolerance;

  // load initial centerlines with duplicate vertices/lines
  GModel * modInitial = new GModel();
  modInitial->load(fileName);
  bool saveAll=false;
  int numVerticesInitial = modInitial->indexMeshVertices(saveAll);
  //int numVerticesInitial = modInitial->getNumMeshVertices();
  //std::cout << "numVerticesInitial " << numVerticesInitial << "\n";
  std::vector<GEntity*> entitiesInitial;
  modInitial->getEntities(entitiesInitial);

  // compute vertices relation between initial mesh and clean mesh
  std::map<int,int> initialVertexIdToCleanVertexId;
  for(unsigned int i = 0; i < entitiesInitial.size(); i++)
    for(unsigned int j = 0; j < entitiesInitial[i]->mesh_vertices.size(); j++)
    {
      MVertex * myvertex = entitiesInitial[i]->mesh_vertices[j];

      MVertex *vertexFind = pos.find(myvertex->x(), myvertex->y(), myvertex->z(), eps);
      if ( !vertexFind )
	continue;
      if ( !vertexFind )
	Msg::Error("vertexFind not find ");

      initialVertexIdToCleanVertexId[myvertex->getIndex()] = vertexFind->getNum();
      //initialVertexIdToCleanVertexId[myvertex->getIndex()] = vertexFind->getIndex();
      //std::cout << "get correspondance "<< myvertex->getIndex() << " vs " << vertexFind->getNum() << " in " << i << "\n";
    }

#if 0
  for(unsigned int i = 0; i < entitiesInitial.size(); i++)
    for(unsigned int j = 0; j < entitiesInitial[i]->mesh_vertices.size(); j++)
    {
      MVertex * myvertexInitial = entitiesInitial[i]->mesh_vertices[j];
      int vertexIdClean = initialVertexIdToCleanVertexId[myvertexInitial->getIndex()];
      int entityIdClean = vertexIdToEntityIdAndLocalId[vertexIdClean].first;
      int localIdClean = vertexIdToEntityIdAndLocalId[vertexIdClean].second;
      MVertex * myvertexClean = entities[entityIdClean]->mesh_vertices[localIdClean];
      std::cout << myvertexInitial->x() << " " << myvertexInitial->y() << " "<< myvertexInitial->z()
		<< " VS "
		<< myvertexClean->x() << " " << myvertexClean->y() << " "<< myvertexClean->z()
		<< "\n";
    }
#endif

  // load vtk Point Data and store
  std::map<std::string,std::vector<std::vector<double> > >  fieldsPointDataInput;
  readVTKPolyDataFields( fileName, fieldsPointDataInput );

  // tranfert into clean centerlines the Point Data Value
  M_centerlinesFieldsPointData.clear();
  for ( auto const& thefield : fieldsPointDataInput )
    {
      M_centerlinesFieldsPointData[thefield.first].resize(numVertices);
      for ( auto & thefieldClean : M_centerlinesFieldsPointData[thefield.first] )
	thefieldClean.resize(thefield.second[0].size());
    }

  for(unsigned int i = 0; i < entitiesInitial.size(); i++)
    for(unsigned int j = 0; j < entitiesInitial[i]->mesh_vertices.size(); j++)
    {
      MVertex * myvertexInitial = entitiesInitial[i]->mesh_vertices[j];
      if ( initialVertexIdToCleanVertexId.find(myvertexInitial->getIndex()) == initialVertexIdToCleanVertexId.end() ) continue;
      int vertexIdClean = initialVertexIdToCleanVertexId[myvertexInitial->getIndex()];

      for ( auto const& thefield : fieldsPointDataInput )
	for ( int comp=0;comp< thefield.second[0].size();++comp )
	{
	  //std::cout << "vertexIdClean " << vertexIdClean << " this->M_mapVertexGmshIdToVtkId[vertexIdClean] " << this->M_mapVertexGmshIdToVtkId[vertexIdClean] << "\n";
	  //std::cout << "add value "<< thefield.first << " : "<< thefield.second[myvertexInitial->getIndex()-1][comp] << "\n";
#if 0
	  // apply a threshold
	  double valTolMin = 1.;
	  if ( thefield.first != "MaximumInscribedSphereRadius" ||thefield.second[myvertexInitial->getIndex()-1][comp] > valTolMin )
	    M_centerlinesFieldsPointData[thefield.first][this->M_mapVertexGmshIdToVtkId[vertexIdClean]][comp] = thefield.second[myvertexInitial->getIndex()-1][comp];
	  else
	    M_centerlinesFieldsPointData[thefield.first][this->M_mapVertexGmshIdToVtkId[vertexIdClean]][comp] = valTolMin;
#else
	  M_centerlinesFieldsPointData[thefield.first][this->M_mapVertexGmshIdToVtkId[vertexIdClean]][comp] = thefield.second[myvertexInitial->getIndex()-1][comp];
#endif

	}
    }
}

void AngioTkCenterline::addFieldBranchIds( std::string const& fieldName )
{
  //if ( M_centerlinesFieldsPointData.find("BranchIds") != M_centerlinesFieldsPointData.end() )
  //  return;
  Msg::Info("AngioTkCenterline: addFieldBranchIds");
  int numVertices = numberOfVertices(edges);
  M_centerlinesFieldsPointData[fieldName].resize( numVertices );

  std::vector<bool> ptDone(numVertices,false);
  int nBranch = edges.size();
  for(int i = 0; i < nBranch; ++i)
  {
    std::vector<MLine*> mylines = edges[i].lines;
    bool firstPtDone = false;
    for(unsigned int k = 0; k < mylines.size(); ++k)
    {
      MLine *l = mylines[k];
      MVertex * myvertex0 = l->getVertex(0);
      MVertex * myvertex1 = l->getVertex(1);
      int v0Id = myvertex0->getIndex(), v1Id = myvertex1->getIndex();
      int v0VtkId = this->M_mapVertexGmshIdToVtkId[v0Id], v1VtkId = this->M_mapVertexGmshIdToVtkId[v1Id];
      if ( !ptDone[v0VtkId] )
	{
	M_centerlinesFieldsPointData[fieldName][v0VtkId] = { (double)i };
	ptDone[v0VtkId] = true;
      }
      if ( !ptDone[v1VtkId] )
      {
	M_centerlinesFieldsPointData[fieldName][v1VtkId] = { (double)i };
	ptDone[v1VtkId] = true;
      }
    }
  }
}

void
AngioTkCenterline::addFieldRadiusMin( std::string const& fieldName )
{
  Msg::Info("AngioTkCenterline: addFieldRadiusMin");

  int numVertices = numberOfVertices(edges);
  if ( colorp.size() != numVertices ) 
    Msg::Error("invalid container size colorp.size()=%d vs numVertices=%d", colorp.size(), numVertices );

  M_centerlinesFieldsPointData[fieldName].clear();
  M_centerlinesFieldsPointData[fieldName].resize( numVertices );

  if( !triangles.empty())
    {
      for ( auto& ptPair : colorp )
	{
	  MVertex* myvertex = ptPair.first;
	  double pt[3] = { myvertex->x(), myvertex->y(), myvertex->z() };
	  ANNidx index[1];
	  ANNdist dist[1];
	  kdtreeR->annkSearch(pt, 1, index, dist);
	  double minRad = sqrt(dist[0]);

	  int vId = myvertex->getIndex();
	  int vVtkId = this->M_mapVertexGmshIdToVtkId[vId];
	  M_centerlinesFieldsPointData[fieldName][vVtkId] = { (double)minRad };
	}
    }
  else
    {
      for ( auto& ptPair : colorp )
	{
	  MVertex* myvertex = ptPair.first;
	  int vId = myvertex->getIndex();
	  int vVtkId = this->M_mapVertexGmshIdToVtkId[vId];
	  M_centerlinesFieldsPointData[fieldName][vVtkId] = { 0. };
	}

    }
}

double
AngioTkCenterline::minRadiusAtVertex( MVertex* myvertex ) const
{
  if ( !kdtreeR )
    Msg::Error("kdtreeR not built");

  double pt[3] = { myvertex->x(), myvertex->y(), myvertex->z() };
  ANNidx index[1];
  ANNdist dist[1];
  kdtreeR->annkSearch(pt, 1, index, dist);
  double minRad = sqrt(dist[0]);
  return minRad;
}

void
AngioTkCenterline::applyFieldThresholdMin(std::string const& fieldName, double value)
{
  this->applyFieldThresholdMin( std::vector<std::string>(1,fieldName),value );
}
void
AngioTkCenterline::applyFieldThresholdMax(std::string const& fieldName, double value)
{
}
void
AngioTkCenterline::applyFieldThresholdMin(std::vector<std::string> const& fieldNames, double value)
{
  int numVertices = numberOfVertices(edges);
  for ( std::string const& fieldName : fieldNames )
    {
      if ( M_centerlinesFieldsPointData.find(fieldName) != M_centerlinesFieldsPointData.end() &&
	   M_centerlinesFieldsPointData.find(fieldName)->second.size() == numVertices )
	{
	  for (int k=0;k<numVertices;++k)
	    {
	      if ( M_centerlinesFieldsPointData[fieldName][k][0] < value )
		M_centerlinesFieldsPointData[fieldName][k][0] = value;
	    }
	}
    }
}
void
AngioTkCenterline::applyFieldThresholdMax(std::vector<std::string> const& fieldNames, double value)
{
  int numVertices = numberOfVertices(edges);
  for ( std::string const& fieldName : fieldNames )
    {
      if ( M_centerlinesFieldsPointData.find(fieldName) != M_centerlinesFieldsPointData.end() &&
	   M_centerlinesFieldsPointData.find(fieldName)->second.size() == numVertices )
	{
	  for (int k=0;k<numVertices;++k)
	    {
	      if ( M_centerlinesFieldsPointData[fieldName][k][0] > value )
		M_centerlinesFieldsPointData[fieldName][k][0] = value;
	    }
	}
    }
}

void
AngioTkCenterline::applyFieldThresholdZoneMin(std::string const& fieldName, double value, std::map<int,std::vector<std::tuple<double,double,double> > > const& mapPointPair )
{
  this->applyFieldThresholdZoneImpl( std::vector<std::string>(1,fieldName),value,mapPointPair,0);
}
void
AngioTkCenterline::applyFieldThresholdZoneMax(std::string const& fieldName, double value, std::map<int,std::vector<std::tuple<double,double,double> > > const& mapPointPair )
{
  this->applyFieldThresholdZoneImpl(std::vector<std::string>(1,fieldName),value,mapPointPair,1);
}
void
AngioTkCenterline::applyFieldThresholdZoneMin(std::vector<std::string> const& fieldNames, double value, std::map<int,std::vector<std::tuple<double,double,double> > > const& mapPointPair )
{
  this->applyFieldThresholdZoneImpl(fieldNames,value,mapPointPair,0);
}
void
AngioTkCenterline::applyFieldThresholdZoneMax(std::vector<std::string> const& fieldNames, double value, std::map<int,std::vector<std::tuple<double,double,double> > > const& mapPointPair )
{
  this->applyFieldThresholdZoneImpl(fieldNames,value,mapPointPair,1);
}
void
AngioTkCenterline::applyFieldThresholdZoneImpl(std::vector<std::string> const& fieldNames, double value, std::map<int,std::vector<std::tuple<double,double,double> > > const& mapPointPair, int type )
{
  std::vector<MLine*> listOfLines;
  for ( auto pointPairBase : mapPointPair )
    {
      if ( pointPairBase.second.size() < 2 )
	Msg::Warning("ignore point pair");
      auto const& pt1 = pointPairBase.second[0];
      auto const& pt2 = pointPairBase.second[1];
      double ptToLocalize1[3] = { std::get<0>( pt1 ),std::get<1>( pt1 ),std::get<2>( pt1 ) };
      double ptToLocalize2[3] = { std::get<0>( pt2 ),std::get<1>( pt2 ),std::get<2>( pt2 ) };
      std::cout << "id " << pointPairBase.first << " ptToLocalize1 : " << ptToLocalize1[0] << " , " << ptToLocalize1[1]  << " , " << ptToLocalize1[2] << "\n";
      std::cout << "ptToLocalize2 : " << ptToLocalize2[0] << " , " << ptToLocalize2[1]  << " , " << ptToLocalize2[2] << "\n";
      MVertex* vertexA = std::get<0>( this->foundClosestPointInCenterlines(ptToLocalize1) );
      MVertex* vertexB = std::get<0>( this->foundClosestPointInCenterlines(ptToLocalize2) );
      //std::cout << "vertexA : " << vertexA->x() << " , " << vertexA->y()  << " , " << vertexA->z() << "\n";
      //std::cout << "vertexB : " << vertexB->x() << " , " << vertexB->y()  << " , " << vertexB->z() << "\n";

      std::vector<MLine*> listOfLinesTemp = std::get<0>(this->pathBetweenVertex( vertexA,vertexB ) );
      listOfLines.insert( listOfLines.end(),listOfLinesTemp.begin(),listOfLinesTemp.end() );
    }

  bool applyMin = (type == 0);
  bool applyMax = (type == 1);
  //bool applyMinMax = (type == 2);

  //if ( type == 0 ) // min
    {
      for ( auto myline : listOfLines )
	{
	  MVertex* v0 = myline->getVertex(0);
	  MVertex* v1 = myline->getVertex(1);
	  int vtkId = M_mapVertexGmshIdToVtkId[v0->getIndex()];
	  for ( std::string const& fieldName : fieldNames )
	    {
	      int nComp = M_centerlinesFieldsPointData[fieldName][vtkId].size();
	      for (int c=0;c<nComp;++c)
		{
		  if ( applyMin && M_centerlinesFieldsPointData[fieldName][vtkId][c] < value ) // min
		    M_centerlinesFieldsPointData[fieldName][vtkId][c]=value;//2.;//100;
		  if ( applyMax && M_centerlinesFieldsPointData[fieldName][vtkId][c] > value ) // max
		    M_centerlinesFieldsPointData[fieldName][vtkId][c]=value;//2.;//100;
		}
	    }
	}
    }
#if 0
  else if ( type == 1 ) // max
    {
      for ( auto myline : listOfLines )
	{
	  MVertex* v0 = myline->getVertex(0);
	  MVertex* v1 = myline->getVertex(1);
	  int vtkId = M_mapVertexGmshIdToVtkId[v0->getIndex()];
	  for ( std::string const& fieldName : fieldNames )
	    {
	      int nComp = M_centerlinesFieldsPointData[fieldName][vtkId].size();
	      for (int c=0;c<nComp;++c)
		{
		  if ( M_centerlinesFieldsPointData[fieldName][vtkId][c] > value )
		    M_centerlinesFieldsPointData[fieldName][vtkId][c]=value;
		}
	    }
	}
    }
#endif

}



void AngioTkCenterline::removeBranchIds( std::set<int> const& _removeBranchIds )
{
  if ( _removeBranchIds.empty() )
    return;

  Msg::Info("AngioTkCenterline: removeBranchIds");

  std::map<int,std::vector<MLine*> > _modEdgesConvert;
  int nBranch = edges.size();
  for(int i = 0; i < nBranch; ++i)
  {
    if ( _removeBranchIds.find( i+1 ) != _removeBranchIds.end() )
      continue;
    std::vector<MLine*> mylines = edges[i].lines;
    for(unsigned int k = 0; k < mylines.size(); ++k)
    {
      MLine *myline = mylines[k];
      _modEdgesConvert[i+1].push_back( myline );
    }
  }

  // copy previous mapping
  std::map<int,int> _previousMapVertexGmshIdToVtkId;
  _previousMapVertexGmshIdToVtkId.insert( this->M_mapVertexGmshIdToVtkId.begin(),this->M_mapVertexGmshIdToVtkId.end() );
  // update for use
  this->updateCenterlinesForUse(_modEdgesConvert);
  // update fields data
  this->updateFieldsDataAfterReduction(_previousMapVertexGmshIdToVtkId);

}


void AngioTkCenterline::removeDuplicateBranch()
{
  Msg::Info("AngioTkCenterline: removeDuplicateBranch");
  Msg::Error("TODO : implemented really this function"); 

  std::map<int,std::vector<MLine*> > _modEdgesConvert;
  int nBranch = edges.size();
  for(int i = 0; i < nBranch; ++i)
  {
    if ( edges[i].length < 1.5*edges[i].maxRad &&
	 ( this->centerlinesExtremities().find( edges[i].vB ) != this->centerlinesExtremities().end() ||
	   this->centerlinesExtremities().find( edges[i].vE ) != this->centerlinesExtremities().end() ) ) { std::cout<<"delete branch id "<<i<<" \n";continue; }
    //if ( edges[i].length < 0.1*edges[i].maxRad ) { std::cout<<"delete branch id "<<i<<" \n";continue; }
    std::vector<MLine*> mylines = edges[i].lines;
    int idNewBranch = _modEdgesConvert.size()+1;
    for(unsigned int k = 0; k < mylines.size(); ++k)
      {
	MLine *myline = mylines[k];
	MVertex *v0 = myline->getVertex(0);
	MVertex *v1 = myline->getVertex(1);
	_modEdgesConvert[idNewBranch].push_back( myline );
	//_modEdgesConvert[i+1].push_back( myline );
      }
  }

  // copy previous mapping
  std::map<int,int> _previousMapVertexGmshIdToVtkId;
  _previousMapVertexGmshIdToVtkId.insert( this->M_mapVertexGmshIdToVtkId.begin(),this->M_mapVertexGmshIdToVtkId.end() );
  // update for use
  this->updateCenterlinesForUse(_modEdgesConvert);
  // update fields data
  this->updateFieldsDataAfterReduction(_previousMapVertexGmshIdToVtkId);

}

void AngioTkCenterline::checkCenterlinesConnectivity()
{
  int nBranch = edges.size();
  for(int i = 0; i < nBranch; ++i)
  {
    auto const& mybranch = edges[i];
    std::vector<MLine*> mylines = mybranch.lines;
    if ( mylines.empty() ) { Msg::Error("branch %d is empty (TODO : clean)",i);continue; }
    MVertex *vLast;
    bool startBy0 = edges[i].vB == mylines.front()->getVertex(0);
    int idLoc0 = 0;// ( startBy0 )? 0 : 1;
    int idLoc1 = 1;//( startBy0 )? 1 : 0;
    if ( edges[i].vB != mylines.front()->getVertex(0) ) Msg::Error("invalid start in lines");
    if ( edges[i].vE != mylines.back()->getVertex(1) ) Msg::Error("invalid end in lines");

    for(unsigned int k = 0; k < mylines.size(); ++k)
      {
	MLine *myline = mylines[k];
	MVertex *v0 = myline->getVertex(idLoc0);
	MVertex *v1 = myline->getVertex(idLoc1);
	if ( k == 0 && edges[i].vB != v0 )
	  Msg::Error("edges[%d].vB invalid %d vs %d",i,mybranch.vB->getIndex(),v0->getIndex());
	if ( k == (mylines.size()-1) && edges[i].vE != v1 )
	  Msg::Error("edges[%d].vE invalid %d vs %d",i,mybranch.vE->getIndex(),v1->getIndex());
	if ( k > 0 )
	  if ( v0 != vLast /*&& v1 != vLast*/ )
	    Msg::Error( "edges[%d] invalid connectivity %d vs %d",i,v0->getIndex(),vLast->getIndex());
	vLast = v1;
      }
  }
}

void AngioTkCenterline::updateFieldsDataAfterReduction( std::map<int,int> const& _previousMapVertexGmshIdToVtkId )
{
  std::map<std::string,std::vector<std::vector<double> > >  newCenterlinesFieldsPointData;

  for (auto const& fieldsPointDataPair : M_centerlinesFieldsPointData )
    {
      newCenterlinesFieldsPointData[fieldsPointDataPair.first].clear();
      newCenterlinesFieldsPointData[fieldsPointDataPair.first].resize(this->M_mapVertexVtkIdToGmshId.size());
	  //std::cout << "fieldsPointDataPair.first " << fieldsPointDataPair.first << " : " << newCenterlinesFieldsPointData[fieldsPointDataPair.first].size() << "\n";
      for (int k=0;k<this->M_mapVertexVtkIdToGmshId.size();++k)
	{
	  int gmshId = this->M_mapVertexVtkIdToGmshId[k];
	  auto itFindGmshIdInPrevious = _previousMapVertexGmshIdToVtkId.find( gmshId );
	  if ( itFindGmshIdInPrevious == _previousMapVertexGmshIdToVtkId.end() )
	    Msg::Error("Error in mapping 2\n");
	  int vtkId = itFindGmshIdInPrevious->second;
	  newCenterlinesFieldsPointData[fieldsPointDataPair.first][k] = this->M_centerlinesFieldsPointData[fieldsPointDataPair.first][vtkId];
	}
    }
  M_centerlinesFieldsPointData.clear();
  M_centerlinesFieldsPointData = newCenterlinesFieldsPointData;
}


void AngioTkCenterline::printSplit() const
{
  FILE * f = Fopen("mySPLIT.pos","w");
  fprintf(f, "View \"\"{\n");
  for(unsigned int i = 0; i < edges.size(); ++i){
    std::vector<MLine*> lines = edges[i].lines;
    for(unsigned int k = 0; k < lines.size(); ++k){
      MLine *l = lines[k];
      fprintf(f, "SL(%g,%g,%g,%g,%g,%g){%g,%g};\n",
	      l->getVertex(0)->x(), l->getVertex(0)->y(), l->getVertex(0)->z(),
	      l->getVertex(1)->x(), l->getVertex(1)->y(), l->getVertex(1)->z(),
	      (double)edges[i].tag, (double)edges[i].tag);
    }
  }
  fprintf(f,"};\n");
  fclose(f);

  // FILE * f3 = Fopen("myJUNCTIONS.pos","w");
  // fprintf(f3, "View \"\"{\n");
  //  std::set<MVertex*>::const_iterator itj = junctions.begin();
  //  while (itj != junctions.end()){
  //    MVertex *v =  *itj;
  //    fprintf(f3, "SP(%g,%g,%g){%g};\n",
  // 	     v->x(),  v->y(), v->z(),
  // 	     (double)v->getNum());
  //    itj++;
  // }
  // fprintf(f3,"};\n");
  // fclose(f3);

  FILE * f4 = Fopen("myRADII.pos","w");
  fprintf(f4, "View \"\"{\n");
  for(unsigned int i = 0; i < lines.size(); ++i){
    MLine *l = lines[i];
    std::map<MLine*,double>::const_iterator itc =  radiusl.find(l);
    fprintf(f4, "SL(%g,%g,%g,%g,%g,%g){%g,%g};\n",
 	    l->getVertex(0)->x(), l->getVertex(0)->y(), l->getVertex(0)->z(),
 	    l->getVertex(1)->x(), l->getVertex(1)->y(), l->getVertex(1)->z(),
 	    itc->second,itc->second);
  }
  fprintf(f4,"};\n");
  fclose(f4);

}

void AngioTkCenterline::updateRelationMapVertex()
{
  M_mapVertexGmshIdToVtkId.clear(); M_mapVertexVtkIdToGmshId.clear();
  this->updateRelationMapVertex(M_mapVertexGmshIdToVtkId,M_mapVertexVtkIdToGmshId);
}
void AngioTkCenterline::updateRelationMapVertex(std::map<int,int> & _mapVertexGmshIdToVtkId,
						std::map<int,int> & _mapVertexVtkIdToGmshId )
{
  int nBranch = edges.size();
  int cptPtId = 0;
  //std::map<int,int> mapInputIdToVtkId;
  //std::map<int,int> mapVtkIdToInputId;
  for(unsigned int i = 0; i < nBranch; ++i)
  {
    std::vector<MLine*> mylines = edges[i].lines;
    for(unsigned int k = 0; k < mylines.size(); ++k)
    {
      MLine *l = mylines[k];
      MVertex *v0 = l->getVertex(0);
      MVertex *v1 = l->getVertex(1);
      int v0Id = v0->getIndex();
      int v1Id = v1->getIndex();
      auto itFind0 = _mapVertexGmshIdToVtkId.find( v0Id );
      if ( itFind0 == _mapVertexGmshIdToVtkId.end() )
      {
	_mapVertexGmshIdToVtkId[v0Id] = cptPtId;
	_mapVertexVtkIdToGmshId[cptPtId] = v0Id;
	++cptPtId;
      }
      auto itFind1 = _mapVertexGmshIdToVtkId.find( v1Id );
      if ( itFind1 == _mapVertexGmshIdToVtkId.end() )
      {
	_mapVertexGmshIdToVtkId[v1Id] = cptPtId;
	_mapVertexVtkIdToGmshId[cptPtId] = v1Id;
	++cptPtId;
      }
    }
  }
}

std::tuple<MVertex*,double>
AngioTkCenterline::foundClosestPointInCenterlines( double ptToLocalize[3] )
{
  if ( !kdtree ) this->buildKdTree();//Msg::Error("kdtree not build");
  ANNidx index[1];

  ANNdist distVec[1];
  kdtree->annkSearch(ptToLocalize, 1, index, distVec);
  ANNpointArray nodes = kdtree->thePoints();
  double dist = distVec[0];
  MVertex* pointFound = NULL;
  if ( index[0] < colorp.size() )
    {
      std::map<MVertex*, int>::iterator itp = colorp.begin();
      for (int k=0;k<index[0];++k)
	++itp;
      pointFound = itp->first;
    }
  else
    {
      int nbPL = 3;  //10 points per line
      int ind = colorp.size();
      for(unsigned int k = 0; k < lines.size(); ++k)
	{
	  MVertex *v0 = lines[k]->getVertex(0);
	  MVertex *v1 = lines[k]->getVertex(1);
	  for (int j = 1; j < nbPL+1; j++)
	    {
	      if ( ind==index[0] )
		{
		  SPoint3 pt(ptToLocalize[0],ptToLocalize[1],ptToLocalize[3]);
		  double dist0 = v0->point().distance( pt );
		  double dist1 = v1->point().distance( pt );
		  if ( dist0 < dist1 )
		    {
		      pointFound = v0;
		      //dist = dist0;
		    }
		  else
		    {
		      pointFound = v1;
		      //dist = dist1;
		    }
		  break;
		}
	      ind++;
	    }
	  if ( pointFound != NULL ) break;
	}
    }
  return std::make_tuple(pointFound,dist);
}

int
AngioTkCenterline::vertexOnSameBranch( MVertex* vertexA, MVertex* vertexB )
{
  int nullBranchId = -1;
  auto itFindVertexA = M_vertexToLinesId.find( vertexA );
  auto itFindVertexB = M_vertexToLinesId.find( vertexB );
  if ( itFindVertexA == M_vertexToLinesId.end() || itFindVertexB == M_vertexToLinesId.end() )
    return nullBranchId;
  for ( auto const& vertexInfoA : itFindVertexA->second )
    {
      int branchIdA =  vertexInfoA.first;
      for ( auto const& vertexInfoB : itFindVertexB->second )
	{
	  int branchIdB = vertexInfoB.first;
	  if ( branchIdA == branchIdB )
	    return branchIdA;
	}
    }
  return nullBranchId;
}

std::tuple<MVertex*, std::vector<int> >
AngioTkCenterline::vertexOnNeighboringBranch( MVertex* vertexA, MVertex* vertexB )
{
  std::vector<int> branchIds;
  MVertex* junc = NULL;
  auto itFindVertexA = M_vertexToLinesId.find( vertexA );
  auto itFindVertexB = M_vertexToLinesId.find( vertexB );
  for ( auto const& vertexInfoA : itFindVertexA->second )
    {
      int branchIdA =  vertexInfoA.first;
      for ( auto const& vertexInfoB : itFindVertexB->second )
	{
	  int branchIdB = vertexInfoB.first;
	  if ( ( edges[branchIdA].vB == edges[branchIdB].vB ) ||
	       ( edges[branchIdA].vB == edges[branchIdB].vE ) ||
	       ( edges[branchIdA].vE == edges[branchIdB].vB ) ||
	       ( edges[branchIdA].vE == edges[branchIdB].vE ) )
	    {
	      branchIds.push_back(branchIdA);
	      branchIds.push_back(branchIdB);
	      if ( ( edges[branchIdA].vB == edges[branchIdB].vB ) ||
		   ( edges[branchIdA].vB == edges[branchIdB].vE ) )
		junc = edges[branchIdA].vB;
	      else
		junc = edges[branchIdA].vE;
	      return std::make_tuple(junc,branchIds);
	    }
	}
    }
  return std::make_tuple(junc,branchIds);
}

bool
AngioTkCenterline::canFindPathBetweenVertex( MVertex* vertexA, MVertex* vertexB )
{
  //return this->vertexOnSameBranch(vertexA, vertexB);

  auto itFindVertexA = M_vertexToLinesId.find( vertexA );
  auto itFindVertexB = M_vertexToLinesId.find( vertexB );
  if ( itFindVertexA == M_vertexToLinesId.end() || itFindVertexB == M_vertexToLinesId.end() )
    return false;

#if 0
  if ( ( M_junctionsVertex.find( vertexA ) != M_junctionsVertex.end() ) ||
       ( M_junctionsVertex.find( vertexB ) != M_junctionsVertex.end() ) )
    return false;
#endif
  for ( auto const& vertexInfoA : itFindVertexA->second )
    {
      int branchIdA =  vertexInfoA.first;
      for ( auto const& vertexInfoB : itFindVertexB->second )
	{
	  int branchIdB = vertexInfoB.first;
	  if ( branchIdB == branchIdA )
	    return true;
	  else if ( ( edges[branchIdA].vB == edges[branchIdB].vB ) ||
		    ( edges[branchIdA].vB == edges[branchIdB].vE ) ||
		    ( edges[branchIdA].vE == edges[branchIdB].vB ) ||
		    ( edges[branchIdA].vE == edges[branchIdB].vE ) )
	    {
	      return true;
	    }
	}
    }
    return false;
}

std::tuple< std::vector<MLine*> , double >
AngioTkCenterline::pathBetweenVertex( MVertex* vertexA, MVertex* vertexB )
{
  //std::cout << "pathBetweenVertex A : " << vertexA->getIndex() << " B : " << vertexB->getIndex() << "\n";

  std::vector<MLine*> listOfLines;
  double lengthPath = 0;

  // case same points
  if ( vertexA == vertexB )
    return std::make_tuple(listOfLines,lengthPath);

  int branchIdCommon = this->vertexOnSameBranch( vertexA,vertexB );
  if ( branchIdCommon >= 0 )
    {
      // case same branch
      std::set<int> ptDone;
      int state=0;
      for ( auto myline : edges[branchIdCommon].lines )
	{
	  std::vector<MVertex*> ptsInLine = { myline->getVertex(0),myline->getVertex(1) };
	  for ( MVertex* myvertex : ptsInLine )
	    {
	      if ( ptDone.find(myvertex->getIndex()) == ptDone.end() || M_junctionsVertex.find( myvertex ) != M_junctionsVertex.end() )
		{
		  if ( vertexA == myvertex || vertexB == myvertex )
		      ++state;
		  ptDone.insert(myvertex->getIndex());
		  if ( state == 2 )
		    break;
		}
	    }
	  if ( state == 1 )
	    {
	      listOfLines.push_back(myline);
	      lengthPath += myline->getLength();
	    }
	  else if ( state == 2 )
	    {
	      listOfLines.push_back(myline);
	      lengthPath += myline->getLength();
	      return std::make_tuple(listOfLines,lengthPath);
	    }
	}
      if ( state != 2 )
	Msg::Error("failure : path computation on same branch");
    } // if ( branchIdCommon >= 0 )
  else
    {
      auto neighborBranchTest = this->vertexOnNeighboringBranch(vertexA,vertexB);
      auto const& neighborBranchIds = std::get<1>( neighborBranchTest );
      if ( !neighborBranchIds.empty() )
	{
	  // case neighboring branch
	  MVertex *junc = std::get<0>( neighborBranchTest );
	  auto path1 = this->pathBetweenVertex(vertexA,junc);
	  auto path2 = this->pathBetweenVertex(vertexB,junc);
	  listOfLines.insert(listOfLines.end(),std::get<0>(path1).begin(),std::get<0>(path1).end());
	  listOfLines.insert(listOfLines.end(),std::get<0>(path2).begin(),std::get<0>(path2).end());
	  lengthPath = std::get<1>(path1) + std::get<1>(path2);
	  return std::make_tuple(listOfLines,lengthPath);
	}
      else
      	Msg::Error("TODO : point pair are not localized in same branch and not have a near connection, must be implemented");

    }

  return std::make_tuple(listOfLines,lengthPath);
}

double AngioTkCenterline::maxScalarValueInPath( std::vector<MLine*> const& path, std::string const& fieldName ) const
{
  double res=0;
  if ( !this->hasField( fieldName ) )
    return res;

  std::set<MVertex*> vertexDone;
  for ( MLine* myline : path )
    {
      std::vector<MVertex*> ptsInLine = { myline->getVertex(0),myline->getVertex(1) };
      for ( MVertex* myvertex : ptsInLine )
	{
	  if ( vertexDone.find( myvertex ) != vertexDone.end() )
	    continue;

	  res = std::max(res, this->centerlinesFieldsPointData(fieldName,myvertex)[0] );
	  vertexDone.insert( myvertex );
	}
    }
  return res;
}

void AngioTkCenterline::applyTubularColisionFix( std::vector<MVertex*> const& vTestedSet )
{
  std::set<int> branchToReApplyTubularColisionFix;
  double spaceMinBetweenTubularStructure = 0.5;//2*0.5;//2*1.2;//0.6;//0.5;
  double radiusMinAllowed = 0.5;

  std::map< MVertex*,std::set<MVertex*> > mapVertexTested;
  for ( MVertex* vTested : vTestedSet )
    {
      double radius = M_centerlinesFieldsPointData["RadiusMin"][M_mapVertexGmshIdToVtkId[vTested->getIndex()]][0];

      std::vector<MVertex*> vTestedSetInSearch;
      std::set<int> ptDoneInSearch;
      for( int i = 0; i < edges.size(); ++i )
	{
#if 1
	  if ( !edges[i].isInsideBox( vTested, spaceMinBetweenTubularStructure+2*radius ) )
	    continue;
#endif

	  for ( auto mylineInSearch : edges[i].lines )
	    {
	      MVertex* v0InSearch = mylineInSearch->getVertex(0);
	      MVertex* v1InSearch = mylineInSearch->getVertex(1);
	      if ( ptDoneInSearch.find(v0InSearch->getIndex()) == ptDoneInSearch.end() )
		{
		  if ( this->canFindPathBetweenVertex( vTested, v0InSearch ) )
		    {
		      auto mypathSearch = this->pathBetweenVertex( vTested, v0InSearch );
		      double lengthPathSearch = std::get<1>( mypathSearch );
		      if ( lengthPathSearch > 10*this->maxScalarValueInPath(std::get<0>(mypathSearch),"RadiusMin") /*12*//*2*5*edges[i].maxRad*/ )
			mapVertexTested[ vTested ].insert( v0InSearch );
		    }
		  else
		    mapVertexTested[ vTested ].insert( v0InSearch );
		  ptDoneInSearch.insert(v0InSearch->getIndex());
		}
	      if ( ptDoneInSearch.find(v1InSearch->getIndex()) == ptDoneInSearch.end() )
		{
		  if ( this->canFindPathBetweenVertex( vTested, v1InSearch ) )
		    {
		      auto mypathSearch = this->pathBetweenVertex( vTested, v1InSearch );
		      double lengthPathSearch = std::get<1>( mypathSearch );
		      if ( lengthPathSearch > 10*this->maxScalarValueInPath(std::get<0>(mypathSearch),"RadiusMin")/*12*//*2*5*edges[i].maxRad*/ )
			mapVertexTested[ vTested ].insert( v1InSearch );
		    }
		  else
		    mapVertexTested[ vTested ].insert( v1InSearch );
		  ptDoneInSearch.insert(v1InSearch->getIndex());
		}
	    }
	} // for(int i = 0; i < edges.size(); ++i)

    } // for ( MVertex* vTested : vTestedSet )

  //this->applyTubularColisionFix( mapVertexTested,1,2 );
  //this->applyTubularColisionFix( mapVertexTested,0,-1 );
  //this->applyTubularColisionFix( mapVertexTested,2,1 );
  this->applyTubularColisionFix( mapVertexTested,1,1 );
  this->applyTubularColisionFix( mapVertexTested,0,-1 );

}


void
AngioTkCenterline::applyTubularColisionFix( std::map< MVertex*,std::set<MVertex*> > const& mapVertexTested, int method, int maxrecurrence, int nrecurrence )
{
  if ( maxrecurrence > 0 && nrecurrence >= maxrecurrence )
    return;
  //if (method == 0)
  //std::cout <<"applyTubularColisionFix start : " << mapVertexTested.size() << " nrecurrence="<<nrecurrence <<"\n";
  double spaceMinBetweenTubularStructure = 2*0.5;//1.5*0.5;//2*0.5;//2*1.2;//0.6;//0.5;
  std::map< MVertex*,std::set<MVertex*> > newMapVertexInColision;
  std::vector<std::pair<double,SPoint3>> newDirs;
  for ( auto const& vertexTestedPair : mapVertexTested )
    {
      MVertex* vTested = vertexTestedPair.first;
      double distTubeColisionMin = 0;
      double radius = M_centerlinesFieldsPointData["RadiusMin"][M_mapVertexGmshIdToVtkId[vTested->getIndex()]][0];
      double newRadius = radius;
      MVertex* vTestedInSearchMin = NULL;
      for ( MVertex* vTestedInSearch : vertexTestedPair.second )
	{
	  double distBetweenPoints = vTested->point().distance( vTestedInSearch->point() );
	  double radiusInSearch = M_centerlinesFieldsPointData["RadiusMin"][M_mapVertexGmshIdToVtkId[vTestedInSearch->getIndex()]][0];
	  double distTubeColision = distBetweenPoints - ( radius + radiusInSearch + spaceMinBetweenTubularStructure );

	  if ( /*false &&*/ distTubeColision < 0 )
	    {
	      newMapVertexInColision[vTestedInSearch].insert( vTested );
	    }
	  if ( method == 2 && distTubeColision < 0 )
	    {
		  SPoint3 dir = vTested->point()-vTestedInSearch->point();
		  double normDir = std::sqrt( std::pow(dir.x(),2) + std::pow(dir.y(),2) + std::pow(dir.z(),2));
		  SPoint3 dirNormalized = ((1./normDir)*dir).point();
		  newDirs.push_back( std::make_pair( -distTubeColision,dirNormalized ) );
	    }


	  if ( distTubeColision < distTubeColisionMin )
	    {
	      distTubeColisionMin = distTubeColision;
	      newRadius = distBetweenPoints - ( radiusInSearch + spaceMinBetweenTubularStructure );
	      vTestedInSearchMin = vTestedInSearch;
	      //std::cout << "newRadius " << newRadius <<"\n";
	      if ( method == 0 && newRadius < 0  )
		Msg::Error("failure : newRadius not usable : %g",newRadius );
	      double radiusMinAllowed = 0.8*radius;
	      if ( newRadius < radiusMinAllowed )
		  newRadius = radiusMinAllowed;
	    }
	}

      if ( distTubeColisionMin < 0 )
	{
	  if ( method == 0 )
	    {
	      //CHECK( newRadius > 1e-5 ) << "new radius must be positive : " << newRadius;
	      M_centerlinesFieldsPointData["RadiusMin"][M_mapVertexGmshIdToVtkId[vTested->getIndex()]][0] = newRadius;
	    }
	  else if ( method == 1 )
	    {
	      if ( vTestedInSearchMin )
		{
		  SPoint3 dir = vTested->point()-vTestedInSearchMin->point();
		  double normDir = std::sqrt( std::pow(dir.x(),2) + std::pow(dir.y(),2) + std::pow(dir.z(),2));
#if 0
		  double distBetweenPoints = vTested->point().distance( vTestedInSearchMin->point() );
		  //SPoint3 dirToUse = (((distBetweenPoints-distTubeColisionMin)/normDir)*dir).point();
		  SPoint3 dirToUse = (((distBetweenPoints+std::abs(newRadius-radius))/normDir)*dir).point();
		  //SPoint3 dirToUse = (-(std::abs(newRadius-radius)/normDir)*dir).point();
		  SPoint3 newPos = vTestedInSearchMin->point() + dirToUse;
#else
		  SPoint3 dirToUse = ((-distTubeColisionMin/normDir)*dir).point();
		  SPoint3 newPos = vTested->point() + dirToUse;
#endif
		  vTested->setXYZ(newPos.x(),newPos.y(),newPos.z());
		}
	    }
	  else if ( method == 2 && !newDirs.empty())
	    {
	      SPoint3 dirToUse(0.,0.,0.);
	      for ( auto const& newDir : newDirs )
		{
		  dirToUse += (newDir.first*newDir.second).point();
		}
	      dirToUse *= (1./newDirs.size());
	      SPoint3 newPos = vTested->point() + dirToUse;
	      vTested->setXYZ(newPos.x(),newPos.y(),newPos.z());
	    }
	}
    } // for ( MVertex* vTested : vTestedSet )

  //std::cout <<"applyTubularColisionFix finish : " << newMapVertexInColision.size() << "\n";

  if ( !newMapVertexInColision.empty() )
    {
	  //if ( nrecurrence < 10 )
	  this->applyTubularColisionFix( newMapVertexInColision, method, maxrecurrence, nrecurrence+1 );
    }


}


void AngioTkCenterline::applyTubularColisionFix( AngioTk::pointpair_data_type const& pointPairData )
{
  Msg::Info("AngioTkCenterline::applyTubularColisionFix");
  double spaceMinBetweenTubularStructure = 0.5;//2*0.5;//2*1.2;//0.6;//0.5;

  this->buildKdTree();

  for ( auto const& pointPair : pointPairData )
    {
      auto const& pt1 = std::get<0>( pointPair.first );
      auto const& pt2 = std::get<0>( pointPair.second );
      double ptToLocalize1[3] = { pt1[0], pt1[1],pt1[2] };
      double ptToLocalize2[3] = { pt2[0], pt2[1],pt2[2] };
      std::cout << "ptToLocalize1 : " << ptToLocalize1[0] << " , " << ptToLocalize1[1]  << " , " << ptToLocalize1[2] << "\n";
      std::cout << "ptToLocalize2 : " << ptToLocalize2[0] << " , " << ptToLocalize2[1]  << " , " << ptToLocalize2[2] << "\n";
      MVertex* vertexA = std::get<0>( this->foundClosestPointInCenterlines(ptToLocalize1) );
      MVertex* vertexB = std::get<0>( this->foundClosestPointInCenterlines(ptToLocalize2) );
      std::cout << "vertexA : " << vertexA->x() << " , " << vertexA->y()  << " , " << vertexA->z() << "\n";
      std::cout << "vertexB : " << vertexB->x() << " , " << vertexB->y()  << " , " << vertexB->z() << "\n";

      auto mypath = this->pathBetweenVertex( vertexA,vertexB );

      std::vector<MLine*> listOfLines = std::get<0>( mypath );
      std::set<int> ptDone;
      for ( auto myline : listOfLines )
	{
	  MVertex* v0 = myline->getVertex(0);
	  MVertex* v1 = myline->getVertex(1);
	  std::vector<MVertex*> vTestedSet;
	  if ( ptDone.find(v0->getIndex()) == ptDone.end() )
	    {
	      vTestedSet.push_back( v0 );
	      ptDone.insert(v0->getIndex());
	    }
	  if ( ptDone.find(v1->getIndex()) == ptDone.end() )
	    {
	      vTestedSet.push_back( v1 );
	      ptDone.insert(v1->getIndex());
	    }

	  this->applyTubularColisionFix( vTestedSet );
	}
    }
}

void AngioTkCenterline::applyTubularColisionFix()
{
  Msg::Info("AngioTkCenterline::applyTubularColisionFix");

  this->buildKdTree();

  for( int i = 0; i < edges.size(); ++i )
    {
      std::set<int> ptDone;
      std::vector<MVertex*> vTestedSet;
      for ( auto myline : edges[i].lines )
	{
	  MVertex* v0 = myline->getVertex(0);
	  MVertex* v1 = myline->getVertex(1);
	  if ( ptDone.find(v0->getIndex()) == ptDone.end() )
	    {
	      vTestedSet.push_back( v0 );
	      ptDone.insert(v0->getIndex());
	    }
	  if ( ptDone.find(v1->getIndex()) == ptDone.end() )
	    {
	      vTestedSet.push_back( v1 );
	      ptDone.insert(v1->getIndex());
	    }
	}
      Msg::Info("AngioTkCenterline::applyTubularColisionFix start for branch id %d",i);
      this->applyTubularColisionFix( vTestedSet );
    }
}



// tubular extension
void AngioTkCenterline::runTubularExtension()
{
    current->createTopologyFromMesh();

    //identify the boundary edges by looping over all discreteFaces
    std::vector<GEdge*> boundEdges;
    int NF = current->getMaxElementaryNumber(2);
    for (int i= 0; i< NF; i++)
    {
        GFace *gf = current->getFaceByTag(i+1);
        std::list<GEdge*> l_edges = gf->edges();
        for(std::list<GEdge*>::iterator it = l_edges.begin(); it != l_edges.end(); it++)
        {
            std::vector<GEdge*>::iterator ite = std::find(boundEdges.begin(),boundEdges.end(), *it);
            if (ite != boundEdges.end())
	      boundEdges.erase(ite);
            else
	      boundEdges.push_back(*it);
        }
    }

#if 0
    // this code below allow to generate the extrusion of inlet edge in the mesh
    if ( !boundEdges.empty() )
    for ( int i=0;i<1;++i)
    {
        std::vector<double> p1( { 0,0,0 } );
        std::vector<double> p2( { 1,1,0 } );
        //current->extrude( boundEdges[i],p1,p2 );
        GEdge * gec = current->getEdgeByTag(boundEdges[i]->tag());
        std::vector<MLine*> lines = gec->lines;
        //std::cout << "lines.size() " << lines.size() << "\n";
        int numEdge = current->getMaxElementaryNumber(1) + 1;
        discreteEdge *discEdge = new discreteEdge(current,numEdge,0,0);
        current->add(discEdge);
        GEntity* discEntity = (GEntity*)discEdge;
        discEdge->setTag( 12345/*itNewTag->second*/ );

        std::map<MVertex*,MVertex*> vertexExtruded;
        for ( int l=0;l< lines.size();++l)
        {
            MVertex * v0 = lines[l]->getVertex(0);
            MVertex * v1 = lines[l]->getVertex(1);
            MVertex * v0New;
            if ( vertexExtruded.find(v0) != vertexExtruded.end() )
                v0New = vertexExtruded.find(v0)->second;
            else
            {
                v0New = new MVertex( v0->x()+0.,v0->y()+0.,v0->z()-1., discEntity/*gec*//*NULL*/ );
                v0New->setIndex( v0New->getNum() );
                vertexExtruded[v0]=v0New;
                discEdge->addMeshVertex(v0New);
            }
            MVertex * v1New;
            if ( vertexExtruded.find(v1) != vertexExtruded.end() )
                v1New = vertexExtruded.find(v1)->second;
            else
            {
                v1New = new MVertex( v1->x()+0.,v1->y()+0.,v1->z()-1., discEntity/*gec*/ );
                v1New->setIndex( v1New->getNum() );
                vertexExtruded[v1]=v1New;
                discEdge->addMeshVertex(v1New);
            }
            //std::cout << "v0New->getIndex() " << v0New->getIndex() << " v0New->getNum() " << v0New->getNum() << "\n";
            //std::cout << "v1New->getIndex() " << v1New->getIndex() << " v1New->getNum() " << v1New->getNum() << "\n";
            CHECK( v0New ) << " not has v0New";
            CHECK( v1New ) << " not has v1New";
            //MLine * lineNew = new MLine(v0New,v1New);
            //discEdge->lines.push_back(lineNew);
            discEdge->lines.push_back(new MLine(v0New,v1New) );
        }
    }
    //current->removeDuplicateMeshVertices(1.e-8);
    current->createTopologyFromMesh();
#endif

    std::vector<std::vector<GEdge *> > myEdgeLoops;
    std::vector<GEdge *> myEdges;
    GEdge * gec = current->getEdgeByTag( /*12345*/boundEdges[0]->tag() );
    myEdges.push_back(gec);
    myEdgeLoops.push_back(myEdges);
    GFace *newFace = current->addPlanarFace(myEdgeLoops);

    int nbElemLayer=3;double hLayer=-0.7;  int dir = 0;
    std::vector<GEntity*> extrudedE = current->extrudeBoundaryLayer(newFace, nbElemLayer,  hLayer, dir, -1/* -5*/);

    current->mesh(2);
    current->remove( (GRegion*) extrudedE[1]);
    current->remove( (GFace*) extrudedE[0]);
    current->remove( (GFace*) newFace);

    //current->exportDiscreteGEOInternals();
    //current->createTopologyFromMesh();
#if 0
    //current->writeSTL( this->outputPath(),  false, false);
    current->writeMSH("myCLIPPARTS.msh", 2.2, false/*binary*/, true/*false*//*saveAll*/);
#endif
}



#endif // ANN
