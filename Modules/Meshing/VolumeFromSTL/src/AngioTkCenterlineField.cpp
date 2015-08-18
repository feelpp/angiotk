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
				std::multimap<MEdge, MTriangle*, Less_Edge> &e2e,
				std::set<MTriangle*> &group,
				std::set<MEdge, Less_Edge> &touched,
				std::set<MEdge, Less_Edge> &theCut)
{
  if (touched.find(e) != touched.end()) return;
  touched.insert(e);
  for (std::multimap <MEdge, MTriangle*, Less_Edge>::iterator it = e2e.lower_bound(e);
       it != e2e.upper_bound(e); ++it){
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

void AngioTkCenterline::importSurfaceFromFile(std::string const& fileName )
{
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

void AngioTkCenterline::importFile(std::string fileName)
{
  current = GModel::current();
  //current->removeDuplicateMeshVertices(1.e-8);
  //if ( !current ) current = new GModel();
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


  //std::set<std::pair<int,int> > registerLinesToRemoveFromPointIdPair;
  //std::set<std::pair<int,int> > M_registerLinesToRemoveFromPointIdPairInModelEdge;
  std::shared_ptr<AngioTkCenterline> newCenterlines;
  // save vertex index in a std::set
  std::set<int > ignoreIndexReplaced;
  int previousMaxVertexIndex=0;
  std::map<GEntity*,std::set<MVertex*> > fusionBranch;

  if ( !modEdges.empty() )
    {
      Msg::Info("AngioTkCenterline: importFile : merge centerlines start");

      //std::shared_ptr<AngioTkCenterline> newCenterlines( new AngioTkCenterline );
      newCenterlines.reset( new AngioTkCenterline );
      //newCenterlines->importFile( fileName );
      //readVTKPolyDataFields( fileName,newCenterlines->centerlinesFieldsPointData );
      newCenterlines->updateCenterlinesFromFile( fileName );
      this->attachAngioTkCenterline(newCenterlines);

      this->buildKdTree();


      std::map<MVertex*,std::pair<std::pair<MLine*,int>, MVertex*> > indexVertexReplaced;
      std::set<std::pair<int,int> > newRegisterLinesToRemove;
      std::set<MLine*> newRegisterLinesToRemove2;
      //if (false)
      for ( auto const& extremityPair : newCenterlines->centerlinesExtremities() )
	{
	  double ptToLocalize[3] = { extremityPair.first->x(),extremityPair.first->y(),extremityPair.first->z() };
	  ANNidx index[1];
	  ANNdist dist[1];
	  kdtree->annkSearch(ptToLocalize, 1, index, dist);
	  ANNpointArray nodes = kdtree->thePoints();


	  int branchId = extremityPair.second.first;
	  int lineIdInBranch = extremityPair.second.second;
	  std::vector<MLine*> mylines = newCenterlines->centerlinesBranch(branchId).lines;
	  MLine* myline = mylines[lineIdInBranch];
	  if ( newCenterlines->centerlinesRadiusl().find(myline) == newCenterlines->centerlinesRadiusl().end() )
	    Msg::Error("radius not find for this line \n");
	  double radius = newCenterlines->centerlinesRadiusl().find(myline)->second;
#if 0
	  std::cout << "extremity ptToLocalize " << ptToLocalize[0] <<","<< ptToLocalize[1] <<","<< ptToLocalize[2]<<"\n";
	  std::cout << "Find KdTree " << nodes[index[0]][0] <<","<< nodes[index[0]][1]<<","<< nodes[index[0]][2] << " dist " << dist[0] << "\n";
	  std::cout << "myRadius " << radius << "\n";
#endif

	  if ( lineIdInBranch != 0 && lineIdInBranch != (mylines.size()-1) )
	    Msg::Error("is not an extremity\n");

	  //if ( 4*radius > dist[0] )
	  //{
	      // get point in colorp (with respect to the initialisation of kdtree)
	      MVertex* pointFound = NULL;
	      if ( index[0] < colorp.size() )
		{
		  std::map<MVertex*, int>::iterator itp = colorp.begin();
		  for (int k=0;k<index[0];++k)
		    ++itp;
		  //std::cout << " itp->first; " << itp->first->x() << "," << itp->first->y() << "," << itp->first->z() <<"\n";
#if 0
		  if ( 4*radius > dist[0] )
		    pointFound = itp->first;
#else
		  if ( 2*(radius+this->minRadiusAtVertex(itp->first) ) > dist[0] )
		    pointFound = itp->first;
		  //if ( 30*(/*radius*/1 ) > dist[0] )
		  //pointFound = itp->first;
#endif
		  //if ( centerlinesFieldsPointData. )
		}
	      else
		{
#if 1
		  int nbPL = 3;  //10 points per line
		  int ind = colorp.size();
		  //bool findPointt = false;
		  for(unsigned int k = 0; k < lines.size(); ++k){
		    MVertex *v0 = lines[k]->getVertex(0);
		    MVertex *v1 = lines[k]->getVertex(1);
		    //SVector3 P0(v0->x(),v0->y(), v0->z());
		    //SVector3 P1(v1->x(),v1->y(), v1->z());
		    for (int j= 1; j < nbPL+1; j++){
		      if ( ind==index[0] ) {
			//if ( 4*radius > dist[0] )
			//  pointFound = v0;
			if ( 2*(radius+this->minRadiusAtVertex(v0) ) > dist[0] )
			  pointFound = v0;
			  //if ( 30*(1/*radius*/ ) > dist[0] )
			//pointFound = v0;

			break;
		      }
			//double inc = (double)j/(double)(nbPL+1);
			//SVector3 Pj = P0+inc*(P1-P0);
		      ind++;
		    }
		    if ( pointFound != NULL ) break;
		  }

#endif
		}

	      if ( false )
		{
		  std::cout << "Find KdTree " << nodes[index[0]][0] <<","<< nodes[index[0]][1]<<","<< nodes[index[0]][2] << " dist " << dist[0] << "\n";
		  //std::cout << "extremity ptToLocalize " << ptToLocalize[0] <<","<< ptToLocalize[1] <<","<< ptToLocalize[2]<<"\n";
		  std::cout << "pointFound  " << pointFound->x() <<","<< pointFound->y() <<","<< pointFound->z()<<"\n";
		  std::cout << "myRadius " << radius << "\n";
		}
#if 0
	      if ( pointFound == NULL )
		std::cout << "point not found with dist " << dist[0] << "\n";
	      else
		std::cout << "point found\n";
#endif
	      if ( pointFound != NULL )
		{

	      MVertex* _ptToReplaced = extremityPair.first;
	      MLine* _lineWhichHasToReplaced = myline;
	      int _idVertexInLine = ( extremityPair.first == myline->getVertex(0) )? 0 : 1;
	      //if ( myline->getVertex(_idVertexInLine) != extremityPair.first ) Msg::Error("Error localToGlobalVertex");

	      MVertex* _mylink = ( extremityPair.first == myline->getVertex(0) )? myline->getVertex(1) : myline->getVertex(0);

	      // search other point closer (near the extremity point) 
	      if ( true )
		{
		  double curMinDist = dist[0];
		  bool isForwardSearch = ( lineIdInBranch == 0 );
		  double lcTotal = 0;
		  int indexSearch = 0;
		  for ( int q=1; ( lcTotal < 4*radius && q < (mylines.size()-1) ) ; ++q )
		    {
		      int newLineId = (isForwardSearch)?lineIdInBranch+q : lineIdInBranch-q;
		      MLine* mylineSearch = mylines[newLineId];
		      int vIdInLine = (isForwardSearch)?0 : 1;
		      MVertex *newV = mylineSearch->getVertex( vIdInLine );
		      if ( _mylink != newV ) Msg::Error("NOT GOOD \n");
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
		  if ( indexSearch > 0 )
		    {
		      for ( int q=0; q<indexSearch ; ++q )
			{
			  int newLineId = (isForwardSearch)?lineIdInBranch+q : lineIdInBranch-q;
			  MLine* mylineSearch = mylines[newLineId];
			  int v0Id = mylineSearch->getVertex(0)->getIndex();
			  int v1Id = mylineSearch->getVertex(1)->getIndex();
			  ignoreIndexReplaced.insert( (isForwardSearch)? v0Id : v1Id );

			  newRegisterLinesToRemove2.insert( mylineSearch );
			  newRegisterLinesToRemove.insert( std::make_pair(v0Id,v1Id) );
			  newRegisterLinesToRemove.insert( std::make_pair(v1Id,v0Id) );
			  //M_registerLinesToRemoveFromPointIdPairInModelEdge.insert( std::make_pair(v0Id,v1Id) );
			  //M_registerLinesToRemoveFromPointIdPairInModelEdge.insert( std::make_pair(v1Id,v0Id) );
			}
		    }
		  //std::cout << "M_registerLinesToRemoveFromPointIdPairInModelEdge.size() " << M_registerLinesToRemoveFromPointIdPairInModelEdge.size() << "\n";
		}

	      ignoreIndexReplaced.insert( _ptToReplaced->getIndex() );
	      indexVertexReplaced[_ptToReplaced] = std::make_pair( std::make_pair(_lineWhichHasToReplaced,_idVertexInLine) ,pointFound );

	      bool ptFoundIsExtremity = false;
	      for ( auto const& previousExtremityPair : this->centerlinesExtremities() )
		{
		  if ( previousExtremityPair.first ==  pointFound )
		    {
		      ptFoundIsExtremity = true;
		      break;
		    }
		}
	      if ( ptFoundIsExtremity )
		{
		  //fusionBranch[_ptToReplaced] = pointFound;
		  //std::map<GEntity*,std::set<GEntity*> > fusionBranch;
		  //std::map<GEntity*,std::set<MVertex*> > fusionBranch;
		  //fusionBranch[ _ptToReplaced->onWhat()/*->model()*/ ].insert( pointFound->onWhat()/*->model()*/ );
		  fusionBranch[ _ptToReplaced->onWhat() ].insert( pointFound );
		}


	    }
	}



      //----------------------------------------------------------------//
      // save vertex index in a std::set
      //std::set<int > ignoreIndexReplaced;
#if 0
      for ( auto idReplacePair : indexVertexReplaced )
	{
	  std::cout << "ignoreIndexReplaced. add " << idReplacePair.first->getIndex() << "\n";
	  ignoreIndexReplaced.insert(idReplacePair.first->getIndex());
	}
#endif
      //std::cout << "ignoreIndexReplaced.size " << ignoreIndexReplaced.size() << "\n";

      //----------------------------------------------------------------//
      // search max tag already register
      int previousMaxTag = 0;
      for (unsigned int i = 0; i < modEdges.size(); i++){
	GEdge *ge = modEdges[i];
	previousMaxTag = std::max(previousMaxTag, ge->tag());
	//std::cout << "ge->tag() " << ge->tag() << "\n";
      }

      previousMaxVertexIndex = maxVerticesIndex( this->edges );
      //std::cout << previousMaxVertexIndex << " VS " << numberOfVertices( this->edges ) << "\n";

      std::set<int> ptDone;
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
	    if ( ptDone.find(v0->getNum()) == ptDone.end() )
	      {
		ptDone.insert(v0->getNum());
		int indexBefore = v0->getIndex();
		v0->setIndex( previousMaxVertexIndex + v0->getIndex() );
		//std::cout << "new vertex0 index " << v0->getIndex() << "\n";
	      }
	    if ( ptDone.find(v1->getNum()) == ptDone.end() )
	      {
		ptDone.insert(v1->getNum());
		v1->setIndex( previousMaxVertexIndex + v1->getIndex() );
		//std::cout << "new vertex1 index " << v1->getIndex() << "\n";
	      }
	  }
      }

      //---------------------------------------------//
      // replace points in line
      for ( auto idReplacePair : indexVertexReplaced )
	{
	  MLine* lineModified = idReplacePair.second.first.first;
	  int ptIdInLine = idReplacePair.second.first.second;
	  MVertex* newPt = idReplacePair.second.second;

	  lineModified->setVertex(ptIdInLine,newPt);
	  if ( false)// lineModified->getLength() > 5 )
	    {
	      std::cout << "lineModified->getLength() "<< lineModified->getLength() << " with ptIdInLine " << ptIdInLine <<"\n";
	      std::cout << " initPt; " << idReplacePair.first->x() << "," << idReplacePair.first->y() << "," << idReplacePair.first->z() <<"with index " << idReplacePair.first->getIndex()<<"\n";
	      std::cout << " newPt; " << newPt->x() << "," << newPt->y() << "," << newPt->z() << "with index " << newPt->getIndex()<<"\n";
	      std::cout << " newline vertex0 " << lineModified->getVertex(0)->x() << "," << lineModified->getVertex(0)->y() << "," << lineModified->getVertex(0)->z() <<"\n"
			<< " newline vertex1 " << lineModified->getVertex(1)->x() << "," << lineModified->getVertex(1)->y() << "," << lineModified->getVertex(1)->z() <<"\n";

	    }
	}
      for ( auto myLine : newRegisterLinesToRemove2 )
	{
	  int v0Id = myLine->getVertex(0)->getIndex(), v1Id=myLine->getVertex(1)->getIndex();
	  M_registerLinesToRemoveFromPointIdPairInModelEdge.insert( std::make_pair(v0Id,v1Id) );
	  M_registerLinesToRemoveFromPointIdPairInModelEdge.insert( std::make_pair(v1Id,v0Id) );
	}
      //---------------------------------------------//

      //modEdges.insert(modEdges.end(),mod->firstEdge(), mod->lastEdge());
      modEdges.insert(modEdges.end(),newCenterlines->mod->firstEdge(), newCenterlines->mod->lastEdge());
      //std::cout << "new distance" << std::distance( newCenterlines->mod->firstEdge(), newCenterlines->mod->lastEdge())<<"\n";
      Msg::Info("AngioTkCenterline: importFile : merge centerlines finish");
    }
  else
    {
      //mod = new GModel();
      mod.reset( new GModel() );
      mod->load(fileName);
      mod->removeDuplicateMeshVertices(1.e-8);
      // call indexMeshVertices allow to assign id on vertices from 1
      bool saveAll = false;
      int numVerticesInCenterlinesGModel = mod->indexMeshVertices(saveAll);
      //mod->writeMSH("myCenterlines.msh", 2.2, false, false);
      //mod->writeVTK("myCenterlines.vtk", false, false);
#if 1
      current->setAsCurrent();
      current->setVisibility(1);
#endif
      //std::cout<< "number of vertice in current GModel "<< current->indexMeshVertices(true) << "\n";
      modEdges.insert(modEdges.end(),mod->firstEdge(), mod->lastEdge());
    }


  // update centerlines : data branch, extremities, kdtree,...
  this->updateCenterlinesForUse();

  //---------------------------------------------------------------------------------------------//

  if ( !newCenterlines )
    {
      this->updateRelationMapVertex();
      this->updateCenterlinesFieldsFromFile(fileName);
    }

  // in case of merging centerline : update centerlinesFieldsPointData
  if ( newCenterlines )
    {
      // get new mapping between GmshId and VtkId
      std::map<int,int>  _mapVertexGmshIdToVtkId,_mapVertexVtkIdToGmshId;
      this->updateRelationMapVertex( _mapVertexGmshIdToVtkId,_mapVertexVtkIdToGmshId);

      // build new storage for fields points data
      std::map<std::string,std::vector<std::vector<double> > >  newCenterlinesFieldsPointData;

      // add info from current centerline merged
      if ( newCenterlines  )
	{
	  for (auto const& fieldsPointDataPair : newCenterlines->centerlinesFieldsPointData )
	    {
	      if ( centerlinesFieldsPointData.find(fieldsPointDataPair.first) == centerlinesFieldsPointData.end() ) continue;

	      newCenterlinesFieldsPointData[fieldsPointDataPair.first].clear();
	      newCenterlinesFieldsPointData[fieldsPointDataPair.first].resize( _mapVertexVtkIdToGmshId.size() );
	      auto const& vecNewPointData = fieldsPointDataPair.second;
	      for (int k=0;k<newCenterlines->M_mapVertexVtkIdToGmshId.size();++k)
		{
		  if ( ignoreIndexReplaced.find( newCenterlines->M_mapVertexVtkIdToGmshId[k] ) != ignoreIndexReplaced.end() ) { ;continue; }

		  auto itFindVtkId =  _mapVertexGmshIdToVtkId.find( (int)(previousMaxVertexIndex+newCenterlines->M_mapVertexVtkIdToGmshId[k]) );
		  if ( itFindVtkId == _mapVertexGmshIdToVtkId.end() )
		    Msg::Error("_mapVertexGmshIdToVtkId not find %d %d and %d %d",
			       previousMaxVertexIndex,newCenterlines->M_mapVertexVtkIdToGmshId[k], _mapVertexGmshIdToVtkId.size(),_mapVertexVtkIdToGmshId.size() ); 
		  int idContainer = itFindVtkId->second;

		  //std::cout << "idContainer " << idContainer << "with " << centerlinesFieldsPointData[fieldsPointDataPair.first].size() << " and " << _mapVertexGmshIdToVtkId.size()<<"\n";
		  newCenterlinesFieldsPointData[fieldsPointDataPair.first][idContainer] = vecNewPointData[k];
		}
	    }
	}

      // delete datas when fieldname is not presence in all centerlines
      std::set<std::string> fieldsToErase;
      for (auto const& fieldsPointDataPair : centerlinesFieldsPointData )
	{
	  if ( newCenterlines->centerlinesFieldsPointData.find(fieldsPointDataPair.first) == newCenterlines->centerlinesFieldsPointData.end() )
	    fieldsToErase.insert( fieldsPointDataPair.first );
	}
      for ( std::string const fieldName : fieldsToErase )
	centerlinesFieldsPointData.erase(fieldName);

      // add info from previous centerlines (already merged)
      for (auto const& fieldsPointDataPair : centerlinesFieldsPointData )
	{
	  auto const& vecNewPointData = fieldsPointDataPair.second;
	  for (int k=0;k<M_mapVertexVtkIdToGmshId.size();++k)
	    {
	      int idContainer = _mapVertexGmshIdToVtkId.find( M_mapVertexVtkIdToGmshId[k] )->second;
	      newCenterlinesFieldsPointData[fieldsPointDataPair.first][idContainer] = vecNewPointData[k];
	      if ( newCenterlinesFieldsPointData[fieldsPointDataPair.first][idContainer].empty() ) Msg::Error("vecNewPointData[k] is empty2");
	    }
	}
      centerlinesFieldsPointData.clear();
      centerlinesFieldsPointData = newCenterlinesFieldsPointData;

      // check component data is not empty
      for (auto const& fieldsPointDataPair : centerlinesFieldsPointData )
	{
	  for (int kk=0;kk<fieldsPointDataPair.second.size();++kk)
	    if ( fieldsPointDataPair.second[kk].empty() )
	      Msg::Error("fieldPointDataPair.second.empty() is empty at %d",kk);
	}
    }
  // update new mapping between GmshId and VtkId
  this->updateRelationMapVertex();

  std::map<int,int>  _previousMapVertexGmshIdToVtkId,_previousMapVertexVtkIdToGmshId;
  this->updateRelationMapVertex(_previousMapVertexGmshIdToVtkId,_previousMapVertexVtkIdToGmshId);

#if 1 // NEW TRUC
  lines.clear();
  colorp.clear();
  colorl.clear();
  int maxN = 0.0;


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
      //std::cout << "branchIdToReverse : " << k << "\n";
      for ( MLine* myline : edges[k].lines )
	{
	  myline->reverse();
	}
    }


  //registerLinesBuildFromPointPair.clear();
  std::set<std::pair<int,int> > registerLinesBuildFromPointPair;
  for(int i = 0; i < edges.size(); i++)
    {
      int geTag = tagMap[i];
      //std::cout << "geTag " << geTag << "\n";
      std::vector<MLine*> mylines = edges[i].lines;
      for (unsigned int j= 0; j < mylines.size(); j++){
	MLine *l = mylines[j];
	MVertex* v0 = l->getVertex(0); MVertex* v1 = l->getVertex(1);
	int v0Id = v0->getIndex(), v1Id = v1->getIndex();

	if ( registerLinesBuildFromPointPair.find( std::make_pair(v0Id,v1Id) ) == registerLinesBuildFromPointPair.end() &&
	     registerLinesBuildFromPointPair.find( std::make_pair(v1Id,v0Id) ) == registerLinesBuildFromPointPair.end() ) {
	  lines.push_back(l);
	  colorl.insert(std::make_pair(l, geTag));
	  maxN = std::max(maxN, geTag);
	  registerLinesBuildFromPointPair.insert( std::make_pair(v1Id,v0Id) );
	  registerLinesBuildFromPointPair.insert( std::make_pair(v0Id,v1Id) );
	}
	std::map<MVertex*, int>::iterator it0 = colorp.find(v0);
	if (it0 == colorp.end()) colorp.insert(std::make_pair(v0, geTag));
	std::map<MVertex*, int>::iterator it1 = colorp.find(v1);
	if (it1 == colorp.end()) colorp.insert(std::make_pair(v1, geTag));

      }

    }
  this->createBranches(maxN);
#endif
  //---------------------------------------------------------------------------------------------//
  // in case of merging centerline : update centerlinesFieldsPointData
    {
      std::map<std::string,std::vector<std::vector<double> > >  newCenterlinesFieldsPointData;
      // get new mapping between GmshId and VtkId
      std::map<int,int>  _mapVertexGmshIdToVtkId,_mapVertexVtkIdToGmshId;
      _mapVertexGmshIdToVtkId.clear(); _mapVertexVtkIdToGmshId.clear();
      this->updateRelationMapVertex( _mapVertexGmshIdToVtkId,_mapVertexVtkIdToGmshId);
      for (auto const& fieldsPointDataPair : centerlinesFieldsPointData )
	{
	  newCenterlinesFieldsPointData[fieldsPointDataPair.first].clear();
	  auto const& vecNewPointData = fieldsPointDataPair.second;
	  newCenterlinesFieldsPointData[fieldsPointDataPair.first].resize(_mapVertexVtkIdToGmshId.size());
	  //std::cout << "fieldsPointDataPair.first " << fieldsPointDataPair.first << " : " << newCenterlinesFieldsPointData[fieldsPointDataPair.first].size() << "\n";
	  for (int k=0;k<_mapVertexVtkIdToGmshId.size();++k)
	    {
	      if ( _previousMapVertexGmshIdToVtkId.find( _mapVertexVtkIdToGmshId[k] ) == _previousMapVertexGmshIdToVtkId.end() )
		Msg::Error("Error in mapping \n");
	      int idContainer = _previousMapVertexGmshIdToVtkId[ _mapVertexVtkIdToGmshId[k] ];

	      //std::cout << "k="<< k << " _mapVertexVtkIdToGmshId[k] "<< _mapVertexVtkIdToGmshId[k]
	      //<< " idContainer " << idContainer << "\n";
	      int compSize = centerlinesFieldsPointData[fieldsPointDataPair.first][idContainer].size();

	      std::vector<double> val(compSize);
	      newCenterlinesFieldsPointData[fieldsPointDataPair.first][k].resize(compSize);
	      //centerlinesFieldsPointData[fieldsPointDataPair.first][k] = vecNewPointData[idContainer];

	      for (int comp=0;comp<compSize;++comp )
		newCenterlinesFieldsPointData[fieldsPointDataPair.first][k][comp] =  vecNewPointData[idContainer][comp];
	    }
	}

	  centerlinesFieldsPointData.clear();
	  centerlinesFieldsPointData = newCenterlinesFieldsPointData;
    }

  // update new mapping between GmshId and VtkId
  this->updateRelationMapVertex();
  //this->updateCenterlinesFieldsFromFile(fileName);

  Msg::Info("AngioTkCenterline: finish importFile");

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
  this->updateCenterlinesForUse();
  this->updateRelationMapVertex();

  Msg::Info("AngioTkCenterline: createFromGeoCenterlinesFile : finish");
}

void AngioTkCenterline::updateCenterlinesForUse()
{
  // create coloring
  lines.clear();
  colorp.clear();
  colorl.clear();
  int maxN = 0.0;

  /*int nbModEdgeAddedFromMerging=0;
    if ( newCenterlines )
    nbModEdgeAddedFromMerging = std::distance( newCenterlines->mod->firstEdge(), newCenterlines->mod->lastEdge() );*/

  // import first build previous centerlines (with maybe some mod) and after add new centerlines
  //std::set<std::pair<MVertex*,MVertex*> > registerLinesBuildFromPointPair;
  std::set<std::pair<int,int> > registerLinesBuildFromPointPair;
  for (unsigned int i = 0; i < modEdges.size(); i++)
    {
      GEdge *ge = modEdges[i];
      int geTag = ge->tag();
      //std::cout<< "init my branch i "<<i<<" has " << ge->lines.size() << " lines \n";
      for(unsigned int j = 0; j < ge->lines.size(); j++)
	{
	  MLine *l = ge->lines[j];
	  MVertex *v0 = l->getVertex(0); MVertex *v1 = l->getVertex(1);
	  int v0Id = v0->getIndex(), v1Id = v1->getIndex();

	  // after remove duplacte vertices, it is possible to have a wrong line with 2 same pts
	  if ( v0 == v1 )
	    {
	      //std::cout << "ignore a wrong line\n";
	      continue;
	    }
	  if ( M_registerLinesToRemoveFromPointIdPairInModelEdge.find( std::make_pair(v0Id,v1Id ) ) != M_registerLinesToRemoveFromPointIdPairInModelEdge.end() ||
	       M_registerLinesToRemoveFromPointIdPairInModelEdge.find( std::make_pair(v1Id,v0Id ) ) != M_registerLinesToRemoveFromPointIdPairInModelEdge.end() )
	    {
	      //std::cout << "ignore a line\n";
	      continue;
	    }

	  //if (it0 == colorp.end() || it1 == colorp.end()){
	  if ( registerLinesBuildFromPointPair.find( std::make_pair(v0Id,v1Id) ) == registerLinesBuildFromPointPair.end() &&
	       registerLinesBuildFromPointPair.find( std::make_pair(v1Id,v0Id) ) == registerLinesBuildFromPointPair.end() ) {
	    lines.push_back(l);
	    colorl.insert(std::make_pair(l, geTag));
	    maxN = std::max(maxN, ge->tag());
	    registerLinesBuildFromPointPair.insert( std::make_pair(v1Id,v0Id) );
	    registerLinesBuildFromPointPair.insert( std::make_pair(v0Id,v1Id) );
	  }
	  std::map<MVertex*, int>::iterator it0 = colorp.find(v0);
	  if (it0 == colorp.end()) colorp.insert(std::make_pair(v0, geTag));
	  std::map<MVertex*, int>::iterator it1 = colorp.find(v1);
	  if (it1 == colorp.end()) colorp.insert(std::make_pair(v1, geTag));
	}
    }
  //std::cout << "lines.size() " << lines.size() << " colorl.size() " << colorl.size() << "\n" << "colorp.size() " << colorp.size() << "\n";

  //---------------------------------------------------------------------------------------------//
  // create centerlines branch
  this->createBranches(maxN);


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
    MLine* l = itl->first;
    allV.insert(l->getVertex(0));
    allV.insert(l->getVertex(1));
    color_edges[color-1].push_back(l);
    itl++;
  }

  Msg::Info("AngioTkCenterline: color_edges built");

  //detect junctions
  junctions.clear();
  std::multiset<MVertex*>::iterator it = allV.begin();
  for ( ; it != allV.end(); ++it){
    if (allV.count(*it) != 2) {
      junctions.insert(*it);
    }
  }
  //std::cout << "junctions.size() " << junctions.size() << "\n";
  Msg::Info("AngioTkCenterline: junction built");

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
      if (vB == vE) {
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
#if 0
  Msg::Info("AngioTkCenterline: in/outlets?? =%d branches =%d ",
            (int)color_edges.size()+1, (int)edges.size());
#endif
  //create children
  for(unsigned int i = 0; i < edges.size(); ++i) {
    MVertex *vE = edges[i].vE;
    std::vector<Branch> myChildren;
    for (std::vector<Branch>::iterator it = edges.begin(); it != edges.end(); ++it){
      Branch myBranch = *it;
      if (myBranch.vB == vE) myChildren.push_back(myBranch);
    }
    edges[i].children = myChildren;
  }

  // compute junction and extremity points
  //std::map<MVertex*,int> numberOfPoints;
  std::map<MVertex*,std::set<std::pair<int,int> > > numberOfPoints;
  for(unsigned int i = 0; i < edges.size(); i++){
    std::vector<MLine*> mylines = edges[i].lines;
      for (unsigned int j= 0; j < mylines.size(); j++){
	if ( mylines[j]->getLength() > 5 )
	  std::cout << "mylines[j]->getLength() "<< mylines[j]->getLength() <<"\n";
	  MVertex *v1 = mylines[j]->getVertex(0);
	  MVertex *v2 = mylines[j]->getVertex(1);
	  //if ( numberOfPoints.find(v1) == numberOfPoints.end() )
	  numberOfPoints[v1].insert( std::make_pair(i,j) );
	    //else
	    //numberOfPoints[v1] = numberOfPoints[v1]+ 1;
	    //if ( numberOfPoints.find(v2) == numberOfPoints.end() )
	  numberOfPoints[v2].insert( std::make_pair(i,j) );
	  //numberOfPoints[v2] = 1;
	  //else
	  //numberOfPoints[v2] = numberOfPoints[v2]+1;
      }
  }
  M_junctionsVertex.clear();
  M_extremityVertex.clear();
  for ( auto const& nPtsPair : numberOfPoints )
    {
      if ( nPtsPair.second.size() == 1 )
	M_extremityVertex[nPtsPair.first] = *(nPtsPair.second.begin());// std::make_pair(0,0);
      else if ( nPtsPair.second.size() > 2 )
	M_junctionsVertex.insert( nPtsPair.first );
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
    //std::cout << "PPPPP faces().size " << pe->faces().size()<< "\n";

    //discreteEdge hola( (const GEdge) *pe);
    //discreteEdge * hola = new discreteEdge( *pe );
    //discreteEdge * hola = new discreteEdge(current, NE+900, pe->getBeginVertex(), pe->getEndVertex());
    //discreteEdge * hola = new discreteEdge(current, NE+900,  pe->getEndVertex(),pe->getBeginVertex());
    //discreteEdge * hola = new discreteEdge(current, NE+900,0,0 );
    //delete hola;

    e_compound.push_back(pe);
    int num_gec = NE+i+1;
    //int num_gec = 2*NE+2*i+2;
    Msg::Info("Create Compound Line (%d)     = %d discrete edge",
              num_gec, pe->tag());
    GEdge *gec =  current->addCompoundEdge(e_compound,num_gec);
    Msg::Info("Create Compound Line (%d) tag = %d discrete edge",
              num_gec, gec->tag());

    if (CTX::instance()->mesh.algo2d != ALGO_2D_BAMG){
      gec->meshAttributes.method = MESH_TRANSFINITE;
      gec->meshAttributes.nbPointsTransfinite = nbPoints+1;
      gec->meshAttributes.typeTransfinite = 0;
      gec->meshAttributes.coeffTransfinite = 1.0;
    }
  }
#if 0
  std::map<int,bool> mapReverseMesh;
  // search inltet to start
  bool find=false;
  for (int i=0; i < NE && !find; i++)
    {
      GEdge *pe = current->getEdgeByTag(i+1);//current edge
      if ( pe->faces().size() == 1 )
	{
	  mapReverseMesh[ pe->faces().front()->tag() ] = false;
	  std::cout << " mapReverseMesh init with " << pe->faces().front()->tag()  << " val " << mapReverseMesh[ pe->faces().front()->tag() ] << "\n";
	  find=true;
	}
      //std::cout << "PPPPP regions().size " << pe->regions().size()<< "\n";
      //std::cout << "PPPPP faces().size " << pe->faces().size()<< "\n";
    }
  if ( !find ) exit(0);
  //int pe->faces().front()->tag();
  //GFace *pf =  current->getFaceByTag( mapReverseMesh.begin()->first );

  while( mapReverseMesh.size() < NF-1 )
    {
      for (int i=0; i < NF; i++)
	{
	  GFace *pf =  current->getFaceByTag( i+1 );
	  if ( mapReverseMesh.find( i+1 ) != mapReverseMesh.end() )
	    {
	      std::list<GEdge*> myEdges = pf->edges();
	      std::list<GEdge*>::iterator itE = myEdges.begin();
	      std::list<GEdge*>::iterator enE = myEdges.end();
	      //for (  pf->edges().size();++kk )
	      for ( ; itE!=enE;++itE ) // foreach edges
		{
		  std::list<GFace*> myFaces = (*itE)->faces();
		  std::list<GFace*>::iterator itF = myFaces.begin();
		  std::list<GFace*>::iterator enF = myFaces.end();
		  for ( ; itF!=enF;++itF ) // foreach edges
		    {
		      if ( (*itF)->tag() == i+1 ) continue;
		      std::cout << "(*itF)->tag()" << (*itF)->tag() << " vs i+1 " << i+1 << "\n";
		      if ( mapReverseMesh.find( (*itF)->tag() ) == mapReverseMesh.end() )
			{
			  mapReverseMesh[ (*itF)->tag() ] = !(mapReverseMesh.find( i+1 )->second);
			  std::cout << " mapReverseMesh add  with " << (*itF)->tag() << " val " << mapReverseMesh[ (*itF)->tag() ] << "\n";
			}
		      else if ( (mapReverseMesh.find( (*itF)->tag() )->second == mapReverseMesh.find( i+1 )->second) )
			{
			  //std::cout << mapReverseMesh[ (*itF)->tag() ]<< " vs " << mapReverseMesh.find( i+1 )->second << "n";
			  exit(0);
			}
		    }

		}
	    }
	}
      std::cout<< " mapReverseMesh.size() " << mapReverseMesh.size() << "\n";
    } // while
#endif
  // Parametrize Compound surfaces
  std::list<GEdge*> U0;
  for (int i=0; i < NF; i++){
    std::vector<GFace*> f_compound;
    GFace *pf =  current->getFaceByTag(i+1);//current face
    f_compound.push_back(pf);
    int num_gfc = NF+i+1;
    Msg::Info("Create Compound Surface (%d)     = %d discrete face",
              num_gfc, pf->tag());

    //1=conf_spectral 4=convex_circle, 7=conf_fe
    GFace *gfc = current->addCompoundFace(f_compound, 7, 0, num_gfc);
    // 4 marche pas mal
    //GFace *gfc = current->addCompoundFace(f_compound, 4, 0, num_gfc);
    //GFace *gfc = current->addCompoundFace(f_compound, 5, 0, num_gfc);

    Msg::Info("Create Compound Surface (%d) tag = %d discrete face",
              num_gfc, gfc->tag());

    //if ( gfc->tag() == 14 )
    //if ( pf->tag() == 1 )
    //if ( mapReverseMesh.find( gfc->tag() ) != mapReverseMesh.end() && mapReverseMesh.find( i+1 )->second )

#if 0
    if ( mapReverseMesh.find( i+1 ) != mapReverseMesh.end() && mapReverseMesh.find( i+1 )->second )
      gfc->revertMeshPATCHAngiotk=true;
#endif
    //gfc->meshAttributes.reverseMesh=true;

    gfc->meshAttributes.recombine = recombine;
    gfc->addPhysicalEntity(1);
    current->setPhysicalName("wall", 2, 1);//tag 1

  }


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
    //std::cout << "before recurConnectByMEdge\n";
    recurConnectByMEdge(me,e2e, group, touched, theCut);
    //std::cout << "touched.size() " << touched.size() << "\n";
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

  int currentTagPhysicalMark = 100;
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
	std::string marker = mapBoundEdgeIdToPhysicalMarkerLumen.find( boundEdges[i]->tag() )->second;
	int thetag = currentTagPhysicalMark;
	if ( mapMarkerToTag.find(marker) != mapMarkerToTag.end() )
	  thetag = mapMarkerToTag.find(marker)->second;
	else
	  {
	    mapMarkerToTag[marker] = thetag;
	    ++currentTagPhysicalMark;
	  }

	newFace->addPhysicalEntity( thetag);
	current->setPhysicalName(marker, 2, thetag);//tag 2
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
  reg->addPhysicalEntity(reg->tag());
  current->setPhysicalName("lumenVolume", 3, reg->tag());

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
    eFace->addPhysicalEntity(5);
    current->setPhysicalName("outerWall", 2, 5);//dim 2 tag 5
    GRegion *eRegion = (GRegion*) extrudedE[1];
    eRegion->addPhysicalEntity(6);
    current->setPhysicalName("wallVolume", 3, 6);//dim 3 tag 6

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


    int currentTagPhysicalMark = 500;
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
	      std::string marker = mapBoundEdgeIdToPhysicalMarkerArerialWall.find( myEdge->tag() )->second;
	      int thetag=currentTagPhysicalMark;
	      if ( mapMarkerToTag.find(marker) != mapMarkerToTag.end() )
		thetag = mapMarkerToTag.find(marker)->second;
	      else
		{
		  mapMarkerToTag[marker] = thetag;
		  ++currentTagPhysicalMark;
		}
	      elFace->addPhysicalEntity( thetag);
	      current->setPhysicalName(marker, 2, thetag);//tag 2
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


  if (is_cut) cutMesh();
  else{
    GFace *gf = current->getFaceByTag(1);
    gf->addPhysicalEntity(4);
    current->setPhysicalName("wall", 2, 4);//tag 1
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
  std::vector<std::tuple<SVector3,SVector3,double> > cutDesc; 

  for(unsigned int i = 0; i < edges.size(); ++i)
    {
      MVertex *vB = edges[i].vB;
      MVertex *vE = edges[i].vE;
      //std::cout << "edge " << i << " vB : " << vB->x() << " " << vB->y() << " " << vB->z() << "\n";
      //std::cout << "edge " << i << " vE : " << vE->x() << " " << vE->y() << " " << vE->z() << "\n";
      bool vBeginConnected = false, vEndConnected = false;

	//for ( Branch const& branchConnected : edges[i].children )
      for(unsigned int i2 = 0; i2 < edges.size(); ++i2)
	{
	  if ( i2 == i ) continue;

	  Branch const& branchConnected =  edges[i2];
	  //std::cout << "branchConnected.vB : " << branchConnected.vB->x() << " " << branchConnected.vB->y() << " " << branchConnected.vB->z() << "\n";
	  //std::cout << "branchConnected.vE : " << branchConnected.vE->x() << " " << branchConnected.vE->y() << " " << branchConnected.vE->z() << "\n";

	  //if ( ( vB->distance( branchConnected.vE ) < 1e-9 ) || (vB->distance( branchConnected.vB ) < 1e-9 ) )
	  if ( vB == branchConnected.vE || vB == branchConnected.vB  )
	       vBeginConnected = true;
	    //if ( ( vE->distance( branchConnected.vB ) < 1e-9 ) || (vE->distance( branchConnected.vE ) < 1e-9 ) )
	  if ( vE == branchConnected.vB || vE == branchConnected.vE )
	    vEndConnected = true;
	}

      if ( !vBeginConnected )
	{
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

	      for ( int k=0 ; k<lines.size() && curLength <= scalingClip*radiusLineB ; ++k )
		{
		  if ( extremityIsPoint0 )
		    {
		      curLine = lines[k];
		      curLength += curLine->getLength();
		      v1 = curLine->getVertex(1);
		      v2 = curLine->getVertex(0);
		    }
		  //v1 = lines.front()->getVertex(1);
		  //v2 = lines.front()->getVertex(0);
		  else
		    {
		      curLine = lines[k];
		      curLength += curLine->getLength();
		      v1 = curLine->getVertex(0);
		      v2 = curLine->getVertex(1);
		    }
		  //v1 = lines.front()->getVertex(0);
		  //v2 = lines.front()->getVertex(1);
		}
	      MVertex *thept = v1;//v2;
	      std::map<MLine*,double>::iterator itr = radiusl.find(curLine);
	      double radius = 1.1*itr->second;

	      cutDesc.push_back( std::make_tuple( SVector3(thept->x(),thept->y(),thept->z()),
						  SVector3(v2->x()-v1->x(),v2->y()-v1->y(),v2->z()-v1->z()),
						  radius ) );
	    }
	}
      if ( !vEndConnected )
	{
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
	      double radiusLineB = 0;
	      int nAvergage = std::min((int)lines.size(),(int)50);
		for ( int ls=0;ls<nAvergage ;++ls)
		{
		  radiusLineB += radiusl.find(lines[lines.size()-1-ls])->second;
		}
		radiusLineB /= nAvergage;
		radiusLineB = std::max( radiusLineB, radiusl.find(lines.back())->second );
#endif


	      MLine* curLine = lines.back();
	      double curLength = 0.;

	      bool extremityIsPoint0 = ( lines.back()->getVertex(0) == vB ) || (lines.back()->getVertex(0) == vE);
	      bool extremityIsPoint1 = ( lines.back()->getVertex(1) == vB ) || (lines.back()->getVertex(1) == vE);
	      if ( !extremityIsPoint0 && !extremityIsPoint1 )
		Msg::Error("error !!");

	      for ( int k=0 ; k<lines.size() && curLength <= scalingClip*radiusLineB ; ++k )
		{
		  if ( extremityIsPoint1 ) 
		    {
		      curLine = lines[lines.size()-1-k];
		      curLength += curLine->getLength();
		      v1 = curLine->getVertex(0);
		      v2 = curLine->getVertex(1);
		    }
		  //v1 = lines.back()->getVertex(0);
		  //v2 = lines.back()->getVertex(1);
		  else
		    {
		      curLine = lines[lines.size()-1-k];
		      curLength += curLine->getLength();
		      v1 = curLine->getVertex(1);
		      v2 = curLine->getVertex(0);
		    }
		  //v1 = lines.back()->getVertex(1);
		  //v2 = lines.back()->getVertex(0);
		}
	      //std::cout << " v1 : " << v1->x() << " " << v1->y() << " " << v1->z() << "\n";
	      //std::cout << " v2 : " << v2->x() << " " << v2->y() << " " << v2->z() << "\n";

	      MVertex *thept = v1;//v2;
	      //std::map<MLine*,double>::iterator itr = radiusl.find(lines[lines.size()-1]);
	      //std::map<MLine*,double>::iterator itr = radiusl.find(lines.back());
	      //double radius = itr->second;//0.8;
	      std::map<MLine*,double>::iterator itr = radiusl.find(curLine);
	      double radius = 1.1*itr->second;

	      cutDesc.push_back( std::make_tuple( SVector3(thept->x(),thept->y(),thept->z()),
						  SVector3(v2->x()-v1->x(),v2->y()-v1->y(),v2->z()-v1->z()),
						  radius ) );
	    }

	}

    }

  Msg::Info("AngioTkCenterline: prepare data for cutting done");


  std::ofstream fileWrited( "sphereremovebranch.data", std::ios::out | std::ios::trunc);
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

    }
  fileWrited.close();
  Msg::Info("AngioTkCenterline: all cuts done");

  //create discreteFaces
  //createFaces();
  this->createFacesFromClip();
  current->createTopologyFromFaces(discFaces);
  current->exportDiscreteGEOInternals();

  //write
  Msg::Info("AngioTkCenterline: writing splitted mesh 'myCLIPPARTS.msh'");
  current->writeMSH("myCLIPPARTS.msh", 2.2, false, false);

#if 0
  //create compounds
  createSplitCompounds();
  Msg::Info("Done splitting mesh by centerlines");
#endif


}



#if 0
void AngioTkCenterline::cutMesh()
{
  Msg::Info("AngioTkCenterline: action (cutMesh) splits surface mesh (%d tris) using %s ",
            triangles.size(), fileName.c_str());

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
#if 1
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

  //create discreteFaces
  createFaces();
  current->createTopologyFromFaces(discFaces);
  current->exportDiscreteGEOInternals();

  //write
  Msg::Info("AngioTkCenterline: writing splitted mesh 'myPARTS.msh'");
  current->writeMSH("myPARTS.msh", 2.2, false, false);
  //exit(0);
  //create compounds
  createSplitCompounds();

  Msg::Info("Done splitting mesh by centerlines");
}
#else

/**
 * NEW
 */
void AngioTkCenterline::cutMesh()
{
  Msg::Info("AngioTkCenterline: action (cutMesh) splits surface mesh (%d tris) using %s ",
            triangles.size(), fileName.c_str());

  // i->j->(pt,radius)
  //std::vector< std::map<int, std::pair<SVector3,double> > > cutDiskToPerform(edges.size());
  std::vector< std::map<int, std::tuple<SVector3,double,double> > > cutDiskToPerform(edges.size());//pt,radius,lengthFromBranchBegin
  // first pass
  for(unsigned int i = 0; i < edges.size(); i++){
    std::vector<MLine*> lines = edges[i].lines;
    double L = edges[i].length;
    //double D = 2.*edges[i].minRad;  //(edges[i].minRad+edges[i].maxRad);
    double D = (edges[i].minRad+edges[i].maxRad);
    double AR = L/D;

    bool vBisJunc = /*junctions*/M_junctionsVertex.find(edges[i].vB) != M_junctionsVertex.end();
    bool vEisJunc = /*junctions*/M_junctionsVertex.find(edges[i].vE) != M_junctionsVertex.end();

    bool vBisInOut = M_extremityVertex.find(edges[i].vB) != M_extremityVertex.end();
    bool vEisInOut = M_extremityVertex.find(edges[i].vE) != M_extremityVertex.end();

    double radiusB = radiusl.find(lines.front())->second;
    double radiusE = radiusl.find(lines.back())->second;


    // printf("*** AngioTkCenterline branch %d (AR=%.1f) \n", edges[i].tag, AR);
#if 0 //VINCENT
    int nbSplit = (int)ceil(AR/2 + 1.1); //AR/2 + 0.9
#else
    //int nbSplit = 4*(int)ceil(AR + 1.1); //AR/2 + 0.9
    int nbSplit = (int)ceil(AR + 1.1); //AR/2 + 0.9
#endif
    //int nbSplit = (int)ceil(AR/4 + 1.1); //AR/2 + 0.9
    //nbSplit = std::min( nbSplit, 8 );
    if( nbSplit > 1 ){
      double li  = L/nbSplit;
      printf("->> cut branch %i in %d parts (L=%f,D=%f,AR=%f,li=%f)\n", i, nbSplit, L,D,AR,li);

      double lc = 0.0, lcTotal = 0.0;
      for (unsigned int j= 0; j < lines.size(); j++){
	lc += lines[j]->getLength();
	MVertex *v1 = lines[j]->getVertex(0);
	MVertex *v2 = lines[j]->getVertex(1);

	if ( lc > li && nbSplit > 1) {
	  SVector3 pt(v1->x(), v1->y(), v1->z());
	  SVector3 dir(v2->x()-v1->x(),v2->y()-v1->y(),v2->z()-v1->z());
	  std::map<MLine*,double>::iterator itr = radiusl.find(lines[j]);
	  double radius= itr->second;
	  //cutByDisk(pt, dir, itr->second);

	  bool applyCut=true;
	  //if ( i!=2 ) applyCut=false;
#if 1
	  if ( vBisInOut && ( lcTotal < 2*radiusB /*radiusl.find(lines.front())->second*/ ) )
	    applyCut=false;
	  if ( vEisInOut && ( lcTotal > (L-2*radiusE/*radiusl.find(lines.front())->second*/ ) ) )
	    applyCut=false;
#endif

	  for(unsigned int ii = 0; ii < edges.size(); ii++){
	    auto/*std::map<int, std::pair<SVector3,double> >::iterator*/ itj = cutDiskToPerform[ii].begin();
	    auto/*std::map<int, std::pair<SVector3,double> >::iterator*/ enj = cutDiskToPerform[ii].end();
	    for (; itj!=enj && applyCut;++itj)
	      {
		int jTest= itj->first;
		SVector3 ptTest = std::get<0>(itj->second);//itj->second.first;
		double radiusTest= std::get<1>(itj->second);//itj->second.second;
		double lengthFromBranchBeginTest= std::get<2>(itj->second);
		double distBetweenCenter = std::sqrt( std::pow( pt.x()-ptTest.x(),2)+std::pow( pt.y()-ptTest.y(),2)+std::pow( pt.z()-ptTest.z(),2) );

		// ignore cut done quite far of branch extremities 
		if ( ii != i &&
		    lengthFromBranchBeginTest > 3*edges[ii].maxRad &&
		    lengthFromBranchBeginTest < (edges[ii].length -3*edges[ii].maxRad) )
		  {
		    continue;
		  }


		MVertex *v1Test = lines[j]->getVertex(0);
		MVertex *v2Test = lines[j]->getVertex(1);
		SVector3 dirTest(v2Test->x()-v1Test->x(),v2Test->y()-v1Test->y(),v2Test->z()-v1Test->z());


		double theScaling = 1./2.;
		theScaling = 1;//1./4;
		if (norm(dir-dirTest) < 1e-2 )
		  theScaling = 1./2;
		if ( lc > 4*li )
		  theScaling = 1./2;

#if 1
		if ( ii == i )
		  {
		//if ( lcTotal > L*(1./3.) && lcTotal < L*(2./3.) )
		double theScalingCoarse =  2;
		if ( !vBisJunc && vEisJunc && lcTotal < (L-4*radiusE)/*edges[i].maxRad*/ )
		  theScaling = theScalingCoarse;
		if ( vBisJunc && !vEisJunc && lcTotal > 4*radiusB/*edges[i].maxRad*/ )
		  theScaling = theScalingCoarse;
		if ( vBisJunc && vEisJunc && lcTotal > 4*radiusB/*edges[i].maxRad*/ && lcTotal < (L-4*radiusB/*edges[i].maxRad*/) )
		  theScaling = theScalingCoarse;
		if ( !vBisJunc && !vEisJunc )
		  theScaling = theScalingCoarse;
		  }
#endif
		if ( distBetweenCenter < (radius+radiusTest)*theScaling ) //if ( distBetweenCenter < (radius+radiusTest) )
		  applyCut=false;
	      }
	  }
	  if ( applyCut )
	    {
	      //cutDiskToPerform[i][j] = std::make_pair(pt,itr->second );
	      cutDiskToPerform[i][j] = std::make_tuple(pt,itr->second,lcTotal );
	      lc = 0.0;
	    }
	  //nbSplit--;
	  //lc = 0.0;
	} // if ( lc > li && nbSplit > 1)
	lcTotal += lines[j]->getLength();
      } // for (unsigned int j= 0; j < lines.size(); j++)
    }
  } // end first for

  for(unsigned int i = 0; i < edges.size(); i++){
    std::vector<MLine*> lines = edges[i].lines;
    auto/*std::map<int, std::pair<SVector3,double> >::iterator*/ itj = cutDiskToPerform[i].begin();
    auto/*std::map<int, std::pair<SVector3,double> >::iterator*/ enj = cutDiskToPerform[i].end();
    //std::cout << "->> cut branch "<< i << "(L="<< L << " ,D="<< D << " ,AR=" << AR << ")"
    //<< " in " << nbSplit << " parts\n";
    std::cout << "->> cut my branch "<< i << " with "<< cutDiskToPerform[i].size() << "\n";
    for (; itj!=enj;++itj)
      {
	int j= itj->first;
	SVector3 pt = std::get<0>(itj->second) ;//itj->second.first;
	double radius= std::get<1>(itj->second);//itj->second.second;

	//if ( i==2 ) std::cout << "pt " << pt.x() << "," << pt.y() << "," << pt.z() << " radius " << radius << "\n";

	MVertex *v1 = lines[j]->getVertex(0);
	MVertex *v2 = lines[j]->getVertex(1);
	SVector3 dir(v2->x()-v1->x(),v2->y()-v1->y(),v2->z()-v1->z());
	cutByDisk(pt, dir, radius);
      }
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
  Msg::Info("AngioTkCenterline: writing splitted mesh 'myPARTS.msh'");
  current->writeMSH("myPARTS.msh", 2.2, false, false);
  //exit(0);
  //create compounds
  createSplitCompounds();

  Msg::Info("Done splitting mesh by centerlines");
}


#endif
bool AngioTkCenterline::cutByDisk(SVector3 &PT, SVector3 &NORM, double &maxRad, int tag)
{
  double a = NORM.x();
  double b = NORM.y();
  double c = NORM.z();
  double d = -a * PT.x() - b * PT.y() - c * PT.z();

  int maxStep = 20;
  //const double EPS = 0.007;
  const double EPS = 1e-5;//0.007;

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
      bool inDisk = ((norm(P1-PT) < rad ) || (norm(P2-PT) < rad)) ? true : false;
      //bool inDisk2 = ((norm(P1-PT) < rad ) && (norm(P2-PT) < rad)) ? true : false;
      double rdist = -V1/(V2-V1);

      if ( std::abs(V1)<=EPS/*1e-4*/ && std::abs(V2)<=EPS/*1e-4*/ && (norm(P1-PT) < rad ) && (norm(P2-PT) < rad) ){
	inters = false;
	//std::cout << "les 2 pts sont dans le plan de coupe\n";
      }
      else if ( std::abs(V1)<=EPS && (norm(P1-PT) < rad ) )
	cutVertices.push_back(me.getVertex(0));
      else if ( std::abs(V2)<=EPS && (norm(P2-PT) < rad ) )
	cutVertices.push_back(me.getVertex(1));
      else if (inters && rdist > EPS && rdist < (1.-EPS) ){
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
  unsigned int nBranch = edges.size();
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
  writeVTKPolyData(mod,edges,fileName,centerlinesFieldsPointData);
}

#include <gmshHeadersMissing/MVertexPositionSet.h>
void AngioTkCenterline::updateCenterlinesFromFile(std::string fileName)
{
  this->importFile( fileName );
  this->updateCenterlinesFieldsFromFile( fileName );
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
  centerlinesFieldsPointData.clear();
  for ( auto const& thefield : fieldsPointDataInput )
    {
      centerlinesFieldsPointData[thefield.first].resize(numVertices);
      for ( auto & thefieldClean : centerlinesFieldsPointData[thefield.first] )
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
	    centerlinesFieldsPointData[thefield.first][this->M_mapVertexGmshIdToVtkId[vertexIdClean]][comp] = thefield.second[myvertexInitial->getIndex()-1][comp];
	  else
	    centerlinesFieldsPointData[thefield.first][this->M_mapVertexGmshIdToVtkId[vertexIdClean]][comp] = valTolMin;
#else
	  centerlinesFieldsPointData[thefield.first][this->M_mapVertexGmshIdToVtkId[vertexIdClean]][comp] = thefield.second[myvertexInitial->getIndex()-1][comp];
#endif

	}
    }
}

void AngioTkCenterline::addFieldBranchIds( std::string const& fieldName )
{
  //if ( centerlinesFieldsPointData.find("BranchIds") != centerlinesFieldsPointData.end() )
  //  return;
  Msg::Info("AngioTkCenterline: addFieldBranchIds");
  int numVertices = numberOfVertices(edges);
  centerlinesFieldsPointData[fieldName].resize( numVertices );

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
	centerlinesFieldsPointData[fieldName][v0VtkId] = { (double)i };
	ptDone[v0VtkId] = true;
      }
      if ( !ptDone[v1VtkId] )
      {
	centerlinesFieldsPointData[fieldName][v1VtkId] = { (double)i };
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

  centerlinesFieldsPointData[fieldName].clear();
  centerlinesFieldsPointData[fieldName].resize( numVertices );
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
      centerlinesFieldsPointData[fieldName][vVtkId] = { (double)minRad };
    }
}

double
AngioTkCenterline::minRadiusAtVertex( MVertex* myvertex )
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
AngioTkCenterline::applyFieldThresholdMin(std::string const& fieldName, double valMinThreshold)
{
  int numVertices = numberOfVertices(edges);
  if ( centerlinesFieldsPointData.find(fieldName) != centerlinesFieldsPointData.end() &&
       centerlinesFieldsPointData.find(fieldName)->second.size() == numVertices )
    {
      for (int k=0;k<numVertices;++k)
	{
	  if ( centerlinesFieldsPointData[fieldName][k][0] < valMinThreshold )
	    centerlinesFieldsPointData[fieldName][k][0] = valMinThreshold;
	}
    }
}



void AngioTkCenterline::removeBranchIds( std::set<int> const& _removeBranchIds )
{
  if ( _removeBranchIds.empty() )
    return;

  lines.clear();
  colorp.clear();
  colorl.clear();
  int maxN = 0.0;

  unsigned int nBranch = edges.size();
  for(unsigned int i = 0; i < nBranch; ++i)
  {
    if ( _removeBranchIds.find( i+1 ) != _removeBranchIds.end() )
      continue;
    std::vector<MLine*> mylines = edges[i].lines;
    for(unsigned int k = 0; k < mylines.size(); ++k)
    {
      MLine *l = mylines[k];
      MVertex *v0 = l->getVertex(0);
      MVertex *v1 = l->getVertex(1);
      std::map<MVertex*, int>::iterator it0 = colorp.find(v0);
      std::map<MVertex*, int>::iterator it1 = colorp.find(v1);
      if (it0 == colorp.end() || it1 == colorp.end())
      {
	lines.push_back(l);
	colorl.insert(std::make_pair(l, i+1/*ge->tag()*/));
	maxN = std::max(maxN, (int)i+1/*ge->tag()*/);
      }
      if (it0 == colorp.end()) colorp.insert(std::make_pair(v0, i/*ge->tag()*/));
      if (it1 == colorp.end()) colorp.insert(std::make_pair(v1, i/*ge->tag()*/));
    }
  }
  //edges.clear();
  createBranches(maxN);

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


#endif
