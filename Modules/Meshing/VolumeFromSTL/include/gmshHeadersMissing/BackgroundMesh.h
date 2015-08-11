// Gmsh - Copyright (C) 1997-2014 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@geuz.org>.

#ifndef _BACKGROUND_MESH_H_
#define _BACKGROUND_MESH_H_

#include <feel/feelconfig.h>
#include "STensor3.h"
#include <vector>
#include <list>
#include "simpleFunction.h"

#if defined(FEELPP_HAVE_ANN_H)
#include <ANN/ANN.h>
class ANNkd_tree;
#endif

class MElementOctree;
class GFace;
class GVertex;
class GEdge;
class MElement;
class MVertex;
class GEntity;

struct crossField2d 
{
  double _angle;
  static void normalizeAngle (double &angle) {
    if (angle < 0) 
      while ( angle <  0 ) angle += (M_PI * .5);
    else if (angle >= M_PI * .5) 
      while ( angle >= M_PI * .5 ) angle -= (M_PI * .5);
  }
  crossField2d (MVertex*, GEdge*);
  crossField2d (double a) : _angle(a){}
  crossField2d & operator += ( const crossField2d & );
};


class backgroundMesh : public simpleFunction<double>
{
  static double sizeFactor;
  MElementOctree *_octree;
  std::vector<MVertex*> _vertices;
  std::vector<MElement*> _triangles;
  std::map<MVertex*,double> _sizes;  
  std::map<MVertex*,MVertex*> _3Dto2D;
  std::map<MVertex*,MVertex*> _2Dto3D;
  std::map<MVertex*,double> _distance;  
  std::map<MVertex*,double> _angles;  
  static backgroundMesh * _current;
  backgroundMesh(GFace *, bool dist = false);
  ~backgroundMesh();
#if defined(FEELPP_HAVE_ANN_H)
   mutable ANNkd_tree *uv_kdtree;
   mutable ANNpointArray nodes;
   ANNidxArray index;
   ANNdistArray dist;
   mutable ANNpointArray angle_nodes;
   mutable ANNkd_tree *angle_kdtree;
   std::vector<double> _cos,_sin;
#endif
 public:
  static void set(GFace *);
  static void setCrossFieldsByDistance(GFace *);
  static void unset();
  static backgroundMesh *current () { return _current; }
  static void setSizeFactor (double s) {sizeFactor = s;}
  void propagate1dMesh(GFace *);
  void propagateCrossField(GFace *, simpleFunction<double> *);
  void propagateCrossFieldHJ(GFace *);
  void propagateCrossField(GFace *);
  void propagateCrossFieldByDistance(GFace *);
  void updateSizes(GFace *);
  double operator () (double u, double v, double w) const; // returns mesh size
  bool inDomain (double u, double v, double w) const; // returns true if in domain
  double getAngle(double u, double v, double w) const ; 
  double getSmoothness(double u, double v, double w)  ; 
  double getSmoothness(MElement*) ; 
  void print(const std::string &filename, GFace *gf, 
	     const std::map<MVertex*, double>&, int smooth = 0) ;
  void print(const std::string &filename, GFace *gf, int choice = 0) 
  {
    switch(choice) {
    case 0 : print(filename, gf, _sizes); return;
    case 2 : print(filename, gf, _sizes, 1); return;
    default : print(filename, gf, _angles); return;
    }
  }
  MElementOctree* get_octree();
  MElement *getMeshElementByCoord(double u, double v, double w, bool strict=true);
  int getNumMeshElements()const{return _triangles.size();}
  std::vector<MVertex*>::iterator begin_vertices(){return _vertices.begin();}
  std::vector<MVertex*>::iterator end_vertices(){return _vertices.end();}
  std::vector<MVertex*>::const_iterator begin_vertices()const{return _vertices.begin();}
  std::vector<MVertex*>::const_iterator end_vertices()const{return _vertices.end();}
  std::vector<MElement*>::iterator begin_triangles(){return _triangles.begin();}
  std::vector<MElement*>::iterator end_triangles(){return _triangles.end();}
  std::vector<MElement*>::const_iterator begin_triangles()const{return _triangles.begin();}
  std::vector<MElement*>::const_iterator end_triangles()const{return _triangles.end();}
};

SMetric3 buildMetricTangentToCurve (SVector3 &t, double l_t, double l_n);
SMetric3 buildMetricTangentToSurface (SVector3 &t1, SVector3 &t2, double l_t1, double l_t2, double l_n);
double BGM_MeshSize(GEntity *ge, double U, double V, double X, double Y, double Z);
SMetric3 BGM_MeshMetric(GEntity *ge, double U, double V, double X, double Y, double Z);
bool Extend1dMeshIn2dSurfaces();
bool Extend2dMeshIn3dVolumes();
SMetric3 max_edge_curvature_metric(const GVertex *gv);
SMetric3 max_edge_curvature_metric(const GEdge *ge, double u, double &l);
SMetric3 metric_based_on_surface_curvature(const GFace *gf, double u, double v, 
					   bool surface_isotropic = false,
					   double d_normal = 1.e12,
					   double d_tangent_max = 1.e12);

#endif
