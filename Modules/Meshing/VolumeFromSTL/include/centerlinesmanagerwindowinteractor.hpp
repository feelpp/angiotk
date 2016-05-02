#ifndef __CENTERLINESMANAGERWINDOWINTERACTOR_H
#define __CENTERLINESMANAGERWINDOWINTERACTOR_H 1

#include <string>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>

#include <vtkWorldPointPicker.h>
#include <vtkPointPicker.h>
#include <vtkCellPicker.h>
#include <vtkPropPicker.h>

#include <vtkSphereSource.h>
#include <vtkLineSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkFollower.h>
#include <vtkVectorText.h>

#include <vtkObjectFactory.h>


#include <vtkPolyData.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include <vtkProperty.h>

#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkTextWidget.h>
#include <vtkTextRepresentation.h>

#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>

#include <vtkCoordinate.h>

#include <cmath>
#include <tuple>

#include <vtkBoxWidget.h>
#include <vtkBoxWidget2.h>
#include <vtkCommand.h>
#include <vtkTransform.h>
#include <vtkBoxRepresentation.h>


#include <feel/feelcore/feel.hpp>

#include <AngioTkCenterlineField.h>
#include <Gmsh.h>
#include <Context.h>

#include <centerlinesmanageriodata.hpp>


class CenterlinesManagerWindowInteractor
{
public :
  CenterlinesManagerWindowInteractor();
  CenterlinesManagerWindowInteractor( CenterlinesManagerWindowInteractor const& o ) = default;
  std::string const& inputSurfacePath() const { return M_inputSurfacePath; }
  std::string const& inputPointSetPath() const { return M_inputPointSetPath; }
  std::string const& inputPointPairPath() const { return M_inputPointPairPath; }
  std::vector<std::string> const& inputCenterlinesPath() const { return M_inputCenterlinesPath; }
  std::string const& inputCenterlinesPath(int k) const { return M_inputCenterlinesPath[k]; }
  int windowWidth() const { return M_windowWidth; }
  int windowHeight() const { return M_windowHeight; }

  void setInputSurfacePath(std::string const& path) { M_inputSurfacePath=path; }
  void setInputPointSetPath(std::string const& path) { M_inputPointSetPath=path; }
  void setInputPointPairPath(std::string const& path) { M_inputPointPairPath=path; }
  void setInputCenterlinesPath(std::string const& path) { M_inputCenterlinesPath={ path }; }
  void setInputCenterlinesPath(std::vector<std::string> const& path) { M_inputCenterlinesPath=path; }
  void setWindowWidth(int i) { M_windowWidth = i; }
  void setWindowHeight(int i) { M_windowHeight = i; }

  void run();


private :
  std::string M_inputSurfacePath, M_inputPointSetPath, M_inputPointPairPath;
  std::vector<std::string> M_inputCenterlinesPath;
  int M_windowWidth, M_windowHeight;

};

#endif // __CENTERLINESMANAGERWINDOWINTERACTOR_H
