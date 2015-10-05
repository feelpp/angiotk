/*=========================================================================

Program:   VMTK
Module:    $RCSfile: vtkvmtkPolyBallLine.h,v $
Language:  C++
Date:      $Date: 2006/04/06 16:46:43 $
Version:   $Revision: 1.4 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
  // .NAME angiotkPolyBallLine - 
  // .SECTION Description
  // ..

#ifndef __angiotkPolyBallLine_h
#define __angiotkPolyBallLine_h

#include <map>
#include <set>

#include <vtkImplicitFunction.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
//#include "vtkvmtkComputationalGeometryWin32Header.h"
#include <vtkvmtkWin32Header.h>

class VTK_VMTK_COMPUTATIONAL_GEOMETRY_EXPORT angiotkPolyBallLine : public vtkImplicitFunction
{
  public:

  static angiotkPolyBallLine *New();
  vtkTypeMacro(angiotkPolyBallLine,vtkImplicitFunction);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description
  // Evaluate polyball.
  double EvaluateFunction(double x[3]);
  double EvaluateFunction(double x, double y, double z)
  {return this->vtkImplicitFunction::EvaluateFunction(x, y, z); } ;

  // Description
  // Evaluate polyball gradient.
  void EvaluateGradient(double x[3], double n[3]);

  // Description:
  // Set / get input poly data.
  vtkSetObjectMacro(Input,vtkPolyData);
  vtkGetObjectMacro(Input,vtkPolyData);

  // Description:
  // Set / get input cell ids used for the function.
  vtkSetObjectMacro(InputCellIds,vtkIdList);
  vtkGetObjectMacro(InputCellIds,vtkIdList);

  // Description:
  // Set / get a single input cell id used for the function.
  vtkSetMacro(InputCellId,vtkIdType);
  vtkGetMacro(InputCellId,vtkIdType);

  // Description:
  // Set / get poly ball radius array name.
  vtkSetStringMacro(PolyBallRadiusArrayName);
  vtkGetStringMacro(PolyBallRadiusArrayName);

  // Description:
  // Get the id of the last nearest poly ball center.
  vtkGetMacro(LastPolyBallCellId,vtkIdType);
  vtkGetMacro(LastPolyBallCellSubId,vtkIdType);
  vtkGetMacro(LastPolyBallCellPCoord,double);
  vtkGetVectorMacro(LastPolyBallCenter,double,3);
  vtkGetMacro(LastPolyBallCenterRadius,double);

  vtkSetMacro(UseRadiusInformation,int);
  vtkGetMacro(UseRadiusInformation,int);
  vtkBooleanMacro(UseRadiusInformation,int);

  void Update();
  std::set<vtkIdType> pointIsOnBranchBounds(double x[3]) const;

  static double ComplexDot(double x[4], double y[4]);

  protected:
  angiotkPolyBallLine();
  ~angiotkPolyBallLine();

  vtkPolyData* Input;
  vtkIdList* InputCellIds;
  vtkIdType InputCellId;

  char* PolyBallRadiusArrayName;
  vtkIdType LastPolyBallCellId;
  vtkIdType LastPolyBallCellSubId;
  double LastPolyBallCellPCoord;
  double LastPolyBallCenter[3];
  double LastPolyBallCenterRadius;

  int UseRadiusInformation;

  vtkDataArray *polyballRadiusArray;
  vtkIdList* cellIds;
  std::map<vtkIdType, std::vector<double> > branchBounds;


  private:
  angiotkPolyBallLine(const angiotkPolyBallLine&);  // Not implemented.
  void operator=(const angiotkPolyBallLine&);  // Not implemented.
};

#endif


