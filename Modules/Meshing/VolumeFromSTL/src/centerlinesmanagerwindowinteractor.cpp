/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <centerlinesmanagerwindowinteractor.hpp>

#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>

#include <vtkWorldPointPicker.h>
#include <vtkPointPicker.h>
#include <vtkPropPicker.h>

#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>

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

class vtkMyCallback : public vtkCommand
{
public:
  //vtkMyCallback()
  //:
    //M_sphereSource(NULL)
  //{}
  static vtkMyCallback *New() 
    { return new vtkMyCallback; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
      vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
      // apply transformation to box and sphere inside
      widget->GetTransform(t);
      widget->GetProp3D()->SetUserTransform(t);
    }

};
 

class SphereSourceObject : public std::tuple< vtkSmartPointer<vtkSphereSource>, vtkSmartPointer<vtkActor>, int >
{
public :
  typedef std::tuple< vtkSmartPointer<vtkSphereSource>, vtkSmartPointer<vtkActor>, int > super_type;

  SphereSourceObject()
    :
    super_type(std::make_tuple( vtkSmartPointer<vtkSphereSource>(),vtkSmartPointer<vtkActor>(),0 ))
  {}
  SphereSourceObject( super_type const& o )
    :
    super_type( o )
  {}
  SphereSourceObject( SphereSourceObject const& o ) = default;

  vtkSmartPointer<vtkSphereSource> const& geometry() const { return std::get<0>(*this); }
  vtkSmartPointer<vtkActor> const& actor() const { return std::get<1>(*this); }
  bool isSourcePoint() const { return (std::get<2>(*this) == 0); }
  bool isTargetPoint() const { return (std::get<2>(*this) == 1); }
  int typePoint() const { return std::get<2>(*this); }
  void setTypePoint(int k) { std::get<2>(*this) = k; }
  void changeTypePoint() { std::get<2>(*this) = (std::get<2>(*this)+1)%2; }

  bool isNull() const { return ( this->geometry() != NULL || this->actor() != NULL); }

  void applyColoring()
  {
    if ( this->actor() == NULL ) return;
    if ( this->isSourcePoint() )
      this->actor()->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
    else if ( this->isTargetPoint() )
      this->actor()->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
  }


};

// Define interaction style
class AngioTkWindowInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static AngioTkWindowInteractorStyle* New();
    vtkTypeMacro(AngioTkWindowInteractorStyle, vtkInteractorStyleTrackballCamera);

    AngioTkWindowInteractorStyle()
        :
        M_activatedMode( ModeType::NO_MODE ),
        M_lastSphereActiveId(-1),
        M_sphereActorSelection(NULL),
        M_sphereActorSelectionId(-1),
        M_lenghtSTL(0),
        M_outputPathPointSetFile("default.data")
    {
        M_LastPickedActor = NULL;
        M_LastPickedProperty = vtkProperty::New();

        M_widgetTextCommandHelp = vtkSmartPointer<vtkTextWidget>::New();
        M_widgetTextCommandHelp->SelectableOff();
        vtkSmartPointer<vtkTextActor> textActorCommandHelpCommon = vtkSmartPointer<vtkTextActor>::New();
        std::ostringstream helpCommandStr;
        helpCommandStr << "Help Commands :\n"
                       << "q : exit  \n"
                       << "h : enable/disable help commands  \n"
                       << "a : enable/disable orientation axis\n"
                       << "s : save points insertion on disk\n"
                       << "0 : no mode (view only)\n"
                       << "1 : mode surface representation\n"
                       << "2 : mode points insertion\n"
                       << "3 : mode centerlines manager\n";
        textActorCommandHelpCommon->SetInput ( helpCommandStr.str().c_str() );
        vtkSmartPointer<vtkTextRepresentation> textRepresentation = vtkSmartPointer<vtkTextRepresentation>::New();
        //textRepresentation->GetPositionCoordinate()->SetValue( .15, .15 );
        //textRepresentation->GetPosition2Coordinate()->SetValue( .7, .2 );
        M_widgetTextCommandHelp->SetRepresentation( textRepresentation );
        textRepresentation->SetTextActor( textActorCommandHelpCommon );
        textRepresentation->SetPosition(0.7,0.2);
        textRepresentation->SetPosition2(0.3,0.8);
        textActorCommandHelpCommon->GetTextProperty()->SetColor ( 0.0,0.0,0.0 );
        textActorCommandHelpCommon->GetTextProperty()->SetJustificationToLeft();
        textActorCommandHelpCommon->GetTextProperty()->SetVerticalJustificationToTop();
        textActorCommandHelpCommon->GetTextProperty()->SetFontSize( 24 );

        M_widgetTextActivatedMode = vtkSmartPointer<vtkTextWidget>::New();
        M_widgetTextActivatedMode->SelectableOff();
        vtkSmartPointer<vtkTextRepresentation> textRepresentation2 = vtkSmartPointer<vtkTextRepresentation>::New();
        M_widgetTextActivatedMode->SetRepresentation( textRepresentation2 );
        vtkSmartPointer<vtkTextActor> textActor2 = vtkSmartPointer<vtkTextActor>::New();
        textRepresentation2->SetTextActor( textActor2 );
        //textRepresentation2->SetPosition(0.01,0.2);
        //textRepresentation2->SetPosition2(0.3,0.8);
        textRepresentation2->SetPosition(0.01,0.95);
        textRepresentation2->SetPosition2(0.3,0.05);
        textActor2->GetTextProperty()->SetColor ( 0.0,0.0,0.0 );
        textActor2->GetTextProperty()->SetJustificationToLeft();
        textActor2->GetTextProperty()->SetVerticalJustificationToTop();
        textActor2->GetTextProperty()->SetFontSize( 24 );



        M_widgetOrientationAxis = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
        M_widgetOrientationAxis->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
        vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
        M_widgetOrientationAxis->SetOrientationMarker( axes );
        M_widgetOrientationAxis->SetViewport( 0.0, 0.0, 0.4, 0.4 );
        /*M_widgetOrientationAxis->SetInteractor( this->Interactor );
          M_widgetOrientationAxis->SetEnabled( 1 );
          M_widgetOrientationAxis->InteractiveOn();*/

        M_widgetBoxAroundSphere = vtkSmartPointer<vtkBoxWidget>::New();
        M_callbackBoxAroundSphere = vtkSmartPointer<vtkMyCallback>::New();
        M_widgetBoxAroundSphere->AddObserver(vtkCommand::InteractionEvent, M_callbackBoxAroundSphere);
        //M_widgetBoxAroundSphere->HandlesOff();
    }


    static
    std::string
    helpTextModeSurfaceRepresentation( bool showHelpCommands )
    {
        std::ostringstream helpStr;
        helpStr << "Mode : Surface Representation :                     \n";
        if ( showHelpCommands )
            helpStr << "w : change representation (surface,wireframe)\n"
                    << "[Up] : increase Opacity\n"
                    << "[Down] : decrease Opacity\n";
        return helpStr.str();
    }

    static
    std::string
    helpTextModePointsInsertion( bool showHelpCommands )
    {
        std::ostringstream helpStr;
        helpStr << "Mode : Points Insertion :                          \n";
        if ( showHelpCommands )
            helpStr << "c : change point type (source,target)\n"
                    << "y : valide point\n"
                    << "r : remove selection\n"
                    << "u : undo insertion\n"
                    << "[Up] : move selection\n"
                    << "[Down] : move selection\n"
                    << "[Left] : move selection\n"
                    << "[Right] : move selection\n"
                    << "o : move selection\n"
                    << "l : move selection\n";
        return helpStr.str();
    }

    static
    std::string
    helpTextModeCenterlinesManger(bool showHelpCommands )
    {
        std::ostringstream helpStr;
        helpStr << "Mode : Centerlines Manager :                          \n";
        if ( showHelpCommands )
            helpStr << "e : points insertion at extremities\n";
        return helpStr.str();
    }
    void saveOnDisk( std::string const& pathFile )
    {
        std::cout << "saveOnDisk " << pathFile << "\n";
        if ( M_vectorSphereSourceObject.empty() ) return;
        if ( M_vectorSphereSourceObject.front().geometry() == NULL ) return;

        std::ofstream fileWrited( pathFile, std::ios::out | std::ios::trunc);
        for (int k = 0; k< M_vectorSphereSourceObject.size() ;++ k)
        {
            if ( M_vectorSphereSourceObject[k].geometry() != NULL )
            {
                double radius = M_vectorSphereSourceObject[k].geometry()->GetRadius();
                double * center = M_vectorSphereSourceObject[k].geometry()->GetCenter();
                int typePt = M_vectorSphereSourceObject[k].typePoint();
                fileWrited << typePt << " " << center[0] << " " << center[1] << " " << center[2] << " " << radius << "\n";
            }
        }
        fileWrited.close();
    }

    void loadFromDisk( std::string const& pathFile )
    {
        std::ifstream fileLoaded( pathFile, std::ios::in);
        while ( !fileLoaded.eof() )
        {
            double * center = new double[3];
            double radius;
            int typePt = -1;
            fileLoaded >> typePt;
            if ( fileLoaded.eof() ) break;
            fileLoaded >> center[0] >> center[1] >> center[2] >> radius;
            //std::cout << "add center " << center[0] <<","<< center[1] <<","<< center[2]<<"\n";
            vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
            sphereSource->SetRadius(radius);
            sphereSource->SetCenter(center);
            delete [] center;
            sphereSource->Update();
            //Create a mapper and actor
            vtkSmartPointer<vtkPolyDataMapper> mapperSphere = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapperSphere->SetInputConnection(sphereSource->GetOutputPort());
            vtkSmartPointer<vtkActor> actorSphere = vtkSmartPointer<vtkActor>::New();
            actorSphere->SetMapper(mapperSphere);
            actorSphere->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)

            M_vectorSphereSourceObject.push_back(std::make_tuple(sphereSource,actorSphere,typePt) );
            M_vectorSphereSourceObject.back().applyColoring();
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actorSphere);
        }

        M_vectorSphereSourceObject.push_back(SphereSourceObject());

        fileLoaded.close();
        this->Interactor->GetRenderWindow()->Render();
    }

    void addSphereAtCenterlinesExtrimities()
    {
        if (!M_angioTkCenterlines) return;

        for ( auto const& extremityPair : M_angioTkCenterlines->centerlinesExtremities() )
        {
            int typePt = 0;//-1;
            double * center = new double[3];
            center[0] = extremityPair.first->x();
            center[1] = extremityPair.first->y();
            center[2] = extremityPair.first->z();
            //std::cout << "add center " << center[0] <<","<< center[1] <<","<< center[2]<<"\n";

            int branchId = extremityPair.second.first;
            int lineIdInBranch = extremityPair.second.second;
            std::vector<MLine*> mylines = M_angioTkCenterlines->centerlinesBranch(branchId).lines;
            MLine* myline = mylines[lineIdInBranch];
            //std::map<MLine*,double>::iterator itr = M_angioTkCenterlines->centerlinesRadiusl().find(myline);
            //auto itr = M_angioTkCenterlines->centerlinesRadiusl().find(myline);
            double radius = -1.;
            //if ( itr != M_angioTkCenterlines->centerlinesRadiusl().end() )
            if ( M_angioTkCenterlines->centerlinesRadiusl().find(myline) != M_angioTkCenterlines->centerlinesRadiusl().end() )
            {
                //radius = itr->second;// not always work, strange!!
                radius = M_angioTkCenterlines->centerlinesRadiusl().find(myline)->second;
                //radius = M_angioTkCenterlines->centerlinesBranch(branchId).minRad;   //itr->second;//0.8;
            }
            //std::cout << "add extremityPair with radius " << radius << " branchId " << branchId << " lineIdInBranch " << lineIdInBranch << "\n";

            vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
            sphereSource->SetRadius(radius);
            sphereSource->SetCenter(center);
            delete [] center;
            sphereSource->Update();
            //Create a mapper and actor
            vtkSmartPointer<vtkPolyDataMapper> mapperSphere = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapperSphere->SetInputConnection(sphereSource->GetOutputPort());
            vtkSmartPointer<vtkActor> actorSphere = vtkSmartPointer<vtkActor>::New();
            actorSphere->SetMapper(mapperSphere);
            actorSphere->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)

            M_vectorSphereSourceObject.push_back(std::make_tuple(sphereSource,actorSphere,typePt) );
            M_vectorSphereSourceObject.back().applyColoring();
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actorSphere);
        }

    }

    int hasActor( vtkActor * _actor ) const
    {
        int k=-1;
        for ( int k = 0 ; k<M_vectorSphereSourceObject.size() ; ++k )
        {
            if ( M_vectorSphereSourceObject[k].actor() == _actor )
                return k;
        }
        return k;
    }

    void setAngioTkCenterlines( boost::shared_ptr<AngioTkCenterline> obj ) { M_angioTkCenterlines = obj; }

    void setActorSTL( vtkSmartPointer<vtkActor> const& actor) { M_actorSTL=actor; }

    void setBoundsSTL(double * bounds ) { M_boundsSTL.insert( M_boundsSTL.begin(),bounds,bounds+6); }
    void setLenghtSTL(double d) { M_lenghtSTL=d; }

    void setOutputPathPointSetFile(std::string const& path) { M_outputPathPointSetFile=path; }


    void activateSphereActorSelection()
    {
        if ( M_sphereActorSelection != NULL )
        {
            // revert actor color
            //M_sphereActorSelection->GetProperty()->SetColor(0.0, 1.0, 0.0); //(R,G,B)
#if 0
            M_sphereActorSelection->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
            M_sphereActorSelection->GetProperty()->SetDiffuse(1.0);
            M_sphereActorSelection->GetProperty()->SetSpecular(0.0);
#endif
            this->Interactor->GetRenderWindow()->Render();
            //M_sphereActorSelection = NULL;


            M_widgetBoxAroundSphere->SetInteractor(this->Interactor);
            //M_widgetBoxAroundSphere->SetPlaceFactor(1.25);
            M_widgetBoxAroundSphere->SetPlaceFactor(1.);
            M_widgetBoxAroundSphere->SetProp3D(M_sphereActorSelection/*coneActor*/);
            M_vectorSphereSourceObject[M_sphereActorSelectionId].actor()->GetProperty()->SetOpacity(0.2);

            M_widgetBoxAroundSphere->PlaceWidget();
            //M_widgetBoxAroundSphere->ScalingEnabledOff();
            //vtkSmartPointer<vtkMyCallback> callback = vtkSmartPointer<vtkMyCallback>::New();
            //M_widgetBoxAroundSphere->AddObserver(vtkCommand::InteractionEvent, callback);

            //M_sphereActorSelection = M_vectorSphereSourceObject[selectId].actor();//LastPickedActor;
            //M_sphereActorSelectionId = selectId;
            //->geometry()
#if 0
            M_widgetBoxAroundSphere->HandlesOn();
            M_widgetBoxAroundSphere->OutlineFaceWiresOff();
            M_widgetBoxAroundSphere->OutlineCursorWiresOff();	
            //M_widgetBoxAroundSphere->TranslationEnabledOff();
            M_widgetBoxAroundSphere->ScalingEnabledOff();
            M_widgetBoxAroundSphere->RotationEnabledOff();
#endif
            double radiusSphere = M_vectorSphereSourceObject[M_sphereActorSelectionId].geometry()->GetRadius();
            M_widgetBoxAroundSphere->SetHandleSize(radiusSphere/10000.);

            M_widgetBoxAroundSphere->On();
            //M_widgetBoxAroundSphere->Print(std::cout);
            //M_callbackBoxAroundSphere->setSphereSource( M_vectorSphereSourceObject[M_sphereActorSelectionId].geometry() );

        }
    }

    void deactivateSphereActorSelection()
    {
        if ( M_sphereActorSelection != NULL )
        {
            // revert actor color
#if 0
            M_sphereActorSelection->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
            M_sphereActorSelection->GetProperty()->SetDiffuse(1.0);
            M_sphereActorSelection->GetProperty()->SetSpecular(0.0);
#endif

#if 1
            double * p=M_widgetBoxAroundSphere->GetProp3D()->GetCenter();
            //std::cout << "Box center: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
            M_vectorSphereSourceObject[M_sphereActorSelectionId].geometry()->SetCenter(p);

            // set identity transformation
            vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
            t->Identity();
            M_widgetBoxAroundSphere->GetProp3D()->SetUserTransform(t);
#endif

            M_vectorSphereSourceObject[M_sphereActorSelectionId].actor()->GetProperty()->SetOpacity(1);


            M_sphereActorSelection = NULL;
            M_sphereActorSelectionId = -1;
            M_widgetBoxAroundSphere->Off();
            //M_callbackBoxAroundSphere->setSphereSource(NULL);

            this->Interactor->GetRenderWindow()->Render();
        }
    }

    virtual void OnLeftButtonDown() 
    {
        if ( M_activatedMode == ModeType::POINTS_INSERTION )
        {

            if ( false )
                std::cout << "Picking pixel: " << this->Interactor->GetEventPosition()[0] << " " << this->Interactor->GetEventPosition()[1] << std::endl;

            //int* clickPos = this->GetInteractor()->GetEventPosition();
            // Pick from this location.
            vtkSmartPointer<vtkPropPicker>  NEWpicker = vtkSmartPointer<vtkPropPicker>::New();
            /*if ( M_lenghtSTL > 50 )
              NEWpicker->SetTolerance(1e-3);
              else
              NEWpicker->SetTolerance(1e-2);*/


            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(M_actorSTL);
            int resultNEW = NEWpicker->Pick(this->Interactor->GetEventPosition()[0], 
                                            this->Interactor->GetEventPosition()[1], 
                                            0,  // always zero.
                                            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(M_actorSTL);
            //if ( LastPickedActor )
            //LastPickedActor->GetProperty()->DeepCopy(LastPickedProperty);
            //if ( M_sphereActorSelectionId >=0 )

            M_LastPickedActor = NEWpicker->GetActor();
            if ( M_LastPickedActor != NULL  && this->hasActor( M_LastPickedActor ) >= 0 )
            {
                int selectId=this->hasActor( M_LastPickedActor );
                if ( M_sphereActorSelectionId < 0 )
                {
                    // Save the property of the picked actor so that we can restore it next time
                    M_sphereActorSelection = M_vectorSphereSourceObject[selectId].actor();//LastPickedActor;
                    M_sphereActorSelectionId = selectId;
                    this->activateSphereActorSelection();
                }
                //else if ( M_sphereActorSelectionId == selectId )
                //{
                //  this->switchModeSphereActorSelection();
                //}
                //this->deactivateSphereActorSelection();
            }
            else if ( M_sphereActorSelectionId < 0 )
            {
                //this->deactivateSphereActorSelection();

                typedef vtkPointPicker picker_type; // vtkPointPicker, vtkWorldPointPicker
                vtkSmartPointer<picker_type> myPicker = vtkSmartPointer<picker_type>::New();
                //myPicker->SetTolerance(0.0005);
                //myPicker->SetTolerance( 1e-4*M_lenghtSTL );
                //myPicker->SetTolerance( (1./M_lenghtSTL) );
                //std::cout << "M_lenghtSTL " << M_lenghtSTL << "\n";
                //myPicker->SetTolerance( 1e-4*M_lenghtSTL );
                if ( M_lenghtSTL > 50 )
                    myPicker->SetTolerance(1e-3);
                else
                    myPicker->SetTolerance(1e-2);
                int result = myPicker->Pick(this->Interactor->GetEventPosition()[0], 
                                            this->Interactor->GetEventPosition()[1], 
                                            0,  // always zero.
                                            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());

                if( result != 0 )
                {
                    double picked[3];
                    //this->Interactor->GetPicker()->GetPickPosition(picked);
                    myPicker->GetPickPosition(picked);
                    if ( false ) std::cout << "Picked value: " << picked[0] << " " << picked[1] << " " << picked[2] << std::endl;

                    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();

                    double lengthX = M_boundsSTL[1]-M_boundsSTL[0];
                    double lengthY = M_boundsSTL[3]-M_boundsSTL[2];
                    double lengthZ = M_boundsSTL[5]-M_boundsSTL[4];
                    double maxLength = std::max( lengthX, std::max(lengthY,lengthZ) );
                    //double myradius = M_lenghtSTL/20.;//maxLength/20.;
                    double myradius = M_lenghtSTL/50.;//maxLength/20.;
                    //std::cout << "myradius " << myradius << "\n";
                    sphereSource->SetRadius(myradius);//20);
                    sphereSource->SetCenter(picked);
                    sphereSource->Update();
                    //Create a mapper and actor
                    vtkSmartPointer<vtkPolyDataMapper> mapperSphere = vtkSmartPointer<vtkPolyDataMapper>::New();
                    mapperSphere->SetInputConnection(sphereSource->GetOutputPort());

                    vtkSmartPointer<vtkActor> actorSphere = vtkSmartPointer<vtkActor>::New();
                    actorSphere->SetMapper(mapperSphere);  
                    //actorSphere->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
                    actorSphere->GetProperty()->SetColor(1.0, 0.5, 0.0); //(R,G,B)

                    if ( !M_vectorSphereSourceObject.empty() )
                    {
                        if ( M_vectorSphereSourceObject.back().geometry() != NULL )
                        {
                            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(M_vectorSphereSourceObject.back().actor());
                        }
                        M_vectorSphereSourceObject.pop_back();
                    }

                    M_vectorSphereSourceObject.push_back(std::make_tuple(sphereSource,actorSphere,0) );
                    //M_vectorActorSphere.push_back(actorSphere);

                    this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actorSphere);
                } // result != 0
            } // else not pick a sphere actor

        }
        // Forward events
        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }

    // allow to remove all preconfigure key (define in vtkInteractorStyle.cxx)
    virtual void OnChar() {}

    virtual void OnKeyPress()
    {
        // Get the keypress
        vtkRenderWindowInteractor *rwi = this->Interactor;
        std::string key = rwi->GetKeySym();

        // Output the key that was pressed
        if ( false )
            std::cout << "Pressed " << key << std::endl;

        if ( key == "q" )
        {
            rwi->ExitCallback();
        }

        // display/remove axis orientation
        if ( key == "a" )
	    {
            static bool isInitOrientationAxisInteractor = false;
            if ( !isInitOrientationAxisInteractor ) {
                M_widgetOrientationAxis->SetInteractor( this->Interactor );
                isInitOrientationAxisInteractor = true;
            }
            static bool isEnableOrientationAxis = false;
            M_widgetOrientationAxis->SetEnabled( !isEnableOrientationAxis );
            isEnableOrientationAxis = !isEnableOrientationAxis;
            if ( isEnableOrientationAxis ) M_widgetOrientationAxis->InteractiveOn();

            this->Interactor->GetRenderWindow()->Render();
        }

        // display/remove help commands
        if ( key == "h" )
        {
            static bool isInitTextCommandHelpInteractor = false;
            if ( !isInitTextCommandHelpInteractor )
            {
                M_widgetTextCommandHelp->SetInteractor(this->Interactor);
                isInitTextCommandHelpInteractor = true;
            }
            M_widgetTextCommandHelp->SetEnabled( (M_widgetTextCommandHelp->GetEnabled()+1)%2 );

            if ( M_activatedMode == ModeType::SURFACE_REPRESENTATION )
            {
                std::string modeMsg = this->helpTextModeSurfaceRepresentation( (bool)M_widgetTextCommandHelp->GetEnabled() );
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->GetTextActor()->SetInput( modeMsg.c_str() );
            }
            else if ( M_activatedMode == ModeType::POINTS_INSERTION )
            {
                std::string modeMsg = this->helpTextModePointsInsertion( (bool)M_widgetTextCommandHelp->GetEnabled() );
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->GetTextActor()->SetInput( modeMsg.c_str() );
            }
            else if ( M_activatedMode == ModeType::CENTERLINES_MANAGER )
            {
                std::string modeMsg = this->helpTextModeCenterlinesManger( (bool)M_widgetTextCommandHelp->GetEnabled() );
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->GetTextActor()->SetInput( modeMsg.c_str() );
            }

            if ( M_widgetTextCommandHelp->GetEnabled() )
            {
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->SetPosition(0.01,0.2);
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->SetPosition2(0.3,0.8);
            }
            else
            {
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->SetPosition(0.01,0.95);
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->SetPosition2(0.3,0.05);
            }

            this->Interactor->GetRenderWindow()->Render();
        }


        if ( key == "0" && M_activatedMode != ModeType::NO_MODE)
        {
            this->deactivateSphereActorSelection();
            M_activatedMode = ModeType::NO_MODE;
            M_widgetTextActivatedMode->SetEnabled( false );
            this->Interactor->GetRenderWindow()->Render();
        }
        if ( key == "1" || key == "2" || key == "3" )
        {
            static bool isInitTextActivatedMode = false;
            if ( !isInitTextActivatedMode )
            {
                M_widgetTextActivatedMode->SetInteractor(this->Interactor);
                isInitTextActivatedMode = true;
            }

            if ( key == "1" && M_activatedMode != ModeType::SURFACE_REPRESENTATION )
            {
                this->deactivateSphereActorSelection();
                M_activatedMode = ModeType::SURFACE_REPRESENTATION;
                std::string modeMsg = this->helpTextModeSurfaceRepresentation( (bool)M_widgetTextCommandHelp->GetEnabled() );
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->GetTextActor()->SetInput( modeMsg.c_str() );
                M_widgetTextActivatedMode->SetEnabled( true );
                this->Interactor->GetRenderWindow()->Render();
            }
            if ( key == "2" && M_activatedMode != ModeType::POINTS_INSERTION )
            {
                M_activatedMode = ModeType::POINTS_INSERTION;
                std::string modeMsg = this->helpTextModePointsInsertion( (bool)M_widgetTextCommandHelp->GetEnabled() );
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->GetTextActor()->SetInput( modeMsg.c_str() );
                M_widgetTextActivatedMode->SetEnabled( true );
                this->Interactor->GetRenderWindow()->Render();
            }
            if ( key == "3" && M_activatedMode != ModeType::CENTERLINES_MANAGER )
            {
                this->deactivateSphereActorSelection();
                M_activatedMode = ModeType::CENTERLINES_MANAGER;
                std::string modeMsg = this->helpTextModeCenterlinesManger( (bool)M_widgetTextCommandHelp->GetEnabled() );
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->GetTextActor()->SetInput( modeMsg.c_str() );
                M_widgetTextActivatedMode->SetEnabled( true );
                this->Interactor->GetRenderWindow()->Render();
            }
        }


      // save current state on disk
      if ( key == "s" )
      {
          this->saveOnDisk( M_outputPathPointSetFile );
      }



      if ( M_activatedMode == ModeType::SURFACE_REPRESENTATION )
      {
          if ( key == "w" )
          {
              int currentRep = M_actorSTL->GetProperty()->GetRepresentation();
              if (currentRep == VTK_SURFACE )
                  M_actorSTL->GetProperty()->SetRepresentation( VTK_WIREFRAME );
              else
                  M_actorSTL->GetProperty()->SetRepresentation( VTK_SURFACE );
              this->Interactor->GetRenderWindow()->Render();
          }
          if( key == "Up" || key == "Down" )
          {
              double changeOpacityStep = ( key == "Up" )? 0.05 : -0.05;
              double previousOpacity = M_actorSTL->GetProperty()->GetOpacity();
              double newOpacity = previousOpacity+changeOpacityStep;
              newOpacity = std::max( 0., newOpacity );
              newOpacity = std::min( 1., newOpacity );
              M_actorSTL->GetProperty()->SetOpacity(newOpacity);
              this->Interactor->GetRenderWindow()->Render();
          }
      }

      if ( M_activatedMode == ModeType::POINTS_INSERTION )
      {
          if ( key == "Escape" )
          {
              this->deactivateSphereActorSelection();
          }

          // increase (p) decrease (m) current sphere radius
          if( key == "p" )
          {
              if ( 	M_sphereActorSelectionId >=0 )
	          {
                  double prevRadius = M_vectorSphereSourceObject[M_sphereActorSelectionId].geometry()->GetRadius();
                  double step = M_lenghtSTL/500.;
                  M_vectorSphereSourceObject[M_sphereActorSelectionId].geometry()->SetRadius(prevRadius+step);
                  this->Interactor->GetRenderWindow()->Render();
              }
          }
          if( key == "m" )
          {
              if ( 	M_sphereActorSelectionId >=0 )
              {
                  double prevRadius = M_vectorSphereSourceObject[M_sphereActorSelectionId].geometry()->GetRadius();
                  double step = M_lenghtSTL/500.;
                  M_vectorSphereSourceObject[M_sphereActorSelectionId].geometry()->SetRadius(prevRadius-step);
                  this->Interactor->GetRenderWindow()->Render();
              }
          }
          // move in x,y,z direction
          if( key == "Left" || key == "Right" || key == "Up" || key == "Down" || key == "o" || key == "l" )
          {
              if ( M_sphereActorSelection != NULL )
              {
                  double lengthX = M_boundsSTL[1]-M_boundsSTL[0];
                  double lengthY = M_boundsSTL[3]-M_boundsSTL[2];
                  double lengthZ = M_boundsSTL[5]-M_boundsSTL[4];
                  double maxLength = std::max( lengthX, std::max(lengthY,lengthZ) );
                  //double movingValue = maxLength/50.;
                  double movingValue = M_lenghtSTL/500.;

                  double * translateVec = new double[3];
                  translateVec[0] = 0.;translateVec[1] = 0.;translateVec[2] = 0.;
                  if(key == "Left")
                      translateVec[0] = -movingValue;
                  if(key == "Right")
                      translateVec[0] = movingValue;
                  if(key == "Up")
                      translateVec[1] = movingValue;
                  if(key == "Down")
                      translateVec[1] = -movingValue;
                  if(key == "o")
                      translateVec[2] = movingValue;
                  if(key == "l")
                      translateVec[2] = -movingValue;

                  vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
                  M_widgetBoxAroundSphere->GetTransform(t);
                  t->Translate(translateVec);
                  M_widgetBoxAroundSphere->SetTransform(t);
                  M_widgetBoxAroundSphere->GetProp3D()->SetUserTransform(t);
                  //M_widgetBoxAroundSphere->GetProp3D()->SetUserTransform(t);
                  //M_widgetBoxAroundSphere->GetProp3D()->Concatenate 
                  delete [] translateVec;
                  this->Interactor->GetRenderWindow()->Render();
              }
          }

          // valid sphere
          if(key == "y")
          {
              this->deactivateSphereActorSelection();
              if ( !M_vectorSphereSourceObject.empty() && M_vectorSphereSourceObject.back().geometry() != NULL )
              {
                  M_vectorSphereSourceObject.back().applyColoring();

                  //M_lastSphereActiveId = M_vectorSphereSourceObject.size()-1;
                  //vtkSmartPointer<vtkSphereSource> sphereSource;
                  //vtkSmartPointer<vtkActor> actorSphere;
                  M_vectorSphereSourceObject.push_back(SphereSourceObject());//std::make_tuple(sphereSource,actorSphere,0));
              }
              this->Interactor->GetRenderWindow()->Render();
          }

          // change propertie of point (used with centerlines algo) : source or target point
          if(key == "c")
          {
              if ( M_sphereActorSelectionId >= 0 )
              {
                  M_vectorSphereSourceObject[M_sphereActorSelectionId].changeTypePoint();
                  if ( M_sphereActorSelectionId < (M_vectorSphereSourceObject.size()-1)  )
                      M_vectorSphereSourceObject[M_sphereActorSelectionId].applyColoring(/*M_sphereActorSelection*/);
#if 0
                  if ( M_vectorSphereSourceObject[M_sphereActorSelectionId].isSourcePoint() )
                      M_sphereActorSelection->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
                  else if ( M_vectorSphereSourceObject[M_sphereActorSelectionId].isTargetPoint() )
                      M_sphereActorSelection->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
#endif
                  this->Interactor->GetRenderWindow()->Render();
              }
          }

          // undo key
          if(key == "u")
          {
              if ( !M_vectorSphereSourceObject.empty() )
              {
                  this->deactivateSphereActorSelection();
                  if ( M_vectorSphereSourceObject.back().geometry() != NULL )
                  {
                      // remove last actor in render
                      this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(M_vectorSphereSourceObject.back().actor());
                      // delete last actor
                      M_vectorSphereSourceObject.pop_back();
                      // add a null pointer to tell that next picker is a test actor
                      //vtkSmartPointer<vtkSphereSource> sphereSource;
                      //vtkSmartPointer<vtkActor> actorSphere;
                      M_vectorSphereSourceObject.push_back(SphereSourceObject());//std::make_tuple(sphereSource,actorSphere,0));
                  }
                  else if ( M_vectorSphereSourceObject.size() > 1 ) // if size==0 else keep null pointer at case 0
                  {
                      // delete null pointer
                      M_vectorSphereSourceObject.pop_back();
                      // remove last actor in render
                      this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(M_vectorSphereSourceObject.back().actor());
                      // delete last actor
                      M_vectorSphereSourceObject.pop_back();
                      //vtkSmartPointer<vtkSphereSource> sphereSource;
                      //vtkSmartPointer<vtkActor> actorSphere;
                      M_vectorSphereSourceObject.push_back(SphereSourceObject());//std::make_tuple(sphereSource,actorSphere,0));
                      --M_lastSphereActiveId;
                  }
                  this->Interactor->GetRenderWindow()->Render();
              }
          }
          if(key == "r")
          {
              if ( M_sphereActorSelectionId >= 0 )
              {
                  int selectId = M_sphereActorSelectionId;
                  this->deactivateSphereActorSelection();
                  this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor( M_vectorSphereSourceObject[selectId].actor() );
                  //M_vectorSphereSourceObject[M_sphereActorSelectionId].changeTypePoint();
                  bool isLastId = ( selectId == (M_vectorSphereSourceObject.size()-1) );
                  M_vectorSphereSourceObject.erase( M_vectorSphereSourceObject.begin()+selectId );
                  if (isLastId )
                  {
                      //vtkSmartPointer<vtkSphereSource> sphereSource;
                      //vtkSmartPointer<vtkActor> actorSphere;
                      M_vectorSphereSourceObject.push_back(SphereSourceObject());//std::make_tuple(sphereSource,actorSphere,0));
                  }

                  //--M_lastSphereActiveId;
                  this->Interactor->GetRenderWindow()->Render();
              }
          }

      }

      if ( M_activatedMode == ModeType::CENTERLINES_MANAGER )
      {
          if(key == "e")
          {
            static bool hasAlreadyDonePointsInsertionAtExtremities = false;
            if ( !hasAlreadyDonePointsInsertionAtExtremities )
            {
                this->addSphereAtCenterlinesExtrimities();
                this->Interactor->GetRenderWindow()->Render();
                hasAlreadyDonePointsInsertionAtExtremities = true;
            }

          }
      }

    }


public :
    enum ModeType { NO_MODE=0, SURFACE_REPRESENTATION = 1, POINTS_INSERTION = 2, CENTERLINES_MANAGER = 3 };


private :
    ModeType M_activatedMode;

    vtkActor * M_LastPickedActor;
    vtkProperty *M_LastPickedProperty;
    int M_lastSphereActiveId;
    vtkActor * M_sphereActorSelection;
    int M_sphereActorSelectionId;


    std::vector<double> M_boundsSTL; double M_lenghtSTL;
    //vtkSphereSource * M_sphereSourceSelection;

    std::vector<SphereSourceObject> M_vectorSphereSourceObject;

    vtkSmartPointer<vtkTextWidget> M_widgetTextCommandHelp;
    vtkSmartPointer<vtkTextWidget> M_widgetTextActivatedMode;
    vtkSmartPointer<vtkOrientationMarkerWidget> M_widgetOrientationAxis; 
    vtkSmartPointer<vtkBoxWidget> M_widgetBoxAroundSphere;
    vtkSmartPointer<vtkMyCallback> M_callbackBoxAroundSphere;

    vtkSmartPointer<vtkActor> M_actorSTL;

    boost::shared_ptr<AngioTkCenterline> M_angioTkCenterlines;

    std::string M_outputPathPointSetFile;

};
vtkStandardNewMacro(AngioTkWindowInteractorStyle);


void
CenterlinesManagerWindowInteractor::run()
{
    typedef vtkPointPicker picker_type; // vtkPointPicker, vtkWorldPointPicker
    //vtkSmartPointer<picker_type> worldPointPicker = vtkSmartPointer<picker_type>::New();
    //worldPointPicker->SetTolerance(0.0005);
    //--------------------------------------------
    // read stl
    if ( this->inputSurfacePath().empty() )
    {
        std::cout << "Required parameters: Filename" << endl;
        return;
    }
    // Create a renderer, render window, and interactor
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    //renderWindowInteractor->SetPicker(worldPointPicker);
    renderWindowInteractor->SetRenderWindow(renderWindow);

#if 0
    //renderer->SetBackground(1,1,1); // Background color white
    renderer->SetBackground(0.5,0.5,0.5); // Background color white
#else
    // Setup the background gradient
    renderer->GradientBackgroundOn();
    renderer->SetBackground(1,1,1);
    renderer->SetBackground2(0,0,1);
#endif
    vtkSmartPointer<AngioTkWindowInteractorStyle> style = vtkSmartPointer<AngioTkWindowInteractorStyle>::New();
    renderWindowInteractor->SetInteractorStyle( style );


    // load surface
    std::string inputFilename = this->inputSurfacePath();
    vtkSmartPointer<vtkSTLReader> readerSTL = vtkSmartPointer<vtkSTLReader>::New();
    readerSTL->SetFileName(inputFilename.c_str());
    readerSTL->Update();
    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapperSTL = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperSTL->SetInputConnection(readerSTL->GetOutputPort());
    vtkSmartPointer<vtkActor> actorSTL = vtkSmartPointer<vtkActor>::New();
    actorSTL->SetMapper(mapperSTL);
    actorSTL->GetProperty()->SetOpacity(0.25);//0.8);//0.7//0.3
    //actorSTL->GetProperty()->SetOpacity(0.8);//0.7//0.3

    // Add the actor to the scene
    renderer->AddActor(actorSTL);

    // update info in style
    style->setActorSTL( actorSTL );
    style->setBoundsSTL(mapperSTL->GetBounds());
    style->setLenghtSTL(mapperSTL->GetLength());

    std::string nameFile = Feel::fs::path(this->inputSurfacePath()).stem().string()+"_pointset.data";
    style->setOutputPathPointSetFile(nameFile);

    boost::shared_ptr<AngioTkCenterline> centerlinesTool;
    for ( int k=0;k<this->inputCenterlinesPath().size();++k)
    {
        if ( this->inputCenterlinesPath(k).empty() || !Feel::fs::exists( this->inputCenterlinesPath(k) ) ) continue;

        vtkSmartPointer<vtkPolyDataReader> readerVTK = vtkSmartPointer<vtkPolyDataReader>::New();
        readerVTK->SetFileName(this->inputCenterlinesPath(k).c_str());
        readerVTK->Update();
        // Create a mapper and actor
        vtkSmartPointer<vtkPolyDataMapper> mapperVTK = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapperVTK->SetInputConnection(readerVTK->GetOutputPort());
        vtkSmartPointer<vtkActor> actorVTK = vtkSmartPointer<vtkActor>::New();
        actorVTK->SetMapper(mapperVTK);

        if (k%8==0)
            actorVTK->GetProperty()->SetColor(0.0, 1.0, 0.0); //(R,G,B)
        if (k%8==1)
            actorVTK->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
        if (k%8==2)
            actorVTK->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
        if (k%8==3)
            actorVTK->GetProperty()->SetColor(0.5, 0.5, 0.0); //(R,G,B)
        if (k%8==5)
            actorVTK->GetProperty()->SetColor(0.5, 0.0, 0.5); //(R,G,B)
        if (k%8==6)
            actorVTK->GetProperty()->SetColor(0.0, 0.5, 0.5); //(R,G,B)
        if (k%8==7)
            actorVTK->GetProperty()->SetColor(0.5, 0.5, 0.5); //(R,G,B)

        actorVTK->GetProperty()->SetLineWidth(3.0);
        // Add the actor to the scene
        renderer->AddActor(actorVTK);
        if ( !centerlinesTool )
        {
            if ( true )
            {
                CTX::instance()->terminal = 1;
                int verbosityLevel = 5;
                Msg::SetVerbosity( verbosityLevel );
            }
            centerlinesTool.reset( new AngioTkCenterline );
            centerlinesTool->importSurfaceFromFile( this->inputSurfacePath() );
            style->setAngioTkCenterlines(centerlinesTool);
            //centerlinesTool->importFile( this->inputCenterlinesPath(k) );
        }
        centerlinesTool->importFile( this->inputCenterlinesPath(k) );
    }
    //if ( centerlinesTool && false )
    //    style->addSphereAtCenterlinesExtrimities();

    if ( !this->inputPointSetPath().empty() )
    {
        std::string previousData = this->inputPointSetPath();
        style->loadFromDisk(previousData);
    }

    renderer->ResetCamera();
    //ren->GetActiveCamera()->Zoom(1.5);

    // Render and interact
    renderWindow->Render();
    renderWindowInteractor->Start();

    //return EXIT_SUCCESS;

}

