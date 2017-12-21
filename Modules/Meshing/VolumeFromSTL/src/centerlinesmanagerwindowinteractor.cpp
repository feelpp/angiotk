/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include "centerlinesmanagerwindowinteractor.hpp"

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

class PointPairConnection;

class vtkCallbackPointPairConnection : public vtkCommand
{
public:
  static vtkCallbackPointPairConnection *New()
    { return new vtkCallbackPointPairConnection; }
    void initialize( PointPairConnection * object,int ptId )
    {
        M_pointPairConnection = object;
        M_ptId = ptId;
    }
    virtual void Execute(vtkObject *caller, unsigned long, void*);

private :
    PointPairConnection * M_pointPairConnection;
    int M_ptId;
};


class SphereSourceObject
{
public :

    SphereSourceObject( vtkSmartPointer<vtkSphereSource> geo, vtkSmartPointer<vtkActor> actor,
                        int typePoint = 0, bool isValidated = false )
        :
        M_geometry( geo ),
        M_actor( actor ),
        M_typePoint( typePoint ),
        M_isValidated( isValidated )
    {}

    SphereSourceObject( SphereSourceObject const& o ) = default;

    vtkSmartPointer<vtkSphereSource> const& geometry() const { return M_geometry; }
    vtkSmartPointer<vtkActor> const& actor() const { return M_actor; }
    int typePoint() const { return M_typePoint; }
    void setTypePoint(int k) { M_typePoint = k; }

    vtkSmartPointer<vtkBoxWidget> const& widgetBoxAroundSphere() const { return M_widgetBoxAroundSphere; }


    bool isValidated() const { return M_isValidated; }
    void setIsValidated(bool b) { M_isValidated = b; }

    bool isSourcePoint() const { return (this->typePoint() == 0); }
    bool isTargetPoint() const { return (this->typePoint() == 1); }
    void changeTypePoint() { this->setTypePoint( (this->typePoint()+1)%2); }

    bool isNull() const { return ( this->geometry() == NULL || this->actor() == NULL); }

    void applyColoring()
    {
        if ( !this->isValidated() ) return;
        if ( this->actor() == NULL ) return;
        if ( this->isSourcePoint() )
            this->actor()->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
        else if ( this->isTargetPoint() )
            this->actor()->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
#if 0
        this->actor()->GetProperty()->SetDiffuse(1.0);
        this->actor()->GetProperty()->SetSpecular(0.0);
#endif
    }

    void
    initializeWidgetBox( vtkRenderWindowInteractor *interactor )
    {
        if ( this->isNull() ) return;

        if ( !M_widgetBoxAroundSphere )
        {
            M_widgetBoxAroundSphere = vtkSmartPointer<vtkBoxWidget>::New();
            M_callbackBoxAroundSphere = vtkSmartPointer<vtkMyCallback>::New();
            M_widgetBoxAroundSphere->AddObserver(vtkCommand::InteractionEvent, M_callbackBoxAroundSphere);
            //M_widgetBoxAroundSphere->HandlesOff();
        }
        M_widgetBoxAroundSphere->SetInteractor( interactor );
        //M_widgetBoxAroundSphere->SetPlaceFactor(1.25);
        M_widgetBoxAroundSphere->SetPlaceFactor(1.);
        M_widgetBoxAroundSphere->SetProp3D( this->actor() );

        M_widgetBoxAroundSphere->Off();
    }


    void
    activateWidgetBox()
    {
        if ( !M_widgetBoxAroundSphere ) return;

        this->actor()->GetProperty()->SetOpacity(0.2);
        M_widgetBoxAroundSphere->PlaceWidget();
        double radiusSphere = this->geometry()->GetRadius();
        M_widgetBoxAroundSphere->SetHandleSize(radiusSphere/10000.);

        M_widgetBoxAroundSphere->On();
    }

    void
    deactivateWidgetBox()
    {
        if (!M_widgetBoxAroundSphere) return;
        double * p=M_widgetBoxAroundSphere->GetProp3D()->GetCenter();
        //std::cout << "Box center: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
        this->geometry()->SetCenter(p);

        // set identity transformation
        vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
        t->Identity();
        M_widgetBoxAroundSphere->GetProp3D()->SetUserTransform(t);

        this->actor()->GetProperty()->SetOpacity(1);

        M_widgetBoxAroundSphere->Off();
    }


    void translateBoxAroundSphere( const double translateVec[3] )
    {
        if (!M_widgetBoxAroundSphere) return;
#if 1
        vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
        M_widgetBoxAroundSphere->GetTransform(t);
        t->Translate(translateVec);
        M_widgetBoxAroundSphere->SetTransform(t);
        M_widgetBoxAroundSphere->GetProp3D()->SetUserTransform(t);
        //M_widgetBoxAroundSphere->GetProp3D()->Concatenate
#else
        double newPos[3] = { M_geometry->GetCenter()[0]+translateVec[0],M_geometry->GetCenter()[1]+translateVec[1],M_geometry->GetCenter()[2]+translateVec[2] };
        M_geometry->SetCenter( newPos/*M_widgetBoxAroundSphere->GetProp3D()->GetCenter()*/ );
        M_geometry->Update();
        //M_widgetBoxAroundSphere->GetProp3D()->SetCenter( newPos );
#endif
    }


private :
    vtkSmartPointer<vtkSphereSource> M_geometry;
    vtkSmartPointer<vtkActor> M_actor;
    int M_typePoint;
    bool M_isValidated;
    vtkSmartPointer<vtkBoxWidget> M_widgetBoxAroundSphere;
    vtkSmartPointer<vtkMyCallback> M_callbackBoxAroundSphere;

};


class PointPairConnection
{
public :
    typedef std::shared_ptr<SphereSourceObject> point_object_ptrtype;
    PointPairConnection() {}
    PointPairConnection( point_object_ptrtype pt1, point_object_ptrtype pt2, int connectionId = -1 )
        :
        M_pointObjet1( pt1 ),
        M_pointObjet2( pt2 ),
        M_connectionId( connectionId )
    {}
    PointPairConnection( PointPairConnection const& o ) = default;

    vtkSmartPointer<vtkActor> const& actor() const { return M_lineActor; }

    point_object_ptrtype
    pointObjet(int k)
    {
        if (k==1)
            return M_pointObjet1;
        else
            return M_pointObjet2;
    }

    bool
    isSameConnection( point_object_ptrtype const& pt1, point_object_ptrtype const& pt2 ) const
    {
        if ( ( pt1.get() == M_pointObjet1.get() && pt2.get() == M_pointObjet2.get() ) ||
             ( pt1.get() == M_pointObjet2.get() && pt2.get() == M_pointObjet1.get() ) )
            return true;
        else
            return false;
    }

    bool
    isConnectedToPoint( point_object_ptrtype const& pt ) const
    {
        if ( pt.get() == M_pointObjet1.get() || pt.get() == M_pointObjet2.get() )
            return true;
        else
            return false;
    }

    void
    addPlot( vtkRenderWindowInteractor *interactor )
    {
        if ( !M_lineSource )
        {
            M_lineSource = vtkSmartPointer<vtkLineSource>::New();
#if 0
            M_lineSource->SetPoint1( M_pointObjet1->geometry()->GetCenter() );
            M_lineSource->SetPoint2( M_pointObjet2->geometry()->GetCenter() );
            M_lineSource->Update();
#else
            this->updatePosition();
#endif

            //Create a mapper and actor
            vtkSmartPointer<vtkPolyDataMapper> lineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            lineMapper->SetInputConnection(M_lineSource->GetOutputPort());
            M_lineActor = vtkSmartPointer<vtkActor>::New();
            M_lineActor->SetMapper(lineMapper);
            M_lineActor->GetProperty()->SetColor(1.0, 0.5, 0.0); //(R,G,B)
            M_lineActor->GetProperty()->SetLineWidth(2.0);//defautl is 1.0

            /* Create an on-screen text indicating the id of the current connection line */
            std::ostringstream oss;
            oss << M_connectionId;
            vtkSmartPointer<vtkVectorText> textSource = vtkSmartPointer<vtkVectorText>::New();
            textSource->SetText( oss.str().c_str() );

            // Create a mapper
            vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection( textSource->GetOutputPort() );

            /* Compute the center of the line */
            double * pos1 = M_lineSource->GetPoint1();
            double * pos2 = M_lineSource->GetPoint2();

            double fpos[3];
            for(int i = 0; i < 3; i++)
            {
                fpos[i] = (pos2[i] + pos1[i]) / 2.0;
            }

            /* Add a follower object for the text: The text will always be facing the camera */
            M_follower = vtkSmartPointer<vtkFollower>::New();
            M_follower->SetMapper( mapper );
            M_follower->GetProperty()->SetColor( 1, 0, 0 );
            M_follower->SetPosition( fpos );
            M_follower->SetScale(8.0);
            M_follower->SetCamera(interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera());

            // add the follower to the scene
            interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(M_follower);

            M_callbackPointPairConnection1 = vtkSmartPointer<vtkCallbackPointPairConnection>::New();
            M_callbackPointPairConnection1->initialize(this,1);
            M_callbackPointPairConnection2 = vtkSmartPointer<vtkCallbackPointPairConnection>::New();
            M_callbackPointPairConnection2->initialize(this,2);
            M_pointObjet1->widgetBoxAroundSphere()->AddObserver(vtkCommand::InteractionEvent, M_callbackPointPairConnection1);
            M_pointObjet2->widgetBoxAroundSphere()->AddObserver(vtkCommand::InteractionEvent, M_callbackPointPairConnection2);

        }
        interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(M_lineActor);
        
    }

    void
    removePlot( vtkRenderWindowInteractor *interactor )
    {
        if ( M_lineActor )
        {
            interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(M_lineActor);
            M_pointObjet1->widgetBoxAroundSphere()->RemoveObserver(M_callbackPointPairConnection1);
            M_pointObjet2->widgetBoxAroundSphere()->RemoveObserver(M_callbackPointPairConnection2);

            interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(M_follower);
        }
    }

    void setPoint( int ptId, /*const*/ double ptCoord[3] )
    {
        if ( !M_lineSource ) return;
        if (ptId == 1 )
            M_lineSource->SetPoint1( ptCoord );
        else if ( ptId == 2 )
            M_lineSource->SetPoint2( ptCoord );
        M_lineSource->Update();
    }

    void updatePosition()
    {
        if ( !M_lineSource ) return;
        if ( M_pointObjet1->widgetBoxAroundSphere() )
            M_lineSource->SetPoint1( M_pointObjet1->widgetBoxAroundSphere()->GetProp3D()->GetCenter() );
        else
            M_lineSource->SetPoint1( M_pointObjet1->geometry()->GetCenter() );

        if ( M_pointObjet2->widgetBoxAroundSphere() )
            M_lineSource->SetPoint2( M_pointObjet2->widgetBoxAroundSphere()->GetProp3D()->GetCenter() );
        else
            M_lineSource->SetPoint2( M_pointObjet2->geometry()->GetCenter() );

        M_lineSource->Update();
    }
private :
    point_object_ptrtype M_pointObjet1, M_pointObjet2;
    vtkSmartPointer<vtkLineSource> M_lineSource;
    vtkSmartPointer<vtkActor> M_lineActor;
    vtkSmartPointer<vtkFollower> M_follower;

    vtkSmartPointer<vtkCallbackPointPairConnection> M_callbackPointPairConnection1,M_callbackPointPairConnection2;

    int M_connectionId;

};

void
vtkCallbackPointPairConnection::Execute(vtkObject *caller, unsigned long, void*)
{
    vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
    vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
    widget->GetTransform(t);

    double pt[3] = { 0.,0.,0. };
    t->TransformPoint( M_pointPairConnection->pointObjet( M_ptId )->geometry()->GetCenter(), pt );

    M_pointPairConnection->setPoint( M_ptId, pt );
}


// Define interaction style
class AngioTkWindowInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static AngioTkWindowInteractorStyle* New();
    vtkTypeMacro(AngioTkWindowInteractorStyle, vtkInteractorStyleTrackballCamera);

    AngioTkWindowInteractorStyle()
        :
        M_windowSizePreviousFullScreen( 2, 0 ), // init correctly when go in full screen (store previous size)
        M_windowFullScreen( false ),
        M_activatedMode( ModeType::NO_MODE ),
        M_sphereActorSelectionId(),
        M_lenghtSTL(0),
        M_outputPathPointSetFile("default_pointset.data"),
        M_outputPathPointPairFile("default_pointpair.data")
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
                       << "f : enable/disable full screen\n"
                       << "a : enable/disable orientation axis\n"
                       << "s : save points insertion on disk\n"
                       << "0 : no mode (view only)\n"
                       << "1 : mode surface representation\n"
                       << "2 : mode points insertion\n"
                       << "3 : mode centerlines manager\n";
        textActorCommandHelpCommon->SetInput ( helpCommandStr.str().c_str() );
        vtkSmartPointer<vtkTextRepresentation> textRepresentation = vtkSmartPointer<vtkTextRepresentation>::New();
        textRepresentation->GetPositionCoordinate()->SetValue( 0.7, 0.75 );
        textRepresentation->GetPosition2Coordinate()->SetValue( 0.3, 0.25 );
        M_widgetTextCommandHelp->SetRepresentation( textRepresentation );
        textRepresentation->SetTextActor( textActorCommandHelpCommon );
        //textRepresentation->SetPosition(0.7,0.5);
        //textRepresentation->SetPosition2(0.3,0.8);
        //textRepresentation->SetProportionalResize(1);
        textActorCommandHelpCommon->GetTextProperty()->SetColor ( 0.0,0.0,0.0 );
        textActorCommandHelpCommon->GetTextProperty()->SetJustificationToLeft();
        textActorCommandHelpCommon->GetTextProperty()->SetVerticalJustificationToTop();
        textActorCommandHelpCommon->GetTextProperty()->SetFontSize( 24 );

        M_widgetTextActivatedMode = vtkSmartPointer<vtkTextWidget>::New();
        M_widgetTextActivatedMode->SelectableOff();
        vtkSmartPointer<vtkTextRepresentation> textRepresentation2 = vtkSmartPointer<vtkTextRepresentation>::New();
        textRepresentation2->GetPositionCoordinate()->SetValue( 0.01, 0.7 );
        textRepresentation2->GetPosition2Coordinate()->SetValue( 0.4, 0.3 );
        M_widgetTextActivatedMode->SetRepresentation( textRepresentation2 );
        vtkSmartPointer<vtkTextActor> textActor2 = vtkSmartPointer<vtkTextActor>::New();
        textRepresentation2->SetTextActor( textActor2 );
        //textRepresentation2->SetPosition(0.01,0.2);
        //textRepresentation2->SetPosition2(0.3,0.8);
        //textRepresentation2->SetPosition(0.01,0.95);
        //textRepresentation2->SetPosition2(0.3,0.05);
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
            helpStr << "c : change point type (source=red,target=blue)\n"
                    << "y : validate point\n"
                    << "r : remove selection\n"
                    //<< "u : undo insertion\n"
                    << "[Up] : move selection (axis Y +)\n"
                    << "[Down] : move selection (axis Y -)\n"
                    << "[Left] : move selection (axis X -)\n"
                    << "[Right] : move selection (axis X +)\n"
                    << "o : move selection (axis Z +)\n"
                    << "l : move selection (axis Z -)\n"
                    << "p/m : increase/decrease sphere size\n"
                    << "d : Connect/disconnect 2 selected points\n";
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
    void savePointSetOnDisk( std::string const& _pathFile )
    {
        if ( M_vectorSphereSourceObject.empty() )
            return;

        namespace fs = boost::filesystem;
        fs::path pathFile = fs::absolute( _pathFile );
        fs::path dir = fs::path(pathFile).parent_path();
        if ( !fs::exists( dir ) )
            fs::create_directories( dir );
        std::string pathFileStr = pathFile.string();
        std::cout << "savePointSetOnDisk : " << pathFileStr << "\n";

        std::ofstream fileWrited( pathFileStr, std::ios::out | std::ios::trunc);
        fileWrited.precision( 16 );
        for (int k = 0; k< M_vectorSphereSourceObject.size() ;++ k)
        {
            double radius = M_vectorSphereSourceObject[k]->geometry()->GetRadius();
            double * center = M_vectorSphereSourceObject[k]->geometry()->GetCenter();
            int typePt = M_vectorSphereSourceObject[k]->typePoint();
            fileWrited << typePt << " " << center[0] << " " << center[1] << " " << center[2] << " " << radius << "\n";
        }
        fileWrited.close();
    }
    void savePointPairOnDisk( std::string const& _pathFile )
    {
        if ( M_vectorPointPairConnection.empty() )
            return;

        namespace fs = boost::filesystem;
        fs::path pathFile = fs::absolute( _pathFile );
        fs::path dir = fs::path(pathFile).parent_path();
        if ( !fs::exists( dir ) )
            fs::create_directories( dir );
        std::string pathFileStr = pathFile.string();
        std::cout << "savePointPairOnDisk : " << pathFileStr << "\n";

        std::ofstream fileWrited( pathFileStr, std::ios::out | std::ios::trunc);
        fileWrited.precision( 16 );
        for (int k = 0; k< M_vectorPointPairConnection.size() ;++ k)
        {
            double radius1 = M_vectorPointPairConnection[k]->pointObjet(1)->geometry()->GetRadius();
            double * center1 = M_vectorPointPairConnection[k]->pointObjet(1)->geometry()->GetCenter();
            double radius2 = M_vectorPointPairConnection[k]->pointObjet(2)->geometry()->GetRadius();
            double * center2 = M_vectorPointPairConnection[k]->pointObjet(2)->geometry()->GetCenter();
            fileWrited << k << " " << center1[0] << " " << center1[1] << " " << center1[2] << " " << radius1 << "\n";
            fileWrited << k << " " << center2[0] << " " << center2[1] << " " << center2[2] << " " << radius2 << "\n";
        }
        fileWrited.close();
    }

    void loadPointSetFromFile( std::string const& pathFile )
    {
        std::cout << "LOADING" << std::endl;
        auto data = AngioTk::loadFromPointSetFile( pathFile );
        for ( auto const& point : data )
        {
            auto const& coord = std::get<0>( point );
            auto const& radius = std::get<1>( point );
            auto const& type = std::get<2>( point );
            int ptId = this->hasPoint( coord );
            if ( ptId < 0 )
                this->addPoint( coord, radius, type );
        }
        this->Interactor->GetRenderWindow()->Render();
    }

    void loadPointPairFromFile( std::string const& pathFile )
    {
        auto data = AngioTk::loadFromPointPairFile( pathFile );
        for ( auto const& pointpair : data )
        {
            auto const& ptdata1 = pointpair.first;
            auto const& ptdata2 = pointpair.second;
            auto const& coord1 = std::get<0>( ptdata1 );
            auto const& coord2 = std::get<0>( ptdata2 );
            double radius1 = std::get<1>( ptdata1 );
            double radius2 = std::get<1>( ptdata2 );
            int ptId1 = this->hasPoint( coord1 );
            if ( ptId1 < 0 )
            {
                this->addPoint( coord1, radius1 );
                ptId1 = M_vectorSphereSourceObject.size()-1;
            }
            int ptId2 = this->hasPoint( coord2 );
            if ( ptId2 < 0 )
            {
                this->addPoint( coord2, radius2 );
                ptId2 = M_vectorSphereSourceObject.size()-1;
            }

            int connectionId = this->hasPointPairConnection( M_vectorSphereSourceObject[ptId1],M_vectorSphereSourceObject[ptId2] );
            if ( connectionId < 0 )
            {
                std::shared_ptr<PointPairConnection> ptPairConnection = std::make_shared<PointPairConnection>( M_vectorSphereSourceObject[ptId1],
                                                                                                               M_vectorSphereSourceObject[ptId2], 
                                                                                                               M_vectorPointPairConnection.size());
                M_vectorPointPairConnection.push_back( ptPairConnection );
                ptPairConnection->addPlot( this->Interactor );
            }

        }
        //std::cout << "size after load " << M_vectorSphereSourceObject.size() << "\n";
        this->Interactor->GetRenderWindow()->Render();
    }

    void addPoint( std::vector<double> const& pt, double radius, int type=0 )
    {
        CHECK( pt.size() == 3 ) << "invalid point size : " << pt.size() << " ( must be 3 )";
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetRadius( radius );
        sphereSource->SetCenter( pt[0],pt[1],pt[2] );
        sphereSource->Update();
        vtkSmartPointer<vtkPolyDataMapper> mapperSphere = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapperSphere->SetInputConnection(sphereSource->GetOutputPort());
        vtkSmartPointer<vtkActor> actorSphere = vtkSmartPointer<vtkActor>::New();
        actorSphere->SetMapper(mapperSphere);
        //actorSphere->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
        M_vectorSphereSourceObject.push_back(std::make_shared<SphereSourceObject>(sphereSource,actorSphere,type,true) );
        M_vectorSphereSourceObject.back()->initializeWidgetBox( this->Interactor );
        M_vectorSphereSourceObject.back()->applyColoring();
        this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actorSphere);

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

            M_vectorSphereSourceObject.push_back( std::make_shared<SphereSourceObject>(sphereSource,actorSphere,typePt) );
            M_vectorSphereSourceObject.back()->initializeWidgetBox( this->Interactor );
            M_vectorSphereSourceObject.back()->applyColoring();
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actorSphere);
        }

    }


    int hasPoint( std::vector<double> const& pt )
    {
        CHECK( pt.size() == 3 ) << "invalid point size : " << pt.size() << " ( must be 3 )";
        return this->hasPoint( pt[0],pt[1],pt[2] );
    }

    int hasPoint( double x,double y,double z )
    {
        for ( int k = 0 ; k<M_vectorSphereSourceObject.size() ; ++k )
        {
            double * center = M_vectorSphereSourceObject[k]->geometry()->GetCenter();
            if ( ( std::abs( center[0] - x ) < 1e-9 ) &&
                 ( std::abs( center[1] - y ) < 1e-9 ) &&
                 ( std::abs( center[2] - z ) < 1e-9 ) )
                return k;
        }
        return -1;
    }

    int hasActor( vtkActor * _actor ) const
    {
        int k=-1;
        for ( int k = 0 ; k<M_vectorSphereSourceObject.size() ; ++k )
        {
            if ( M_vectorSphereSourceObject[k]->actor() == _actor )
                return k;
        }
        return k;
    }

    int hasPointPairConnection( std::shared_ptr<SphereSourceObject> const& pt1, std::shared_ptr<SphereSourceObject> const& pt2 ) const
    {
        int k=-1;
        for ( int k = 0 ; k<M_vectorPointPairConnection.size() ; ++k )
            if ( M_vectorPointPairConnection[k]->isSameConnection( pt1,pt2 ) )
                return k;
        return k;
    }
    std::set<int> hasPointPairConnection( std::shared_ptr<SphereSourceObject> const& pt ) const
    {
        std::set<int> res;
        for ( int k = 0 ; k<M_vectorPointPairConnection.size() ; ++k )
            if ( M_vectorPointPairConnection[k]->isConnectedToPoint( pt ) )
                res.insert( k );
        return res;
    }


    void setAngioTkCenterlines( boost::shared_ptr<AngioTkCenterline> obj ) { M_angioTkCenterlines = obj; }

    void setActorSTL( vtkSmartPointer<vtkActor> const& actor) { M_actorSTL=actor; }

    void setBoundsSTL(double * bounds ) { M_boundsSTL.insert( M_boundsSTL.begin(),bounds,bounds+6); }
    void setLenghtSTL(double d) { M_lenghtSTL=d; }

    void setOutputPathPointSetFile(std::string const& path) { M_outputPathPointSetFile=path; }
    void setOutputPathPointPairFile(std::string const& path) { M_outputPathPointPairFile=path; }


    void activateSphereActorSelection( int selectId )
    {
        if ( selectId >= 0 && M_sphereActorSelectionId.find(selectId) == M_sphereActorSelectionId.end() )
        {
            M_vectorSphereSourceObject[selectId]->activateWidgetBox();
            M_sphereActorSelectionId.insert( selectId );
        }
    }

    void deactivateSphereActorSelection()
    {
        for ( int selectId : M_sphereActorSelectionId )
        {
            M_vectorSphereSourceObject[selectId]->deactivateWidgetBox();
        }
        M_sphereActorSelectionId.clear();
    }

    void removeSphereActorNonValidated()
    {
        if ( !M_vectorSphereSourceObject.empty() && !M_vectorSphereSourceObject.back()->isValidated() )
        {
            if ( M_vectorSphereSourceObject.back()->geometry() != NULL )
                this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(M_vectorSphereSourceObject.back()->actor());
            M_vectorSphereSourceObject.pop_back();
        }
    }

    void updatePointPairConnection()
    {
        for ( int k=0;k<M_vectorPointPairConnection.size();++k)
            M_vectorPointPairConnection[k]->updatePosition();
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
                //std::cout << "activate sphere " << selectId << "\n";
                this->activateSphereActorSelection( selectId );
            }
            else if ( M_sphereActorSelectionId.empty() )
            {
                //typedef vtkPointPicker picker_type; // vtkPointPicker, vtkWorldPointPicker
                typedef vtkCellPicker picker_type; // vtkPointPicker, vtkWorldPointPicker
                vtkSmartPointer<picker_type> myPicker = vtkSmartPointer<picker_type>::New();
                //myPicker->SetTolerance(0.0005);
                //myPicker->SetTolerance( 1e-4*M_lenghtSTL );
                //myPicker->SetTolerance( (1./M_lenghtSTL) );
                //std::cout << "M_lenghtSTL " << M_lenghtSTL << "\n";
                //myPicker->SetTolerance( 1e-4*M_lenghtSTL );
#if 0
                if ( M_lenghtSTL > 50 )
                    myPicker->SetTolerance(1e-3);
                else
                    myPicker->SetTolerance(1e-2);
#else
                myPicker->SetTolerance(5e-3);
#endif
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

                    if ( !M_vectorSphereSourceObject.empty() && !M_vectorSphereSourceObject.back()->isValidated() )
                    {
                        if ( M_vectorSphereSourceObject.back()->geometry() != NULL )
                            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(M_vectorSphereSourceObject.back()->actor());
                        M_vectorSphereSourceObject.pop_back();
                    }

                    M_vectorSphereSourceObject.push_back( std::make_shared<SphereSourceObject>(sphereSource,actorSphere,0,false) );
                    M_vectorSphereSourceObject.back()->initializeWidgetBox( this->Interactor );

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

        if ( key == "f" )
        {
#if 0
            // NOT WORKS !
            int fullScreenState = this->Interactor->GetRenderWindow()->GetFullScreen();
            this->Interactor->GetRenderWindow()->SetFullScreen( (fullScreenState+1)%2 );
            this->Interactor->GetRenderWindow()->Render();
#else
            if ( !M_windowFullScreen )
            {
                this->Interactor->GetRenderWindow()->Render();
                int * windowsize = this->Interactor->GetRenderWindow()->GetSize();
                M_windowSizePreviousFullScreen[0] = windowsize[0];
                M_windowSizePreviousFullScreen[1] = windowsize[1];
                int * fullscreensize = this->Interactor->GetRenderWindow()->GetScreenSize();
                this->Interactor->GetRenderWindow()->SetSize(fullscreensize[0],fullscreensize[1]);
                M_windowFullScreen = true;
            }
            else
            {
                this->Interactor->GetRenderWindow()->SetSize(M_windowSizePreviousFullScreen[0],M_windowSizePreviousFullScreen[1]);
                M_windowFullScreen = false;
            }
            this->Interactor->GetRenderWindow()->Render();
            //this->Interactor->GetRenderWindow()->Modified();
            //this->Interactor->GetRenderWindow()->Render();
#endif
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
                //vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->SetPosition(0.01,0.2);
                //vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->SetPosition2(0.3,0.8);

                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->GetPositionCoordinate()->SetValue( 0.01, 0.7 );
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->GetPosition2Coordinate()->SetValue( 0.4, 0.3 );
            }
            else
            {
                //vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->SetPosition(0.01,0.95);
                //vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->SetPosition2(0.3,0.05);

                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->GetPositionCoordinate()->SetValue( 0.01, 0.7 );
                vtkTextRepresentation::SafeDownCast( M_widgetTextActivatedMode->GetRepresentation() )->GetPosition2Coordinate()->SetValue( 0.4, 0.3 );
            }

            this->Interactor->GetRenderWindow()->Render();
        }


        if ( key == "0" && M_activatedMode != ModeType::NO_MODE)
        {
            this->deactivateSphereActorSelection();
            this->removeSphereActorNonValidated();
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
                this->removeSphereActorNonValidated();
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
                this->removeSphereActorNonValidated();
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
            this->savePointSetOnDisk( M_outputPathPointSetFile );
            this->savePointPairOnDisk( M_outputPathPointPairFile );
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
              for ( int selectId : M_sphereActorSelectionId )
	          {
                  double prevRadius = M_vectorSphereSourceObject[selectId]->geometry()->GetRadius();
                  double step = M_lenghtSTL/500.;
                  M_vectorSphereSourceObject[selectId]->geometry()->SetRadius(prevRadius+step);
                  this->Interactor->GetRenderWindow()->Render();
              }
          }
          if( key == "m" )
          {
              for ( int selectId : M_sphereActorSelectionId )
              {
                  double prevRadius = M_vectorSphereSourceObject[selectId]->geometry()->GetRadius();
                  double step = M_lenghtSTL/500.;
                  M_vectorSphereSourceObject[selectId]->geometry()->SetRadius(prevRadius-step);
                  this->Interactor->GetRenderWindow()->Render();
              }
          }
          // move in x,y,z direction
          if( key == "Left" || key == "Right" || key == "Up" || key == "Down" || key == "o" || key == "l" )
          {
              //if ( M_sphereActorSelectionId.size() == 1 )
              for ( int selectId : M_sphereActorSelectionId )
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

                  M_vectorSphereSourceObject[selectId]->translateBoxAroundSphere( translateVec );
                  delete [] translateVec;

                  this->updatePointPairConnection();

                  this->Interactor->GetRenderWindow()->Render();
              }
          }

          // valid sphere
          if(key == "y")
          {
              this->deactivateSphereActorSelection();
              if ( !M_vectorSphereSourceObject.empty() && M_vectorSphereSourceObject.back()->geometry() != NULL &&
                   !M_vectorSphereSourceObject.back()->isValidated() )
              {
                  M_vectorSphereSourceObject.back()->setIsValidated(true);
                  M_vectorSphereSourceObject.back()->applyColoring();
              }
              this->Interactor->GetRenderWindow()->Render();
          }

          // change propertie of point (used with centerlines algo) : source or target point
          if(key == "c")
          {
              for ( int selectId : M_sphereActorSelectionId )
              {
                  M_vectorSphereSourceObject[selectId]->changeTypePoint();
                  M_vectorSphereSourceObject[selectId]->applyColoring();
                  this->Interactor->GetRenderWindow()->Render();
              }
          }

          if(key == "d")
          {
              if ( M_sphereActorSelectionId.size() == 2 )
              {
                  int ptId1 = *M_sphereActorSelectionId.begin();
                  int ptId2 = *M_sphereActorSelectionId.rbegin();
                  int connectionId = this->hasPointPairConnection( M_vectorSphereSourceObject[ptId1],M_vectorSphereSourceObject[ptId2] );
                  if ( connectionId < 0 )
                  {
                      std::shared_ptr<PointPairConnection> ptPairConnection = std::make_shared<PointPairConnection>( M_vectorSphereSourceObject[ptId1],
                                                                                                                     M_vectorSphereSourceObject[ptId2],
                                                                                                                     M_vectorPointPairConnection.size() );
                      M_vectorPointPairConnection.push_back( ptPairConnection );
                      ptPairConnection->addPlot( this->Interactor );
                  }
                  else
                  {
                      M_vectorPointPairConnection[connectionId]->removePlot( this->Interactor );
                      M_vectorPointPairConnection.erase( M_vectorPointPairConnection.begin()+connectionId );
                  }

                  this->Interactor->GetRenderWindow()->Render();
              }
          }



#if 0
          // undo key
          if(key == "u")
          {
              if ( !M_vectorSphereSourceObject.empty() )
              {
                  this->deactivateSphereActorSelection();
                  if ( M_vectorSphereSourceObject.back()->geometry() != NULL )
                  {
                      // remove last actor in render
                      this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(M_vectorSphereSourceObject.back()->actor());
                      // delete last actor
                      M_vectorSphereSourceObject.pop_back();
                  }
                  else if ( M_vectorSphereSourceObject.size() > 1 ) // if size==0 else keep null pointer at case 0
                  {
                      // delete null pointer
                      M_vectorSphereSourceObject.pop_back();
                      // remove last actor in render
                      this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(M_vectorSphereSourceObject.back()->actor());
                      // delete last actor
                      M_vectorSphereSourceObject.pop_back();
                  }
                  this->Interactor->GetRenderWindow()->Render();
              }
          }
#endif
          if(key == "r")
          {
              if ( !M_sphereActorSelectionId.empty() )
              {
                  std::set<int> idsRemove( M_sphereActorSelectionId.begin(),M_sphereActorSelectionId.end() );
                  this->deactivateSphereActorSelection();
                  std::vector<std::shared_ptr<SphereSourceObject> > newVectorSphereSourceObject;
                  // if ( M_vectorSphereSourceObject.back()->geometry() != NULL )

                  for (int objectId=0 ; objectId<M_vectorSphereSourceObject.size() ; ++objectId)
                  {
                      if ( idsRemove.find( objectId ) != idsRemove.end() )
                          this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor( M_vectorSphereSourceObject[objectId]->actor() );
                      else
                          newVectorSphereSourceObject.push_back( M_vectorSphereSourceObject[objectId] );
                  }

                  std::set<int> idsConnectionToRemove;
                  for (int k : idsRemove )
                  {
                      std::set<int> idsConnection = this->hasPointPairConnection( M_vectorSphereSourceObject[k] );
                      idsConnectionToRemove.insert( idsConnection.begin(), idsConnection.end() );
                  }
                  if ( !idsConnectionToRemove.empty() )
                  {
                      std::vector<std::shared_ptr<PointPairConnection> > newVectorPointPairConnection;
                      for (int objectId=0 ; objectId<M_vectorPointPairConnection.size() ; ++objectId)
                      {
                          if ( idsConnectionToRemove.find( objectId ) != idsConnectionToRemove.end() )
                              M_vectorPointPairConnection[objectId]->removePlot( this->Interactor );
                          else
                              newVectorPointPairConnection.push_back( M_vectorPointPairConnection[objectId] );
                      }
                      // update new connection objects
                      M_vectorPointPairConnection.clear();
                      M_vectorPointPairConnection = newVectorPointPairConnection;
                  }

                  // update new point objects
                  M_vectorSphereSourceObject.clear();
                  M_vectorSphereSourceObject = newVectorSphereSourceObject;

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

    std::vector<int> M_windowSizePreviousFullScreen;
    bool M_windowFullScreen;

    ModeType M_activatedMode;

    vtkActor * M_LastPickedActor;
    vtkProperty *M_LastPickedProperty;
    std::set<int> M_sphereActorSelectionId;

    std::vector<double> M_boundsSTL; double M_lenghtSTL;

    std::vector<std::shared_ptr<SphereSourceObject> > M_vectorSphereSourceObject;
    std::vector<std::shared_ptr<PointPairConnection> > M_vectorPointPairConnection;

    vtkSmartPointer<vtkTextWidget> M_widgetTextCommandHelp;
    vtkSmartPointer<vtkTextWidget> M_widgetTextActivatedMode;
    vtkSmartPointer<vtkOrientationMarkerWidget> M_widgetOrientationAxis;

    vtkSmartPointer<vtkActor> M_actorSTL;

    boost::shared_ptr<AngioTkCenterline> M_angioTkCenterlines;

    std::string M_outputPathPointSetFile, M_outputPathPointPairFile;

};
vtkStandardNewMacro(AngioTkWindowInteractorStyle);

CenterlinesManagerWindowInteractor::CenterlinesManagerWindowInteractor()
    :
    M_windowWidth( 1024 ),
    M_windowHeight (768)
{}

void
CenterlinesManagerWindowInteractor::run()//bool fullscreen, int windowWidth, int windowHeight)
{
    typedef vtkPointPicker picker_type; // vtkPointPicker, vtkWorldPointPicker
    //vtkSmartPointer<picker_type> worldPointPicker = vtkSmartPointer<picker_type>::New();
    //worldPointPicker->SetTolerance(0.0005);
    //--------------------------------------------
    // Create a renderer, render window, and interactor
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetSize(this->windowWidth(), this->windowHeight());

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
    if ( !this->inputSurfacePath().empty() )
    {
        std::string inputFilename = this->inputSurfacePath();
        vtkSmartPointer<vtkSTLReader> readerSTL = vtkSmartPointer<vtkSTLReader>::New();
        readerSTL->SetFileName(inputFilename.c_str());
        readerSTL->Update();

        vtkSmartPointer<vtkPolyData> polyData = readerSTL->GetOutput();
        std::cout << "Read STL file " << this->inputSurfacePath() << std::endl;
        std::cout << "- Number of Points: " << polyData->GetNumberOfPoints() << std::endl;
        std::cout << "- Number of Cells: " << polyData->GetNumberOfCells() << std::endl;

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

        std::string namePointSetFile = Feel::fs::path(this->inputSurfacePath()).stem().string()+"_pointset.data";
        std::string namePointPairFile = Feel::fs::path(this->inputSurfacePath()).stem().string()+"_pointpair.data";
        style->setOutputPathPointSetFile(namePointSetFile);
        style->setOutputPathPointPairFile(namePointPairFile);
    }
    else
    {
        std::cout << "The program won't load surface data as it is empty. Specify it with --input.surface.path" << std::endl;
    }

    boost::shared_ptr<AngioTkCenterline> centerlinesTool;
    for ( int k=0;k<this->inputCenterlinesPath().size();++k)
    {
        if ( this->inputCenterlinesPath(k).empty() || !Feel::fs::exists( this->inputCenterlinesPath(k) ) ) continue;

        if ( true )
            {
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
            }
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

#if 1
    if ( false )
        {
            //centerlinesTool->addFieldBranchIds();
            centerlinesTool->writeCenterlinesVTK( "toto.vtk" );
            vtkSmartPointer<vtkPolyDataReader> readerVTK = vtkSmartPointer<vtkPolyDataReader>::New();
            readerVTK->SetFileName("toto.vtk");
            readerVTK->Update();
            // Create a mapper and actor
            vtkSmartPointer<vtkPolyDataMapper> mapperVTK = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapperVTK->SetInputConnection(readerVTK->GetOutputPort());
            vtkSmartPointer<vtkActor> actorVTK = vtkSmartPointer<vtkActor>::New();
            actorVTK->SetMapper(mapperVTK);
            actorVTK->GetProperty()->SetLineWidth(3.0);
            // Add the actor to the scene
            renderer->AddActor(actorVTK);
        }
#endif

    if ( !this->inputPointSetPath().empty() )
    {
        std::string previousData = this->inputPointSetPath();
        style->loadPointSetFromFile(previousData);
    }
    if ( !this->inputPointPairPath().empty() )
        style->loadPointPairFromFile( this->inputPointPairPath() );

    renderer->ResetCamera();
    //ren->GetActiveCamera()->Zoom(1.5);

    // Render and interact
    renderWindow->Render();
    renderWindowInteractor->Start();

    //return EXIT_SUCCESS;

}

