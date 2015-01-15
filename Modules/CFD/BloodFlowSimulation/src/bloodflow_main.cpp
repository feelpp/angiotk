/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feel.hpp>


namespace Feel
{

enum BCType { BC_INLET=0,BC_OUTLET=1 };

class InletOutletData : public boost::tuple<BCType,std::string,std::string,std::vector<double> >
{
    typedef boost::tuple<BCType,std::string,std::string,std::vector<double> > super_type;
public :
    InletOutletData( BCType bctype,std::string markerLumen,std::string markerArterialWall,double x, double y, double z)
        :
        super_type(bctype,markerLumen,markerArterialWall,{x,y,z} )
    {}
    InletOutletData( InletOutletData const& e )
        :
        super_type( e )
    {}
public :
    BCType bcType() const { return this->get<0>(); }
    std::string markerLumen() const { return this->get<1>(); }
    std::string markerArterialWall() const { return this->get<2>(); }
    std::vector<double> const& node() const { return this->get<3>(); }
    double nodeX() const { return this->node()[0]; }
    double nodeY() const { return this->node()[1]; }
    double nodeZ() const { return this->node()[2]; }

};
class InletOutletDesc : public std::vector<InletOutletData>
{
    typedef std::vector<InletOutletData> super_type;
public :

    InletOutletDesc( std::string const& path )
    {
        std::ifstream fileloaded(path.c_str(), std::ios::in);  // load file
        if( fileloaded ) // if open sucess
        {
            std::string bctype,markerLumen,markerArterialWall;
            double ptx,pty,ptz;
            while ( !fileloaded.eof() )
            {
                //fileloaded >> markerLumen >> markerArterialWall >> ptx >> pty >> ptz;
                fileloaded >> bctype;
                if ( fileloaded.eof() ) break;
                fileloaded >> markerLumen >> markerArterialWall >> ptx >> pty >> ptz;
                //if( fileloaded.eof() break;
                CHECK( bctype == "INLET" || bctype == "OUTLET" ) << "invalid type " << bctype;
                BCType mybctype = ( bctype == "INLET" )? BCType::BC_INLET : BCType::BC_OUTLET;
                //std::cout << "add " << bctype  << " " << markerLumen << " " << markerArterialWall<<"\n";
                this->add( InletOutletData( mybctype,markerLumen,markerArterialWall,ptx,pty,ptz ) );
            }
            fileloaded.close();
        }
    }

    InletOutletDesc( InletOutletDesc const& e )
        :
        super_type( e ),
        M_path( e.M_path )
    {}

    void add( InletOutletData const& data )
    {
        this->push_back( data );
    }

 private :
    std::string M_path;

};

} // namespace Feel

int main(int argc, char**argv )
{
    //! [marker1]
    using namespace Feel;

    po::options_description stokesoptions( "Stokes options" );
    stokesoptions.add_options()
        ( "lumen.filename", po::value<std::string>()->default_value( "" ), "lumen filename" )
        ( "desc.filename", po::value<std::string>()->default_value( "" ), "desc.filename" )
        ( "blood-model", po::value<std::string>()->default_value( "stokes" ), "blood model" )
        ( "rho", po::value<double>()->default_value( 1.0 ), "density" )
        ( "mu", po::value<double>()->default_value( 1.0 ), "viscosity" )
        ( "inlet.pressure", po::value<std::string>()->default_value( "1" ), "inlet pressure" )
        ( "output.directory", Feel::po::value<std::string>()->default_value(""), "(string) output directory")
        ;
    Environment env( _argc=argc, _argv=argv,
                     _desc=stokesoptions,
                     _about=about(_name="cfd_bloodflow",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    double rho = doption(_name="rho");
    double mu = doption(_name="mu");
    std::string inputLumenMesh = soption(_name="lumen.filename");
    CHECK( fs::exists( inputLumenMesh ) ) << "inputLumenMesh not exist " << inputLumenMesh;
    auto inletPressureExpr = expr( soption(_name="inlet.pressure") );

    std::list<std::string> markersInlet;//={ "inlet0" };
    CHECK( fs::exists(soption(_name="desc.filename") ) ) << "desc.filename not exist";
    InletOutletDesc iodesc( soption(_name="desc.filename") );
    for ( auto const& data : iodesc )
    {
        if ( data.bcType() == BCType::BC_INLET )
            markersInlet.push_back( data.markerLumen() );
    }

    /*if (Environment::isMasterRank() )
        for ( std::string mark : markersInlet )
        std::cout << "mark " << mark << "\n";*/

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>,
                         _filename=inputLumenMesh );

    auto Vh = THch<1>( mesh );
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto q = V.element<1>();

    //auto deft = gradt( u );
    //auto def = grad( v );
    auto deft = sym(gradt(u));
    auto def = sym(grad(v));

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate( _range=elements( mesh ), _expr=mu*inner( deft,def ) );
    a +=integrate( _range=elements( mesh ), _expr=-div( v )*idt( p ) + divt( u )*id( q ) );

    auto l = form1( _test=Vh );
    l +=integrate( _range=markedfaces( mesh, markersInlet ),
                   _expr=-inletPressureExpr*inner( vf::N(),id(v) ) );

    a+=on(_range=markedfaces(mesh,"wall"), _rhs=l, _element=u,
          _expr=zero<3,1>() ) ;

    a.solve(_rhs=l,_solution=U);


    std::string outputDirectory = soption(_name="output.directory");
    if ( outputDirectory.empty() )
        outputDirectory=fs::current_path().string();
    else
        outputDirectory = (fs::path(Environment::rootRepository())/fs::path(outputDirectory)).string();
    auto e = exporter( _mesh=mesh,
                       _name="Export",
                       _path= outputDirectory );
    e->add( "u", u );
    e->add( "p", p );
    e->save();

    return 0;
}
