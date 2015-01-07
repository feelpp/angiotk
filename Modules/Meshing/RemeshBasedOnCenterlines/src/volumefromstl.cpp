/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <volumefromstl.hpp>


namespace Feel
{

CenterlinesFromSTL::CenterlinesFromSTL( std::string prefix )
    :
    M_prefix( prefix ),
    M_inputPath( soption(_name="input.filename",_prefix=this->prefix()) ),
    M_outputDirectory( soption(_name="output.directory",_prefix=this->prefix()) ),
    M_forceRebuild( boption(_name="force-rebuild",_prefix=this->prefix() ) )
{
    std::vector<int> sids,tids;
    if ( Environment::vm().count(prefixvm(this->prefix(),"source-ids").c_str()) )
        sids = Environment::vm()[prefixvm(this->prefix(),"source-ids").c_str()].as<std::vector<int> >();
    if ( Environment::vm().count(prefixvm(this->prefix(),"target-ids").c_str() ) )
        tids = Environment::vm()[prefixvm(this->prefix(),"target-ids").c_str()].as<std::vector<int> >();
    for ( int id : sids )
        M_sourceids.insert( id );
    for ( int id : tids )
        M_targetids.insert( id );

    if ( !M_inputPath.empty() && M_outputPath.empty() )
    {
        this->updateOutputPathFromInputFileName();
    }
}

CenterlinesFromSTL::CenterlinesFromSTL( CenterlinesFromSTL const& e )
    :
    M_prefix( e.M_prefix ),
    M_inputPath( e.M_inputPath ),
    M_outputPath( e.M_outputPath ),
    M_outputDirectory( e.M_outputDirectory ),
    M_targetids( e.M_targetids ),
    M_sourceids( e.M_sourceids ),
    M_forceRebuild( e.M_forceRebuild )
{}


} // namespace Feel
