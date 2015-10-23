/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#ifndef __CENTERLINESMANAGERIODATA_H
#define __CENTERLINESMANAGERIODATA_H 1


#include <utility>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>
#include <feel/feelcore/feel.hpp>


namespace AngioTk
{
    // vector of ( coord,radius,type )
    using pointset_data_type = std::vector< std::tuple< std::vector<double>,double,int > >;

    pointset_data_type
    loadFromPointSetFile( std::string const& pathFile );

    using pointpair_data_type = std::vector< std::pair< std::tuple< std::vector<double>, double >, std::tuple< std::vector<double>, double >  > > ;

    pointpair_data_type
    loadFromPointPairFile( std::string const& pathFile );

} // namespace AngioTk

#endif // __CENTERLINESMANAGERIODATA_H
