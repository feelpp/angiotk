/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <centerlinesmanageriodata.hpp>

namespace AngioTk
{

// vector of ( coord,radius,type )
pointset_data_type
loadFromPointSetFile( std::string const& pathFile )
{
    namespace fs = boost::filesystem;

    std::vector< std::tuple< std::vector<double>,double,int > > res;
    if ( !fs::exists( pathFile ) )
        return res;

    std::ifstream fileLoaded( pathFile, std::ios::in);
    while ( !fileLoaded.eof() )
      {
          std::vector<double> pt(3);
          double radius;
          int typePt = -1;
          fileLoaded >> typePt;
          if ( fileLoaded.eof() ) break;
          fileLoaded >> pt[0] >> pt[1] >> pt[2] >> radius;

          res.push_back(std::make_tuple(pt,radius,typePt));
      }
    fileLoaded.close();

    return res;
}

pointpair_data_type
loadFromPointPairFile( std::string const& pathFile )
{
  namespace fs = boost::filesystem;

  pointpair_data_type res;
  if ( !fs::exists( pathFile ) )
    return res;

    std::ifstream fileLoaded( pathFile, std::ios::in);
    while ( !fileLoaded.eof() )
    {
        std::vector<double> pt1(3);
        std::vector<double> pt2(3);
        double radius1=0,radius2=0;
        int typePt1 = -1,typePt2 = -1;
        fileLoaded >> typePt1;
        if ( fileLoaded.eof() ) break;
        fileLoaded >> pt1[0] >> pt1[1] >> pt1[2] >> radius1;

        fileLoaded >> typePt2;
        if ( fileLoaded.eof() ) break;
        fileLoaded >> pt2[0] >> pt2[1] >> pt2[2] >> radius2;

        res.push_back( std::make_pair(std::make_tuple(pt1,radius1),std::make_tuple(pt2,radius2)) );
    }
    fileLoaded.close();
    return res;
}

} // namespace AngioTk
