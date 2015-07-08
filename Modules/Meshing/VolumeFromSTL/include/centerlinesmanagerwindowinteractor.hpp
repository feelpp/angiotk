
#include <string>


class CenterlinesManagerWindowInteractor
{
public :
  std::string const& inputSurfacePath() const { return M_inputSurfacePath; }
  std::string const& inputPointSetPath() const { return M_inputPointSetPath; }
  std::string const& inputCenterlinesPath() const { return M_inputCenterlinesPath; }
  void setInputSurfacePath(std::string const& path) { M_inputSurfacePath=path; }
  void setInputPointSetPath(std::string const& path) { M_inputPointSetPath=path; }
  void setInputCenterlinesPath(std::string const& path) { M_inputCenterlinesPath=path; }

  void run();

private :
  std::string M_inputSurfacePath, M_inputPointSetPath, M_inputCenterlinesPath;

};

