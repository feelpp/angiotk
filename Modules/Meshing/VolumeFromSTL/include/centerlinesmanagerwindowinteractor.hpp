
#include <string>
#include <vector>


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

