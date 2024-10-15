#include <pybind11/pybind11.h>

void parameterize(std::string meshFilename, std::string curvaturesFilename,
                  std::string scaleFactorsFilename, std::string ffieldFilename,
                  double greedyConeMaxU, std::string outputMeshFilename,
                  std::string outputLinearTextureFilename,
                  std::string outputMatrixFilename,
                  std::string outputVertexMapFilename,
                  std::string outputLogFilename, bool useExactCones,
                  bool noFreeBoundary, bool viz, bool verbose, bool version,
                  bool help);

PYBIND11_MODULE(PyCEPS, m) {
    m.def("parameterize", &parameterize, "CEPS Parameterization ");
}