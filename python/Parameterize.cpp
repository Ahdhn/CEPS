#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/SurfaceProjectiveParameterizationQuantity.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "ConeFlattening.h"
#include "IO.h"
#include "Logger.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace CEPS;

//== Input data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::vector<size_t> imaginaryFaceIndices;
std::vector<std::pair<Vertex, double>> prescribedCurvatures;
std::vector<std::pair<Vertex, double>> prescribedScaleFactors;
double uTol  = 5; // maximum allowed area distortion in greedy cone placement
bool verbose = false;

//== Result data
ParameterizationResult result;

//== Visualization data
polyscope::SurfaceMesh* psMesh;

const int IM_STR_LEN = 128;
static char meshSaveName[IM_STR_LEN];
static char indexMapSaveName[IM_STR_LEN];

enum class InterpolationType { PROJECTIVE, LINEAR };
int currentInterpolationType                       = 0;
const char* prettyInterpolationTypeOptions[]       = {"Homogeneous", "Linear"};
const InterpolationType interpolationTypeOptions[] = {
    InterpolationType::PROJECTIVE, InterpolationType::LINEAR};

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    if (ImGui::Button("Place greedy cones")) {
        prescribedScaleFactors.clear();
        prescribedCurvatures.clear();

        std::vector<Vertex> coneVertices =
            placeCones(*mesh, *geometry, uTol, 4, 4, verbose);

        // list of cone indices, with a dummy second variable since polyscope
        // wants to render a function on the cones
        std::vector<std::pair<size_t, int>> coneIndices;
        for (Vertex v : coneVertices) {
            coneIndices.emplace_back(v.getIndex(), 0);
            prescribedScaleFactors.emplace_back(v, 0);
        }

        // Minimal area distortion boundary conditions.
        for (BoundaryLoop b : mesh->boundaryLoops()) {
            for (Vertex v : b.adjacentVertices()) {
                prescribedScaleFactors.emplace_back(v, 0);
            }
        }

        psMesh->addVertexCountQuantity("cones", coneIndices);
    }

    if (ImGui::Button("Uniformize with greedy cones")) {
        bool viz              = true;
        bool checkInjectivity = true;
        result = parameterizeWithGreedyCones(*mesh, *geometry, viz,
                                             checkInjectivity, uTol, verbose);
        psMesh->setEnabled(false);

        std::cout << "nInvertedTriangles: " << result.nFlippedTriangles
                  << "\tnZeroAreaTriangles: " << result.nZeroAreaTriangles
                  << std::endl;

        psMesh->addVertexScalarQuantity("parentMap", result.parentMap);
    }

    if (ImGui::Button("Uniformize")) {
        bool viz              = true;
        bool checkInjectivity = true;
        result = parameterize(*mesh, *geometry, prescribedScaleFactors,
                              prescribedCurvatures, imaginaryFaceIndices, viz,
                              checkInjectivity, verbose);
        psMesh->setEnabled(false);

        std::cout << "nInvertedTriangles: " << result.nFlippedTriangles
                  << "\tnZeroAreaTriangles: " << result.nZeroAreaTriangles
                  << std::endl;
    }

    ImGui::Separator();
    ImGui::InputText("###PtexturedMeshSaveName", meshSaveName,
                     IM_ARRAYSIZE(meshSaveName));
    ImGui::SameLine();
    if (ImGui::Button("Save Texture")) {
        switch (interpolationTypeOptions[currentInterpolationType]) {
        case InterpolationType::PROJECTIVE:
            writeMeshWithProjectiveTextureCoords(
                result.mesh, result.param, std::string(meshSaveName) + ".obj");
            break;
        case InterpolationType::LINEAR:
            writeMeshWithOrdinaryTextureCoords(
                result.mesh, result.param, std::string(meshSaveName) + ".obj");
            break;
        }
        saveMatrix(result.interpolationMatrix,
                   std::string(meshSaveName) + ".spmat");
    }
    ImGui::Combo("Saved Texture Type", &currentInterpolationType,
                 prettyInterpolationTypeOptions,
                 IM_ARRAYSIZE(interpolationTypeOptions));
    ImGui::Separator();
    ImGui::InputText("###IndexMapSaveName", indexMapSaveName,
                     IM_ARRAYSIZE(indexMapSaveName));
    ImGui::SameLine();
    if (ImGui::Button("Save Vertex Map")) {
        writeVertexMap(*mesh, result.parentMap,
                       std::string(indexMapSaveName) + ".txt");
    }
}

void parameterize(std::string meshFilename, std::string curvaturesFilename,
                  std::string scaleFactorsFilename, std::string ffieldFilename,
                  double greedyConeMaxU, std::string outputMeshFilename,
                  std::string outputLinearTextureFilename,
                  std::string outputMatrixFilename,
                  std::string outputVertexMapFilename,
                  std::string outputLogFilename, bool useExactCones,
                  bool noFreeBoundary, bool viz, bool verbose, bool version,
                  bool help) {

    if (version) {
        std::cout << "parameterize version 1.2" << std::endl;
        return;
    }

    if (meshFilename.empty()) {
        std::cout << "Please provide a mesh file as input." << std::endl;
        return;
    }

    // Load mesh
    std::tie(mesh, geometry) = loadMesh(meshFilename);
    verbose_assert(mesh->nConnectedComponents() == 1, "mesh must be connected");

    std::string nicename = polyscope::guessNiceNameFromPath(meshFilename);

    std::string saveNameGuess  = nicename + "_ceps";
    std::string indexNameGuess = nicename + "_ceps_vtx_map";
    // Initialize ImGui's C-string to our niceName
    // https://stackoverflow.com/a/347959
    // truncate saveNameGuess to fit in ImGui's string
    if (saveNameGuess.length() + 1 > IM_STR_LEN)
        saveNameGuess.resize(IM_STR_LEN - 1);
    if (indexNameGuess.length() + 1 > IM_STR_LEN)
        indexNameGuess.resize(IM_STR_LEN - 1);
    // copy over string contents
    std::copy(saveNameGuess.begin(), saveNameGuess.end(), meshSaveName);
    std::copy(indexNameGuess.begin(), indexNameGuess.end(), indexMapSaveName);
    // null-terminate string
    meshSaveName[saveNameGuess.size()]      = '\0';
    indexMapSaveName[indexNameGuess.size()] = '\0';

    bool loadedCones = false;
    // Read prescribed curvatures and scale factors if present
    std::vector<std::pair<size_t, double>> prescribedCurvatureInput,
        prescribedScaleFactorInput, ffieldConeInput;

    if (!curvaturesFilename.empty()) {
        loadedCones              = true;
        prescribedCurvatureInput = readPrescribedConeAngles(curvaturesFilename);
        for (const std::pair<size_t, double>& c : prescribedCurvatureInput)
            prescribedCurvatures.emplace_back(mesh->vertex(c.first), c.second);
    }
    if (!scaleFactorsFilename.empty()) {
        loadedCones = true;
        prescribedScaleFactorInput =
            readPrescribedScaleFactors(scaleFactorsFilename);
        for (const std::pair<size_t, double>& c : prescribedScaleFactorInput)
            prescribedScaleFactors.emplace_back(mesh->vertex(c.first),
                                                c.second);
    }
    if (!ffieldFilename.empty()) {
        loadedCones = true;
        std::vector<std::pair<Vertex, double>> ffieldCones =
            readFFieldCones(ffieldFilename, *mesh, *geometry);

        auto processedCones =
            useExactCones ? ffieldCones : lumpCones(*mesh, ffieldCones);

        for (const std::pair<Vertex, double>& c : processedCones)
            ffieldConeInput.emplace_back(c.first.getIndex(), c.second);

        prescribedCurvatures.insert(std::end(prescribedCurvatures),
                                    std::begin(processedCones),
                                    std::end(processedCones));
    }

    // Set u=0 at any remaining boundary vertices unless instructed not to
    if (!noFreeBoundary) {
        // Identify which vertices are already fixed
        VertexData<bool> prescribed(*mesh);
        for (const std::pair<Vertex, double> cone : prescribedCurvatures)
            prescribed[cone.first] = true;
        for (const std::pair<Vertex, double> cone : prescribedScaleFactors)
            prescribed[cone.first] = true;

        // Set u=0 at all other boundary vertices
        for (BoundaryLoop b : mesh->boundaryLoops()) {
            for (Vertex v : b.adjacentVertices()) {
                if (!prescribed[v]) prescribedScaleFactors.emplace_back(v, 0);
            }
        }
    }


    uTol = greedyConeMaxU;


    if (viz) {
        // Initialize polyscope
        polyscope::init();

        // Set the callback function
        polyscope::state::userCallback = myCallback;

        // Register the mesh with polyscope
        psMesh = polyscope::registerSurfaceMesh(
            "input_mesh", geometry->inputVertexPositions,
            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        if (viz) {
            if (!prescribedCurvatureInput.empty()) {
                auto q = psMesh->addVertexIsolatedScalarQuantity(
                    "prescribed curvatures", prescribedCurvatureInput);
                symmetrizeVizRange(*q);
                q->setEnabled(true);
            }

            if (!prescribedScaleFactorInput.empty()) {
                auto q = psMesh->addVertexIsolatedScalarQuantity(
                    "prescribed scaleFactors", prescribedScaleFactorInput);
                symmetrizeVizRange(*q);
                q->setEnabled(true);
            }

            if (!ffieldConeInput.empty()) {
                auto q = psMesh->addVertexIsolatedScalarQuantity(
                    "ffield cones", ffieldConeInput);
                symmetrizeVizRange(*q);
                q->setEnabled(true);
            }
        }

        std::cout << "Loaded mesh " << meshFilename << std::endl;
        std::cout << "nBoundaryLoops: " << mesh->nBoundaryLoops() << std::endl;
        std::cout << "Genus: " << mesh->genus() << std::endl;
        std::cout << "Euler characteristic: " << mesh->eulerCharacteristic()
                  << std::endl;

        // Give control to the polyscope gui
        polyscope::show();
    } else {
        if (!outputLogFilename.empty()) {
            std::ofstream out;

            // std::ios::trunc ensures that we overwrite old versions
            out.open(outputLogFilename, std::ios::trunc);
            if (out.is_open()) {
                // output a temporary symbol so we can tell if the program
                // crashes before writing the real log
                out << ":'(" << std::endl;
                out.close();
            } else {
                std::cout << "Error: failed to write to " << outputLogFilename
                          << vendl;
            }
        }

        Logger logger;
        bool logStats = false;
        if (!outputLogFilename.empty()) logStats = true;

        if (logStats) {
            logger.log("name", nicename);
            logger.log("nVertices", mesh->nVertices());
        }

        double duration;
        if (loadedCones) {
            std::clock_t start = std::clock();

            bool viz              = false;
            bool checkInjectivity = logStats;
            result = parameterize(*mesh, *geometry, prescribedScaleFactors,
                                  prescribedCurvatures, imaginaryFaceIndices,
                                  viz, checkInjectivity, verbose);

            duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        } else {
            std::clock_t start = std::clock();

            bool viz                  = false;
            bool checkInjectivity     = logStats;
            std::vector<Vertex> cones = placeCones(*mesh, *geometry, uTol);
            result = parameterizeWithGivenCones(*mesh, *geometry, cones, viz,
                                                checkInjectivity, verbose);

            duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        }

        if (logStats) {
            logger.log("duration", duration);
            logger.log("nFlippedTriangles", result.nFlippedTriangles);
            logger.log("nZeroAreaTriangles", result.nZeroAreaTriangles);
        }

        if (!outputMeshFilename.empty()) {
            writeMeshWithProjectiveTextureCoords(result.mesh, result.param,
                                                 outputMeshFilename);
        }

        if (!outputLinearTextureFilename.empty()) {
            writeMeshWithOrdinaryTextureCoords(result.mesh, result.param,
                                               outputLinearTextureFilename);
        }

        if (!outputVertexMapFilename.empty()) {
            writeVertexMap(*mesh, result.parentMap, outputVertexMapFilename);
        }

        if (!outputMatrixFilename.empty()) {
            saveMatrix(result.interpolationMatrix, outputMatrixFilename);
        }

        if (!outputLogFilename.empty()) {
            logger.writeLog(outputLogFilename);
        }
    }
}