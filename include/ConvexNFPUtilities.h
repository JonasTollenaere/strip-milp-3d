//
// Created by Jonas Tollenaere on 17/03/2025.
//

#ifndef EXTENDEDMESHCORE_CONVEXNFPUTILITIES_H
#define EXTENDEDMESHCORE_CONVEXNFPUTILITIES_H

#include <meshcore/optimization/StripPackingSolution.h>
#include <meshcore/core/Plane.h>
#include <meshcore/geometric/Intersection.h>

#include "glm/gtx/component_wise.hpp"

class ConvexNFPUtilities {
public:
    struct Slice {
        std::vector<Plane> boundingPlanes;
    };

    static std::shared_ptr<ModelSpaceMesh> meshFromPlanes(const std::vector<Plane>& planes);

    static std::shared_ptr<ModelSpaceMesh> computeSliceMesh(const std::shared_ptr<const StripPackingSolution>& solution, int a, int b, const Slice& slice, int sliceIndex);

    static std::shared_ptr<ModelSpaceMesh> computeNFP(const std::shared_ptr<ModelSpaceMesh>& meshA, const std::shared_ptr<ModelSpaceMesh>& meshB);

    template<bool SYMMETRY_BREAKING=true, bool USE_DISJUNCTIVE_SLICES=true, bool FILTER_EMPTY_SLICES=true>
    static std::vector<Slice> computeSlices(const std::shared_ptr<const StripPackingSolution>& solution, const std::shared_ptr<ModelSpaceMesh>& meshA, const std::shared_ptr<ModelSpaceMesh>& meshB, size_t a, size_t b) {
        const auto& NFP = computeNFP(meshA, meshB);
        const auto& NFPVertices = NFP->getVertices();

        // One slice for every face of the NFP
        std::vector<Slice> slices;
        for (const auto &face: NFP->getFaces()) {

            // Newell's Method to calculate the face normal
            const auto& faceIndices = face.vertexIndices;
            glm::vec3 faceNormal(0.0f);
            for (int fi = 0; fi < faceIndices.size(); ++fi) {
                auto indexA = faceIndices[fi];
                auto indexB = faceIndices[(fi + 1) % faceIndices.size()];

                Vertex vertexA = NFPVertices[indexA];
                Vertex vertexB = NFPVertices[indexB];

                faceNormal.x += (vertexA.y - vertexB.y) * (vertexA.z + vertexB.z);
                faceNormal.y += (vertexA.z - vertexB.z) * (vertexA.x + vertexB.x);
                faceNormal.z += (vertexA.x - vertexB.x) * (vertexA.y + vertexB.y);
            }

            auto n = glm::normalize(faceNormal);
            auto d = -glm::dot(n, NFPVertices[faceIndices[0]]);

            // The face itself is the first plane that bounds the slice
            Slice newSlice;
            newSlice.boundingPlanes.emplace_back(n, d);

            // Symmetry breaking: ordering around the z-axis
            if(SYMMETRY_BREAKING && meshA == meshB){
                newSlice.boundingPlanes.emplace_back(glm::vec3(0, 0, 1), 0);
            }

            slices.emplace_back(newSlice);
        }

        // If disjunct slices are desired, generate bounding planes between neighbouring faces
        if(USE_DISJUNCTIVE_SLICES) {
            for (int faceIndex = 0; faceIndex < NFP->getFaces().size(); ++faceIndex) {
                const auto& face = NFP->getFaces()[faceIndex];

                for (int otherFaceIndex = 0; otherFaceIndex < NFP->getFaces().size(); ++otherFaceIndex){

                    if(faceIndex == otherFaceIndex) continue;

                    const auto& otherFace = NFP->getFaces()[otherFaceIndex];

                    // Find edges shared by the faces and add bounding planes for them
                    for (int fi = 0; fi < face.vertexIndices.size(); ++fi){

                        // fi-th edge of the face
                        size_t edgeIndexA = face.vertexIndices[fi];
                        size_t edgeIndexB = face.vertexIndices[(fi + 1) % face.vertexIndices.size()];

                        for(int ofi = 0; ofi < otherFace.vertexIndices.size(); ++ofi){

                            size_t otherEdgeIndexA = otherFace.vertexIndices[ofi];
                            size_t otherEdgeIndexB = otherFace.vertexIndices[(ofi + 1) % otherFace.vertexIndices.size()];

                            // Opposite winding
                            if(edgeIndexA == otherEdgeIndexB && edgeIndexB == otherEdgeIndexA){

                                // The faces share an edge
                                // Triangles don't share the same facet, define an additional cut for Ta - Tb
                                auto nA = glm::normalize(slices[faceIndex].boundingPlanes[0].getNormal());
                                auto nB = glm::normalize(slices[otherFaceIndex].boundingPlanes[0].getNormal());

                                // Normalization of tangent and cotangent is required for robustness
                                auto tangent = glm::normalize(nA + nB); // Tangent to the plane that cuts for Ta - Tb
                                // Just use the edge itself as cotangent: much more robust compared to the cross product of the normals if they are nearly parallel
                                auto cotangent = glm::normalize(NFP->getVertices()[edgeIndexA] - NFP->getVertices()[edgeIndexB]);

                                auto normal = glm::normalize(glm::cross(tangent, cotangent)); // Normal to the plane that cuts for Ta - Tb

                                auto orientedNormalA = glm::dot(normal, nA) > 0 ? normal : -normal; // Should not be oriented away from nA
                                auto dA = -glm::dot(orientedNormalA, NFP->getVertices()[edgeIndexA]);
                                auto& sliceA = slices[faceIndex];
                                sliceA.boundingPlanes.emplace_back(orientedNormalA, dA);

                                auto orientedNormalB = glm::dot(normal, nB) > 0 ? normal : -normal; // Should not be oriented away from nB
                                auto dB = -glm::dot(orientedNormalB, NFP->getVertices()[edgeIndexA]);
                                auto& sliceB = slices[otherFaceIndex];
                                sliceB.boundingPlanes.emplace_back(orientedNormalB, dB);

                            }
                        }
                    }
                }
            }
        }

        // Filter out empty slices
        if(FILTER_EMPTY_SLICES){
            assert(slices.size() == NFP->getFaces().size());
            slices.erase(std::remove_if(slices.begin(), slices.end(), [&](const Slice& slice){
                // When no vertices remain when computing the slice mesh, a nullptr is returned
                return computeSliceMesh(solution, a, b, slice, 0)==nullptr;
            }), slices.end());
        }

        return slices;
    }

    template<bool SYMMETRY_BREAKING=true, bool USE_DISJUNCTIVE_SLICES=true, bool FILTER_EMPTY_SLICES=true>
    static std::map<std::pair<size_t, size_t>, std::vector<Slice>> generateSlices(const std::shared_ptr<const StripPackingSolution>& solution){

        std::map<std::pair<size_t, size_t>, std::vector<Slice>> result;

        for (int a = 0; a < solution->getItems().size(); ++a){
            for (int b = a + 1; b < solution->getItems().size(); ++b) {

                const auto& meshA = solution->getItem(a)->getModelSpaceMesh();
                const auto& meshB = solution->getItem(b)->getModelSpaceMesh();
                result[std::make_pair(a, b)] = computeSlices<SYMMETRY_BREAKING, USE_DISJUNCTIVE_SLICES, FILTER_EMPTY_SLICES>(solution, meshA, meshB, a, b);
            }
        }

        return result;
    }
};


#endif //EXTENDEDMESHCORE_CONVEXNFPUTILITIES_H
