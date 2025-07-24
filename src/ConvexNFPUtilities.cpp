//
// Created by Jonas Tollenaere on 04/07/2025.
//

#include "ConvexNFPUtilities.h"

std::shared_ptr<ModelSpaceMesh> ConvexNFPUtilities::meshFromPlanes(const std::vector<Plane> &planes) {

    std::vector<glm::vec3> intersections;

    for (int i = 0; i < planes.size(); ++i) {
        for (int j = i + 1; j < planes.size(); ++j) {
            for (int k = j + 1; k < planes.size(); ++k) {
                std::optional<glm::vec3> point = Intersection::intersect(planes[i], planes[j], planes[k]);
                if(point.has_value()){
                    intersections.push_back(point.value());
                }
            }
        }
    }

    // Filter out intersections at negative side of one of the planes
    std::vector<glm::vec3> filteredVertices;
    for (const auto& vertex: intersections){
        bool inside = true;
        for (const auto& plane: planes){
            if(glm::dot(plane.getNormal(), vertex) + plane.getD() < -1e-2){
                inside = false;
                break;
            }

        }
        if(inside){
            filteredVertices.push_back(vertex);
        }
    }

    if(filteredVertices.empty()){
        return nullptr;
    }

    return ModelSpaceMesh(filteredVertices).getConvexHull();
}

std::shared_ptr<ModelSpaceMesh> ConvexNFPUtilities::computeSliceMesh(
    const std::shared_ptr<const StripPackingSolution> &solution, int a, int b, const Slice &slice, int sliceIndex) {
    const auto& meshA = solution->getItems()[a]->getModelSpaceMesh();
    const auto& meshB = solution->getItems()[b]->getModelSpaceMesh();

    // Determine the extreme relative positions possible for the meshes
    auto container = solution->getProblem()->getContainer();
    auto sensiblePositionsBMaximum = container.getMaximum() - meshB->getBounds().getMaximum();
    auto sensiblePositionsBMinimum = container.getMinimum() - meshB->getBounds().getMinimum();
    auto sensiblePositionsAMinimum = container.getMinimum() - meshA->getBounds().getMinimum();
    auto sensiblePositionsAMaximum = container.getMaximum() - meshA->getBounds().getMaximum();
    auto BMinAMax = sensiblePositionsBMaximum - sensiblePositionsAMinimum;
    auto BMinAMin = sensiblePositionsBMinimum - sensiblePositionsAMaximum;

    // Gather the planes that bound these extreme positions
    std::vector<Plane> boundingPlanes;
    boundingPlanes.emplace_back(glm::vec3(1, 0, 0), BMinAMin);
    boundingPlanes.emplace_back(glm::vec3(-1, 0, 0), BMinAMax);
    boundingPlanes.emplace_back(glm::vec3(0, 1, 0), BMinAMin);
    boundingPlanes.emplace_back(glm::vec3(0, -1, 0), BMinAMax);
    boundingPlanes.emplace_back(glm::vec3(0, 0, 1), BMinAMin);
    boundingPlanes.emplace_back(glm::vec3(0, 0, -1), BMinAMax);

    // Add the bounding planes of the slice
    boundingPlanes.insert(boundingPlanes.end(), slice.boundingPlanes.begin(), slice.boundingPlanes.end());

    // Create a mesh from the planes
    auto sliceMesh = meshFromPlanes(boundingPlanes);

    if(sliceMesh == nullptr){
        return nullptr;
    }

    if(!slice.boundingPlanes.empty()){
        auto disconnected = true;
        auto& basePlane = slice.boundingPlanes[0];
        for (const auto &item: sliceMesh->getVertices()){
            auto distance = basePlane.distance(item);
            if(distance <= 1e-4f){
                disconnected = false;
                break;
            }
        }

        if(disconnected){
            return nullptr;
        }
    }

    sliceMesh->setName("Slice_" + std::to_string(a) + "_" + std::to_string(b) + "_" + std::to_string(sliceIndex));
    return sliceMesh;
}

std::shared_ptr<ModelSpaceMesh> ConvexNFPUtilities::computeNFP(const std::shared_ptr<ModelSpaceMesh> &meshA,
    const std::shared_ptr<ModelSpaceMesh> &meshB) {

    assert(meshA->isConvex());
    assert(meshB->isConvex());

    // Naive NFP calculation: Minkowski difference of all vertices, compute their convex hull
    std::vector<Vertex> nfpVertices;
    nfpVertices.reserve(meshA->getVertices().size() * meshB->getVertices().size());
    for (const auto &va: meshA->getVertices()){
        for (const auto &vb: meshB->getVertices()){
            nfpVertices.push_back(va - vb);
        }
    }
    ModelSpaceMesh nfpMesh(nfpVertices);
    return nfpMesh.getConvexHull();
}
