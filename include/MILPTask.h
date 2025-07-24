//
// Created by Jonas Tollenaere on 04/07/2025.
//

#ifndef MIPTASK_H
#define MIPTASK_H

#include "gurobi_c++.h"

#include <fstream>
#include <utility>

#include <meshcore/optimization/StripPackingSolution.h>
#include <meshcore/tasks/AbstractTask.h>

#include "BranchStrategy.h"
#include "ConvexNFPUtilities.h"

#define FILTER_EMPTY_SLICES true

template<bool USE_SYMMETRY_BREAKING, bool USE_DISJUNCTIVE_SLICES, bool USE_LIFTING, BranchStrategy branchStrategy, bool EXPORT_ONLY=false>
class MILPTask: public AbstractTask {

    std::string instancePath;

public:
    explicit MILPTask(const std::string& instancePath): instancePath(instancePath) {}

    class Callback: public GRBCallback {
    public:

        std::shared_ptr<StripPackingSolution> solution;

        GRBVar HVar;

        std::vector<GRBVar> xVariables;
        std::vector<GRBVar> yVariables;
        std::vector<GRBVar> zVariables;

        MILPTask* task;

        Callback(MILPTask* task, const std::shared_ptr<StripPackingSolution>& callbackSolution,
                 const GRBVar& HVar,
                 const std::vector<GRBVar>& xVariables,
                 const std::vector<GRBVar>& yVariables,
                 const std::vector<GRBVar>& zVariables):
                HVar(HVar),
                xVariables(xVariables),
                yVariables(yVariables),
                zVariables(zVariables),
                solution(callbackSolution), task(task) {}

    protected:
        void callback() override {
            if (where == GRB_CB_MIPSOL) {

                for (int i = 0; i < solution->getItems().size(); ++i){

                    const auto& item = solution->getItems()[i];

                    Transformation transformation;
                    transformation.setPositionX(getSolution(xVariables[i]));
                    transformation.setPositionY(getSolution(yVariables[i]));
                    transformation.setPositionZ(getSolution(zVariables[i]));
                    solution->setItemTransformation(i, transformation);
                }

                task->notifyObserversSolution(solution);

                const auto& height = getSolution(HVar);

                task->notifyObserversStatus("Found solution with height " + std::to_string(height));
            }

            if (this->task->stopCalled) {
                abort();
            }
        }
    };

    void run() override {

        this->notifyObserversStatus("Initialising");

        // Load the problem
        auto start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        const std::shared_ptr<StripPackingProblem> problem = StripPackingProblem::fromInstancePath(instancePath, ObjectOrigin::AlignToMinimum);
        auto end = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        std::cout << "Loading problem took " << (end-start) << "ms" << std::endl;

        // Create and notify the solution
        notifyObserversStatus("Constructing solution");
        auto startms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        auto solution = std::make_shared<StripPackingSolution>(problem);
        auto endms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        std::cout << "Creating solution took " << (endms-startms) << "ms" << std::endl;

        // Derive a lower bound based on the total volume of the items
        const auto& container = solution->getProblem()->getContainer();
        float totalVolume = problem->getTotalItemVolume();
        auto containerDimensions = container.getMaximum() - container.getMinimum();
        double minimumHeight = totalVolume/(containerDimensions.x * containerDimensions.y);
        printf("Lower bound based off item volumes: %f\n", minimumHeight);

        // Check if an item is higher than this lower bound, if so, update the lower bound
        {
            bool updated = false;
            for (size_t itemIndex = 0; itemIndex < solution->getItems().size(); ++itemIndex) {

                auto aabb = solution->getItemAABB(itemIndex);
                auto height = aabb.getMaximum().z - aabb.getMinimum().z;

                if(height >= minimumHeight){
                    minimumHeight = height; // The final solution will be at least as high as the highest item
                    updated = true;
                }
            }
            if(updated) printf("Lower bound updated with largest item height: %f\n", minimumHeight);
        }

        try {
            // Create an environment
            GRBEnv env = GRBEnv(true);
            env.set("LogFile", "mip.log");
            env.start();

            // Create an empty model
            GRBModel model = GRBModel(env);
            model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

            // Set 1 starting solution
            model.set(GRB_IntParam_PoolSearchMode, 2);

            // Variable height as objective, other container dimensions fixed
            GRBVar HVar = model.addVar(container.getMinimum().z + minimumHeight, container.getMaximum().z, 1.0, GRB_CONTINUOUS, "H");
            HVar.set(GRB_DoubleAttr_Start, container.getMaximum().z);

            const auto HMax = problem->getContainer().getMaximum().z;
            const auto L = problem->getContainer().getMaximum().x;
            const auto W = problem->getContainer().getMaximum().y;

            // Define decision variables and constraints per item
            std::vector<GRBVar> xVariables;
            std::vector<GRBVar> yVariables;
            std::vector<GRBVar> zVariables;

            auto initialStartingHeight = 0.0; // Incremented z-coordinates for initial solution
            for (int itemIndex = 0; itemIndex < solution->getItems().size(); ++itemIndex) {
                const auto& item = solution->getItems()[itemIndex];
                const auto& aabb = solution->getItem(itemIndex)->getModelSpaceMesh()->getBounds();

                // Notation used in paper
                const auto li = aabb.getMaximum().x;
                const auto wi = aabb.getMaximum().y;
                const auto hi = aabb.getMaximum().z;

                AABB sensiblePositions(container.getMinimum() - aabb.getMinimum(), container.getMaximum() - aabb.getMaximum());

                GRBVar xi = model.addVar(0, L - li, 0, GRB_CONTINUOUS, "x" + std::to_string(itemIndex));
                GRBVar yi = model.addVar(0, W - wi, 0, GRB_CONTINUOUS, "y" + std::to_string(itemIndex));
                GRBVar zi = model.addVar(0, HMax - hi, 0, GRB_CONTINUOUS, "z" + std::to_string(itemIndex));

                // Set starting positions
                xi.set(GRB_DoubleAttr_Start, 0);
                yi.set(GRB_DoubleAttr_Start, 0);
                zi.set(GRB_DoubleAttr_Start, initialStartingHeight);
                initialStartingHeight += hi;

                xVariables.push_back(xi);
                yVariables.push_back(yi);
                zVariables.push_back(zi);

                // Constrain vertices to be inside the container
                model.addConstr(hi + zi <= HVar, "Z_upper_bound_" + std::to_string(itemIndex));
            }

            // Symmetry breaking constraints for identical items
            if (USE_SYMMETRY_BREAKING) {
                size_t numberOfItems = 0;
                for (unsigned long long count : problem->getRequiredItemCounts()){
                    for(auto a = 0; a < count; ++a){
                        for(auto b = a + 1; b < count; ++b){
                            model.addConstr(zVariables[numberOfItems + a] <= zVariables[numberOfItems + b], "Symmetry breaking");
                        }
                    }
                    numberOfItems += count;
                }
            }

            // Keep track of item index when sorted
            std::vector<int> sortedItemOrder(solution->getItems().size());
            {
                std::vector<size_t> sortedItemIndices(solution->getItems().size());
                std::iota(sortedItemIndices.begin(), sortedItemIndices.end(), 0);

                if(branchStrategy == SORTED_LARGEST_ITEM_VOLUME){
                    // Sort largest volume items first
                    std::sort(sortedItemIndices.begin(), sortedItemIndices.end(), [&](size_t a, size_t b){
                        const auto& meshA = solution->getItems()[a]->getModelSpaceMesh();
                        const auto& meshB = solution->getItems()[b]->getModelSpaceMesh();
                        return meshA->getVolume() > meshB->getVolume();
                    });
                }
                else if(branchStrategy == SORTED_LARGEST_ITEM_AABB_VOLUME) {
                    // Sort largest bounds volume items first
                    std::sort(sortedItemIndices.begin(), sortedItemIndices.end(), [&](size_t a, size_t b){
                        const auto& meshA = solution->getItems()[a]->getModelSpaceMesh();
                        const auto& meshB = solution->getItems()[b]->getModelSpaceMesh();
                        return meshA->getBounds().getVolume() > meshB->getBounds().getVolume();
                    });
                }
                else if(branchStrategy == SORTED_LARGEST_ITEM_HEIGHT){
                    // Sort tallest items first
                    std::sort(sortedItemIndices.begin(), sortedItemIndices.end(), [&](size_t a, size_t b){
                        const auto& meshA = solution->getItems()[a]->getModelSpaceMesh();
                        const auto& meshB = solution->getItems()[b]->getModelSpaceMesh();

                        const auto heightA = meshA->getBounds().getMaximum().z - meshA->getBounds().getMinimum().z;
                        const auto heightB = meshB->getBounds().getMaximum().z - meshB->getBounds().getMinimum().z;

                        return heightA > heightB;
                    });
                }
                else if(branchStrategy == SORTED_LARGEST_ITEM_SURFACE_XY){
                    // Sort largest xy surface area items first
                    std::sort(sortedItemIndices.begin(), sortedItemIndices.end(), [&](size_t a, size_t b){
                        const auto& meshA = solution->getItems()[a]->getModelSpaceMesh();
                        const auto& meshB = solution->getItems()[b]->getModelSpaceMesh();

                        const auto boundsXLengthA = meshA->getBounds().getMaximum().x - meshA->getBounds().getMinimum().x;
                        const auto boundsYLengthA = meshA->getBounds().getMaximum().y - meshA->getBounds().getMinimum().y;
                        const auto boundsXLengthB = meshB->getBounds().getMaximum().x - meshB->getBounds().getMinimum().x;
                        const auto boundsYLengthB = meshB->getBounds().getMaximum().y - meshB->getBounds().getMinimum().y;

                        return boundsXLengthA * boundsYLengthA > boundsXLengthB * boundsYLengthB;
                    });
                }

                // Keep track of where each item is in the sorted list
                std::cout << "Sorted item indices: ";
                for (int i = 0; i < sortedItemIndices.size(); ++i){
                    std::cout << sortedItemIndices[i] << " ";
                    sortedItemOrder[sortedItemIndices[i]] = i;
                }
                std::cout << std::endl;
            }

            // Enforce separation between each pair of items
            auto sliceMap = ConvexNFPUtilities::generateSlices<USE_SYMMETRY_BREAKING, USE_DISJUNCTIVE_SLICES, FILTER_EMPTY_SLICES>(solution);
            for (int a = 0; a < solution->getItems().size(); ++a){
                for (int b = a + 1; b < solution->getItems().size(); ++b) {

                    const auto& meshA = solution->getItems()[a]->getModelSpaceMesh();
                    const auto& meshB = solution->getItems()[b]->getModelSpaceMesh();

                    int branchPriority = 0;

                    if(branchStrategy == SORTED_LARGEST_ITEM_SURFACE_XY || branchStrategy == SORTED_LARGEST_ITEM_HEIGHT || branchStrategy == SORTED_LARGEST_ITEM_VOLUME || branchStrategy == SORTED_LARGEST_ITEM_AABB_VOLUME){
                        auto maxIndex = std::max(sortedItemOrder[a],sortedItemOrder[b]);
                        auto minIndex = std::min(sortedItemOrder[a],sortedItemOrder[b]);
                        branchPriority = - maxIndex * int(solution->getItems().size()) - minIndex;
                    }
                    else if(branchStrategy == LARGEST_NFP_AABB_VOLUME){
                        branchPriority = 100 * int(ConvexNFPUtilities::computeNFP(meshA, meshB)->getBounds().getVolume()) + b * solution->getItems().size() + a;
                    }
                    else if(branchStrategy == LEAST_NFP_FACES){
                        branchPriority = -ConvexNFPUtilities::computeNFP(meshA, meshB)->getFaces().size();
                    }

                    // Sum of slices that are enforced
                    GRBLinExpr sum;

                    // Translation variables for both models
                    const auto& xa = xVariables[a];
                    const auto& ya = yVariables[a];
                    const auto& za = zVariables[a];

                    const auto& xb = xVariables[b];
                    const auto& yb = yVariables[b];
                    const auto& zb = zVariables[b];

                    const auto la = meshA->getBounds().getMaximum().x;
                    const auto wa = meshA->getBounds().getMaximum().y;
                    const auto ha = meshA->getBounds().getMaximum().z;

                    const auto lb = meshB->getBounds().getMaximum().x;
                    const auto wb = meshB->getBounds().getMaximum().y;
                    const auto hb = meshB->getBounds().getMaximum().z;

                    // Support variables used for tighter big M constraints
                    // The extreme values that the relative translation of a and b can take
                    auto BMinAMax = glm::vec3(L-lb, W-wb, HMax-hb); // When A is in (0,0,0) and B in the maximum position
                    auto BMinAMin = - glm::vec3(L-la, W-wa, HMax-ha); // When B is in (0,0,0) and A in the minimum position

                    // Precompute the slice meshes
                    const auto& slices = sliceMap[std::make_pair(a, b)];
                    auto sliceMeshes = std::vector<std::shared_ptr<ModelSpaceMesh>>();
                    sliceMeshes.reserve(slices.size());
                    for (int k = 0; k < slices.size(); ++k){
                        sliceMeshes.push_back(ConvexNFPUtilities::computeSliceMesh(solution, a, b, slices[k], k));
                    }

                    std::vector<GRBVar> enforcedVariables;
                    for (int k = 0; k < slices.size(); ++k){

                        auto& sliceK = slices[k];

                        GRBVar enforced = model.addVar(0, 1, 0, GRB_BINARY, "Enforced_" + std::to_string(a) + "_" + std::to_string(b) + "_" + std::to_string(k));
                        enforcedVariables.push_back(enforced);
                        sum += enforced;

                        if(branchStrategy!=DEFAULT) enforced.set(GRB_IntAttr_BranchPriority, branchPriority); // Set branch priority determined above

                        for (const auto &plane: sliceK.boundingPlanes){

                            const auto& n = plane.getNormal();
                            const auto& d = plane.getD();

                            // tightM should be the min value that n * (TB - TA) + d can take so >=-tightM always holds when !enforced
                            auto minXDot = n.x <= 0 ? n.x * BMinAMax.x : n.x * BMinAMin.x;
                            auto minYDot = n.y <= 0 ? n.y * BMinAMax.y : n.y * BMinAMin.y;
                            auto minZDot = n.z <= 0 ? n.z * BMinAMax.z : n.z * BMinAMin.z;
                            const auto tightM = minXDot + minYDot + minZDot + d;

                            model.addConstr(n.x * (xb - xa) + n.y * (yb - ya) + n.z * (zb - za) + d >= tightM * (1 - enforced), "Separation_Constraint_" + std::to_string(a) + "_" + std::to_string(b) + "_" + std::to_string(k) + "_" + std::to_string(&plane - &sliceK.boundingPlanes[0]));
                        }
                    }
                    model.addConstr(sum == 1, "One slice per NFP should be enforced");

                    if(USE_LIFTING){

                        GRBLinExpr sumXMinEnforced;
                        GRBLinExpr sumYMinEnforced;
                        GRBLinExpr sumZMinEnforced;
                        GRBLinExpr sumXMaxEnforced;
                        GRBLinExpr sumYMaxEnforced;
                        GRBLinExpr sumZMaxEnforced;
                        for (int h = 0; h < slices.size(); ++h) {

                            auto& sliceHEnforced = enforcedVariables[h];

                            const auto& sliceHMesh = sliceMeshes[h];
                            assert(!sliceHMesh->getVertices().empty());

                            auto bounds = sliceHMesh->getBounds();

                            sumXMinEnforced += bounds.getMinimum().x * sliceHEnforced;
                            sumXMaxEnforced += bounds.getMaximum().x * sliceHEnforced;
                            sumYMinEnforced += bounds.getMinimum().y * sliceHEnforced;
                            sumYMaxEnforced += bounds.getMaximum().y * sliceHEnforced;
                            sumZMinEnforced += bounds.getMinimum().z * sliceHEnforced;
                            sumZMaxEnforced += bounds.getMaximum().z * sliceHEnforced;
                        }
                        model.addConstr(xb - xa >= sumXMinEnforced, "X lifting per slice, min");
                        model.addConstr(xb - xa <= sumXMaxEnforced, "X lifting per slice, max");
                        model.addConstr(yb - ya >= sumYMinEnforced, "Y lifting per slice, min");
                        model.addConstr(yb - ya <= sumYMaxEnforced, "Y lifting per slice, max");
                        model.addConstr(zb - za >= sumZMinEnforced, "Z lifting per slice, min");
                        model.addConstr(zb - za <= sumZMaxEnforced, "Z lifting per slice, max");
                    }
                }
            }

            if (EXPORT_ONLY){
                std::string name = problem->getName();
                if(USE_SYMMETRY_BREAKING) name += "_SYMMETRY-BREAKING";
                if(USE_DISJUNCTIVE_SLICES) name +=  "_DISJUNCT";
                if(USE_LIFTING) name += "_LIFTED";
                name += "_" + branchStrategyToString(branchStrategy);
                if(branchStrategy != DEFAULT){
                    model.write(MODEL_DIR + name + ".ord");
                    std::cout << "Priorities written to " << name << ".ord" << std::endl;
                }
                model.write(MODEL_DIR + name + ".mps");
                std::cout << "Model written to " << name << ".mps" << std::endl;

                // Write the starting point to an .mst file
                model.write(MODEL_DIR + name + ".mst");
            }
            else {
                // Configure callback
                Callback cb(this, solution, HVar, xVariables, yVariables, zVariables);
                model.setCallback(&cb);

                // Optimize model
                notifyObserversStatus("Optimizing");
                model.optimize();

                for (int i = 0; i < solution->getItems().size(); ++i){
                    const auto& item = solution->getItems()[i];
                    Transformation transformation;
                    transformation.setPositionX(xVariables[i].get(GRB_DoubleAttr_X));
                    transformation.setPositionY(yVariables[i].get(GRB_DoubleAttr_X));
                    transformation.setPositionZ(zVariables[i].get(GRB_DoubleAttr_X));
                    solution->setItemTransformation(i, transformation);
                }
                notifyObserversSolution(solution);
                notifyObserversProgress(1.0);

                // Export solution
                auto solutionJSON = solution->toJson();
                solutionJSON["totalHeight"] = HVar.get(GRB_DoubleAttr_X);
                solutionJSON["gap"] = model.get(GRB_DoubleAttr_MIPGap);
                std::ofstream jsonOutputFile(SOLUTION_DIR + problem->getName() + ".json");
                jsonOutputFile << solutionJSON.dump(4);
                jsonOutputFile.close();
            }

        } catch(GRBException& e) {
            std::cout << "Error code = " << e.getErrorCode() << " (" << e.getMessage() << ")" << std::endl;
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }
    }
};

#endif //MIPTASK_H
