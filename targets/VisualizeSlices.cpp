//
// Created by Jonas Tollenaere on 04/07/2025.
//

#include <meshcore/optimization/StripPackingSolution.h>
#include <meshcore/rendering/ApplicationWindow.h>

#include "ConvexNFPUtilities.h"

#define RENDER_NFP true
#define RENDER_NFP_SLICE true
#define RENDER_NFP_SLICE_PLANES false

#define USE_SYMMETRY_BREAKING true
#define USE_DISJUNCT_SLICES true
#define FILTER_EMPTY_SLICES true

class LoadSolutionTask: public AbstractTask {
public:
    void run() override {

        // Create and notify the solution
        notifyObserversStatus("Loading solution");
        auto start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

        std::string path = std::string(SOLUTION_DIR) + "STOYAN_2005_EXAMPLE_1.json";

        auto solution = StripPackingSolution::fromJson(path);
        auto end = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        std::cout << "Creating solution took " << (end-start) << "ms" << std::endl;

        notifyObserversSolution(solution);

        auto result = solution->toJson();

        std::cout << "Feasible: " << solution->isFeasible() << std::endl; // The mathematical solution is not feasible due to constraint tolerances
    }
};

int main(int argc, char *argv[]){

    LoadSolutionTask task;

// How to render the solutions
    std::function<void(RenderWidget* renderWidget, std::shared_ptr<const AbstractSolution> solution)> onSolutionNotified = [&](RenderWidget* renderWidget, const std::shared_ptr<const AbstractSolution>& sol){

        auto solution = std::dynamic_pointer_cast<const StripPackingSolution>(sol);

        renderWidget->renderBox("Container", "Bounds", solution->getProblem()->getContainer());
        renderWidget->clearGroup("MinimalContainer");
        float maximumHeight = 0.0f;

        for (int itemIndex = 0; itemIndex < solution->getItems().size(); ++itemIndex){
            auto item =   solution->getItem(itemIndex);
            auto maximumItemHeight = solution->getItemAABB(itemIndex).getMaximum().z;
            if(maximumItemHeight > maximumHeight){
                maximumHeight = maximumItemHeight;
            }

            renderWidget->renderWorldSpaceMesh("Items", item,  StripPackingProblem::getItemColor(solution->getItemName(itemIndex)));
        }
        auto min = solution->getProblem()->getContainer().getMinimum();
        auto max = solution->getProblem()->getContainer().getMaximum();

        renderWidget->renderBox("MinimalContainer", "AABB", {min, {max.x, max.y, maximumHeight}});

        // Draw NFVs and their potential bounding planes
        renderWidget->clearGroup("NFP");
        if(RENDER_NFP) { // Whether to render the NFPs

            for (int a = 0; a < solution->getItems().size(); ++a){
                for(int b = a + 1; b < solution->getItems().size(); ++b){
                    const auto& itemA = solution->getItems()[a];
                    const auto& itemB = solution->getItems()[b];
                    const auto& meshA = solution->getItem(a)->getModelSpaceMesh();
                    const auto& meshB = solution->getItem(b)->getModelSpaceMesh();

                    // Average color
                    auto color = 0.5f * StripPackingProblem::getItemColor(solution->getItemName(a)) +
                                            0.5f * StripPackingProblem::getItemColor(solution->getItemName(b));

                    auto NFP = ConvexNFPUtilities::computeNFP(meshA, meshB);
                    NFP->setName("NFP_" + std::to_string(a) + "_" + std::to_string(b));
                    renderWidget->renderWorldSpaceMesh("NFP", std::make_shared<WorldSpaceMesh>(NFP), Color(color));
                }
            }

            if(RENDER_NFP_SLICE){ // Whether to render the slices
                auto slicesMap = ConvexNFPUtilities::generateSlices<USE_SYMMETRY_BREAKING, USE_DISJUNCT_SLICES, FILTER_EMPTY_SLICES>(solution);

                for (int a = 0; a < solution->getItems().size(); ++a){
                    for(int b = a + 1; b < solution->getItems().size(); ++b){

                        // Name
                        std::string name = "Slices_" + std::to_string(a) + "_" + std::to_string(b);
                        renderWidget->clearGroup(name);

                        // Average color
                        const auto& itemA = solution->getItems()[a];
                        const auto& itemB = solution->getItems()[b];
                        auto color = 0.25f * StripPackingProblem::getItemColor(solution->getItemName(a))
                                     + 0.25f * StripPackingProblem::getItemColor(solution->getItemName(b))
                                     + 0.5f * Color(0.8, 0.8, 0.8, 0.6);

                        const auto slices = slicesMap[std::make_pair(a, b)];
                        for (int i = 0; i < slices.size(); ++i){
                            const auto& slice = slices[i];
                            const auto& mesh = ConvexNFPUtilities::computeSliceMesh(solution, a, b, slice, i);
                            renderWidget->renderWorldSpaceMesh(name, std::make_shared<WorldSpaceMesh>(mesh), Color(color));
                        }
                    }
                }

                if(RENDER_NFP_SLICE_PLANES){ // Whether to render the slice's bounding planes
                    for (int a = 0; a < solution->getItems().size(); ++a){
                        for(int b = a + 1; b < solution->getItems().size(); ++b){

                            // Average color
                            const auto& itemA = solution->getItems()[a];
                            const auto& itemB = solution->getItems()[b];
                            auto color = 0.25f * StripPackingProblem::getItemColor(solution->getItemName(a))
                                         + 0.25f * StripPackingProblem::getItemColor(solution->getItemName(b))
                                         + 0.5f * Color(0.8, 0.8, 0.8, 0.6);

                            const auto& slices = slicesMap[std::make_pair(a, b)];
                            for(int k = 0; k < slices.size(); ++k){
                                const auto& slice = slices[k];
                                for(int l = 0; l < slice.boundingPlanes.size(); ++l){
                                    const auto& plane = slice.boundingPlanes[l];
                                    renderWidget->renderPlane("Planes", "Plane_" + std::to_string(a) + "_" + std::to_string(b) + "_" + std::to_string(k) + "_" + std::to_string(l),
                                                              plane, Color(color));
                                }
                            }
                        }
                    }

                }
            }
        }
    };


    QApplication app(argc, argv);
    ApplicationWindow window;
    window.show();
    RenderWidget* renderWidget = window.getRenderWidget();
    renderWidget->observeTask(&task, onSolutionNotified);
    task.start();
    int returnCode = QApplication::exec();
    task.stop();
    task.join();
    return returnCode;
}

