//
// Created by Jonas on 30/08/2024.
//

#include <meshcore/rendering/ApplicationWindow.h>
#include <../datasets/StripPackingInstances.h>

#include "../include/MILPTask.h"

#define USE_SYMMETRY_BREAKING true
#define USE_DISJUNCT_SLICES true
#define USE_LIFTING true
#define BRANCHING_STRATEGY SORTED_LARGEST_ITEM_SURFACE_XY

int main(int argc, char *argv[]){

    MILPTask<USE_SYMMETRY_BREAKING,USE_DISJUNCT_SLICES,USE_LIFTING,BRANCHING_STRATEGY> task(TOLLENAERE_2025_STOYAN_2005_09_ITEMS);

    QApplication app(argc, argv);
    ApplicationWindow window;
    window.show();
    RenderWidget* renderWidget = window.getRenderWidget();
    renderWidget->observeTask(&task);
    task.start();
    const int returnCode = QApplication::exec();
    task.stop();
    task.join();
    return returnCode;
}