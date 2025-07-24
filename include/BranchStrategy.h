//
// Created by Jonas Tollenaere on 04/07/2025.
//

#pragma once

#ifndef BRANCHSTRATEGY_H
#define BRANCHSTRATEGY_H
#include <string>

enum BranchStrategy {
    DEFAULT,                        // Gurobi default
    SORTED_LARGEST_ITEM_VOLUME,     // Increment build, starting with largest volume items
    SORTED_LARGEST_ITEM_AABB_VOLUME,// Increment build, starting with largest volume items
    SORTED_LARGEST_ITEM_HEIGHT,     // Increment build, starting with tallest items
    SORTED_LARGEST_ITEM_SURFACE_XY, // Increment build, starting with largest xy surface area items
    LARGEST_NFP_AABB_VOLUME,        // Prioritize NFPs with largest AABB volume
    LEAST_NFP_FACES,                // Prioritize NFPs with the least faces
};

inline std::string branchStrategyToString(BranchStrategy strategy){
    switch (strategy) {
        case DEFAULT:
            return "DEFAULT";
        case SORTED_LARGEST_ITEM_VOLUME:
            return "SORTED-LARGEST-ITEM-VOLUME";
        case SORTED_LARGEST_ITEM_AABB_VOLUME:
            return "SORTED-LARGEST-ITEM-AABB-VOLUME";
        case SORTED_LARGEST_ITEM_HEIGHT:
            return "SORTED-LARGEST-ITEM-HEIGHT";
        case SORTED_LARGEST_ITEM_SURFACE_XY:
            return "SORTED-LARGEST-ITEM-SURFACE-XY";
        case LARGEST_NFP_AABB_VOLUME:
            return "LARGEST-NFP-AABB-VOLUME";
        case LEAST_NFP_FACES:
            return "LEAST-NFP-FACES";
        default:
            return "UNKNOWN";
    }
}

#endif //BRANCHSTRATEGY_H
