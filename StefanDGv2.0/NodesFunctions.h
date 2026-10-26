#pragma once
#include <cstdint>
#include "GeometryConstants.h"

namespace NodesFunctions
{
    void extractFaceNodesTags(const uint8_t localIndex, const size_t tetrahedronNodes[constants::tetrahedron::N_NODES], size_t faceNodes[constants::triangle::N_NODES]);

    uint8_t computeChangingIndexes(const size_t original[constants::triangle::N_NODES], const size_t comparison[constants::triangle::N_NODES]);

}

