#include "NodesFunctions.h"

namespace NodesFunctions
{
    void extractFaceNodesTags(const uint8_t localIndex, const size_t tetrahedronNodes[constants::tetrahedron::N_NODES], size_t faceNodes[constants::triangle::N_NODES])
    {
        switch (localIndex)
        {
        case 0:
        {
            faceNodes[0] = tetrahedronNodes[0];
            faceNodes[1] = tetrahedronNodes[2];
            faceNodes[2] = tetrahedronNodes[1];
            break;
        }
        case 1:
        {
            faceNodes[0] = tetrahedronNodes[0];
            faceNodes[1] = tetrahedronNodes[1];
            faceNodes[2] = tetrahedronNodes[3];
            break;
        }
        case 2:
        {
            faceNodes[0] = tetrahedronNodes[0];
            faceNodes[1] = tetrahedronNodes[3];
            faceNodes[2] = tetrahedronNodes[2];
            break;
        }
        case 3:
        {
            faceNodes[0] = tetrahedronNodes[3];
            faceNodes[1] = tetrahedronNodes[1];
            faceNodes[2] = tetrahedronNodes[2];
            break;
        }
        }
    }

    uint8_t computeChangingIndexes(const size_t original[constants::triangle::N_NODES], const size_t comparison[constants::triangle::N_NODES])
    {
        uint8_t changingIndexes = 0;
        if (*original == *comparison)
        {
            ++original;
            ++comparison;
            if (*original == *comparison)
            {
                changingIndexes |= 1 << 2;
                changingIndexes |= 2 << 4;

                return changingIndexes;
            }

            changingIndexes |= 2 << 2;
            changingIndexes |= 1 << 4;
            return changingIndexes;
        }

        ++comparison;
        if (*original == *comparison)
        {
            changingIndexes |= 1;
            ++original;
            ++comparison;

            if (*original == *comparison)
            {
                changingIndexes |= 2 << 2;
                return changingIndexes;
            }

            changingIndexes |= 2 << 4;
            return changingIndexes;

        }

        changingIndexes |= 2;
        ++original;
        if (*original == *comparison)
        {
            changingIndexes |= 1 << 2;
            return changingIndexes;
        }

        changingIndexes |= 1 << 4;
        return changingIndexes;
    }
}
