// ============================================================================
// Test quantile normalization
// ============================================================================

#include <iostream>

#include "../../src/quantile_normalization.h"

int main()
{
    std::vector<std::vector<std::pair<float, unsigned> > > input{{std::make_pair(10,0), std::make_pair(8,0), std::make_pair(4,0), std::make_pair(4,0)},
                                                             {std::make_pair(20,1), std::make_pair(6,1), std::make_pair(0,1), std::make_pair(4,1)},
                                                             {std::make_pair(50,2), std::make_pair(4,2), std::make_pair(8,2), std::make_pair(10,2)},
                                                             {std::make_pair(30,3), std::make_pair(2,3), std::make_pair(6,3), std::make_pair(16,3)},
                                                             {std::make_pair(20,4), std::make_pair(1,4), std::make_pair(2,4), std::make_pair(20,4)}};

    std::vector<std::vector<float> > control{{3.75, 21.5, 9.5, 3.75}, {7, 14.5, 3.75, 7}, {21.5, 9.5, 21.5, 9.5}, {14.5, 7, 14.5, 14.5}, {9.5, 3.75, 7, 21.5}};

    MatrixPair inputMatrix;
    inputMatrix.matrix = input;

    Matrix outputMatrix;
    quantileNormalization(outputMatrix, inputMatrix);

    for (unsigned i = 0; i < outputMatrix.matrix.size(); i++)
    {
        for (unsigned j = 0; j < outputMatrix.matrix[i].size(); j++)
        {
            if (outputMatrix.matrix[i][j] != control[i][j])
            {
                std::cerr << "ERROR: Control and output matrix not the same." << std::endl;
                return 1;
            }
        }
    }

    std::cout << "Test quantile normalization successful." << std::endl;
    return 0;
}
