// ============================================================================
// Test quantile normalization
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#include <gtest/gtest.h>
#include <iostream>

#include "../src/quantile_normalization.h"

TEST(quantile_normalization_test, quantile_normalization)
{
    std::vector<std::vector<std::pair<float, unsigned> > > input{{std::make_pair(10,0), std::make_pair(8,0), std::make_pair(4,0), std::make_pair(4,0)},
                                                             {std::make_pair(20,1), std::make_pair(6,1), std::make_pair(0,1), std::make_pair(4,1)},
                                                             {std::make_pair(50,2), std::make_pair(4,2), std::make_pair(8,2), std::make_pair(10,2)},
                                                             {std::make_pair(30,3), std::make_pair(2,3), std::make_pair(6,3), std::make_pair(16,3)},
                                                             {std::make_pair(20,4), std::make_pair(1,4), std::make_pair(2,4), std::make_pair(20,4)}};

    std::vector<std::vector<float> > control{{3.75, 21.5, 9.5, 3.75}, {7, 14.5, 3.75, 7}, {21.5, 9.5, 21.5, 9.5}, {14.5, 7, 14.5, 14.5}, {9.5, 3.75, 7, 21.5}};

    MatrixPair inputMatrix;
    inputMatrix.matrix = input;

    MatrixSingle outputMatrix;
    quantileNormalization(outputMatrix, inputMatrix);

    EXPECT_EQ(outputMatrix.matrix, control);
}
