// ============================================================================
// Test functions from create dataset
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#include <gtest/gtest.h>
#include <iostream>

#include <seqan/sequence.h>

#include "../src/create_dataset.h"

TEST(create_dataset, cpg_content)
{
    Dna5String seq = "ACCGTGGCACCGACTAACAGATACAGATACCAGAGCGCGCA";
    double cpgContent = calculateCpGContent(seq);
    double control = 1.2400756;

    EXPECT_NEAR(cpgContent, control, 0.0001);

    Dna5String seq2 = "ACCGTGGCACCGACTAACAGATACAGATACCAGAGCGCGCG";
    double cpgContent2 = calculateCpGContent(seq2);
    double control2 = 1.42361;

    EXPECT_NEAR(cpgContent2, control2, 0.0001);

    Dna5String seq3 = "ATATATATATATATATATATATATATATATATATATATAT";
    double cpgContent3 = calculateCpGContent(seq3);

    EXPECT_EQ(cpgContent3, 0);

    Dna5String seq4 = "";
    double cpgContent4 = calculateCpGContent(seq4);

    EXPECT_EQ(cpgContent4, 0);
}
