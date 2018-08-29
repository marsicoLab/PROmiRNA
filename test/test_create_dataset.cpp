// ============================================================================
// Test functions from create dataset
// ============================================================================

#include <iostream>

#include <seqan/sequence.h>

#include "../../src/create_dataset.h"

int main()
{
    Dna5String seq = "ACCGTGGCACCGACTAACAGATACAGATACCAGAGCGCGCA";
    double cpgContent = calculateCpGContent(seq);
    double control = 1.2400756;

    if (abs(cpgContent - control) > 0.0001)
    {
        std::cerr << "ERROR: Test CpG content calculation failed for regular sequence." << std::endl;
        return 1;
    }

    Dna5String seq2 = "ACCGTGGCACCGACTAACAGATACAGATACCAGAGCGCGCG";
    double cpgContent2 = calculateCpGContent(seq2);
    double control2 = 1.42361;

    if (abs(cpgContent2 - control2) > 0.0001)
    {
        std::cerr << "ERROR: Test CpG content calculation failed for regular sequence with CG at end." << std::endl;
        return 1;
    }

    Dna5String seq3 = "ATATATATATATATATATATATATATATATATATATATAT";
    double cpgContent3 = calculateCpGContent(seq3);

    if (cpgContent3 != 0)
    {
        std::cerr << "ERROR: Test CpG content calculation failed for sequence without GC." << std::endl;
        return 1;
    }

    Dna5String seq4 = "";
    double cpgContent4 = calculateCpGContent(seq4);

    if (cpgContent4 != 0)
    {
        std::cerr << "ERROR: Test CpG content calculation failed for empty string." << std::endl;
        return 1;
    }

    std::cout << "Test CpG content calculation successful." << std::endl;

    return 0;
}
