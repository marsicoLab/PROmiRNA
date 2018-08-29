// ============================================================================
// Test merging of bg sequences
// ============================================================================

#include <iostream>

#include "../../src/create_background.h"

template <typename TMap>
bool setCompare(TSet const & lhs, TSet const & rhs)
{
    return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin(), [] (auto a, auto b) { return a == b; });
}


int main()
{
    GenomicLocation gl1;
    gl1.chr = "chr1";
    gl1.start = 10000;
    gl1.end = 11000;
    gl1.strand = '+';

    GenomicLocation gl2;
    gl2.chr = "chr1";
    gl2.start = 100;
    gl2.end = 2200;
    gl2.strand = '+';

    GenomicLocation gl3;
    gl3.chr = "chr1";
    gl3.start = 10500;
    gl3.end = 13000;
    gl3.strand = '-';

    GenomicLocation gl4;
    gl4.chr = "chr2";
    gl4.start = 50000;
    gl4.end = 54000;
    gl4.strand = '-';

    GenomicLocation gl5;
    gl5.chr = "chr1";
    gl5.start = 10900;
    gl5.end = 14000;
    gl5.strand = '+';

    GenomicLocation gl6;
    gl6.chr = "chr1";
    gl6.start = 11000;
    gl6.end = 15000;
    gl6.strand = '+';

    GenomicLocation gl7;
    gl7.chr = "chr2";
    gl7.start = 54000;
    gl7.end = 55000;
    gl7.strand = '-';

    std::set<GenomicLocation> backgroundSet;
    backgroundSet.insert(gl1);
    backgroundSet.insert(gl2);
    backgroundSet.insert(gl3);
    backgroundSet.insert(gl4);
    backgroundSet.insert(gl5);
    backgroundSet.insert(gl6);
    backgroundSet.insert(gl7);

    mergeBackgroundRegions(backgroundSet);

    GenomicLocation merged1;
    merged1.chr = "chr1";
    merged1.start = 10000;
    merged1.end = 15000;
    merged1.strand = '+';

    std::set<GenomicLocation> controlSet;
    controlSet.insert(merged1));
    controlSet.insert(gl2);
    controlSet.insert(gl3);
    controlSet.insert(gl4);
    controlSet.insert(gl7);

    if (!setCompare(backgroundSet, controlSet))
    {
        std::cerr << "ERROR: Differing elements in test and control map." << std::endl;
        return 1;
    }

    std::cout << "Test merge background successful." << std::endl;

    return 0;
}
