// ============================================================================
// Test single overlap bg and TSS
// ============================================================================

#include <iostream>

#include <seqan/sequence.h>

#include "../../src/overlap.h"

template <typename TMap>
bool mapCompare (TMap const & lhs, TMap const & rhs)
{
    return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

int main()
{
    GenomicLocation gl1;
    gl1.chr = "chr1";
    gl1.start = 10000;
    gl1.end = 10050;
    gl1.strand = '+';

    GenomicLocation gl2;
    gl2.chr = "chr1";
    gl2.start = 10000;
    gl2.end = 10050;
    gl2.strand = '-';

    GenomicLocation gl3;
    gl3.chr = "chrX";
    gl3.start = 10000;
    gl3.end = 10050;
    gl3.strand = '+';

    GenomicLocation gl4;
    gl4.chr = "chr2";
    gl4.start = 10000;
    gl4.end = 10050;
    gl4.strand = '+';

    GenomicLocation gl5;
    gl5.chr = "chr2";
    gl5.start = 20000;
    gl5.end = 20050;
    gl5.strand = '+';

    GenomicLocation gl6;
    gl6.chr = "chr2";
    gl6.start = 30000;
    gl6.end = 30050;
    gl6.strand = '-';

    // miRNA regions
    miRNAIntTrees originalRegions;

    MiRNA mirna1;
    mirna1.name = "mirna1";
    mirna1.classification.insert("intergenic");
    MiRNA mirna2;
    mirna2.name = "mirna2";
    MiRNA mirna3;
    mirna3.name = "mirna3";
    mirna3.classification.insert("intergenic");
    MiRNA mirna4;
    mirna4.name = "mirna4";
    MiRNA mirna5;
    mirna5.name = "mirna5";
    mirna5.classification.insert("intergenic");
    MiRNA mirna6;
    mirna6.name = "mirna6";

    std::vector<std::pair<GenomicLocation, MiRNA> > originalRegionsBase;
    originalRegionsBase.push_back(std::make_pair(gl1, mirna1));
    originalRegionsBase.push_back(std::make_pair(gl2, mirna2));
    originalRegionsBase.push_back(std::make_pair(gl3, mirna3));
    originalRegionsBase.push_back(std::make_pair(gl4, mirna4));
    originalRegionsBase.push_back(std::make_pair(gl5, mirna5));
    originalRegionsBase.push_back(std::make_pair(gl6, mirna6));

    std::random_shuffle(originalRegionsBase.begin(), originalRegionsBase.end());

    // Insert elements into tree
    for (unsigned i = 0; i < originalRegionsBase.size(); i++)
    {
        std::pair<CharString, char> currentKey = std::make_pair(originalRegionsBase[i].first.chr, originalRegionsBase[i].first.strand);
        originalRegions[currentKey].insert(originalRegionsBase[i]);
    }

    // Test overlap (left, right, included inner/outer) -> also test for interval tree functionality
    // Test no overlap
    // Test overlap by one base
    // Test no overlap by one base
    // Test no overlap due to strand

    // Inner overlap mirna1
    GenomicLocation test1;
    test1.chr = "chr1";
    test1.start = 10045;
    test1.end = 10046;
    test1.strand = '+';

    // Outer overlap mirna1
    GenomicLocation test2;
    test2.chr = "chr1";
    test2.start = 9000;
    test2.end = 20000;
    test2.strand = '+';

    // Left overlap mirna4
    GenomicLocation test3;
    test3.chr = "chr2";
    test3.start = 9000;
    test3.end = 10010;
    test3.strand = '+';

    // Right overlap mirna5
    GenomicLocation test4;
    test4.chr = "chr2";
    test4.start = 20020;
    test4.end = 20060;
    test4.strand = '+';

    // Same overlap mirna1
    GenomicLocation test5;
    test5.chr = "chr1";
    test5.start = 10000;
    test5.end = 10050;
    test5.strand = '+';

    // Overlap by one base mirna2
    GenomicLocation test6;
    test6.chr = "chr1";
    test6.start = 9900;
    test6.end = 10001;
    test6.strand = '-';

    // No overlap
    GenomicLocation test7;
    test7.chr = "chrX";
    test7.start = 200000;
    test7.end = 300000;
    test7.strand = '+';

    // No overlap by one base
    GenomicLocation test8;
    test8.chr = "chr2";
    test8.start = 10050;
    test8.end = 20000;
    test8.strand = '+';

    // No overlap due to strand
    GenomicLocation test9;
    test9.chr = "chr2";
    test9.start = 30045;
    test9.end = 40000;
    test9.strand = '+';

    std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > overlapRegions;
    overlapRecord(overlapRegions, test1, 1, originalRegions);
    overlapRecord(overlapRegions, test2, 2, originalRegions);
    overlapRecord(overlapRegions, test2, 2, originalRegions);
    overlapRecord(overlapRegions, test3, 3, originalRegions);
    overlapRecord(overlapRegions, test4, 4, originalRegions);
    overlapRecord(overlapRegions, test5, 5, originalRegions);
    overlapRecord(overlapRegions, test6, 6, originalRegions);
    overlapRecord(overlapRegions, test7, 7, originalRegions);
    overlapRecord(overlapRegions, test8, 8, originalRegions);
    overlapRecord(overlapRegions, test9, 9, originalRegions);

    std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > controlRegions;
    controlRegions[std::make_pair("mirna1", test1)] = std::make_pair(1, true);
    controlRegions[std::make_pair("mirna1", test2)] = std::make_pair(4, true);
    controlRegions[std::make_pair("mirna4", test3)] = std::make_pair(3, false);
    controlRegions[std::make_pair("mirna5", test4)] = std::make_pair(4, true);
    controlRegions[std::make_pair("mirna1", test5)] = std::make_pair(5, true);
    controlRegions[std::make_pair("mirna2", test6)] = std::make_pair(6, false);

    if (!mapCompare(overlapRegions, controlRegions))
    {
        std::cerr << "ERROR: Differing elements in test and control miRNA overlap regions." << std::endl;
        return 1;
    }

    // Background regions
    backgroundIntTrees originalRegionsBG;
    std::vector<std::pair<GenomicLocation, Background> > originalRegionsBaseBG;

    Background testBackground;
    testBackground.overlaps = 0;

    originalRegionsBaseBG.push_back(std::make_pair(gl1, testBackground));
    originalRegionsBaseBG.push_back(std::make_pair(gl2, testBackground));
    originalRegionsBaseBG.push_back(std::make_pair(gl3, testBackground));
    originalRegionsBaseBG.push_back(std::make_pair(gl4, testBackground));
    originalRegionsBaseBG.push_back(std::make_pair(gl5, testBackground));
    originalRegionsBaseBG.push_back(std::make_pair(gl6, testBackground));

    std::random_shuffle(originalRegionsBaseBG.begin(), originalRegionsBaseBG.end());

    // Insert elements into tree
    for (unsigned i = 0; i < originalRegionsBaseBG.size(); i++)
    {
        std::pair<CharString, char> currentKey = std::make_pair(originalRegionsBaseBG[i].first.chr, originalRegionsBaseBG[i].first.strand);
        originalRegionsBG[currentKey].insert(originalRegionsBaseBG[i]);
    }

    std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > overlapRegionsBG;
    overlapRecord(overlapRegionsBG, test1, 1, originalRegionsBG);
    overlapRecord(overlapRegionsBG, test2, 2, originalRegionsBG);
    overlapRecord(overlapRegionsBG, test2, 2, originalRegionsBG);
    overlapRecord(overlapRegionsBG, test3, 3, originalRegionsBG);
    overlapRecord(overlapRegionsBG, test4, 4, originalRegionsBG);
    overlapRecord(overlapRegionsBG, test5, 5, originalRegionsBG);
    overlapRecord(overlapRegionsBG, test6, 6, originalRegionsBG);
    overlapRecord(overlapRegionsBG, test7, 7, originalRegionsBG);
    overlapRecord(overlapRegionsBG, test8, 8, originalRegionsBG);
    overlapRecord(overlapRegionsBG, test9, 9, originalRegionsBG);

    std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > controlRegionsBG;
    controlRegionsBG[std::make_pair("", gl1)] = std::make_pair(10, true);
    controlRegionsBG[std::make_pair("", gl4)] = std::make_pair(3, true);
    controlRegionsBG[std::make_pair("", gl5)] = std::make_pair(4, true);
    controlRegionsBG[std::make_pair("", gl2)] = std::make_pair(6, true);

    if (!mapCompare(overlapRegionsBG, controlRegionsBG))
    {
        std::cerr << "ERROR: Differing elements in test and control background overlap regions." << std::endl;
        return 1;
    }

    std::cout << "Test overlap successful." << std::endl;

    return 0;
}
