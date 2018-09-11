// ============================================================================
// Test build matrix function
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#include <gtest/gtest.h>
#include <iostream>

#include "../src/build_matrix.h"

TEST(build_matrix_test, no_overlap_background)
{
    std::vector<std::string> files{"fileA.bed", "fileB.bed", "fileC.bed"};

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

    std::pair<unsigned, bool> testValue = std::make_pair(1, true);

    // Test all files have different regions (background)
    std::vector<std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > > overlapRegions1;
    overlapRegions1.resize(3);
    overlapRegions1[0].insert(std::make_pair(std::make_pair("", gl1), testValue));
    overlapRegions1[0].insert(std::make_pair(std::make_pair("", gl2), testValue));
    overlapRegions1[1].insert(std::make_pair(std::make_pair("", gl3), testValue));
    overlapRegions1[1].insert(std::make_pair(std::make_pair("", gl4), testValue));
    overlapRegions1[2].insert(std::make_pair(std::make_pair("", gl5), testValue));
    overlapRegions1[2].insert(std::make_pair(std::make_pair("", gl6), testValue));

    MatrixSingle normMatrix1;
    MatrixPair countMatrix1;
    std::map<std::pair<CharString, GenomicLocation>, unsigned> uniqueRegions1;
    buildMatrix(normMatrix1, countMatrix1, uniqueRegions1, overlapRegions1, files);

    // NOTE Order of matrix is in order of overlap vectors/maps
    MatrixPair controlMatrix1;
    controlMatrix1.matrix = {{std::make_pair(1, 0), std::make_pair(0, 0), std::make_pair(0, 0)},
                             {std::make_pair(1, 1), std::make_pair(0, 1), std::make_pair(0, 1)},
                             {std::make_pair(0, 2), std::make_pair(1, 2), std::make_pair(0, 2)},
                             {std::make_pair(0, 3), std::make_pair(1, 3), std::make_pair(0, 3)},
                             {std::make_pair(0, 4), std::make_pair(0, 4), std::make_pair(1, 4)},
                             {std::make_pair(0, 5), std::make_pair(0, 5), std::make_pair(1, 5)}};

    controlMatrix1.rownames = {"chr1-10000-10050-+", "chr1-10000-10050--", "chr2-10000-10050-+", "chrX-10000-10050-+", "chr2-20000-20050-+", "chr2-30000-30050--"};
    controlMatrix1.colnames = {"fileA.bed", "fileB.bed", "fileC.bed"};

    EXPECT_EQ(countMatrix1, controlMatrix1);
}

TEST(build_matrix_test, no_overlap_tss)
{
    std::vector<std::string> files{"fileA.bed", "fileB.bed", "fileC.bed"};

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

    std::pair<unsigned, bool> testValue = std::make_pair(1, true);

    // Test all files have different regions (TSS)
    std::vector<std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > > overlapRegions2;
    overlapRegions2.resize(3);
    overlapRegions2[0].insert(std::make_pair(std::make_pair("mirna1", gl1), testValue));
    overlapRegions2[0].insert(std::make_pair(std::make_pair("mirna1", gl2), testValue));
    overlapRegions2[1].insert(std::make_pair(std::make_pair("mirna2", gl1), testValue));
    overlapRegions2[1].insert(std::make_pair(std::make_pair("mirna2", gl4), testValue));
    overlapRegions2[2].insert(std::make_pair(std::make_pair("mirna3", gl2), testValue));
    overlapRegions2[2].insert(std::make_pair(std::make_pair("mirna3", gl6), testValue));

    MatrixSingle normMatrix2;
    MatrixPair countMatrix2;
    std::map<std::pair<CharString, GenomicLocation>, unsigned> uniqueRegions2;
    buildMatrix(normMatrix2, countMatrix2, uniqueRegions2, overlapRegions2, files);

    // NOTE Order of matrix is in order of overlap vectors/maps
    MatrixPair controlMatrix2;
    controlMatrix2.matrix = {{std::make_pair(1, 0), std::make_pair(0, 0), std::make_pair(0, 0)},
                             {std::make_pair(1, 1), std::make_pair(0, 1), std::make_pair(0, 1)},
                             {std::make_pair(0, 2), std::make_pair(1, 2), std::make_pair(0, 2)},
                             {std::make_pair(0, 3), std::make_pair(1, 3), std::make_pair(0, 3)},
                             {std::make_pair(0, 4), std::make_pair(0, 4), std::make_pair(1, 4)},
                             {std::make_pair(0, 5), std::make_pair(0, 5), std::make_pair(1, 5)}};
    controlMatrix2.rownames = {"chr1-10000-10050-mirna1", "chr1-10000-10050-mirna1", "chr1-10000-10050-mirna2", "chr2-10000-10050-mirna2", "chr1-10000-10050-mirna3", "chr2-30000-30050-mirna3"};
    controlMatrix2.colnames = {"fileA.bed", "fileB.bed", "fileC.bed"};

    EXPECT_EQ(countMatrix2, controlMatrix2);
}

TEST(build_matrix_test, overlap_background)
{
    std::vector<std::string> files{"fileA.bed", "fileB.bed", "fileC.bed"};

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

    std::pair<unsigned, bool> testValue = std::make_pair(1, true);

    // Test duplicated regions in files (background)
    std::vector<std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > > overlapRegions3;
    overlapRegions3.resize(3);
    overlapRegions3[0].insert(std::make_pair(std::make_pair("", gl1), testValue));
    overlapRegions3[0].insert(std::make_pair(std::make_pair("", gl4), testValue));
    overlapRegions3[0].insert(std::make_pair(std::make_pair("", gl2), testValue));
    overlapRegions3[1].insert(std::make_pair(std::make_pair("", gl3), testValue));
    overlapRegions3[2].insert(std::make_pair(std::make_pair("", gl1), testValue));
    overlapRegions3[2].insert(std::make_pair(std::make_pair("", gl4), testValue));

    MatrixSingle normMatrix3;
    MatrixPair countMatrix3;
    std::map<std::pair<CharString, GenomicLocation>, unsigned> uniqueRegions3;
    buildMatrix(normMatrix3, countMatrix3, uniqueRegions3, overlapRegions3, files);

    // NOTE Order of matrix is in order of overlap vectors/maps
    MatrixPair controlMatrix3;
    controlMatrix3.matrix = {{std::make_pair(1, 0), std::make_pair(0, 0), std::make_pair(1, 0)},
                             {std::make_pair(1, 1), std::make_pair(0, 1), std::make_pair(1, 1)},
                             {std::make_pair(1, 2), std::make_pair(0, 2), std::make_pair(0, 2)},
                             {std::make_pair(0, 3), std::make_pair(1, 3), std::make_pair(0, 3)}};
    controlMatrix3.rownames = {"chr1-10000-10050-+", "chr2-10000-10050-+", "chr1-10000-10050--", "chrX-10000-10050-+"};
    controlMatrix3.colnames = {"fileA.bed", "fileB.bed", "fileC.bed"};

    EXPECT_EQ(countMatrix3, controlMatrix3);
}

TEST(build_matrix_test, overlap_tss)
{
    std::vector<std::string> files{"fileA.bed", "fileB.bed", "fileC.bed"};

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

    std::pair<unsigned, bool> testValue = std::make_pair(1, true);

    // Test duplicated regions in files (TSS)
    std::vector<std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > > overlapRegions4;
    overlapRegions4.resize(3);
    overlapRegions4[0].insert(std::make_pair(std::make_pair("mirna1", gl1), testValue));
    overlapRegions4[0].insert(std::make_pair(std::make_pair("mirna2", gl2), testValue));
    overlapRegions4[0].insert(std::make_pair(std::make_pair("mirna3", gl3), testValue));
    overlapRegions4[1].insert(std::make_pair(std::make_pair("mirna1", gl1), testValue));
    overlapRegions4[2].insert(std::make_pair(std::make_pair("mirna2", gl2), testValue));
    overlapRegions4[2].insert(std::make_pair(std::make_pair("mirna4", gl3), testValue));

    MatrixSingle normMatrix4;
    MatrixPair countMatrix4;
    std::map<std::pair<CharString, GenomicLocation>, unsigned> uniqueRegions4;
    buildMatrix(normMatrix4, countMatrix4, uniqueRegions4, overlapRegions4, files);

    // NOTE Order of matrix is in order of overlap vectors/maps
    MatrixPair controlMatrix4;
    controlMatrix4.matrix = {{std::make_pair(1, 0), std::make_pair(1, 0), std::make_pair(0, 0)},
                             {std::make_pair(1, 1), std::make_pair(0, 1), std::make_pair(1, 1)},
                             {std::make_pair(1, 2), std::make_pair(0, 2), std::make_pair(0, 2)},
                             {std::make_pair(0, 3), std::make_pair(0, 3), std::make_pair(1, 3)}};
    controlMatrix4.rownames = {"chr1-10000-10050-mirna1", "chr1-10000-10050-mirna2", "chrX-10000-10050-mirna3", "chrX-10000-10050-mirna4"};
    controlMatrix4.colnames = {"fileA.bed", "fileB.bed", "fileC.bed"};

    EXPECT_EQ(countMatrix4, controlMatrix4);
}

TEST(build_matrix_test, empty_file)
{
    std::vector<std::string> files{"fileA.bed", "fileB.bed", "fileC.bed"};

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

    std::pair<unsigned, bool> testValue = std::make_pair(1, true);

    // Test that files are empty
    std::vector<std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > > overlapRegions5;
    overlapRegions5.resize(3);
    overlapRegions5[0].insert(std::make_pair(std::make_pair("", gl1), testValue));
    overlapRegions5[0].insert(std::make_pair(std::make_pair("", gl2), testValue));
    overlapRegions5[1].insert(std::make_pair(std::make_pair("", gl1), testValue));
    overlapRegions5[1].insert(std::make_pair(std::make_pair("", gl4), testValue));

    MatrixSingle normMatrix5;
    MatrixPair countMatrix5;
    std::map<std::pair<CharString, GenomicLocation>, unsigned> uniqueRegions5;
    buildMatrix(normMatrix5, countMatrix5, uniqueRegions5, overlapRegions5, files);

    std::cout << countMatrix5.matrix[0][0].first << std::endl;

    // NOTE Order of matrix is in order of overlap vectors/maps
    MatrixPair controlMatrix5;
    controlMatrix5.matrix = {{std::make_pair(1, 0), std::make_pair(1, 0), std::make_pair(0, 0)},
                             {std::make_pair(1, 1), std::make_pair(0, 1), std::make_pair(0, 1)},
                             {std::make_pair(0, 2), std::make_pair(1, 2), std::make_pair(0, 2)}};
    controlMatrix5.rownames = {"chr1-10000-10050-+", "chr1-10000-10050--", "chr2-10000-10050-+"};
    controlMatrix5.colnames = {"fileA.bed", "fileB.bed", "fileC.bed"};

    EXPECT_EQ(countMatrix5, controlMatrix5);
}
