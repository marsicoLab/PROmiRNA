// ============================================================================
// Test merging of bg sequences
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#include <gtest/gtest.h>
#include <iostream>

#include "../src/mirna_classification.h"

void traverseIntTreeCompareMirnas(IntTreeNode<std::pair<GenomicLocation, MiRNA> > *root, std::map<GenomicLocation, MiRNA> const & controlRecords)
{
    if (root == NULL) return;
    traverseIntTreeCompareMirnas(root->left, controlRecords);

    auto mit = controlRecords.find((root->nodeInterval).first);
    EXPECT_NE(mit, controlRecords.end());

    EXPECT_EQ((root->nodeInterval).second, mit->second);

    traverseIntTreeCompareMirnas(root->right, controlRecords);
}

TEST(mirna_classification, classification)
{
    MiRNA mirna1;
    mirna1.classification.insert("intron");
    mirna1.id = "86506";
    mirna1.accession = "MI0022705";
    mirna1.name = "hsa-mir-6859-1";

    GenomicLocation gl1;
    gl1.chr = "chr1";
    gl1.start = 17435;
    gl1.end = 67435;
    gl1.strand = '-';

    MiRNA mirna2;
    mirna2.classification.insert("intron");
    mirna2.classification.insert("exon");
    mirna2.id = "70011";
    mirna2.accession = "MI0006363";
    mirna2.name = "hsa-mir-1302-2";

    GenomicLocation gl2;
    gl2.chr = "chr1";
    gl2.start = 0;
    gl2.end = 30365;
    gl2.strand = '+';

    MiRNA mirna3;
    mirna3.classification.insert("intron");
    mirna3.id = "84818";
    mirna3.accession = "MI0022558";
    mirna3.name = "hsa-mir-6723";

    GenomicLocation gl3;
    gl3.chr = "chr1";
    gl3.start = 567792;
    gl3.end = 617792;
    gl3.strand = '-';

    MiRNA mirna4;
    mirna4.classification.insert("intergenic");
    mirna4.id = "80708";
    mirna4.accession = "MI0017330";
    mirna4.name = "hsa-mir-4697";

    GenomicLocation gl4;
    gl4.chr = "chr11";
    gl4.start = 133768475;
    gl4.end = 133818475;
    gl4.strand = '-';

    MiRNA mirna5;
    mirna5.classification.insert("intron");
    mirna5.id = "79434";
    mirna5.accession = "MI0016049";
    mirna5.name = "hsa-mir-3649";

    GenomicLocation gl5;
    gl5.chr = "chr12";
    gl5.start = 1769545;
    gl5.end = 1819545;
    gl5.strand = '-';

    MiRNA mirna6;
    mirna6.classification.insert("intron");
    mirna6.id = "68198";
    mirna6.accession = "MI0003675";
    mirna6.name = "hsa-mir-411";

    GenomicLocation gl6;
    gl6.chr = "chr14";
    gl6.start = 101439661;
    gl6.end = 101489661;
    gl6.strand = '+';

    MiRNA mirna7;
    mirna7.classification.insert("intron");
    mirna7.id = "65390";
    mirna7.accession = "MI0000744";
    mirna7.name = "hsa-mir-299";

    GenomicLocation gl7;
    gl7.chr = "chr14";
    gl7.start = 101440130;
    gl7.end = 101490130;
    gl7.strand = '+';

    MiRNA mirna8;
    mirna8.classification.insert("intron");
    mirna8.id = "65431";
    mirna8.accession = "MI0000788";
    mirna8.name = "hsa-mir-380";

    GenomicLocation gl8;
    gl8.chr = "chr14";
    gl8.start = 101441353;
    gl8.end = 101491353;
    gl8.strand = '+';

    MiRNA mirna9;
    mirna9.classification.insert("intergenic");
    mirna9.id = "86495";
    mirna9.accession = "MI0022694";
    mirna9.name = "hsa-mir-6848";

    GenomicLocation gl9;
    gl9.chr = "chr8";
    gl9.start = 145540977;
    gl9.end = 145590977;
    gl9.strand = '-';

    MiRNA mirna10;
    mirna10.classification.insert("intergenic");
    mirna10.id = "69208";
    mirna10.accession = "MI0005530";
    mirna10.name = "hsa-mir-509-2";

    GenomicLocation gl10;
    gl10.chr = "chrX";
    gl10.start = 146340367;
    gl10.end = 146390367;
    gl10.strand = '-';

    std::map<GenomicLocation, MiRNA> controlRecords;
    controlRecords[gl1] = mirna1;
    controlRecords[gl2] = mirna2;
    controlRecords[gl3] = mirna3;
    controlRecords[gl4] = mirna4;
    controlRecords[gl5] = mirna5;
    controlRecords[gl6] = mirna6;
    controlRecords[gl7] = mirna7;
    controlRecords[gl8] = mirna8;
    controlRecords[gl9] = mirna9;
    controlRecords[gl10] = mirna10;

    miRNAIntTrees classifiedMirnas;
    std::map<CharString, unsigned> mirnaStarts;
    std::map<CharString, CharString> mirnaExonChecks;
    CharString const gffPath = "../test/test_files/test_hsa.gff";
    CharString const mirnaPath = "../external_data/mirna.txt";
    CharString const mirnaContextPath = "../external_data/mirna_context.txt";

    mirnaClassification(classifiedMirnas, mirnaStarts, mirnaExonChecks, gffPath, mirnaPath, mirnaContextPath);

    for (auto it = classifiedMirnas.begin(); it != classifiedMirnas.end(); it++)
    {
        traverseIntTreeCompareMirnas((it->second).root, controlRecords);
    }
}
