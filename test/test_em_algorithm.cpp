// ============================================================================
// Test EM algorithm
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#include <gtest/gtest.h>
#include <iostream>
#include <algorithm>
#include <vector>

#include <seqan/find.h>
#include <seqan/sequence.h>

#include "../src/create_dataset.h"
#include "../src/em_algorithm.h"
#include "../src/structs.h"

using namespace seqan;

TEST(em_algorithm_test, em_algorithm)
{
    // Read in beta parameters
    // Control values have been computed with the python version, default values and external_data from website
    std::vector<double> betaParams{-10.3280678842784, 15.7354397363585, 26.9956132389546, 6.09039873451161};

    // Read in dataset
    std::ifstream input1;
    std::string line;
    std::vector<DatasetRecord> dataset; // Read from file

    input1.open("../test/test_files/data_set.txt");
    while(getline(input1, line))
    {
        StringSet<CharString> currentRecord;
        strSplit(currentRecord, line, EqualsChar<'\t'>());
        DatasetRecord datasetRecord;
        datasetRecord.location.chr = std::string("chr") + toCString(currentRecord[2]);
        datasetRecord.location.start = std::stoi(toCString(currentRecord[3])) - 1;
        datasetRecord.location.end = std::stoi(toCString(currentRecord[4]));
        datasetRecord.location.strand = currentRecord[5][0];
        datasetRecord.tagCount = std::stod(toCString(currentRecord[6]));
        datasetRecord.cpgContent = std::stod(toCString(currentRecord[7]));
        datasetRecord.conservation = std::stod(toCString(currentRecord[8]));
        datasetRecord.affinity = std::stod(toCString(currentRecord[9]));
        datasetRecord.distScore = std::stod(toCString(currentRecord[10]));

        if (currentRecord[11] == "background")
            datasetRecord.bg = true;
        else
            datasetRecord.bg = false;
        dataset.push_back(datasetRecord);
    }

    std::ifstream input2;
    input2.open("../test/test_files/lib_number.txt");
    unsigned count = 0;
    while(getline(input2, line))
    {
        dataset[count].libs = std::stoi(line);
        count++;
    }

    // Initial parameters are now randomized and not taken from the front as in python version
    // In order to test the EM algorithm parameters are therefore taken from the python version
    double mu1 = 13.9653333333;
    double mu2 = 0.639542275906;
    double lambda1 = 1.93249880225;
    double lambda2 = 0.295475492659;

    std::vector<unsigned> indexValidSamplesTSS;
    std::vector<unsigned> indexValidSamplesBG;

    std::vector<double> outputParams;
    std::vector<std::vector<double> > outputForMStep;
    emAlgorithm(outputParams, outputForMStep, dataset, mu1, mu2, lambda1, lambda2, betaParams);

    // We compare if output parameters are approximately the same between python version and new version
    // Bugs in the EM fixed make it impossible to compare exact output predictions
    double mu1FinalControl = 1.47043371445;
    double mu2FinalControl = 0.933621533169;
    double lambda1FinalControl = 6.2091265777;
    double lambda2FinalControl = 16.5356541586;

    // Make max difference depending on order of magnitude
    std::vector<double> comparisonMaxDifference;
    for (auto out : outputParams)
    {
        if (round(log10(out)) < 1)
            comparisonMaxDifference.push_back(0.05);
        else
            comparisonMaxDifference.push_back(0.5);
    }

    EXPECT_NEAR(outputParams[0], mu1FinalControl, comparisonMaxDifference[0]);
    EXPECT_NEAR(outputParams[1], mu2FinalControl, comparisonMaxDifference[1]);
    EXPECT_NEAR(outputParams[2], lambda1FinalControl, comparisonMaxDifference[2]);
    EXPECT_NEAR(outputParams[3], lambda2FinalControl, comparisonMaxDifference[3]);
}
