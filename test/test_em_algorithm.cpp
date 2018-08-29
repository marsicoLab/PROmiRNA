// ============================================================================
// Test EM algorithm
// ============================================================================
//
// Compiler flags needed: -std=c++17 -lstdc++fs -fopenmp -I /path/to/seqan/include/

#include <iostream>
#include <algorithm>
#include <vector>

#include <seqan/find.h>
#include <seqan/sequence.h>

#include "../../src/em_algorithm.h"
#include "../../src/structs.h"

using namespace seqan;

int main()
{
    // Read in beta parameters
    std::ifstream input;
    std::string line;
    std::vector<double> betaParams; // Read from file

    input.open("test_files/beta_parameters.tmp");
    if (input.is_open())
    {
        while(getline(input, line))
        {
            betaParams.push_back(stod(line));
        }
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open beta parameter file.");
        return 1;
    }
    betaParams.push_back(std::min(std::min(betaParams[1], betaParams[2]), betaParams[3]));

    // Read in dataset
    std::ifstream input2;
    std::vector<DatasetRecord> dataset; // Read from file

    input2.open("test_files/data_set.txt");
    if (input2.is_open())
    {
        while(getline(input2, line))
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
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open dataset file.");
        return 1;
    }

    std::ifstream input3;

    input3.open("test_files/lib_number.txt");
    if (input3.is_open())
    {
        unsigned count = 0;
        while(getline(input3, line))
        {
            dataset[count].libs = std::stoi(line);
            count++;
        }
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open lib number file.");
        return 1;
    }

    double mu1, mu2, lambda1, lambda2;
    unsigned const minLib = 1;
    unsigned const minReadCount = 0; // For comparison with old version that does not have this requirement

    estimateInitialParameters(mu1, mu2, lambda1, lambda2, dataset, minLib, minReadCount);

    std::vector<double> outputParams;
    std::vector<std::vector<double> > outputForMStep;
    emAlgorithm(outputParams, outputForMStep, dataset, mu1, mu2, lambda1, lambda2, betaParams);

    outputParams.insert(outputParams.end(), betaParams.begin(), betaParams.end()); // Concatenate for comparison

    const char * pythonOutputParams = "test_files/final_parameters.txt";
    const char * pythonOutput = "test_files/output_EM_complete.txt";

    std::ifstream controlParams;

    controlParams.open(pythonOutputParams);
    if (controlParams.is_open())
    {
        unsigned count = 0;
        while(getline(controlParams, line))
        {
            StringSet<CharString> currentRecord;
            strSplit(currentRecord, line, EqualsChar<':'>());
            if (round(outputParams[count] * 10000) / 10000 != round(std::stod(toCString(currentRecord[1])) * 10000) / 10000)
            {
                std::cerr << "ERROR: Final parameters are not the same between old and new PROmiRNA code."  << std::endl;
                return 1;
            }
            count++;
        }
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open parameter control file.");
        return 1;
    }

    std::ifstream control;
    control.open(pythonOutput);
    if (control.is_open())
    {
        unsigned count = 0;
        while(getline(control, line))
        {
            if (count != 0)
            {
                StringSet<CharString> currentRecord;
                strSplit(currentRecord, line, EqualsChar<'\t'>());
                if (round(outputForMStep[count - 1][4] * 10000) / 10000 != round(std::stod(toCString(currentRecord[5])) * 10000) / 10000 &&
                    round(outputForMStep[count - 1][0] * 10000) / 10000 != round(std::stod(toCString(currentRecord[6])) * 10000) / 10000 &&
                    round(outputForMStep[count - 1][1] * 10000) / 10000 != round(std::stod(toCString(currentRecord[7])) * 10000) / 10000)
                {
                    std::cerr << "ERROR: Final candidates are not the same between old and new PROmiRNA code."  << std::endl;
                    return 1;
                }
            }
            count++;
        }
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open EM output control file.");
        return 1;
    }

    std::cout << "Test of EM algorithm successful." << std::endl;

    return 0;
}
