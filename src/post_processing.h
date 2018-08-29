// ============================================================================
// Functions and routines for post-processing the EM algorithm output
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#pragma once

#include <iostream>
#include <fstream>

#include "CD-HIT/cdhit-common.h"
#include "CD-HIT/cdhit-utility.h"
#include "CD-HIT/cdhit-est.h"

#include "create_background.h"

/*!
 * @fn getFastaSeqs
 * @brief Extract FASTA sequences and IDs for final candidate sequences
 *
 * @signature void getFastaSeqs(fastaIdentifier, fastaLines, promoterIndex, faiIndex, dataset)
 *
 * @param[in,out]   fastaIdentifier             Identifiers of fasta sequences storing their position
 * @param[in,out]   fastaLines                  Extracted fasta sequences
 * @param[in]       promoterIndex               Index of regions classified as promoters
 * @param[in]       faiIndex                    FAI index of the reference genome
 * @param[in]       dataset                     Dataset for EM algorithm
*/
void getFastaSeqs(std::map<std::string, unsigned> & fastaIdentifier, std::vector<std::string> & fastaLines,
std::vector<unsigned> const & promoterIndex, FaiIndex & faiIndex, std::vector<DatasetRecord> const & dataset)
{
    fastaLines.reserve(2 * promoterIndex.size());
    for (auto i : promoterIndex)
    {
        Dna5String currentSeq;
        CharString indexChr;
        if (sequenceName(faiIndex, 0)[0] == 'c') // If reference genome has chromosome names starting with "chr"
            indexChr = dataset[i].location.chr;
        else
            indexChr = suffix(dataset[i].location.chr, 3);

        std::string currentId = std::string(">") + toCString(indexChr) + std::string(":") + std::to_string(dataset[i].location.start) + std::string("-") + std::to_string(dataset[i].location.end); // 0-based position, half-open interval
        extractSequenceFromIndex(currentSeq, faiIndex, indexChr, dataset[i].location.start, dataset[i].location.end);

        fastaLines.push_back(currentId);
        fastaLines.push_back(toCString((CharString) currentSeq));
        fastaIdentifier[currentId] = i;
    }
}

/*!
 * @fn reduceSequenceRedundancy
 * @brief Function wrapping CD-HIT in order to reduce sequence redundancy and find representative sequences for final output
 *
 * @signature void reduceSequenceRedundancy(fastaLines, threads)
 *
 * @param[in]       fastaLines                  Extracted fasta sequences
 * @param[in]       threads                     Number of threads to use
*/
void reduceSequenceRedundancy(std::vector<std::string> & fastaLines, unsigned const threads)
{
    const char * outputFile = "cdhit_output";
    cdhit::cdhitEst(fastaLines, outputFile, threads);
}

/*!
 * @fn readClusterOutput
 * @brief Parse the CD-HIT output file for final promoter prediction output
 *
 * @signature void readClusterOutput(selectedIdentifier, fastaIdentifier, dataset)
 *
 * @param[in,out]   selectedIdentifier          Identifiers of representative sequences selected by CD-HIT
 * @param[in,out]   fastaIdentifier             Identifiers of fasta sequences storing their position
 * @param[in]       dataset                     Dataset for EM algorithm
*/
void readClusterOutput(std::vector<std::string> & selectedIdentifier, std::map<std::string, unsigned> & fastaIdentifier, std::vector<DatasetRecord> const & dataset)
{
    const char * clusterPath = "cdhit_output.clstr";

    std::ifstream clusterFile;
    std::string line;
    clusterFile.open(clusterPath);

    // The representative promoter for each cluster is the one with the maximum tag count
    std::string representativeId;
    double maxTagCount = -1;

    if (clusterFile.is_open())
    {
        while(std::getline(clusterFile, line))
        {
            if (line.substr(0, 1) == ">")
            {
                if (maxTagCount != -1)
                {
                    selectedIdentifier.push_back(representativeId);
                    maxTagCount = -1;
                }
            }
            else
            {
                unsigned first = line.find(">");
                unsigned last = line.find("...");
                std::string currentId = line.substr(first, last - first);
                if (dataset[fastaIdentifier[currentId]].tagCount > maxTagCount)
                {
                    maxTagCount = dataset[fastaIdentifier[currentId]].tagCount;
                    representativeId = currentId;
                }
            }
        }
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open cluster file.");
    }
}

/*!
 * @fn postProcessing
 * @brief Post-processing of the EM output in order to obtain final list of predicted promoters
 *
 * @signature void postProcessing(outputPath, promoterProbs, dataset, oldPositionsTSSDataset, faiIndex, threads)
 *
 * @param[in]       outputPath                  Output path for final list of predicted promoters
 * @param[in]       promoterProbs               Predicted probabilities of being a promoter from EM output
 * @param[in]       dataset                     Dataset for EM algorithm
 * @param[in]       oldPositionsTSSDataset      Real positions of TSS (not 1000 bp core promoter)
 * @param[in]       faiIndex                    FAI index of the reference genome
 * @param[in]       threads                     Number of threads to use
*/
void postProcessing(CharString const & outputPath, std::vector<std::vector<double> > const & promoterProbs, std::vector<DatasetRecord> & dataset,
std::vector<std::pair<unsigned, unsigned> > const & oldPositionsTSSDataset, FaiIndex & faiIndex, unsigned const threads)
{
    // Get index of promoter entries
    std::vector<unsigned> promoterIndex;
    promoterIndex.reserve(promoterProbs.size());

    for (unsigned i = 0; i < promoterProbs.size(); i++)
    {
        if (promoterProbs[i][0] > promoterProbs[i][1] && !dataset[i].bg)
        {
            promoterIndex.push_back(i);
        }
    }

    std::cout << "Promoter for post-processing: " << promoterIndex.size() << std::endl;

    std::map<std::string, unsigned> fastaIdentifier;
    std::vector<std::string> fastaLines;
    getFastaSeqs(fastaIdentifier, fastaLines, promoterIndex, faiIndex, dataset);

    reduceSequenceRedundancy(fastaLines, threads);

    std::vector<std::string> selectedIdentifier;
    readClusterOutput(selectedIdentifier, fastaIdentifier, dataset);

    std::ofstream output;
    output.open(toCString(outputPath));
    if (output.is_open())
    {
        // First 6 columns can be read in as BED6 format for downstream analyses
        output << "#chr\t" << "start\t" << "stop\t" << "miRNA\t" << "tags\t" << "strand\t" << "distance\n";
        for (auto id : selectedIdentifier)
        {
            output << dataset[fastaIdentifier[id]].location.chr << "\t"
                   << oldPositionsTSSDataset[fastaIdentifier[id]].first << "\t"
                   << oldPositionsTSSDataset[fastaIdentifier[id]].second << "\t"
                   << dataset[fastaIdentifier[id]].miRNA << "\t"
                   << dataset[fastaIdentifier[id]].tagCount << "\t"
                   << dataset[fastaIdentifier[id]].location.strand << "\t"
                   << dataset[fastaIdentifier[id]].distance << "\n";
        }
        output.close();
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open output file.");
    }
}
