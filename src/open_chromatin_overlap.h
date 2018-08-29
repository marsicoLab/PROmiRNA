// ============================================================================
// Determining high confidence promoters based on overlap with DHSs
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#pragma once

#include <iostream>

#include <seqan/find.h>

#include "interval_tree.h"

using namespace seqan;

// A map of tree structures that stores information about DNase Hypersensitive sites (DHSs)
// Key is the chromosome because the interval trees can only process simple intervals [start, stop)
// Value is a tree of genomic locations because no other information is needed
typedef std::map<CharString, IntTree<GenomicLocation> > dhsTree;

/*!
 * @fn overlapDHS
 * @brief Overlap selected promoters for output with DHSs such that high confidence promoters can be determined
 *
 * @signature void postProcessing(outputPath, promoterProbs, dataset, oldPositionsTSSDataset, faiIndex, threads)
 *
 * @param[in]       selectedIdentifier          Identifiers of representative sequences selected by CD-HIT
 * @param[in]       fastaIdentifier             Identifiers of fasta sequences storing their position
 * @param[in]       oldPositionsTSSDataset      Real positions of TSS (not 1000 bp core promoter)
 * @param[in]       dataset                     Dataset for EM algorithm
 * @param[in]       bedFile                     The BED file containing positional information of DHSs
 *
 * Overlap with DHSs means open chromatin and therefore a higher confidence for actually being a promoter.
*/
void overlapDHS(std::vector<std::string> & selectedIdentifier, std::map<std::string, unsigned> & fastaIdentifier,
std::vector<std::pair<unsigned, unsigned> > const & oldPositionsTSSDataset, std::vector<DatasetRecord> & dataset, CharString const & bedFile)
{
    std::vector<GenomicLocation> dhsVector;
    dhsTree dhs;

    std::ifstream input;
    std::string line;

    input.open(toCString(bedFile));
	if (input.is_open())
    {
        while(getline(input, line))
        {
            StringSet<CharString> currentRecord;
            strSplit(currentRecord, line, EqualsChar<'\t'>());

            GenomicLocation gl;
            gl.chr = currentRecord[0];
            gl.start = std::atoi(toCString(currentRecord[1]));
            gl.end = std::atoi(toCString(currentRecord[2]));

            dhsVector.push_back(gl);
        }
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open DHS file.");
        return;
    }

    std::random_shuffle(dhsVector.begin(), dhsVector.end());

    for (unsigned i = 0; i < dhsVector.size(); i++)
	{
		dhs[dhsVector[i].chr].insert(dhsVector[i]);
	}

    std::vector<std::string> highConfidenceIdentifier;
    highConfidenceIdentifier.reserve(selectedIdentifier.size());

    for (auto id : selectedIdentifier)
    {
        GenomicLocation * overlap = dhs[dataset[fastaIdentifier[id]].location.chr].overlapSearchSingle(dataset[fastaIdentifier[id]].location);
        if (overlap != NULL)
        {
            highConfidenceIdentifier.push_back(id);
        }
    }

    std::ofstream output;
    output.open("high_confidence_predicted_promoters.txt");
    if (output.is_open())
    {
        // First 6 columns can be read in as BED6 format for downstream analyses
        output << "#chr\t" << "start\t" << "stop\t" << "miRNA\t" << "tags\t" << "strand\t" << "distance\n";
        for (auto id : highConfidenceIdentifier)
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
