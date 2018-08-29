// ============================================================================
// Functions for building matrix and merging CAGE tags
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)
#pragma once

#include <iterator>
#include <string>
#include <vector>

#include <seqan/basic.h>

#include "overlap.h"
#include "structs.h"
#include "quantile_normalization.h"

using namespace seqan;

/*!
 * @fn buildMatrix
 * @brief Build overlap matrix for miRNAs/background and BED files and perform quantile normalization
 *
 * @signature void buildMatrix(normMatrix, countMatrix, uniqueRegions, overlapRegions, files)
 *
 * @param[in,out]   normMatrix          Resulting normalized overlap matrix
 * @param[in,out]   countMatrix         Resulting overlap matrix
 * @param[in,out]   uniqueRegions       Unique regions across CAGE files from overlap
 * @param[in]       overlapRegions      Vector of miRNA overlap with CAGE data
 * @param[in]       files               Vector of all BED files from a directory
 */
void buildMatrix(MatrixSingle & normMatrix, MatrixPair & countMatrix,
std::map<std::pair<CharString, GenomicLocation>, unsigned> & uniqueRegions,
std::vector<std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > > const & overlapRegions,
std::vector<std::string> const & files)
{
    // Hash map with all unique miRNA+regions from all CAGE BED files
    unsigned count = 0;
    for (unsigned i = 0; i < overlapRegions.size(); i++)
    {
        for (auto const & [key, val] : overlapRegions[i])
        {
            auto ret = uniqueRegions.insert(std::make_pair(key, count));

            if (ret.second)
            {
                count++;
            }
        }
    }

    if (uniqueRegions.size() == 0)
    {
        throw std::runtime_error("ERROR: No overlap found between CAGE tags and miRNAs/background. Rerun with more/other CAGE data.");
        return;
    }

    // Construct count matrix
    countMatrix.colnames = files;
    countMatrix.rownames.resize(uniqueRegions.size());
    countMatrix.matrix.resize(uniqueRegions.size());
    for (unsigned i = 0; i < countMatrix.matrix.size(); i++)
        countMatrix.matrix[i].resize(files.size());

    for (auto const & [key, val] : uniqueRegions)
    {
        if (key.first == "")
        {
            // Background
            countMatrix.rownames[val] = toCString(key.second.chr) + std::string("-") + std::to_string(key.second.start) + std::string("-") + std::to_string(key.second.end) + std::string("-") + key.second.strand;
        }
        else
        {
            // miRNA
            countMatrix.rownames[val] = toCString(key.second.chr) + std::string("-") + std::to_string(key.second.start) + std::string("-") + std::to_string(key.second.end) + std::string("-") + toCString(key.first);
        }

        for (unsigned i = 0; i < overlapRegions.size(); i++)
        {
            auto mit = overlapRegions[i].find(key);
            if (mit != overlapRegions[i].end())
            {
                countMatrix.matrix[val][i].first = mit->second.first;
            }
            else
            {
                countMatrix.matrix[val][i].first = 0;
            }
            countMatrix.matrix[val][i].second = val;
        }
    }

    // Quantile normalization
    quantileNormalization(normMatrix, countMatrix);
}
