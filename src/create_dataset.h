// ============================================================================
// Functions for creating dataset
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <time.h>

#include <cmath>

#include <seqan/find.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include "create_background.h"
#include "structs.h"
#include "TRAP/annotate.h"
#include "conservation.h"

using namespace seqan;

/*!
 * @fn calculateCpGContent
 * @brief Calculates CpG content of a single sequence
 *
 * @signature void calculateCpGContent(sequence)
 *
 * @param[in]       sequence        The DNA sequence that should be scanned
 */
double calculateCpGContent(Dna5String & sequence)
{
    unsigned doubleCount = 0;
    unsigned singleCount = 0;

    for (unsigned i = 0; i < length(sequence); i++)
    {
        if (sequence[i] == 'C')
        {
            if (i < (length(sequence) - 1) && sequence[i + 1] == 'G')
            {
                doubleCount++;
            }
            singleCount++;
        }
        else if (sequence[i] == 'G')
        {
            singleCount++;
        }
    }

    // Calculate the CpG content
    // Avoid divison by zero errors and set cpgContent to zero manually if dividers are not greater than zero
    if (singleCount > 0 && length(sequence) > 0)
    {
        return ((double) doubleCount / (double) length(sequence)) /  pow(((double) singleCount / ((double) 2 * length(sequence))), 2);
    }
    else
    {
        return 0;
    }
}

/*!
 * @fn gffToFasta
 * @brief Get FASTA sequences from a GFF file, write sequences and compute CpG content
 *
 * @signature void gffToFasta(dataset, fastaPath, indexPath)
 *
 * @param[in,out]   dataset         Vector of dataset records
 * @param[in]       fastaPath       Path for output FASTA
 * @param[in]       indexPath       Path to the genome that should be indexed or where FAI index
 *                                  with the same name lies in the same directory
 *
 * @throw Exception if GFF file or index cannot be opened.
 */
void gffToFasta(std::vector<DatasetRecord> & dataset, const char * fastaPath, FaiIndex & faiIndex)
{
    // Open output file
	SeqFileOut seqOut;

    if (!open(seqOut, fastaPath))
    {
        throw std::runtime_error("ERROR: Could not open FASTA output file.");
    }

    // Save FASTA IDs and sequences that are extracted from the index in StringSets that are written altogether to the output file
    // NOTE: SeqAn translates 1-based positions and closed intervals internally into zero-based, half-open intervals
    //       This is also expected by the index, so positions don't need to be adapted
    // In additon for each sequence the CpG content is calculated at the same time, so that sequences don't have to be read in again
    StringSet<CharString> fastaIds;
    resize(fastaIds, dataset.size());
    StringSet<Dna5String> fastaSeqs;
    resize(fastaSeqs, dataset.size());

    for (unsigned i = 0; i < dataset.size(); i++)
    {
        CharString currentId = toCString(dataset[i].location.chr) + std::string(":") + std::to_string(dataset[i].location.start + 1) + std::string("-") + std::to_string(dataset[i].location.end); // Same as in bedtools, change 0-based to 1-based again
        assignValue(fastaIds, i, currentId);

        Dna5String currentSeq;
        CharString indexChr;
        if (sequenceName(faiIndex, 0)[0] == 'c') // If reference genome has chromosome names starting with "chr"
            indexChr = dataset[i].location.chr;
        else
            indexChr = suffix(dataset[i].location.chr, 3);

        extractSequenceFromIndex(currentSeq, faiIndex, indexChr, dataset[i].location.start, dataset[i].location.end);
        if (dataset[i].location.strand == '-')
            reverseComplement(currentSeq);

        assignValue(fastaSeqs, i, currentSeq);
        dataset[i].cpgContent = calculateCpGContent(currentSeq);
    }

    try
    {
        writeRecords(seqOut, fastaIds, fastaSeqs);
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return;
    }
}

/*!
 * @fn createDataset
 * @brief Create the dataset for the EM algorithm
 *
 * @signature void createDataset(dataset, psemPath, wigPath, faiIndex)
 *
 * @param[in,out]   dataset         Dataset for EM algorithm
 * @param[in]       psemPath        Path to output FASTA file
 * @param[in]       wigPath         Path to WIG file
 * @param[in]       faiIndex        FAI index of the reference genome to extract sequences for TRAP and CpG content calculation
 *
 * Records of dataset have positional information (GenomicLocation), tag count, CpG content, conservation scores, affinity, lib number and whether
 * they are background regions or not. Conservation is extracted and averaged by hgWiggle. Affinity is computed by TRAP.
 */
void createDataset(std::vector<DatasetRecord> & dataset, CharString const psemPath, CharString wigPath, FaiIndex & faiIndex)
{
    // Get FASTA sequences and calculate CpG content
    // FASTA sequences needed for TRAP - can we avoid writing them to a file?
    const char * fastaPath = "temp_seqs.fa";
    gffToFasta(dataset, fastaPath, faiIndex);

    std::string psemPathForTrap = toCString(psemPath);
    std::vector<double> affinities;

    CharString pref = prefix(wigPath, length(wigPath) - 4);
    CharString wibPath = toCString(pref) + std::string(".wib");
     std::vector<double> conservationScores(dataset.size(), 0.0);

    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {
            trap::annotate(affinities, const_cast<char*>(fastaPath), const_cast<char*>(psemPathForTrap.c_str()), 0.5);
        }

        #pragma omp section
        {
            getAverageConservation(conservationScores, dataset, wigPath, wibPath);
        }
    }

    for (unsigned i = 0; i < affinities.size(); i++)
    {
        dataset[i].affinity = affinities[i];
        dataset[i].conservation = conservationScores[i];
    }
}

/*!
 * @fn getValidSamplesForParameterEstimation
 * @brief Get all indices that are valid for both beta parameter and initial parameter estimation
 *
 * @signature void getValidSamplesForParameterEstimation(indexValidSamplesTSS, indexValidSamplesBG, dataset, backgroundStart, minLib)
 *
 * @param[in,out]   indexValidSamplesTSS        Index of samples in dataset that pass strongly indicated TSS requirements
 * @param[in,out]   indexValidSamplesBG         Index of samples in dataset that pass strongly indicated background requirements
 * @param[in]       dataset                     Dataset for EM algorithm
 * @param[in]       backgroundStart             Index of start of background regions in dataset vector
 * @param[in]       minLib                      Minimum number of libraries required
 */
void getValidSamplesForParameterEstimation(std::vector<unsigned> & indexValidSamplesTSS, std::vector<unsigned> & indexValidSamplesBG,
std::vector<DatasetRecord> const & dataset, unsigned const backgroundStart, unsigned const minLib)
{
    indexValidSamplesTSS.reserve(backgroundStart);
    indexValidSamplesBG.reserve(dataset.size() - backgroundStart);

    for (unsigned i = 0; i < dataset.size(); i++)
    {
        if (!dataset[i].bg && (dataset[i].libs >= minLib))
        {
            indexValidSamplesTSS.push_back(i);
        }
        else if (dataset[i].bg && (dataset[i].cpgContent <= 0.5) && (dataset[i].conservation <= 0.2) && (dataset[i].affinity <= 0.1))
        {
            indexValidSamplesBG.push_back(i);
        }
    }

    if (indexValidSamplesTSS.size() < 2000 || indexValidSamplesBG.size() < 2000)
    {
        throw std::runtime_error("ERROR: Not enough characteristic regions for classification present. Please provide more CAGE data.");
    }

    std::random_shuffle(indexValidSamplesTSS.begin(), indexValidSamplesTSS.end());
    std::random_shuffle(indexValidSamplesBG.begin(), indexValidSamplesBG.end());
}

/*!
 * @fn getExamplesForBetaParams
 * @brief Extract examples for beta parameter estimation
 *
 * @signature void getExamplesForBetaParams(classificationExamples, variablesExamples, dataset, indexValidSamplesTSS, indexValidSamplesBG)
 *
 * @param[in,out]   classificationExamples      Class of selected examples for beta parameter estimation
 * @param[in,out]   variablesExamples           Variables (X matrix) of selected examples for beta parameter estimation
 * @param[in]       dataset                     Dataset for EM algorithm
 * @param[in]       indexValidSamplesTSS        Index of samples in dataset that pass strongly indicated TSS requirements
 * @param[in]       indexValidSamplesBG         Index of samples in dataset that pass strongly indicated background requirements
*/
void getExamplesForBetaParams(std::vector<unsigned> & classificationExamples, std::vector<std::vector<double> > & variablesExamples,
std::vector<DatasetRecord> const & dataset, std::vector<unsigned> const & indexValidSamplesTSS, std::vector<unsigned> const & indexValidSamplesBG)
{
    classificationExamples.reserve(2000);
    variablesExamples.resize(3);
    for (auto var : variablesExamples)
        var.reserve(2000);

    for (unsigned i = 0; i < 1000; i++)
    {
        classificationExamples.push_back(1);
        variablesExamples[0].push_back(dataset[indexValidSamplesTSS[i]].cpgContent);
        variablesExamples[1].push_back(dataset[indexValidSamplesTSS[i]].conservation);
        variablesExamples[2].push_back(dataset[indexValidSamplesTSS[i]].affinity);
    }

    for (unsigned i = 0; i < 1000; i++)
    {
        classificationExamples.push_back(0);
        variablesExamples[0].push_back(dataset[indexValidSamplesBG[i]].cpgContent);
        variablesExamples[1].push_back(dataset[indexValidSamplesBG[i]].conservation);
        variablesExamples[2].push_back(dataset[indexValidSamplesBG[i]].affinity);
    }
}
