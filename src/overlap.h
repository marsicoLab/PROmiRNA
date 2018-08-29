// ============================================================================
// Functions for overlapping miRNAs with CAGE data
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#pragma once

#include <algorithm>
#include <experimental/filesystem>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <seqan/bed_io.h>

#include "interval_tree.h"
#include "structs.h"

using namespace seqan;

// A map of tree structures that stores information about miRNAs
// Key is a pair that defines chromosome and strand because the interval trees can only process simple intervals [start, stop)
// Value is a tree of pairs of genomic locations and miRNA class/name/IDs
typedef std::map<std::pair<CharString, char>, IntTree<std::pair<GenomicLocation, MiRNA> > > miRNAIntTrees;

// A map of tree structures that stores information about background regions
// Key is a pair that defines chromosome and strand because the interval trees can only process simple intervals [start, stop)
// Value is a tree of genomic locations
typedef std::map<std::pair<CharString, char>, IntTree<GenomicLocation> > backgroundIntTrees;

/*!
 * @fn isOverlap
 * @brief Check whether to regions overlap (start and stop on the same chromosome)
 *
 * @signature bool isOverlap(start1, end1, start2, end2)
 *
 * @param[in]   start1      Start of first region
 * @param[in]   end1        End of first region
 * @param[in]   start2      Start of second region
 * @param[in]   end2        End of second region
 */
inline bool isOverlap(unsigned start1, unsigned end1, unsigned start2, unsigned end2)
{
       return std::max(start1, start2) < std::min(end1, end2);
}

/*!
 * @fn listFiles
 * @brief List BED files in a directory (not recurrent)
 *
 * @signature void listFiles(files, path)
 *
 * @param[in,out]   files       Vector of BED files in the directory
 * @param[in]       path        Path to check for files
 */
void listFiles(std::vector<std::string> & files, CharString & path)
{
    for (auto & p : std::experimental::filesystem::directory_iterator(toCString(path)))
    {
        if (std::experimental::filesystem::path(p).extension() == ".bed")
        {
            files.push_back(p.path().string());
        }
    }
}

/*!
 * @fn overlapRecord
 * @brief Overlap single BED record (CAGE tag) with all miRNA regions
 *
 * @signature void overlapRecord(overlapRegions, currentCage, currentScore, originalRegions)
 *
 * @param[in,out]   overlapRegions      Vector of background regions overlapping with CAGE tags (value: unsigned count,
 *                                      bool intergenic)
 * @param[in]       currentCage         Single CAGE tag that is checked for overlap
 * @param[in]       currentScore        Current score if count already present in the BED file
 * @param[in]       originalRegions     Vector of original background regions
 */
void overlapRecord(std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > & overlapRegions,
GenomicLocation const & currentCage, int const currentScore, miRNAIntTrees & originalRegions)
{
    MiRNA dummyMirna;
    std::pair<GenomicLocation, MiRNA> currentKey = std::make_pair(currentCage, dummyMirna);
    std::vector<std::pair<GenomicLocation, MiRNA>* > allOverlaps;
    originalRegions[std::make_pair(currentCage.chr, currentCage.strand)].overlapSearchAll(currentKey, allOverlaps);

    for (unsigned i = 0; i < allOverlaps.size(); i++)
    {
        std::pair<CharString, GenomicLocation> currentOverlapKey = std::make_pair((*allOverlaps[i]).second.name, currentCage);
        std::pair<unsigned, bool> currentOverlapValue;
        if ((*(*allOverlaps[i]).second.classification.begin()) == "intergenic")
        {
            currentOverlapValue.second = true;
        }
        else
        {
            currentOverlapValue.second = false;
        }
        currentOverlapValue.first = currentScore;

        auto ret = overlapRegions.insert(std::make_pair(currentOverlapKey, currentOverlapValue));
        if (!ret.second)
        {
            ret.first->second.first += currentScore;
        }
    }
}

/*!
 * @fn overlapRecord
 * @brief Overlap single BED record (CAGE tag) with all background regions
 *
 * @signature void overlapRecord(overlapRegions, currentCage, currentScore, originalRegions)
 *
 * @param[in,out]   overlapRegions      Vector of background regions overlapping with CAGE tags (value: unsigned count,
 *                                      bool intergenic)
 * @param[in]       currentCage         Single CAGE tag that is checked for overlap
 * @param[in]       currentScore        Current score if count already present in the BED file
 * @param[in]       originalRegions     Vector of original background regions
 */
void overlapRecord(std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > & overlapRegions, GenomicLocation & currentCage,
int const currentScore, backgroundIntTrees & originalRegions)
{
    std::vector<GenomicLocation* > allOverlaps;
    originalRegions[std::make_pair(currentCage.chr, currentCage.strand)].overlapSearchAll(currentCage, allOverlaps);

    for (unsigned i = 0; i < allOverlaps.size(); i++)
    {
        std::pair<CharString, GenomicLocation> currentOverlapKey = std::make_pair("", (*allOverlaps[i]));
        std::pair<unsigned, bool> currentOverlapValue;
        currentOverlapValue = std::make_pair(currentScore, true);

        auto ret = overlapRegions.insert(std::make_pair(currentOverlapKey, currentOverlapValue));
        if (!ret.second)
        {
            ret.first->second.first += currentScore;
        }
    }
}

/*!
 * @fn overlapSingleBedFile
 * @brief Overlap a single BED file of CAGE tags with miRNAs or background regions
 *
 * @signature void overlapSingleBedFile(overlapRegions, bedPath, originalRegions)
 *
 * @param[in,out]   overlapRegions      Vector of miRNAs overlapping with CAGE tags
 * @param[in]       bedPath             Path to BED file that should be overlapped with miRNAs or background regions
 * @param[in]       originalRegions     Vector of classified miRNAs or original background regions
 *
 * @throw Exception if input BED file cannot be opened
 */
template<typename TSeq>
void overlapSingleBedFile(std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > & overlapRegions, const char * bedPath,
std::map<std::pair<CharString, char>, IntTree<TSeq> > & originalRegions)
{
    BedFileIn bedIn;
    if (!open(bedIn, bedPath))
    {
        throw std::runtime_error("ERROR: Could not open BED file.");
        return;
    }

    BedRecord<Bed6> bedRecord;

    try
    {
        while (!atEnd(bedIn))
        {
            readRecord(bedRecord, bedIn);
            GenomicLocation currentCage;
            currentCage.chr = bedRecord.ref;
            currentCage.start = bedRecord.beginPos;
            currentCage.end = bedRecord.endPos;
            currentCage.strand = bedRecord.strand;

            // Distinguish ENCODE (score column does not refer to tag counts) from FANTOM data (score column does refer to tag counts) based on name column
            bool fantom = false;
            if (prefix(bedRecord.name, 3) == "chr")
            {
                fantom = true;
            }

            // The score needs to be stored in case the file format stores collapsed duplicate number in this column (FANTOM5)
            int currentScore;
            if (fantom && isdigit(bedRecord.score[0])) // TODO maybe change to check whole string
            {
                currentScore = std::stoi(toCString(bedRecord.score));
            }
            else
            {
                currentScore = 1;
            }

            overlapRecord(overlapRegions, currentCage, currentScore, originalRegions);
        }
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return;
    }
}

/*!
 * @fn overlapAllBedFiles
 * @brief Overlap all BED files in a directory with classified miRNAs or background regions
 *
 * @signature void overlapAllBedFiles(overlapRegions, bedFiles, originalRegions, threads)
 *
 * @param[in,out]   overlapRegions      Vector of miRNAs overlapping with CAGE tags
 * @param[in,out]   bedFiles            Vector of all BED files from a directory
 * @param[in]       originalRegions     Vector of classified miRNAs or background regions
 * @param[in]       threads             Number of threads
 *
 * @throw Exception if input files cannot be opened
 */
template<typename TSeq>
void overlapAllBedFiles(std::vector<std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > > & overlapRegions,
std::vector<std::string> & bedFiles, std::map<std::pair<CharString, char>, IntTree<TSeq> > & originalRegions, unsigned threads)
{
    overlapRegions.resize(bedFiles.size());

    // Read BED files
    omp_set_num_threads(threads);

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < bedFiles.size(); i++)
    {
        overlapSingleBedFile(overlapRegions[i], bedFiles[i].c_str(), originalRegions);
        std::cout << bedFiles[i] << "\t" << overlapRegions[i].size() << std::endl;
    }
}
