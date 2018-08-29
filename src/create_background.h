// ============================================================================
// Functions for creating background
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#pragma once

#include <iostream>
#include <vector>
#include <time.h>

#include <seqan/find.h>
#include <seqan/bed_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include "interval_tree.h"
#include "overlap.h"
#include "structs.h"

using namespace seqan;

// A map of tree structures that stores information about background regions
// Key is a pair that defines chromosome and strand because the interval trees can only process simple intervals [start, stop)
// Value is a tree of genomic locations
typedef std::map<std::pair<CharString, char>, IntTree<GenomicLocation> > backgroundIntTrees;

/*!
 * @fn extractSequenceFromIndex
 * @brief Extracts a sequence from an fai index based on chromosome, start and end.
 *
 * @signature void extractSequenceFromIndex(sequence, faiIndex, id, start, end)
 *
 * @param[in,out]   sequence        Extracted sequence
 * @param[in]       faiIndex        The fai index to be searched
 * @param[in]       id              The fasta ID were the sequence is located (chromosome in this case)
 * @param[in]       start           Start position of the sequence (0-based)
 * @param[in]       end             End position of the sequence (not included, half open intervals)
 *
 * @throw Exception if sequence cannot be read in
 */
void extractSequenceFromIndex(Dna5String & sequence, FaiIndex & faiIndex, CharString & id, unsigned start, unsigned end)
{
    // Translate sequence name to index.
    unsigned idx = 0;
    if (!getIdByName(idx, faiIndex, id))
    {
        throw std::out_of_range("ERROR: Index out of range.");
    }

    // Make sure start <= end <= sequenceLength
    if (start > sequenceLength(faiIndex, idx))
        start = sequenceLength(faiIndex, idx);
    if (end > sequenceLength(faiIndex, idx))
        end = sequenceLength(faiIndex, idx);
    if (start > end)
        end = start;

    try
    {
        readRegion(sequence, faiIndex, idx, start, end);
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return;
    }
}

/*!
 * @fn chromosomeList
 * @brief Returns a list of chromosomes present in the index excluding chrUn, chrM, random and hap
 *
 * @signature void chromosomeList(chromosomes, faiIndex)
 *
 * @param[in,out]   chromosomes     Set of chromosomes that can be used for background generation
 * @param[in]       faiIndex        The fai index to be searched
 *
 * The function lists all sequences/chromosomes present in the index and scans via pattern matching for names that should not be
 * used for background generation. All other sequence names are returned.
 *
 * NOTE Function not used in current version but kept for potential future usage.
 */
void chromosomeList(std::set<CharString> & chromosomes, FaiIndex & faiIndex)
{
    StringSet<CharString> keywords;
    appendValue(keywords, "random");
    appendValue(keywords, "Un");
    appendValue(keywords, "M");
    appendValue(keywords, "hap");

    for (unsigned i = 0; i < numSeqs(faiIndex); i++)
    {
        CharString currentName = sequenceName(faiIndex, i);
        Finder<CharString> finder(currentName);
        Pattern<StringSet<CharString>, WuManber> pattern(keywords);
        if (!find(finder, pattern))
        {
            chromosomes.insert(currentName);
        }
    }
}

/*!
 * @fn containsN
 * @brief Checks whether a DNA5 sequence contains one or more Ns.
 *
 * @signature bool containsN(seq)
 *
 * @param[in]       seq        The DNA5 sequence to be searched.
 *
 * Returns true if at least one N found.
 */
bool containsN(Dna5String & seq)
{
    for (unsigned i = 0; i < length(seq); i++)
    {
        if (seq[i] == 'N')
        {
            return true;
        }
    }
    return false;
}

/*!
 * @fn mergeBackgroundRegions
 * @brief Merge overlapping background regions.
 *
 * @signature void mergeBackgroundRegions(bgSequencesSet)
 *
 * @param[in,out]        bgSequencesSet         Set containing background sequences
 */
void mergeBackgroundRegions(std::set<GenomicLocation> & bgSequencesSet)
{
    for (auto it = bgSequencesSet.begin(); it != bgSequencesSet.end(); )
    {
        if (it != bgSequencesSet.begin())
        {
            if (((*std::prev(it)).chr == (*it).chr) &&
               ((*std::prev(it)).strand == (*it).strand) &&
               (isOverlap((*std::prev(it)).start, (*std::prev(it)).end, (*it).start, (*it).end)))
            {
                GenomicLocation mergedLocation;
                mergedLocation.chr = (*std::prev(it)).chr;
                mergedLocation.strand = (*std::prev(it)).strand;
                mergedLocation.start = (*std::prev(it)).start;

                if((*std::prev(it)).end < (*it).end)
                {
                    mergedLocation.end = (*it).end;
                }
                else
                {
                    mergedLocation.end = (*std::prev(it)).end;
                }
                bgSequencesSet.erase(std::prev(it));
                bgSequencesSet.erase(it);
                it = bgSequencesSet.insert(mergedLocation).first;
                it++;

            }
            else
            {
                ++it;
            }
        }
        else
        {
            ++it;
        }
    }
}

/*!
 * @fn createBackground
 * @brief Creates background regions.
 *
 * @signature void createBackground(bgSequences, chromosomes, totalNumber, minLength, maxLength, codingPath, repeatPath, faiIndex)
 *
 * @param[in,out]        bgSequences        The generated background sequences
 * @param[in]            chromosomes        Set of chromosomes that can be used for background generation
 * @param[in]            totalNumber        Total number of sequences to be generated
 * @param[in]            minLength          Minimum length of sequences
 * @param[in]            maxLength          Maximum length of sequences
 * @param[in]            repeatPath         Path to repeat regions
 * @param[in]            faiIndex           The fai index to be searched
 * @param[in]            regionsVector      Vector containing already gene regions for filtering of background sequences
 *
 * Background regions are generated randomly from available chromosomes, corresponding positions and lengths.
 * Sequences must not overlap with repeat or coding regions. If generated sequences overlap by chance, they are merged.
 *
 * @throw Exception if file(s) cannot be read in
 */
void createBackground(backgroundIntTrees & bgSequences, std::set<CharString> const & chromosomes,
unsigned const totalNumber, unsigned const minLength, unsigned const maxLength, CharString const & repeatPath,
FaiIndex & faiIndex, std::vector<GenomicLocation> & regionsVector)
{
    std::set<GenomicLocation> bgSequencesSet;

    // Read in repeat regions and store the positional informaton in an interval tree
    // Coding regions already present from reading in genes file in merge and filter step
    std::map<CharString, IntTree<GenomicLocation> > regions;

    // Repeat file
    BedFileIn bedIn;
    if (!open(bedIn, toCString(repeatPath)))
    {
        throw std::runtime_error("ERROR: Could not open BED file.");
        return;
    }

    BedRecord<Bed5> bedRecord;

    try
    {
        while (!atEnd(bedIn))
        {
            readRecord(bedRecord, bedIn);
            if (chromosomes.find(suffix(bedRecord.ref, 3)) == chromosomes.end())
                continue;
            GenomicLocation currentLocation;
            currentLocation.chr = bedRecord.ref;
            currentLocation.start = bedRecord.beginPos;
            currentLocation.end = bedRecord.endPos;

            regionsVector.push_back(currentLocation);
        }
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return;
    }

    // Random shuffle elements in order to avoid having a sorted input for tree which will heavily increase runtime
    std::random_shuffle(regionsVector.begin(), regionsVector.end());

    // Insert elements into tree
    for (unsigned i = 0; i < regionsVector.size(); i++)
    {
        regions[regionsVector[i].chr].insert(regionsVector[i]);
    }

    // Create background sequences
    while (bgSequencesSet.size() < totalNumber)
    {
        unsigned randomChrNum = rand() % chromosomes.size();                // Choose random chromosome
        std::set<CharString>::const_iterator mit = chromosomes.begin();
        std::advance(mit, randomChrNum);
        CharString randomChr = *mit;
        unsigned randomSize = minLength + rand() % (maxLength - minLength); // Choose random sequence length
        unsigned randomChrId;
        getIdByName(randomChrId, faiIndex, randomChr);                      // Cannot use randomChrNumber since this does not match FAI index that has additional chromosomes we do not use
        unsigned seqLength = sequenceLength(faiIndex, randomChrId);         // Size of chromosome with ID of chromosome
        unsigned randomPos = rand() % seqLength;                            // Choose random start position given the chromosome length

        if ((randomPos + randomSize - 1) < seqLength)
        {
            unsigned strand = rand() % 2; // Choose random strand

            GenomicLocation currentLocation;
            if (sequenceName(faiIndex, 0)[0] == 'c')                        // If reference genome has chromosome names starting with "chr"
                currentLocation.chr = randomChr;
            else
                currentLocation.chr = std::string("chr") + toCString(randomChr);

            currentLocation.start = randomPos;
            currentLocation.end = randomPos + randomSize;
            if (strand)
            {
                currentLocation.strand = '-';
            }
            else
            {
                currentLocation.strand = '+';
            }

            GenomicLocation *res = regions[currentLocation.chr].overlapSearchSingle(currentLocation);
            if (res == NULL)
            {
                Dna5String seq;
                extractSequenceFromIndex(seq, faiIndex, randomChr, randomPos, randomPos + randomSize);
                if (!containsN(seq))
                {
                    bgSequencesSet.insert(currentLocation);
                }
            }
        }
    }
    // Merge overlapping background regions
    // Cannot merge them earlier otherwise more background sequences than required would be generated
    mergeBackgroundRegions(bgSequencesSet);

    // Insert into tree
    for (auto bg : bgSequencesSet)
    {
        bgSequences[std::make_pair(bg.chr, bg.strand)].insert(bg);
    }
}
