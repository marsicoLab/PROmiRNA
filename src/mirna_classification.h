// ============================================================================
// Functions for miRNA classification
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <seqan/gff_io.h>

#include "interval_tree.h"
#include "structs.h"

using namespace seqan;

// A map of tree structures that stores information about miRNAs
// Key is a pair that defines chromosome and strand because the interval trees can only process simple intervals [start, stop)
// Value is a tree of pairs of genomic locations and miRNA class/name/IDs
// GenomicLocation needs to be separate for comparison with CAGE tags later
typedef std::map<std::pair<CharString, char>, IntTree<std::pair<GenomicLocation, MiRNA> > > miRNAIntTrees;

/*!
 * @fn mirnaClassification
 * @brief Classify miRNAs
 *
 * @signature void mirnaClassification(classifiedMirnas, mirnaStarts, mirnaExonChecks, gffPath, mirnaPath, mirnaContextPath)
 *
 * @param[in,out]   classifiedMirnas    Classified miRNAs
 * @param[in,out]   mirnaStarts         Store start of real miRNAs for later calculation of distance between promoter and miRNA
 * @param[in,out]   mirnaExonChecks     Detect miRNAs that are also protein coding - needed for filter step later
 * @param[in]       gffPath             Path of input GFF file
 * @param[in]       mirnaPath           Path of miRNA txt file
 * @param[in]       mirnaContextPath    Path of miRNA context txt file
 *
 * @throw Exception if input files cannot be opened
 *
 * This function overlaps three files from miRBase in order to obtain miRNAs that have been annotated, where positional
 * information is present and where the class can be determined (exon, intron, intergenic, etc.). Those are used for
 * promoter prediction/all downstream analysis. Only primary transcripts are used. From miRNAs a region 50kb upstream
 * is extracted which is in the next step overlapped by CAGE tags.
 */
void mirnaClassification(miRNAIntTrees & classifiedMirnas, std::map<CharString, unsigned> & mirnaStarts,
std::map<CharString, CharString> & mirnaExonChecks, CharString const & gffPath, CharString const & mirnaPath,
CharString const & mirnaContextPath)
{
    std::vector<std::pair<GenomicLocation, MiRNA> > tempMirnas;

    // Read in mirBase files
	std::ifstream input;
    std::string line;

    std::map<CharString, CharString> mirnaMap;
    input.open(toCString(mirnaPath));
	if (input.is_open())
    {
        while(getline(input, line))
        {
            StringSet<CharString> currentRecord;
            strSplit(currentRecord, line, EqualsChar<'\t'>()); // No need to check for '#', no header for this file
            mirnaMap[currentRecord[1]] = currentRecord[0];
        }
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open mirBase file.");
        return;
    }

    // miRNA context
    std::map<CharString, std::set<CharString> > mirnaContextMap;

    std::ifstream input2;
    input2.open(toCString(mirnaContextPath));
    if (input2.is_open())
    {
        while(getline(input2, line))
        {
            StringSet<CharString> currentRecord;
            strSplit(currentRecord, line, EqualsChar<'\t'>());

            std::pair<std::map<CharString, std::set<CharString> >::iterator, bool> ret;
            std::set<CharString> currentClass;
            ret = mirnaContextMap.insert(std::make_pair(currentRecord[0], currentClass));
            ret.first->second.insert(currentRecord[3]);
			mirnaExonChecks[currentRecord[1]] = currentRecord[0];
        }
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open mirBase file.");
        return;
    }

    // Read in GFF file
    // Require GFF3 file extension as in current miRBase version 22
    // For debug purposes old version 18 (from old PROmiRNA) can be enabled when building with cmake

    std::ifstream input3;
    input3.open(toCString(gffPath));

    if (input3.is_open())
    {
        while(getline(input3, line))
        {
            if (line[0] != '#')
            {
                StringSet<CharString> currentRecord;
                strSplit(currentRecord, line, EqualsChar<'\t'>());

                // Iterate over attributes
                StringSet<CharString> attributes;
                strSplit(attributes, currentRecord[8], EqualsChar<';'>());

                MiRNA currentMirna;
                GenomicLocation currentMirnaLocation;

                #ifndef OLD_MIRBASE_VERSION
                    // New/release mode: gff3 from miRBase 22
                    for (unsigned i = 0; i < length(attributes); i++)
                    {
                        if (length(attributes[i]) >= 2 && prefix(attributes[i], 2) == "ID")
                        {
                            currentMirna.accession = suffix(attributes[i], 3);
                        }
                        else if (length(attributes[i]) >= 4 && prefix(attributes[i], 4) == "Name")
                        {
                            currentMirna.name = suffix(attributes[i], 5);
                        }
                    }
                #else
                    // Old/debug mode: gff2 from miRBase 18
                    for (unsigned i = 0; i < length(attributes); i++)
                    {
                        if (length(attributes[i]) >= 3 && prefix(attributes[i], 3) == "ACC")
                        {
                            currentMirna.accession = infix(attributes[i], 5, length(attributes[i]) - 1);
                        }
                        else if (length(attributes[i]) >= 3 && prefix(attributes[i], 3) == " ID")
                        {
                            currentMirna.name = infix(attributes[i], 5, length(attributes[i]) - 1);
                        }
                    }
                #endif

                unsigned beginPos = std::atoi(toCString(currentRecord[3])) - 1;
                unsigned endPos = std::atoi(toCString(currentRecord[4]));
                char strand = currentRecord[6][0];
                mirnaStarts[currentMirna.name] = beginPos;

                std::map<CharString, CharString>::iterator mit = mirnaMap.find(currentMirna.accession);
                if (mit != mirnaMap.end())
                {
                    currentMirnaLocation.chr = currentRecord[0];
                    // Extract 50 kb upstream
                    if (strand == '+')
                    {
                        if (beginPos < 50000)
                        {
                            currentMirnaLocation.start = 0;
                        }
                        else
                        {
                            currentMirnaLocation.start = beginPos - 50000; // 0-based, half open intervals
                        }
                        currentMirnaLocation.end = beginPos; // 0-based, half open intervals
                    }
                    else
                    {
                        currentMirnaLocation.start = endPos - 1; // 0-based, half open intervals
                        currentMirnaLocation.end = endPos + 50000 - 1; // 0-based, half open intervals
                    }
                    currentMirnaLocation.strand = strand;
                    currentMirna.id = mit->second;

                    std::map<CharString, std::set<CharString> >::iterator mcit = mirnaContextMap.find(mit->second);
                    if (mcit != mirnaContextMap.end())
                    {
                        currentMirna.classification = mcit->second;
                    }
                    else
                    {
                        std::set<CharString> currentClass;
                        currentClass.insert("intergenic");
                        currentMirna.classification = currentClass;
                    }
                }
                else
                {
                    continue;
                }

                // Store miRNAs in temporary vector
                // Tree is not balanced therefore need to random shuffle vector and insert into tree afterwards
                tempMirnas.push_back(std::make_pair(currentMirnaLocation, currentMirna));
            }
        }
    }
    else
    {
        throw std::runtime_error("ERROR: Could not open miRBase file.");
        return;
    }

    std::cout << "Number of miRNAs for analysis: " << tempMirnas.size() << std::endl;

    // Random shuffle elements in order to avoid having a sorted input for tree which will heavily increase runtime
    std::random_shuffle(tempMirnas.begin(), tempMirnas.end());

    // Insert elements into tree
    for (unsigned i = 0; i < tempMirnas.size(); i++)
    {
        std::pair<CharString, char> currentKey = std::make_pair(tempMirnas[i].first.chr, tempMirnas[i].first.strand);
        classifiedMirnas[currentKey].insert(tempMirnas[i]);
    }
}
