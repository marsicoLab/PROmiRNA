// ============================================================================
// Functions for merging and filtering of miRNAs
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#pragma once

#include <algorithm>
#include <vector>

#include <seqan/find.h>
#include <seqan/gff_io.h>
#include <seqan/sequence.h>

#include "interval_tree.h"
#include "overlap.h"
#include "structs.h"

using namespace seqan;

// A map of tree structures that stores information about gene starts
// Key is a pair that defines chromosome and strand because the interval trees can only process simple intervals [start, stop)
// Value is a tree of genomic locations because no other information is needed for the filtering step from gene starts
typedef std::map<std::pair<CharString, char>, IntTree<GenomicLocation> > geneStartTrees;

// A map of tree structures that stores information about genes
// Key is a pair that defines chromosome and strand because the interval trees can only process simple intervals [start, stop)
// Value is a tree of pairs of additional information and genomic locations
// The additional information consists of a pair that stores first the type (exon. CDS) and second the gene name in order to filter for miRNA names
typedef std::map<std::pair<CharString, char>, IntTree<std::pair<std::pair<CharString, CharString>, GenomicLocation> > > geneTrees;

/*!
 * @fn isMirna
 * @brief Tests wether a gene name is a miRNA
 *
 * @signature bool isMirna(geneName)
 *
 * @param[in]   geneName             The gene name to check for via online pattern matching
 */
bool isMirna(CharString & geneName)
{
	// Patterns
	String<CharString> mirnaIdentifier;
	appendValue(mirnaIdentifier, "-mir-");
	appendValue(mirnaIdentifier, "-miR-");
	appendValue(mirnaIdentifier, "MIR");

	Finder<CharString> finder(geneName);

	Pattern<String<CharString>, WuManber> pattern(mirnaIdentifier);
	if (find(finder, pattern))
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*!
 * @fn isValidTag
 * @brief Checks whether a merged regions passes certain filter criteria
 *
 * @signature bool isValidTag(tag, intergenic, geneStarts, genes, mirnaExonChecks, geneToTrans)
 *
 * @param[in]       tag                The current merged region to check
 * @param[in]       intergenic         Bool defining whether tag is an intergenic region or not
 * @param[in]       geneStarts         Trees containing information about gene starts
 * @param[in]       genes              Trees containing information about genes
 * @param[in]		mirnaExonChecks	   Map of transcript IDs to miR IDs (mirBase) from mirna_context
 * @param[in]		genetoTrans		   Map of gene name -> vector of transcript IDs (for filtering)
 *
 * If the tag is an intergenic region, the overlap with genes is checked. If the merged region overlaps one or more
 * coding region and the type is equal to CDS the region is not valid, except the coding region corresponds to a miRNA.
 * If the tag is an intragenic region, first the overlap with gene starts is checked. Any overlap with a gene start
 * makes the region valid and the function returns true. If not, overlaps with genes are tested and the tag is declared
 * not valid if an overlap is found where the type is an exon and the coding region corresponds not to a miRNA.
 */
bool isValidTag(DatasetRecord & tag, bool intergenic, geneStartTrees & geneStarts, geneTrees & genes,
std::map<CharString, CharString> const & mirnaExonChecks, std::map<CharString, std::vector<CharString> > const & geneToTrans)
{
	if (intergenic)
	{
		std::pair<std::pair<CharString, CharString>, GenomicLocation> currentKey = std::make_pair(std::make_pair("", ""), tag.location);
		std::vector<std::pair<std::pair<CharString, CharString>, GenomicLocation>* > allOverlaps;
		genes[std::make_pair(tag.location.chr, tag.location.strand)].overlapSearchAll(currentKey, allOverlaps);

		// Filter out intergenic regions with CDS if they don't have -mir-, -miR- or MIR (filtered out before)
		for (unsigned i = 0; i < allOverlaps.size(); i++)
		{
			if ((*allOverlaps[i]).first.first == "CDS" && !isMirna((*allOverlaps[i]).first.second))
			{
				return false;
			}
		}
	}
	else
	{
		GenomicLocation adaptedTag = tag.location;
		if (adaptedTag.strand == '+')
		{
			adaptedTag.start -= 100;
			adaptedTag.end++;
		}
		else
		{
			adaptedTag.start++;
			adaptedTag.end += 100;
		}

		GenomicLocation * overlap = geneStarts[std::make_pair(adaptedTag.chr, adaptedTag.strand)].overlapSearchSingle(adaptedTag);
		if (overlap != NULL)
        {
			return true;
		}

		std::pair<std::pair<CharString, CharString>, GenomicLocation> currentKey = std::make_pair(std::make_pair("", ""), tag.location);
		std::vector<std::pair<std::pair<CharString, CharString>, GenomicLocation>* > allOverlaps;
		genes[std::make_pair(tag.location.chr, tag.location.strand)].overlapSearchAll(currentKey, allOverlaps);


		// Filter out intergenic regions with exon if they don't have -mir-, -miR- or MIR (filtered out before)
		for (unsigned i = 0; i < allOverlaps.size(); i++)
		{
			if ((*allOverlaps[i]).first.first == "exon" && !isMirna((*allOverlaps[i]).first.second))
			{
				auto mit = geneToTrans.find((*allOverlaps[i]).first.second);
				if (mit != geneToTrans.end())
                {
					bool truth = false;
					for (int k = 0; k < mit->second.size(); ++k)
                    {
						auto pit = mirnaExonChecks.find(mit->second[k]);
						if (pit != mirnaExonChecks.end())
                        {
							truth = true;
						}
					}
					if (!truth)
                    {
						return false;
					}
				}
			}
		}
	}

	return true;
}

void nonZeroElementComparison(std::vector<unsigned> & nonZeroElements, std::vector<unsigned> const & v1, std::vector<unsigned> const & v2)
{
	if (v1.size() == v2.size() && v1.size() == nonZeroElements.size())
	{
		for (unsigned i = 0; i < v1.size(); i++)
		{
			if (v1[i] > 0 && v2[i] > 0)
			{
				nonZeroElements[i] = 1;
			}
			else
			{
				nonZeroElements[i] = 0;
			}
		}
	}
	else
	{
		throw std::runtime_error("ERROR: Vectors have different size, cannot compare element-wise.");
	}
}

/*!
 * @fn mergeTags
 * @brief Merges tags that overlap with miRNAs found in a single file if overlapping each other
 *
 * @signature void mergeTags(dataset, countMatrix, normMatrix, overlapRegions, uniqueRegions, mirnaStarts, geneStarts, genes, index, mirnaExonChecks, geneToTrans)
 *
 * @param[in, out]  dataset            The resulting dataset for the EM algorithm
 * @param[in]       countMatrix        The count matrix for all overlap regions and BED files
 * @param[in]       normMatrix         The normalized count matrix
 * @param[in]       overlapRegions     All tags found to overlap miRNAs for all BED files
 * @param[in]       uniqueRegions      Unique tags found to overlap miRNAs in all BED files
 * @param[in]       mirnaStarts        A map storing the name (as key) and start of miRNAs (as value) needed for
 *                                     distance score calculation
 * @param[in]       geneStarts         Trees containing information about gene starts
 * @param[in]       genes              Trees containing information about genes
 * @param[in]       index              Index in countMatrix and normMatrix depending on the column/BED file that is
 *                                     currently processed
 * @param[in]		mirnaExonChecks	   Map of transcript IDs to miR IDs (mirBase) from mirna_context
 * @param[in]		geneToTrans		   Map of gene name -> vector of transcript IDs (for filtering)
 *
 * The function merges overlapping regions of each BED file that were already shown to overlap a miRNA. If the merged
 * regions pass a filter step, a dataset record is created with information about tag count (added up from normalized
 * matrix), lib number (max from countMatrix), distance score and positional information.
 */
void mergeTags(std::vector<DatasetRecord> & dataset, MatrixPair const & countMatrix, MatrixSingle const & normMatrix,
std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > const & overlapRegions,
std::map<std::pair<CharString, GenomicLocation>, unsigned> const & uniqueRegions,
std::map<CharString, unsigned> & mirnaStarts, geneStartTrees & geneStarts, geneTrees & genes,
unsigned const index, std::map<CharString, CharString> const & mirnaExonChecks,
std::map<CharString, std::vector<CharString> > const & geneToTrans)
{
	DatasetRecord currentRecord;
	std::vector<unsigned> countVectorMerged(countMatrix.colnames.size(), 1); // Used to calculate merged tag counts
	for (auto it = overlapRegions.begin(); it != overlapRegions.end(); it++)
	{
		StringSet<std::string> identifier;
		strSplit(identifier, normMatrix.rownames[uniqueRegions.find(it->first)->second], EqualsChar<'-'>());

		std::string currentMirna;
		for (unsigned i = 3; i < (length(identifier) - 1); i++)
			currentMirna = currentMirna + identifier[i] + std::string("-");
		currentMirna = currentMirna + identifier[length(identifier) - 1];

		if (it != overlapRegions.begin())
		{
			if ((std::prev(it)->first.first == it->first.first) &&
				(currentRecord.location.chr == it->first.second.chr) &&
				(currentRecord.location.strand == it->first.second.strand) &&
				(isOverlap(currentRecord.location.start - 20, currentRecord.location.end + 20+1, it->first.second.start, it->first.second.end)))
			{
				if(currentRecord.location.end < it->first.second.end)
				{
					currentRecord.location.end = it->first.second.end;
				}
				currentRecord.tagCount += normMatrix.matrix[uniqueRegions.find(it->first)->second][index];

				auto & currentLibVector = countMatrix.matrix[uniqueRegions.find(it->first)->second];
				unsigned c = 0;
				for (auto i : currentLibVector)
				{
					if (i.first > 0)
						c++;
				}
				currentRecord.libs = std::min(currentRecord.libs, c);
			}
			else
			{
				if (isValidTag(currentRecord, std::prev(it)->second.second, geneStarts, genes, mirnaExonChecks, geneToTrans))
				{
					dataset.push_back(currentRecord);
				}
				currentRecord.location = it->first.second;
				currentRecord.tagCount = normMatrix.matrix[uniqueRegions.find(it->first)->second][index];

				auto & currentLibVector = countMatrix.matrix[uniqueRegions.find(it->first)->second];
				unsigned c = 0;
				for (auto i : currentLibVector)
				{
					if (i.first > 0)
						c++;
				}
				currentRecord.libs = c;

				if (currentRecord.location.strand == '+')
				{
					currentRecord.distance = mirnaStarts[currentMirna] - currentRecord.location.start;
				}
				else
				{
					currentRecord.distance = currentRecord.location.start - mirnaStarts[currentMirna];
				}
				currentRecord.distScore = pow((1 + currentRecord.distance / 1000), -1);

				currentRecord.bg = false;
				currentRecord.miRNA = currentMirna;
			}
		}
		else
		{
			currentRecord.location = it->first.second;
			currentRecord.tagCount = normMatrix.matrix[uniqueRegions.find(it->first)->second][index];

			auto & currentLibVector = countMatrix.matrix[uniqueRegions.find(it->first)->second];
			unsigned c = 0;
			for (auto i : currentLibVector)
			{
				if (i.first > 0)
					c++;
			}
			currentRecord.libs = c;

			if (currentRecord.location.strand == '+')
			{
				currentRecord.distance = mirnaStarts[currentMirna] - currentRecord.location.start;
			}
			else
			{
				currentRecord.distance = currentRecord.location.start - mirnaStarts[currentMirna];
			}
			currentRecord.distScore = pow((1 + currentRecord.distance / 1000), -1);

			currentRecord.bg = false;
			currentRecord.miRNA = currentMirna;
		}
	}

	if (isValidTag(currentRecord, std::prev(overlapRegions.end())->second.second, geneStarts, genes, mirnaExonChecks, geneToTrans))
	{
		dataset.push_back(currentRecord);
	}
}

/*!
 * @fn extractCorePromoter
 * @brief Replaces original positions of regions with 1000 bp core promoter region
 *
 * @signature void extractCorePromoter(dataset)
 *
 * @param[in, out]  oldPositionsTSSDataset      Old positions of the regions
 * @param[in, out]  dataset                     The dataset for the EM algorithm
 */
void extractCorePromoter(std::vector<std::pair<unsigned, unsigned> > & oldPositionsTSSDataset,
std::vector<DatasetRecord> & dataset)
{
	// Old positions needed in the end for the output file
	oldPositionsTSSDataset.resize(dataset.size());

	for (unsigned i = 0; i < dataset.size(); i++)
	{
		oldPositionsTSSDataset[i] = std::make_pair(dataset[i].location.start, dataset[i].location.end);

		if ((dataset[i].location.end - dataset[i].location.start) < 1000)
		{
			unsigned dist = dataset[i].location.end - dataset[i].location.start;
			dataset[i].location.start -= ((double) (1000 - dist) / (double) 2);
			dataset[i].location.end += ((double) (1000 - dist) / (double) 2);
		}
	}
}

/*!
 * @fn mergeAndFilter
 * @brief Merges tags that overlap with miRNAs found in for all BED files if overlapping each other
 *
 * @signature void mergeAndFilter(dataset, countMatrix, normMatrix, overlapRegions, uniqueRegions, mirnaStarts, geneStartsPath, genesPath)
 *
 * @param[in, out]  dataset            The resulting dataset for the EM algorithm
 * @param[in]       countMatrix        The count matrix for all overlap regions and BED files
 * @param[in]       normMatrix         The normalized count matrix
 * @param[in]       overlapRegions     All tags found to overlap miRNAs for all BED files
 * @param[in]       uniqueRegions      Unique tags found to overlap miRNAs in all BED files
 * @param[in]       mirnaStarts        A map storing the name (as key) and start of miRNAs (as value) needed for distance score calculation
 * @param[in]       geneStarts         Path to gene starts
 * @param[in]       genes              Path to genes
 * @param[in]		mirnaExonChecks	   Map of transcript IDs to miR IDs (mirBase) from mirna_context
 * @param[in]		regionsVectorForBG Store gene regions already for background
 * @param[in]		threads            Number of threads to use
 *
 * @throw Exception if input files cannot be opened
 *
 * Creates trees with gene starts and genes information and calls merge function for each BED file.
 */
void mergeAndFilter(std::vector<DatasetRecord> & dataset, MatrixPair const & countMatrix, MatrixSingle const & normMatrix,
std::vector<std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > > const & overlapRegions,
std::map<std::pair<CharString, GenomicLocation>, unsigned> const & uniqueRegions,
std::map<CharString, unsigned> & mirnaStarts, CharString const & geneStartsPath, CharString const & genesPath,
std::map<CharString, CharString> const & mirnaExonChecks, std::vector<GenomicLocation> & regionsVectorForBG, unsigned const threads)
{
	// Read in GFF file with gene starts
	GffFileIn gffIn;
	if (!open(gffIn, toCString(geneStartsPath)))
	{
		throw std::runtime_error("ERROR: Could not open GFF file.");
		return;
	}

	GffRecord gffRecord;
	std::vector<GenomicLocation> geneStartsVector;
	geneStartTrees geneStarts;

	try
	{
		while (!atEnd(gffIn))
		{
			readRecord(gffRecord, gffIn);
			GenomicLocation currentLocation;
			currentLocation.chr = std::string("chr") + toCString(gffRecord.ref);
			currentLocation.start = gffRecord.beginPos;
			currentLocation.end = gffRecord.endPos;
			currentLocation.strand = gffRecord.strand;

			geneStartsVector.push_back(currentLocation);
		}
	}
	catch (Exception const & e)
	{
		std::cout << e.what() << std::endl;
		return;
	}

	GffFileIn gffIn2;
	if (!open(gffIn2, toCString(genesPath)))
	{
		throw std::runtime_error("ERROR: Could not open GTF file.");
		return;
	}

	std::random_shuffle(geneStartsVector.begin(), geneStartsVector.end());

	// Insert elements into tree
	for (unsigned i = 0; i < geneStartsVector.size(); i++)
	{
		geneStarts[std::make_pair(geneStartsVector[i].chr, geneStartsVector[i].strand)].insert(geneStartsVector[i]);
	}

	// Read in GFF file with genes
	std::vector<std::pair<std::pair<CharString, CharString>, GenomicLocation> > genesVector;
	geneTrees genes;

	unsigned geneNameIndex;
	unsigned transNameIndex;
	unsigned count = 0;
	std::map<CharString, std::vector<CharString> > geneToTrans;

	try
	{
		while (!atEnd(gffIn2))
		{
			readRecord(gffRecord, gffIn2);
			if(count == 0)
			{
				for (unsigned i = 0; i < length(gffRecord.tagNames); i++)
				{
					if (gffRecord.tagNames[i] == "gene_name")
					{
						geneNameIndex = i;
					}
                    else if (gffRecord.tagNames[i] == "transcript_id")
                    {
						transNameIndex = i;
					}
				}
				geneToTrans[gffRecord.tagValues[geneNameIndex]].push_back(gffRecord.tagValues[transNameIndex]);
				count++;
			}

			GenomicLocation currentLocation;
			currentLocation.chr = std::string("chr") + toCString(gffRecord.ref);
			currentLocation.start = gffRecord.beginPos;
			currentLocation.end = gffRecord.endPos;
			currentLocation.strand = gffRecord.strand;

			genesVector.push_back(std::make_pair(std::make_pair(gffRecord.type, gffRecord.tagValues[transNameIndex]), currentLocation));
            regionsVectorForBG.push_back(currentLocation);
		}
	}
	catch (Exception const & e)
	{
		std::cout << e.what() << std::endl;
		return;
	}
	std::random_shuffle(genesVector.begin(), genesVector.end());

	// Insert elements into tree
	for (unsigned i = 0; i < genesVector.size(); i++)
	{
		genes[std::make_pair(genesVector[i].second.chr, genesVector[i].second.strand)].insert(genesVector[i]);
	}

    // Parallelize over BED files
    // Each thread gets own dataset
    // Lock for concatenating datasets
    omp_set_num_threads(threads);

    omp_lock_t writelock;
    omp_init_lock(&writelock);

    #pragma omp parallel for schedule(dynamic)
	for (unsigned i = 0; i < normMatrix.colnames.size(); i++)
	{
        std::vector<DatasetRecord> subset;
		mergeTags(subset, countMatrix, normMatrix, overlapRegions[i], uniqueRegions, mirnaStarts, geneStarts, genes, i, mirnaExonChecks, geneToTrans);

        omp_set_lock(&writelock);
        dataset.insert(dataset.end(), subset.begin(), subset.end());
        omp_unset_lock(&writelock);
	}
    omp_destroy_lock(&writelock);
}

/*!
 * @fn prepareBackground
 * @brief Create dataset records for background with count and normalized matrix
 *
 * @signature void buildMatrix(dataset, countMatrix, normMatrix)
 *
 * @param[in,out]   dataset             The resulting dataset with background sequence information
 * @param[in]       countMatrix         Count matrix of background sequence overlap
 * @param[in]       uniqueRegions       Normalized count matrix
 */
void prepareBackground(std::vector<DatasetRecord> & dataset, MatrixPair const & countMatrix, MatrixSingle const & normMatrix,
unsigned const totalNumber, unsigned const maxTagValue)
{
	// Make all pairs of column and row indices that are valid
	std::vector<std::pair<unsigned, unsigned> > allValidIndices;
	for (unsigned i = 0; i < normMatrix.colnames.size(); i++)
	{
		for (unsigned j = 0; j < normMatrix.rownames.size(); j++)
		{
			if (countMatrix.matrix[j][i].first != 0 && countMatrix.matrix[j][i].first <= maxTagValue)
			{
				allValidIndices.push_back(std::make_pair(j, i));
			}
		}
	}

	std::random_shuffle(allValidIndices.begin(), allValidIndices.end());
	unsigned upperLimit = (allValidIndices.size() <= totalNumber) ? (countMatrix.rownames.size()) : (totalNumber);

	for (unsigned i = 0; i < upperLimit; i++)
	{
		StringSet<std::string> currentRegion;
		strSplit(currentRegion, normMatrix.rownames[allValidIndices[i].first], EqualsChar<'-'>());

		DatasetRecord currentRecord;
		currentRecord.location.chr = currentRegion[0];
		currentRecord.location.start = stoi(currentRegion[1]);
		currentRecord.location.end = stoi(currentRegion[2]);
		currentRecord.location.strand = currentRegion[3][0];
		currentRecord.tagCount = normMatrix.matrix[allValidIndices[i].first][allValidIndices[i].second];
		currentRecord.distance = 0; // No distance to miRNA possible, set to zero
		currentRecord.distScore = 0.01;
		currentRecord.libs = 0; // No requirement for examples, value does not matter
		currentRecord.bg = true;

		dataset.push_back(currentRecord);
	}
}
