// ============================================================================
// This is the main function of PROmiRNA
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#include <iostream>
#include <vector>

#include <RcppArmadillo.h>      	// For Armadillo as well as Rcpp
#include <RInside.h>                // For the embedded R via RInside

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include "build_matrix.h"
#include "create_background.h"
#include "create_dataset.h"
#include "em_algorithm.h"
#include "filter_mirna.h"
#include "mirna_classification.h"
#include "overlap.h"
#include "post_processing.h"
#include "structs.h"

using namespace seqan;

int main(int argc, char *argv[])
{
    // Argument parser
    ArgumentParser parser("PROmiRNA");
    addDescription(parser, "Prediction of miRNA promoters from CAGE data.");

    addOption(parser, ArgParseOption("g", "genome", "Path to genome FASTA file", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "genome", "fa fasta");
    setRequired(parser, "genome");

    addOption(parser, ArgParseOption("c", "coding", "Path to coding/gene regions GTF file", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "coding", "gtf");
    setRequired(parser, "coding");

    addOption(parser, ArgParseOption("s", "starts", "Path to gene start regions GFF file", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "starts", "gff");
    setRequired(parser, "starts");

    addOption(parser, ArgParseOption("r", "repeat", "Path to repeat regions BED file", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "repeat", "bed");
    setRequired(parser, "repeat");

    addOption(parser, ArgParseOption("m", "mirna", "Path to miRNA txt file", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "mirna", "txt");
    setRequired(parser, "mirna");

    addOption(parser, ArgParseOption("n", "context", "Path to miRNA context txt file", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "context", "txt");
    setRequired(parser, "context");

    addOption(parser, ArgParseOption("a", "annotation", "Path to miRNA annotation GFF file", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "annotation", "gff gff2 gff3");
    setRequired(parser, "annotation");

    addOption(parser, ArgParseOption("p", "psem", "Path to PSEM matrix file", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "psem", "psem");
    setRequired(parser, "psem");

    addOption(parser, ArgParseOption("w", "wig", "Path to wig file (corresponding wib file with same name needs to present in the same folder)", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "wig", "wig txt");
    setRequired(parser, "wig");

    addOption(parser, ArgParseOption("i", "cage", "Path to folder with CAGE tag BED files", ArgParseArgument::INPUT_DIRECTORY, "IN"));
    setRequired(parser, "cage");

    addOption(parser, ArgParseOption("t", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", "1");

    addOption(parser, ArgParseOption("l", "minlib", "Minimum number of libraries required", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "minlib", "0");

    addOption(parser, ArgParseOption("o", "output", "Path to output txt file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setValidValues(parser, "output", "txt");
    setDefaultValue(parser, "output", "output_promoter_mirnas.txt");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    unsigned threads, minLib;
    CharString genomePath, genesPath, geneStartsPath, repeatPath, gffPath, mirnaPath, mirnaContextPath, psemPath, wigPath, bedDirPath, outputPath;

    getOptionValue(genomePath, parser, "genome");
    getOptionValue(genesPath, parser, "coding");
    getOptionValue(geneStartsPath, parser, "starts");
    getOptionValue(repeatPath, parser, "repeat");
    getOptionValue(gffPath, parser, "annotation");
    getOptionValue(mirnaPath, parser, "mirna");
    getOptionValue(mirnaContextPath, parser, "context");
    getOptionValue(psemPath, parser, "psem");
    getOptionValue(wigPath, parser, "wig");
    getOptionValue(bedDirPath, parser, "cage");
    getOptionValue(threads, parser, "threads");
    getOptionValue(minLib, parser, "minlib");
    getOptionValue(outputPath, parser, "output");

    std::cout << "Starting miRNA promoter prediction" << std::endl;

    try
    {
        // Set random seed for all following randomized functions
        srand(time(NULL));

        // miRNA classification
        miRNAIntTrees classifiedMirnas;
        std::map<CharString, unsigned> mirnaStarts;
        std::map<CharString, CharString> mirnaExonChecks;
        mirnaClassification(classifiedMirnas, mirnaStarts, mirnaExonChecks, gffPath, mirnaPath, mirnaContextPath);


        // Overlap miRNAs with CAGE tags
        std::vector<std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > > overlapMirnas;
        std::vector<std::string> bedFiles;
        listFiles(bedFiles, bedDirPath);
        std::cout << "Number of overlaps between tags and miRNAs:" << std::endl;
        overlapAllBedFiles(overlapMirnas, bedFiles, classifiedMirnas, threads);

        // Build count and quantile normalized matrix for miRNAs
        MatrixSingle normMatrix;
        MatrixPair countMatrix;
        std::map<std::pair<CharString, GenomicLocation>, unsigned> uniqueRegions;
        buildMatrix(normMatrix, countMatrix, uniqueRegions, overlapMirnas, bedFiles);
        std::cout << "Unique regions TSS: " << uniqueRegions.size() << std::endl;

        // Merge and filter miRNAs, extract 1000 bp
        std::vector<DatasetRecord> dataset;
        dataset.reserve(normMatrix.rownames.size() * normMatrix.colnames.size() + (normMatrix.rownames.size() * 0.6) * normMatrix.colnames.size());
        std::vector<GenomicLocation> regionsVectorForBG;
        mergeAndFilter(dataset, countMatrix, normMatrix, overlapMirnas, uniqueRegions, mirnaStarts, geneStartsPath, genesPath, mirnaExonChecks, regionsVectorForBG, threads);

        std::vector<std::pair<unsigned, unsigned> > oldPositionsTSSDataset;
        extractCorePromoter(oldPositionsTSSDataset, dataset);

        unsigned backgroundStart = dataset.size();
        std::cout << "Number of TSS in dataset for EM algorithm: " << dataset.size() << std::endl;

        // Background creation
        backgroundIntTrees bgSequences;
        unsigned minLength = 1000;
        unsigned maxLength = 5000;
        unsigned totalNumber = dataset.size() * 0.6; // We need 20% in the end but since some can get sorted out due to no overlaps this is a compromise; still not ideal

        FaiIndex faiIndex;
        if (!open(faiIndex, toCString(genomePath)))
        {
            if (!build(faiIndex, toCString(genomePath)))
            {
                throw std::runtime_error("ERROR: Index could not be loaded or built.");
            }
            if (!save(faiIndex))
            {
                throw std::runtime_error("ERROR: Index could not be written do disk.");
            }
        }

        std::set<CharString> chromosomes;
        if (sequenceName(faiIndex, 0)[0] == 'c')
        {
            for (auto i : uniqueRegions)
            {
                chromosomes.insert(i.first.second.chr);
            }
        }
        else
        {
            for (auto i : uniqueRegions)
            {
                chromosomes.insert(suffix(i.first.second.chr, 3));
            }
        }

        createBackground(bgSequences, chromosomes, totalNumber, minLength, maxLength, repeatPath, faiIndex, regionsVectorForBG);

        // Overlap background regions with CAGE tags
        std::vector<std::map<std::pair<CharString, GenomicLocation>, std::pair<unsigned, bool> > > overlapBGRegions;
        std::cout << "Number of overlaps between tags and background regions:" << std::endl;
        overlapAllBedFiles(overlapBGRegions, bedFiles, bgSequences, threads);

        // Build count and quantile normalized matrix for background regions
        MatrixSingle normMatrixBG;
        MatrixPair countMatrixBG;
        std::map<std::pair<CharString, GenomicLocation>, unsigned> uniqueRegionsBG;
        buildMatrix(normMatrixBG, countMatrixBG, uniqueRegionsBG, overlapBGRegions, bedFiles);
        std::cout << "Unique regions background: " << uniqueRegionsBG.size() << std::endl;

        // Extract 1000 bp
        unsigned const maxReadCount = 50; // Measured via exploring distribution of tag counts for bg and miRNA
        prepareBackground(dataset, countMatrixBG, normMatrixBG, dataset.size() * 0.2, maxReadCount);
        std::cout << "Number of BG in dataset for EM algorithm: " << dataset.size() - backgroundStart << std::endl;
        std::cout << "Total number of records in dataset for EM algorithm: " << dataset.size() << std::endl;

        // TODO Verify exact threshold
        // Need sufficient number of samples for beta parameters and initial parameters estimation
        if ((dataset.size() - backgroundStart) < 6000)
        {
            throw std::runtime_error("ERROR: The data provided is too sparse to run a sensible prediction. Please provide more CAGE data as input.");
        }

        if (minLib < 2)
        {
            std::cerr << "WARNING: Minimum number of libraries either not user-defined or too low. The default settings will be used." << std::endl;
            minLib = max(2, (1/3) * bedFiles.size());
        }

        // Create dataset
        createDataset(dataset, psemPath, wigPath, faiIndex);

        // Get all indices of samples that can be used for parameter estimation
        // Random shuffle so that they can be used for parameters without sampling
        std::vector<unsigned> indexValidSamplesTSS, indexValidSamplesBG;
        getValidSamplesForParameterEstimation(indexValidSamplesTSS, indexValidSamplesBG, dataset, backgroundStart, minLib);

        // Get examples for beta parameter estimation - takes the first 1000 random indices
        std::vector<unsigned> classificationExamples;
        std::vector<std::vector<double> > variablesExamples;
        getExamplesForBetaParams(classificationExamples, variablesExamples, dataset, indexValidSamplesTSS, indexValidSamplesBG);


        // Estimate beta parameters
        // Features: Tag count, CpG content and conservation from extracted samples
        // Labels: 1 for promoter, 0 for background
        std::vector<double> betaParams;
        RInside R(argc, argv);
        estimateBetaParameters(betaParams, R, classificationExamples, variablesExamples);
        betaParams.push_back(std::min<double>(std::min<double>(betaParams[1], betaParams[2]), betaParams[3]));

        std::cout << "Beta params: ";
        for (auto i : betaParams)
            std::cout << i << "\t";
        std::cout << std::endl;

        std::cout << "Starting EM algorithm" << std::endl;
        // Estimate initial parameters
        double mu1, mu2, lambda1, lambda2;
        estimateInitialParameters(mu1, mu2, lambda1, lambda2, dataset, indexValidSamplesTSS, indexValidSamplesBG);

        // EM algorithm
        std::vector<double> outputParams;
        std::vector<std::vector<double> > promoterProbs;
        emAlgorithm(outputParams, promoterProbs, dataset, mu1, mu2, lambda1, lambda2, betaParams);

        // Post-processing
        postProcessing(outputPath, promoterProbs, dataset, oldPositionsTSSDataset, faiIndex, threads);
    }
    catch (Exception const & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    std::cout << "Finished promoter prediction" << std::endl;
	return 0;
}
