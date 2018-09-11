// ============================================================================
// Functions for the EM algorithm
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#pragma once

#include <iostream>
#include <vector>

// The tests do not need to embed R
#ifndef TESTCASE
    #include <RcppArmadillo.h>      	// For Armadillo as well as Rcpp
    #include <RInside.h>                // For the embedded R via RInside
#endif

#include <cmath>

#include "structs.h"

using namespace seqan;

/*!
 * @fn estimateBetaParameters
 * @brief Estimate the beta parameters for the EM algorithm
 *
 * @signature void estimateBetaParameters(betaParams, R, classification, variables)
 *
 * @param[in,out]   betaParams       Beta parameters calculated by the generalized linear model
 * @param[in]       R                Embedded R instance (need to be initiated in the main function with argc and argv)
 * @param[in]       classification   Classification of samples (1 for promoter, 0 for background)
 * @param[in]       variables        X (variables) matrix for regression model (tag count, CpG content and conservation)
 *
 * Calls the glm function with binomial family in R with samples from EM algorithm dataset.
 *
 * NOTE: This function is not built for test cases.
 */
#ifndef TESTCASE
void estimateBetaParameters(std::vector<double> & betaParams, RInside & R, std::vector<unsigned> const & classification,
std::vector<std::vector<double> > const & variables)
{
    std::string txt = "suppressMessages(require(stats))";
    R.parseEvalQ(txt);

    R["y"] = classification;
    R["x1"] = variables[0];
    R["x2"] = variables[1];
    R["x3"] = variables[2];

    std::string evalstr = "glm(y ~ x1 + x2 + x3, family = 'binomial')$coefficients";
    betaParams = Rcpp::as<std::vector<double> >(R.parseEval(evalstr));
}
#endif

/*!
 * @fn average
 * @brief Computes the mean of a numerical vector.
 *
 * @signature double average(inputVector)
 *
 * @param[in]       inputVector      Vector with numerical values to be averaged
 *
 * @throws Exception if vector is empty.
 */
template<typename TNum>
double average(std::vector<TNum> const & inputVector)
{
    if (inputVector.size() > 0)
    {
        double sumVector = 0;
        for (auto & n : inputVector)
            sumVector += n;
        return sumVector / (double) inputVector.size();
    }
    else
    {
        throw std::runtime_error("ERROR: Cannot compute the mean of an empty vector.");
    }
}

/*!
 * @fn probDensityFunction
 * @brief Probability density function.
 *
 * @signature double probDensityFunction(count, mu, lambda)
 *
 * @param[in]       count       The logarithm of tag count minus two
 * @param[in]       mu          Mean estimator
 * @param[in]       lambda      Standard deviation estimator
 */
double probDensityFunction(double const count, double const mu, double const lambda)
{
    return pow((lambda / (2 * M_PI * pow(count, 3))), 0.5) * exp(-lambda * pow((count - mu), 2) / (2 * pow(mu, 2) * count));
}

/*!
 * @fn estimateInitialParameters
 * @brief Estimate initial parameters for the EM algorithm
 *
 * @signature void estimateInitialParameters(mu1, mu2, lambda1, lambda2, dataset, indexValidSamplesTSS, indexValidSamplesBG)
 *
 * @param[in,out]       mu1                         Mean of TSS examples (first estimate of mean of first distribution)
 * @param[in,out]       mu2                         Mean of background examples (first estimate of mean of second distribution)
 * @param[in,out]       lambda1                     SD of TSS examples (first estimate of SD of first distribution)
 * @param[in,out]       lambda2                     SD of background examples (first estimate of SD of second distribution)
 * @param[in]           dataset                     Dataset for the EM algorithm
 * @param[in]           indexValidSamplesTSS        Index of samples in dataset that pass strongly indicated TSS requirements
 * @param[in]           indexValidSamplesBG         Index of samples in dataset that pass strongly indicated background requirements
 */
void estimateInitialParameters(double & mu1, double & mu2, double & lambda1, double & lambda2, std::vector<DatasetRecord> const & dataset,
std::vector<unsigned> const & indexValidSamplesTSS, std::vector<unsigned> const & indexValidSamplesBG)
{
    std::vector<double> samplesTSS, samplesBG;
    samplesTSS.reserve(1000);
    samplesBG.reserve(1000);

    for (unsigned i = 1000; i < 2000; i++) // First 1000 random samples used for beta parameter estimation, ensured that number of indices is at least 2000 previously
        samplesTSS.push_back(dataset[indexValidSamplesTSS[i]].tagCount);

    for (unsigned i = 1000; i < 2000; i++) // First 1000 random samples used for beta parameter estimation, ensured that number of indices is at least 2000 previously
        samplesBG.push_back(dataset[indexValidSamplesBG[i]].tagCount);

    mu1 = average(samplesTSS);
    mu2 = average(samplesBG);

    double sum1 = 0;

    for (auto i : samplesTSS)
    {
        sum1 += ((double) 1 / i - (double) 1 / mu1);
    }

    lambda1 = (double) 1 / (sum1 / (double) samplesTSS.size());

    double  sum2 = 0;

    for (auto i : samplesBG)
    {
        sum2 += ((double) 1 / i - (double) 1 / mu2);
    }
    lambda2 = (double) 1 / (sum2 / (double) samplesBG.size());
}

/*!
 * @fn eStep
 * @brief Expectation step of the EM algorithm
 *
 * @signature void eStep(outputForMStep, dataset, mu1, mu2, lambda1, lambda2, betaParams)
 *
 * @param[in,out]       outputForMStep              Current probabilities for regions of being promoter or background and parameters
 * @param[in]           dataset                     Dataset for the EM algorithm
 * @param[in]           mu1                         Current estimate of mean of first distribution (promoters)
 * @param[in]           mu2                         Current estimate of mean of second distribution (background)
 * @param[in]           lambda1                     Current estimate of SD of first distribution (promoters)
 * @param[in]           lambda2                     Current estimate of SD of second distribution (background)
 * @param[in]           betaParams                  Estimated beta parameters
 */
void eStep(std::vector<std::vector<double> > & outputForMStep, std::vector<DatasetRecord> const & dataset,
double const mu1, double const mu2, double const lambda1, double const lambda2, std::vector<double> const & betaParams)
{
    outputForMStep.resize(dataset.size());

    for (unsigned i = 0; i < dataset.size(); i++)
    {
        outputForMStep[i].resize(6); // z1, z2, f1, f2, pi1, pi2

        double f1 = probDensityFunction(log(dataset[i].tagCount + 2), mu1, lambda1);
        double f2 = probDensityFunction(log(dataset[i].tagCount + 2), mu2, lambda2);

        // Discretizing the CpG content
        // CpG content set to 0.1 or 0.9 for background or TSS respectively
        double pi1;

        if (dataset[i].cpgContent < 0.45)
        {
            pi1 = (double) 1 / (1 + exp(-(betaParams[0] + betaParams[1] * 0.1 + betaParams[2] * dataset[i].conservation + betaParams[3] * dataset[i].affinity + betaParams[4] * dataset[i].distScore)));
        }
        else
        {
            pi1 = (double) 1 / (1 + exp(-(betaParams[0] + betaParams[1] * 0.9 + betaParams[2] * dataset[i].conservation + betaParams[3] * dataset[i].affinity + betaParams[4] * dataset[i].distScore)));
        }

        double pi2 = 1 - pi1;

        if (dataset[i].bg)
        {
            outputForMStep[i][0] = 0; // Probability of being a promoter
            outputForMStep[i][1] = 1; // Probability of being a background region
        }
        else
        {
            outputForMStep[i][0] = ((double) f1 * pi1) / ((double) f1 * pi1 + (double) f2 * pi2); // Probability of being a promoter
            outputForMStep[i][1] = ((double) f2 * pi2) / ((double) f1 * pi1 + (double) f2 * pi2); // Probability of being a background region
        }

        outputForMStep[i][2] = f1;
        outputForMStep[i][3] = f2;
        outputForMStep[i][4] = pi1;
        outputForMStep[i][5] = pi2;
    }
}

/*!
 * @fn mStep
 * @brief Maximization step of the EM algorithm
 *
 * @signature void mStep(outputParams, outputForMStep, dataset)
 *
 * @param[in,out]       outputParams                New estimates of mu1, mu2, lambda1 and lambda2
 * @param[in]           outputForMStep              Current probabilities for regions of being promoter or background and parameters (output of E-Step)
 * @param[in]           dataset                     Dataset for the EM algorithm
 *
 * Re-estimate the parameters based on the marginal probability distributions.
 */
void mStep(std::vector<double> & outputParams, std::vector<std::vector<double> > const & outputForMStep, std::vector<DatasetRecord> const & dataset)
{
    outputParams.resize(4); // mu1, mu2, lambda1, lambda2

    // Re-estimate mu1
    double muNum1 = 0;
    double muDen1 = 0;
    for (unsigned i = 0; i < dataset.size(); i++)
    {
        muNum1 += (log(dataset[i].tagCount + 2) * outputForMStep[i][0]);
        muDen1 += outputForMStep[i][0];
    }

    outputParams[0] = muNum1 / muDen1;

    // Re-estimate mu2
    double muNum2 = 0;
    double muDen2 = 0;
    for (unsigned i = 0; i < dataset.size(); i++)
    {
        muNum2 += (log(dataset[i].tagCount + 2) * outputForMStep[i][1]);
        muDen2 += outputForMStep[i][1];
    }

    outputParams[1] = muNum2 / muDen2;

    // Re-estimate lambda1
    double lambdaNum1 = 0;
    double lambdaDen1 = 0;

    for (unsigned i = 0; i < dataset.size(); i++)
    {
        lambdaNum1 += outputForMStep[i][0];
        lambdaDen1 += (outputForMStep[i][0] * (pow((log(dataset[i].tagCount + 2) - outputParams[0]), 2) / log(dataset[i].tagCount + 2))); // TODO: Avoid recalculation of log
    }

    outputParams[2] = (lambdaNum1 / lambdaDen1) * pow(outputParams[0], 2);

    // Re-estimate lambda2
    double lambdaNum2 = 0;
    double lambdaDen2 = 0;

    for (unsigned i = 0; i < dataset.size(); i++)
    {
        lambdaNum2 += outputForMStep[i][1];
        lambdaDen2 += (outputForMStep[i][1] * (pow((log(dataset[i].tagCount + 2) - outputParams[1]), 2) / log(dataset[i].tagCount + 2))); // TODO: Avoid recalculation of log
    }

    outputParams[3] = (lambdaNum2 / lambdaDen2) * pow(outputParams[1], 2);
}

/*!
 * @fn calculateLikelihood
 * @brief Calculate the likelihood
 *
 * @signature double calculateLikelihood(dataset, outputForMStep, outputParams)
 *
 * @param[in]           dataset                     Dataset for the EM algorithm
 * @param[in]           outputForMStep              Current probabilities for regions of being promoter or background and parameters (output of E-Step)
 * @param[in]           outputParams                Current estimates of mu1, mu2, lambda1 and lambda2
 */
double calculateLikelihood(std::vector<DatasetRecord> const & dataset,
std::vector<std::vector<double> > const & outputForMStep, std::vector<double> const & outputParams)
{
    double likelihood = 0;
    for (unsigned i = 0; i < dataset.size(); i++)
    {
        double f1 = probDensityFunction(log(dataset[i].tagCount + 2), outputParams[0], outputParams[2]);
        double f2 = probDensityFunction(log(dataset[i].tagCount + 2), outputParams[1], outputParams[3]);

        likelihood += log((outputForMStep[i][4] * f1) + (outputForMStep[i][5] * f2));
    }
    return likelihood;
}

/*!
 * @fn countPromoters
 * @brief Count number of regions identified as promoters
 *
 * @signature void countPromoters(counts, outputForMStep)
 *
 * @param[in,out]       counts              Counts of promoters and background regions
 * @param[in]           outputForMStep      Final probabilities of being background or promoter for all dataset regions
 */
void countPromoters(std::pair<unsigned, unsigned> & counts, std::vector<std::vector<double> > const & outputForMStep)
{
    counts.first = 0; // Promoters
    counts.second = 0; // Background

    for (unsigned i = 0; i < outputForMStep.size(); i++)
    {
        if (outputForMStep[i][0] > outputForMStep[i][1])
            counts.first++;
        if (outputForMStep[i][0] <= outputForMStep[i][1])
            counts.second++;
    }
}

/*!
 * @fn emAlgorithm
 * @brief The EM algorithm that classifies promoters and background sequences
 *
 * @signature void emAlgorithm(outputParams, outputForMStep, dataset, mu1, mu2, lambda1, lambda2, betaParams)
 *
 * @param[in,out]       outputParams                Final estimates for mu1, mu2, lambda1 and lambda2 ans calculated by the EM
 * @param[in,out]       outputForMStep              Probabilities for regions of being promoter or background and parameters
 * @param[in]           dataset                     Dataset for the EM algorithm
 * @param[in]           mu1                         First estimate of mean of first distribution (promoters)
 * @param[in]           mu2                         First estimate of mean of second distribution (background)
 * @param[in]           lambda1                     First estimate of SD of first distribution (promoters)
 * @param[in]           lambda2                     First estimate of SD of second distribution (background)
 * @param[in]           betaParams                  Estimated beta parameters
 */
void emAlgorithm(std::vector<double> & outputParams, std::vector<std::vector<double> > & outputForMStep,
std::vector<DatasetRecord> const & dataset, double const mu1, double const mu2, double const lambda1,
double const lambda2, std::vector<double> const & betaParams)
{
    eStep(outputForMStep, dataset, mu1, mu2, lambda1, lambda2, betaParams);
    mStep(outputParams, outputForMStep, dataset);

    double likelihood = -DBL_MAX;
    double newLikelihood = calculateLikelihood(dataset, outputForMStep, outputParams);
    std::cout.precision(2);
    std::cout << std::fixed << newLikelihood << std::endl;

    while ((newLikelihood - likelihood) > 0.1)
    {
        likelihood = newLikelihood;
        eStep(outputForMStep, dataset, outputParams[0], outputParams[1], outputParams[2], outputParams[3], betaParams);
        mStep(outputParams, outputForMStep, dataset);
        newLikelihood = calculateLikelihood(dataset, outputForMStep, outputParams);
        std::cout.precision(2);
        std::cout << std::fixed << newLikelihood << std::endl;
    }

    std::pair<unsigned, unsigned> counts;
    countPromoters(counts, outputForMStep);
    std::cout << "Regions classified as promoter: " << counts.first << std::endl;
    std::cout << "Regions classified as background: " << counts.second << std::endl;
}
