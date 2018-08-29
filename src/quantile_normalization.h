// ============================================================================
// Functions for quantile normalization of a matrix
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#pragma once

#include <vector>

#include "structs.h"

using namespace seqan;

/*!
 * @fn transpose
 * @brief Transpose a matrix pair, return type matrix pair
 *
 * @signature void transpose(outputMatrix, inputMatrix)
 *
 * @param[in,out]   outputMatrix        Resulting transposed matrix
 * @param[in]       inputMatrix         Input matrix
 */
void transpose(MatrixPair & outputMatrix, MatrixPair const & inputMatrix)
{
    outputMatrix.matrix.resize(inputMatrix.matrix[0].size());
    for (unsigned i = 0; i < outputMatrix.matrix.size(); i++)
        outputMatrix.matrix[i].resize(inputMatrix.matrix.size());

    outputMatrix.colnames = inputMatrix.rownames;
	outputMatrix.rownames = inputMatrix.colnames;

	for (unsigned i = 0; i < inputMatrix.matrix[0].size(); i++)
    {
		for (unsigned j = 0; j < inputMatrix.matrix.size(); j++)
        {
			outputMatrix.matrix[i][j].first = inputMatrix.matrix[j][i].first;
			outputMatrix.matrix[i][j].second = j;
		}
	}
}

/*!
 * @fn transpose
 * @brief Transpose a matrix pair, return type matrix
 *
 * @signature void transpose(outputMatrix, inputMatrix)
 *
 * @param[in,out]   outputMatrix        Resulting transposed matrix
 * @param[in]       inputMatrix         Input matrix
 */
void transpose(MatrixSingle & outputMatrix, MatrixPair const & inputMatrix)
{
    outputMatrix.matrix.resize(inputMatrix.matrix[0].size());
    for (unsigned i = 0; i < outputMatrix.matrix.size(); i++)
        outputMatrix.matrix[i].resize(inputMatrix.matrix.size());

    outputMatrix.colnames = inputMatrix.rownames;
	outputMatrix.rownames = inputMatrix.colnames;

	for (unsigned i = 0; i < inputMatrix.matrix[0].size(); i++)
    {
		for (unsigned j = 0; j < inputMatrix.matrix.size(); j++)
        {
			outputMatrix.matrix[i][j] = inputMatrix.matrix[j][i].first;
		}
	}
}

/*!
 * @fn quantileNormalization
 * @brief Perform quantile normalization of a matrix
 *
 * @signature void quantileNormalization(normMatrix, inputMatrix)
 *
 * @param[in,out]   normMatrix        Resulting matrix
 * @param[in]       inputMatrix       Input matrix
 */
void quantileNormalization(MatrixSingle & normMatrix, MatrixPair const & inputMatrix)
{
    // Initialize output matrix
	MatrixPair transMatrix;
    transpose(transMatrix, inputMatrix);

    // Sort columns
	for (unsigned i = 0; i < transMatrix.matrix.size(); i++)
    {
		std::sort(transMatrix.matrix[i].begin(), transMatrix.matrix[i].end());
	}

	// Calculate row means
    std::vector<float> rowMeans;
	for (unsigned i = 0; i < transMatrix.matrix[0].size(); i++)
    {
		float currMean = 0;
		for (unsigned j = 0; j < transMatrix.matrix.size(); j++)
        {
			currMean += transMatrix.matrix[j][i].first;
		}
		currMean /= (float) transMatrix.matrix.size();
		rowMeans.push_back(currMean);
	}

    // Replace entries with mean
    for (unsigned i = 0; i < rowMeans.size(); i++)
    {
        for (unsigned j = 0; j < transMatrix.matrix.size(); j++)
        {
            transMatrix.matrix[j][i].first = rowMeans[i];
        }
    }

    // Sort according to original order
    for (unsigned i = 0; i < transMatrix.matrix.size(); i++)
    {
        std::sort(transMatrix.matrix[i].begin(), transMatrix.matrix[i].end(), [](auto &n1, auto &n2) { return n1.second < n2.second; });
    }

    transpose(normMatrix, transMatrix);
}
