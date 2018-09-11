// ============================================================================
// Struct definitions used by various functions in PROmiRNA
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#pragma once

#include <string>
#include <vector>

#include <seqan/sequence.h>

using namespace seqan;

/*!
 * @struct GenomicLocation
 * @brief Stores information of a genomic location including chromosome, start, end and strand.
 *
 * @signature struct GenomicLocation
 *
 * Start and end position are 0-based with half open intervals: [start,stop)
 */
struct GenomicLocation
{
    unsigned start;
    unsigned end;
    CharString chr;
    char strand;
};

// Equality operator for GenomicLocation
inline bool operator==(GenomicLocation const & gl1, GenomicLocation const & gl2)
{
    return gl1.chr == gl2.chr &&
           gl1.strand == gl2.strand &&
           gl1.start == gl2.start &&
           gl1.end == gl2.end;
}

// Less operator for GenomicLocation
inline bool operator<(GenomicLocation const & gl1, GenomicLocation const & gl2)
{
    return gl1.strand < gl2.strand ||
           (gl1.strand == gl2.strand && gl1.chr < gl2.chr) ||
           (gl1.strand == gl2.strand && gl1.chr == gl2.chr && gl1.start < gl2.start) ||
           (gl1.strand == gl2.strand && gl1.chr == gl2.chr && gl1.start == gl2.start && gl1.end < gl2.end);
}

/*!
 * @struct MiRNA
 * @brief Stores information about a miRNA except their genomic location.
 *
 * @signature struct MiRNA
 *
 * Stores information about classification which includes all found classifications (could be multiple).
 * In addition, ID, accession and name are stored.
 * This struct does not contain a genomic location since for the interval trees a genomic location is needed as key.
 * Also for comparisons with other tags it is important that the genomic location is stores separately.
 */
struct MiRNA
{
    std::set<CharString> classification; // Can have multiple classifications if intragenic
	CharString id;
	CharString accession;
	CharString name;
};

// Equality operator for MiRNA
inline bool operator==(MiRNA const & m1, MiRNA const & m2)
{
    return m1.classification == m2.classification &&
           m1.id == m2.id &&
           m1.accession == m2.accession &&
           m1.name == m2.name;
}

inline bool operator!=(MiRNA const & m1, MiRNA const & m2)
{
    return m1.classification != m2.classification ||
           m1.id != m2.id ||
           m1.accession != m2.accession ||
           m1.name != m2.name;
}

/*!
 * @class MatrixSingle
 * @brief Matrix with one value per cell, row and column names.
 *
 * @signature class MatrixSingle
 */
class MatrixSingle
{
public:
    std::vector<std::string> rownames;
    std::vector<std::string> colnames;
    std::vector<std::vector<float> > matrix;
    MatrixSingle(){};
    MatrixSingle(unsigned n1, unsigned n2)
    {
        rownames.resize(n1);
        colnames.resize(n2);
        matrix.resize(n1);
        for (unsigned i = 0; i < matrix.size(); i++)
            matrix[i].resize(n2);
    };
    MatrixSingle(MatrixSingle const &) = default;
    MatrixSingle(MatrixSingle &&) = default;
    MatrixSingle & operator =(MatrixSingle const &) = default;
    MatrixSingle & operator =(MatrixSingle &&) = default;
    ~MatrixSingle() = default;
    void ZERO();
};

/*!
 * @class MatrixPair
 * @brief Matrix with two values (pair) per cell, row and column names.
 *
 * @signature class MatrixPair
 *
 * Used for quantile normalization where the matrix needs to store the order of initial rows as well.
 */
class MatrixPair
{
public:
    std::vector<std::string> rownames;
    std::vector<std::string> colnames;
    std::vector<std::vector<std::pair<float, unsigned> > > matrix;
    MatrixPair(){};
    MatrixPair(unsigned n1, unsigned n2)
    {
        rownames.resize(n1);
        colnames.resize(n2);
        matrix.resize(n1);
        for (unsigned i = 0; i < matrix.size(); i++)
            matrix[i].resize(n2);
    };
    MatrixPair(MatrixPair const &) = default;
    MatrixPair(MatrixPair &&) = default;
    MatrixPair & operator =(MatrixPair const &) = default;
    MatrixPair & operator =(MatrixPair &&) = default;
    ~MatrixPair() = default;
    void ZERO();
};

// Equality operator for MatrixPair
inline bool operator==(MatrixPair const & mp1, MatrixPair const & mp2)
{
    return mp1.matrix == mp2.matrix &&
           mp1.rownames == mp2.rownames &&
           mp1.colnames == mp2.colnames;
}

inline bool operator!=(MatrixPair const & mp1, MatrixPair const & mp2)
{
    return mp1.matrix != mp2.matrix ||
           mp1.rownames != mp2.rownames ||
           mp1.colnames != mp2.colnames;
}

/*!
 * @struct DatasetRecord
 * @brief Stores information of a dataset record later used for the EM algorithm.
 *
 * @signature struct DatasetRecord
 *
 * Contains information about genomic location, the tag count (from normalized count matrix), CpG content, conservation, affinity, distance score,
 * number of libraries the sequence appeared and whether the sequence is a background region or not.
 */
struct DatasetRecord
{
    GenomicLocation location;
    double tagCount;    // These are the added values from the normalized count matrix
    double cpgContent;
    double conservation;
    double affinity;
    double distScore;   // Needed for create dataset
    unsigned distance;  // Actual distance between TSS and miRNA
    unsigned libs;
	unsigned numCons=0;
    bool bg;
    std::string miRNA;
};
