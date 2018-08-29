// ============================================================================
// Functions and objects needed for an interval tree
// ============================================================================
// Author: Sara Hetzel (hetzel@molgen.mpg.de)

#pragma once

#include <iostream>

#include "structs.h"

// Comparators needed for equality of chromosome and strand
// Comparators for GenomicLocation only
// Needed for background overlap with genes and repeats
// Also needed for CAGE overlap with gene starts
inline unsigned low(GenomicLocation & gl)
{
    return gl.start;
}

inline unsigned high(GenomicLocation & gl)
{
    return gl.end;
}

// Comparators for GenomicLocation and miRNA or background region
// Needed for initial overlap of miRNA or background region with CAGE tags
template<typename TSeq>
inline unsigned low(std::pair<GenomicLocation, TSeq> & gl)
{
    return low(gl.first);
}

template<typename TSeq>
inline unsigned high(std::pair<GenomicLocation, TSeq> & gl)
{
    return high(gl.first);
}

template<typename TAddInfo>
inline unsigned low(std::pair<TAddInfo, GenomicLocation> & gl)
{
    return low(gl.second);
}

template<typename TAddInfo>
inline unsigned high(std::pair<TAddInfo, GenomicLocation> & gl)
{
    return high(gl.second);
}

// Node in interval tree
template<typename TInterval>
struct IntTreeNode
{
    TInterval nodeInterval;
    unsigned max;
    IntTreeNode *left, *right;
};

// Create new node
template<typename TInterval>
IntTreeNode<TInterval> * newNode(TInterval & currentInterval)
{
    IntTreeNode<TInterval> *temp = new IntTreeNode<TInterval>;
    temp->nodeInterval = currentInterval;
    temp->max = high(temp->nodeInterval);
    temp->left = temp->right = NULL;
    return temp;
};


// Insert a new node into an interval tree
template<typename TInterval>
IntTreeNode<TInterval> * insertImpl(IntTreeNode<TInterval> *root, TInterval & currentInterval)
{
    // If tree is empty, new node becomes root
    if (root == NULL)
    {
        return newNode(currentInterval);
    }
    // Get low value of interval at root
    unsigned currentLowVal = low(root->nodeInterval);

    // If low value of root is smaller, go to left subtree
    if (low(currentInterval) < currentLowVal)
        root->left = insertImpl(root->left, currentInterval);

    // Else, go to right subtree
    else
        root->right = insertImpl(root->right, currentInterval);

    // Update the max value of the ancestor if needed
    if (root->max < high(currentInterval))
        root->max = high(currentInterval);

    return root;
}

// Overlap of two intervals
template<typename TInterval>
bool overlapInterval(TInterval & interval1, TInterval &interval2)
{
    if (low(interval1) < high(interval2) && low(interval2) < high(interval1))
    {
        return true;
    }

    return false;
}

// Search interval in tree
template<typename TInterval>
TInterval * overlapSearchSingleImpl(IntTreeNode<TInterval> *root, TInterval & currentInterval)
{
    // Tree is empty
    if (root == NULL) return NULL;

    // If given interval overlaps with root
    if (overlapInterval(root->nodeInterval, currentInterval))
        return &(root->nodeInterval);

    // If left child of root is present and max of left child is
    // greater than or equal to given interval, then currentInterval may
    // overlap with an interval in left subtree
    if (root->left != NULL && root->left->max >= low(currentInterval))
        return overlapSearchSingleImpl(root->left, currentInterval);

    // Else interval can only overlap with right subtree
    return overlapSearchSingleImpl(root->right, currentInterval);
}

template<typename TInterval>
void overlapSearchAllImpl(IntTreeNode<TInterval> *root, TInterval & currentInterval, std::vector<TInterval*> & allOverlaps)
{
    // Tree is empty
    if (root == NULL) return;

    // If given interval overlaps with root, must search both subtrees
    if (overlapInterval(root->nodeInterval, currentInterval))
    {
        allOverlaps.push_back(&(root->nodeInterval));
        overlapSearchAllImpl(root->left, currentInterval, allOverlaps);
        overlapSearchAllImpl(root->right, currentInterval, allOverlaps);
        return;
    }

    // If left child of root is present and max of left child is
    // greater than or equal to given interval, then currentInterval may
    // overlap with an interval in left subtree
    if (root->left != NULL && root->left->max >= low(currentInterval))
    {
        overlapSearchAllImpl(root->left, currentInterval, allOverlaps);
    }

    // Else interval can only overlap with right subtree
    overlapSearchAllImpl(root->right, currentInterval, allOverlaps);
}

// Delete tree
template<typename TInterval>
void deleteTree(IntTreeNode<TInterval> *root)
{
  if(root != 0) {
        deleteTree(root->left);
        deleteTree(root->right);
        delete root;
     }
}

template<typename TInterval>
struct IntTree
{
    IntTreeNode<TInterval>* root;
    IntTree() : root(NULL) {}
    ~IntTree() { deleteTree(root); }

private:
    IntTree(IntTree const &);
    IntTree & operator=(const IntTree &);

public:
    void insert(TInterval & val) {root = insertImpl(root, val);};
    TInterval * overlapSearchSingle(TInterval & val) {return overlapSearchSingleImpl(root, val);}
    void overlapSearchAll(TInterval & val, std::vector<TInterval*> & allOverlaps) {overlapSearchAllImpl(root, val, allOverlaps);}
    void traverseIntTree() {traverseIntTreeImpl(root);};
};

template<typename TInterval>
void traverseIntTreeImpl(IntTreeNode<TInterval> *root)
{
    if (root == NULL) return;

    traverseIntTreeImpl(root->left);

    std::cout << "[" << low(root->nodeInterval) << ", " << high(root->nodeInterval) << "]" << std::endl;
    traverseIntTreeImpl(root->right);
}
