/**
 * @file binary_tree.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Binary tree for rsLQR algorithm
 * @version 0.1
 * @date 2022-02-01
 * 
 * @copyright Copyright (c) 2022
 * 
 * @addtogroup rsLQR 
 * @{
 */
#pragma once
#include <stdbool.h>

/**
 * @brief Represents a range of consecutive integers 
 * 
 * Abstract representation of the set of consecutive integers from
 * UnitRange::start to UnitRange::stop.
 * 
 */
typedef struct {
  int start;  // inclusive
  int stop;   // exclusive
} UnitRange;

/**
 * @brief One of the nodes in the binary tree for the rsLQR solver
 * 
 */
typedef struct BinaryNode_s BinaryNode;
struct BinaryNode_s {
  int idx;               ///< knot point index
  int level;             ///< level in the tree
  int levelidx;          ///< leaf index at the current level
  UnitRange left_inds;   ///< range of knot point indices of all left children 
  UnitRange right_inds;  ///< range of knot point indices of all right children

  BinaryNode* parent;       ///< parent node
  BinaryNode* left_child;   ///< left child
  BinaryNode* right_child;  ///< right child
};

/**
 * @brief The binary tree for the rsLQR solver
 * 
 * Caches useful information to speed up some of the computations during the 
 * rsLQR solver, mostly around converting form linear knot point indices to 
 * the hierarchical indexing of the algorithm.
 * 
 * ## Construction and destruction
 * A new tree can be constructed solely given the length of the time horizon, which must 
 * be a power of two. Use ndlqr_BuildTree(N) to build a new tree, which can be 
 * de-allocated using ndlqr_FreeTree().
 * 
 * ## Methods
 * - ndlqr_BuildTree()
 * - ndlqr_FreeTree()
 * - ndlqr_GetIndexFromLeaf()
 * - ndlqr_GetIndexLevel()
 * - ndlqr_GetIndexAtLevel()
 */
typedef struct {
  BinaryNode* root;       ///< root of the tree. Corresponds to the "middle" knot point.
  BinaryNode* node_list;  ///< a list of all the nodes, ordered by their knot point index
  int num_elements;       ///< length of the OrderedBinaryTree::node_list
  int depth;              ///< total depth of the tree
} OrderedBinaryTree;

/**
 * @brief Construct a new binary tree for a horizon of length @p N
 * 
 * Must be paired with a corresponding call to ndlqr_FreeTree().
 * 
 * @param  N horizon length. Must be a power of 2.
 * @return A new binary tree
 */
OrderedBinaryTree ndlqr_BuildTree(int N);

/**
 * @brief Frees the data in @p tree
 * 
 * @param tree An initialized binary tree
 * @return 0 if successful
 */
int ndlqr_FreeTree(OrderedBinaryTree* tree);

/**
 * @brief Get the knot point index given the leaf index at a given level.
 * 
 * @param tree  An initialized binary tree for the problem horizon
 * @param leaf  Leaf index
 * @param level Level of tree from which to get the index
 * @return
 */
int ndlqr_GetIndexFromLeaf(const OrderedBinaryTree* tree, int leaf, int level);

/**
 * @brief Get the level for a given knot point index
 * 
 * @param tree  An initialized binary tree for the problem horizon
 * @param index Knot point index
 * @return      The level for the given index
 */
int ndlqr_GetIndexLevel(const OrderedBinaryTree* tree, int index);

/**
 * @brief Get the index in `level` that corresponds to `index`.
 * 
 * If the level is higher than the level of the given index, it's simply the parent 
 * at that level. If it's lower, then it's the index that's closest to the given one, with
 * ties broken by choosing the left (or smaller) of the two.
 * 
 * @param tree  Precomputed binary tree
 * @param index Start index of the search. The result will be the index closest to this index.
 * @param level The level in which the returned index should belong to.
 * @return int  The index closest to the provided one, in the given level. -1 if unsucessful.
 */
int ndlqr_GetIndexAtLevel(const OrderedBinaryTree* tree, int index, int level);

/**@} */
