#pragma once
#include <stdbool.h>

typedef struct {
  int start;  // inclusive
  int stop;   // exclusive
} UnitRange;

typedef struct BinaryNode_s BinaryNode;
struct BinaryNode_s {
  int idx;
  int level;
  int levelidx;
  UnitRange left_inds; 
  UnitRange right_inds;

  BinaryNode* parent;
  BinaryNode* left_child;
  BinaryNode* right_child;
};

typedef struct {
  BinaryNode* root;
  BinaryNode* node_list;
  int num_elements;
  int depth;
} OrderedBinaryTree;

OrderedBinaryTree ndlqr_BuildTree(int N);
int ndlqr_FreeTree(OrderedBinaryTree* tree);

int ndlqr_GetIndexFromLeaf(const OrderedBinaryTree* tree, int leaf, int level);
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
