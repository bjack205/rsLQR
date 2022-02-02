#include "binary_tree.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

static BinaryNode* BuildSubTree(BinaryNode* start, int len) {
  int mid = (len + 1) / 2;
  int new_len = mid - 1;
  if (len > 1) {
    BinaryNode* root = start + new_len;
    BinaryNode* left_root = BuildSubTree(start, new_len);
    BinaryNode* right_root = BuildSubTree(start + mid, new_len);
    root->left_child = left_root;
    root->right_child = right_root;
    left_root->parent = root;
    right_root->parent = root;
    root->left_inds.start = left_root->left_inds.start;
    root->left_inds.stop = left_root->right_inds.stop;
    root->right_inds.start = right_root->left_inds.start;
    root->right_inds.stop = right_root->right_inds.stop;
    root->level = left_root->level + 1;
    return root;
  } else {
    int k = start->idx;
    start->left_inds.start = k;
    start->left_inds.stop = k;
    start->right_inds.start = k + 1;
    start->right_inds.stop = k + 1;
    start->level = 0;
    start->left_child = NULL;
    start->right_child = NULL;
    return start;
  }
}

OrderedBinaryTree ndlqr_BuildTree(int nhorizon) {
  assert(IsPowerOfTwo(nhorizon));

  BinaryNode* node_list = (BinaryNode*)malloc(nhorizon * sizeof(*node_list));
  if (!node_list) {
  }
  OrderedBinaryTree tree;
  for (int i = 0; i < nhorizon; ++i) {
    node_list[i].idx = i;
  }
  tree.num_elements = nhorizon;
  tree.node_list = node_list;
  tree.depth = LogOfTwo(nhorizon);

  // Build the tree
  tree.root = BuildSubTree(node_list, nhorizon - 1);

  return tree;
}

int ndlqr_FreeTree(OrderedBinaryTree* tree) {
  if (!tree) return -1;
  free(tree->node_list);
  return 0;
}

int ndlqr_GetIndexFromLeaf(const OrderedBinaryTree* tree, int leaf, int level) {
  (void)tree;
  int linear_index = PowerOfTwo(level) * (2 * leaf + 1) - 1;
  return linear_index;
}

int ndlqr_GetIndexLevel(const OrderedBinaryTree* tree, int index) {
  return tree->node_list[index].level;
}

const BinaryNode* GetNodeAtLevel(const BinaryNode* node, int index, int level) {
  if (node->level == level) {
    return node;
  } else if (node->level > level) {
    if (index <= node->idx) {
      return GetNodeAtLevel(node->left_child, index, level);
    } else {
      return GetNodeAtLevel(node->right_child, index, level);
    }
  } else {
    return GetNodeAtLevel(node->parent, index, level);
  }
}

int ndlqr_GetIndexAtLevel(const OrderedBinaryTree* tree, int index, int level) {
  if (!tree) return -1;
  if (index < 0 || index >= tree->num_elements) {
    fprintf(stderr, "ERROR: Invalid index (%d). Should be between %d and %d.\n", index, 0,
            tree->num_elements - 1);
  }
  if (level < 0 || level >= tree->depth) {
    fprintf(stderr, "ERROR: Invalid level (%d). Should be between %d and %d.\n", level, 0,
            tree->depth - 1);
  }
  BinaryNode* node = tree->node_list + index;
  if (index == tree->num_elements - 1) {
    node = tree->node_list + index - 1;
  }
  const BinaryNode* base_node = GetNodeAtLevel(node, index, level);
  return base_node->idx;
  // return base_node->idx;
}
