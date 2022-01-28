#include "binary_tree.h"
#include "test/minunit.h"

int TestBuildTree() {
  OrderedBinaryTree tree = ndlqr_BuildTree(8);
  mu_assert(tree.num_elements == 8);
  mu_assert(tree.root->idx == 3);
  mu_assert(tree.root->level == 2);
  mu_assert(tree.root->left_inds.start == 0);
  mu_assert(tree.root->left_inds.stop == 3);
  mu_assert(tree.root->right_inds.start == 4);
  mu_assert(tree.root->right_inds.stop == 7);
  mu_assert(tree.root->left_child->idx == 1);
  mu_assert(tree.root->left_child->level == 1);
  mu_assert(tree.root->left_child->left_child->idx == 0);
  mu_assert(tree.root->left_child->left_child->level == 0);
  mu_assert(tree.root->right_child->idx == 5);
  mu_assert(tree.root->right_child->right_child->idx == 6);
  ndlqr_FreeTree(&tree);
  return 1;
}

int GetIndexLevel() {
  OrderedBinaryTree tree = ndlqr_BuildTree(8);
  mu_assert(ndlqr_GetIndexLevel(&tree, 0) == 0); 
  mu_assert(ndlqr_GetIndexLevel(&tree, 1) == 1);
  mu_assert(ndlqr_GetIndexLevel(&tree, 2) == 0);
  mu_assert(ndlqr_GetIndexLevel(&tree, 3) == 2);
  mu_assert(ndlqr_GetIndexLevel(&tree, 4) == 0);
  mu_assert(ndlqr_GetIndexLevel(&tree, 5) == 1);
  mu_assert(ndlqr_GetIndexLevel(&tree, 6) == 0);
  ndlqr_FreeTree(&tree);
  return 1;
}

int GetIndexAtLevel() {
  OrderedBinaryTree tree = ndlqr_BuildTree(8);
  int index = 5;
  int level = 0;
  // int base_index = tree.depth;
  int base_index = ndlqr_GetIndexAtLevel(&tree, index, level);
  mu_assert(base_index == 4);

  index = 3;
  base_index = ndlqr_GetIndexAtLevel(&tree, index, level);
  mu_assert(base_index == 2);

  index = 2;
  level = 2;
  base_index = ndlqr_GetIndexAtLevel(&tree, index, level);
  mu_assert(base_index == 3);

  index = 7;
  base_index = ndlqr_GetIndexAtLevel(&tree, index, level);
  mu_assert(base_index == 3);

  level = 0;
  base_index = ndlqr_GetIndexAtLevel(&tree, index, level);
  mu_assert(base_index == 6);
  return 1;
}
 
void AllTests() {
  mu_run_test(TestBuildTree);
  mu_run_test(GetIndexLevel);
  mu_run_test(GetIndexAtLevel);
}

mu_test_main
