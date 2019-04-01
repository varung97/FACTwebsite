#ifndef RMQ_H_
#define RMQ_H_

class RMQ {
public:
    // FreqDiffTree must be ordered such that for any node, the first child is the one with largest leafset
    RMQ(FreqDiffTree* tree, lca_t* lca_prep);

    // Return max weight along path from v to w, excluding w
    // w must be a proper ancestor of v
    int max_weight_path(FreqDiffTree::Node* v, FreqDiffTree::Node* w);

    // Pointer to lca structure for tree
    lca_t* lca_prep;

private:
    // Pointer from each node to the root of its centroid path
    std::vector<FreqDiffTree::Node*> cp_roots;

    // Pointers to rmq for each of the centroid paths
    std::vector<gen_rmq_t*> cp_rmqs;

    // Index of the centroid path each node belongs to
    std::vector<int> cp_indices;

    // Depth of each centroid path (defined as number of centroid path roots that are proper ancestors of the root of this path)
    std::vector<int> cp_depths;

    // Depth of each node in the tree
    std::vector<int> depths;

    // Pointers from each node to the rmq for some leaf in its leafset
    std::vector<gen_rmq_t*> leaf_rmq_for_nodes;

    // Numbering leaves from left to right, pointer to the leftmost leaf in the subtree of each node
    std::vector<FreqDiffTree::Node*> left_most_leaf;

    // Numbering leaves from left to right, index of each leaf
    std::vector<int> taxa_l2r_indices;

    // For each node, for each leaf, a pointer to the child of the node that is an ancestor of the leaf
    std::vector<std::vector<FreqDiffTree::Node*>> child_for_leaf_vectors;

    // Index of smallest leaf in side trees of node
    std::vector<int> child_for_leaf_offset;

    // Pointer to the heaviest child of each node
    std::vector<FreqDiffTree::Node*> heaviest_child;

    void preprocess_depths(FreqDiffTree* tree);
    void preprocess_centroid_paths(FreqDiffTree* tree);
    void preprocess_leaf_rmqs(FreqDiffTree* tree);
    int preprocess_child_for_leaf_helper(FreqDiffTree::Node* node, int first_available_index);
    void preprocess_child_for_leaf(FreqDiffTree* tree);

    int rmq_q1(FreqDiffTree::Node* v, FreqDiffTree::Node* w);
    int rmq_q2_qg1(FreqDiffTree::Node* v, FreqDiffTree::Node*w);
    int rmq_qg(FreqDiffTree::Node* v, FreqDiffTree::Node* w);
    int rmq(FreqDiffTree::Node* v, FreqDiffTree::Node* w);
    FreqDiffTree::Node* child_for_descendant(FreqDiffTree::Node* v, FreqDiffTree::Node* w);
};

#endif /* RMQ_hpp */
