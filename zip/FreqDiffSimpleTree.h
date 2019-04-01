#ifndef FREQ_DIFF_SIMPLETREE_H_
#define FREQ_DIFF_SIMPLETREE_H_

#include <vector>

class FreqDiffSimpleTree {
public:
	FreqDiffSimpleTree(size_t nodes_num_hint = 0);
	virtual ~FreqDiffSimpleTree();

	class SimpleNode;

	SimpleNode* get_node(int i);
	size_t get_nodes_num();
	SimpleNode* add_node(int taxa = -1, int label = 1);

	SimpleNode* root;

private:
	std::vector<SimpleNode*> nodes;
};

class FreqDiffSimpleTree::SimpleNode {
public:
	static const int NONE = -1;

	std::vector<SimpleNode*> children;

	int id;
	int taxa;
    int label;

	SimpleNode(int id, int taxa, int label);

	void add_child(SimpleNode* child);
	bool is_leaf();
};

#endif
