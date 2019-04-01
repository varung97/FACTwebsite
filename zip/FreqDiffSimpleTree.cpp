#include <vector>
#include <cassert>

#include "FreqDiffSimpleTree.h"

FreqDiffSimpleTree::FreqDiffSimpleTree(size_t nodes_num_hint) {
	if (nodes_num_hint > 0) {
		nodes.reserve(nodes_num_hint);
	}
}

FreqDiffSimpleTree::~FreqDiffSimpleTree() {
	for (SimpleNode* node : nodes) {
		delete node;
	}
}

size_t FreqDiffSimpleTree::get_nodes_num() {
	return nodes.size();
}

FreqDiffSimpleTree::SimpleNode* FreqDiffSimpleTree::get_node(int i) {
	return nodes[i];
}

FreqDiffSimpleTree::SimpleNode* FreqDiffSimpleTree::add_node(int taxa, int label) {
	FreqDiffSimpleTree::SimpleNode* newnode = new FreqDiffSimpleTree::SimpleNode(get_nodes_num(), taxa, label);
	nodes.push_back(newnode);
	return newnode;
}

FreqDiffSimpleTree::SimpleNode::SimpleNode(int id, int taxa, int label) : id(id), taxa(taxa), label(label) {}

void FreqDiffSimpleTree::SimpleNode::add_child(SimpleNode* child) {
	children.push_back(child);
}

bool FreqDiffSimpleTree::SimpleNode::is_leaf() {
	return taxa != NONE;
}
