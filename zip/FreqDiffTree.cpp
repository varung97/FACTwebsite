#include <vector>
#include <unordered_map>
#include <cassert>
#include <sstream>
#include <iostream>
#include <queue>
#include <stack>

#include "FreqDiffTree.h"

std::unordered_map<std::string,int> FreqDiffTree::taxa_ids;
std::unordered_map<int,std::string> FreqDiffTree::taxa_names;

FreqDiffTree::FreqDiffTree(size_t nodes_num_hint) : leaves_num(0) {
	if (nodes_num_hint > 0) {
		nodes.reserve(nodes_num_hint);
	}
}

FreqDiffTree::FreqDiffTree(std::string& newick_str) : leaves_num(0) {
	const char* str = newick_str.c_str();
	while (*str != '(') str++;
	this->build_tree(str);

	//taxa_to_leaf.resize(FreqDiffTree::get_taxas_num());
	for (Node* node : nodes) {
		if (node->is_leaf()) {
			// std::cout << "Leaf: " << node->id << std::endl;
			//taxa_to_leaf[node->taxa] = node;
			taxa_to_leaf_map[node->taxa] = node;
			leaves_num++;
		} else {
			// std::cout << "Internal: " << node->id << " " << node->children[0]->id << " " << node->children[1]->id << std::endl;
		}

		// std::cout << node->to_string() << std::endl;

		if (!node->is_root()) {
			node->depth = node->parent->depth + 1;
		}
	}

	// calc sizes
	for (int i = nodes.size()-1; i > 0 ; i--) {
		if (nodes[i]->is_leaf()) {
			nodes[i]->size = 1;
		}
		nodes[i]->parent->size += nodes[i]->size;
	}
}

FreqDiffTree::FreqDiffTree(FreqDiffTree* other) : leaves_num(other->leaves_num) {
	nodes.reserve(other->get_nodes_num());
	//taxa_to_leaf.resize(FreqDiffTree::get_taxas_num());
	for (Node* node : other->nodes) {
		FreqDiffTree::Node* newnode = new Node(node->id, node->taxa);
		newnode->weight = node->weight;
		newnode->size = node->size;
		newnode->depth = node->depth;
		nodes.push_back(newnode);

		if (!node->is_root()) {
			nodes[node->parent->id]->add_child(newnode);
		}
		if (newnode->is_leaf()) {
			//taxa_to_leaf[newnode->taxa] = newnode;
			taxa_to_leaf_map[newnode->taxa] = newnode;
		}
	}
}

FreqDiffTree::~FreqDiffTree() {
	for (Node* node : nodes) {
		delete node;
	}
}

std::string FreqDiffTree::to_string() {
	std::stringstream ss;
	for (Node* node : nodes) {
		ss << node->to_string() << std::endl;
	}
	return ss.str();
}

std::string FreqDiffTree::to_newick(Node* node) {
	if (node == NULL) {
		node = get_root();
	}
	return node->to_newick();
}

size_t FreqDiffTree::get_taxas_num() {
	return taxa_ids.size();
}

size_t FreqDiffTree::get_nodes_num() {
	return nodes.size();
}

FreqDiffTree::Node* FreqDiffTree::get_node(int i) {
	return nodes[i];
}

FreqDiffTree::Node* FreqDiffTree::get_root() {
	return get_node(0);
}

FreqDiffTree::Node* FreqDiffTree::get_leaf(int taxa) {
	//return taxa_to_leaf[taxa];
	return taxa_to_leaf_map[taxa];
}

size_t FreqDiffTree::get_leaves_num() {
	return leaves_num;
}

FreqDiffTree::Node* FreqDiffTree::add_node(int taxa, int label) {
	FreqDiffTree::Node* newnode = new FreqDiffTree::Node(get_nodes_num(), taxa, label);
	nodes.push_back(newnode);
	if (taxa >= 0) {
		// taxa_to_leaf.insert(taxa_to_leaf.begin()+taxa, newnode); TODO
		leaves_num++;
		taxa_to_leaf_map[taxa] = newnode;
	}
	return newnode;
}

void FreqDiffTree::delete_nodes(bool* to_delete) {
	for (size_t i = 1; i < get_nodes_num(); i++) {
		if (to_delete[i]) {
			for (FreqDiffTree::Node* child : nodes[i]->children) {
				nodes[i]->parent->add_child(child);
			}
			nodes[i]->parent->null_child(nodes[i]->pos_in_parent);
			delete nodes[i];
			nodes[i] = NULL;
		}
	}
	fix_tree();
}

void FreqDiffTree::fix_tree(FreqDiffTree::Node* root) {
	if (root == NULL)
		root = get_root();
	nodes.clear();
	fix_tree_supp(root);

	// recalc sizes
	nodes[0]->size = 0;
	nodes[0]->depth = 0;
	for (size_t i = 1; i < nodes.size(); i++) {
		nodes[i]->size = 0;
		nodes[i]->depth = nodes[i]->parent->depth + 1;
	}
	for (int i = nodes.size()-1; i > 0 ; i--) {
		if (nodes[i]->is_leaf()) {
			nodes[i]->size = 1;
		}
		nodes[i]->parent->size += nodes[i]->size;
	}
}

void FreqDiffTree::fix_tree_supp(FreqDiffTree::Node* curr) {
	curr->id = nodes.size();
	nodes.push_back(curr);
	curr->fix_children();
	for (FreqDiffTree::Node* child : curr->children) {
		fix_tree_supp(child);
	}
}

void FreqDiffTree::reorder() {
	for (FreqDiffTree::Node* node : nodes) {
		if (node->is_leaf()) continue;

		int heaviest = 0;
		for (size_t i = 1; i < node->get_children_num(); i++) {
			if (node->children[i]->size > node->children[heaviest]->size) {
				heaviest = i;
			}
		}
		FreqDiffTree::Node* heaviest_node = node->children[heaviest];
		node->set_child(node->children[0], heaviest);
		node->set_child(heaviest_node, 0);
	}
}

FreqDiffTree::Node* FreqDiffTree::build_tree(const char*& str) {
	assert(*str == '(');
	str++;

	int vecpos = nodes.size();
	nodes.push_back(new FreqDiffTree::Node(vecpos));

	assert(*str == '(' || isdigit(*str));
	while (true) {
		FreqDiffTree::Node* subtree;
		if (*str == '(') {
			subtree = build_tree(str);
		} else {
			std::string taxa;
			while (*str != ':' && *str != ',' && *str != ')') {
				taxa += *str;
				str++;
			}
			subtree = new Node(nodes.size(), get_taxa_id(taxa));
			nodes.push_back(subtree);
		}
		nodes[vecpos]->add_child(subtree);
		while (*str != ',' && *str != ')') str++;
		if (*str == ',') {
			str++;
		} else if (*str == ')') {
			str++;
			break;
		}
	}

	return nodes[vecpos];
}

int FreqDiffTree::get_taxa_id(std::string& taxa) {
	auto taxa_entry = FreqDiffTree::taxa_ids.find(taxa);
	if (taxa_entry == FreqDiffTree::taxa_ids.end()) {
		int taxa_id = taxa_ids.size();
		taxa_ids.insert({taxa,taxa_id});
		taxa_names.insert({taxa_id,taxa});
		return taxa_id;
	} else {
		return taxa_entry->second;
	}
}

FreqDiffTree::Node::Node(int id) : parent(NULL), pos_in_parent(NONE), id(id), secondary_id(id), taxa(NONE), label(NONE), weight(0), size(0), depth(0), counter(0) {}
FreqDiffTree::Node::Node(int id, int taxa) : parent(NULL), pos_in_parent(NONE), id(id), secondary_id(id), taxa(taxa), label(NONE), weight(0), size(0), depth(0), counter(0) {}
FreqDiffTree::Node::Node(int id, int taxa, int label) : parent(NULL), pos_in_parent(NONE), id(id), secondary_id(id), taxa(taxa), label(label), weight(0), size(0), depth(0), counter(0) {}

void FreqDiffTree::Node::fix_children() {
	size_t curr_pos = 0;
	for (size_t i = 0; i < children.size(); i++) {
		if (children[i] != NULL) {
			children[curr_pos] = children[i];
			children[curr_pos]->pos_in_parent = curr_pos;
			curr_pos++;
		}
	}
	children.resize(curr_pos);
}


void FreqDiffTree::Node::add_child(Node* child) {
	child->parent = this;
	child->pos_in_parent = children.size();
	children.push_back(child);
}
void FreqDiffTree::Node::set_child(Node* child, size_t pos) {
	child->parent = this;
	child->pos_in_parent = pos;
	children[pos] = child;
}

bool FreqDiffTree::Node::is_leaf() {
	return taxa != NONE;
}

size_t FreqDiffTree::Node::get_children_num() {
	return children.size();
}

bool FreqDiffTree::Node::is_root() {
	return parent == NULL;
}

void FreqDiffTree::Node::null_child(size_t pos) {
	children[pos] = NULL;
}

void FreqDiffTree::Node::clear_children() {
	children.clear();
}

std::string FreqDiffTree::Node::to_string() {
	std::stringstream ss;
	ss << id << " ";
	// std::cout << "B" << std::endl;
	if (taxa == NONE) {
		ss << "[";
		for (size_t i = 0; i < children.size(); i++) {
			// std::cout << i << std::endl;
			ss << (i > 0 ? "," : "") << children[i]->id;
		}
		ss << "]";
	} else {
		ss << taxa_names[taxa];
	}
//	if (!is_root()) {
	ss << ", weight: " << weight << ", sec id: " << secondary_id << ", label: " << label;
//	}
	return ss.str();
}

std::string FreqDiffTree::Node::to_newick() {
	if (is_leaf()) {
		return taxa_names[taxa];
	}
	std::stringstream ss;
	ss << "(";
	for (size_t i = 0; i < get_children_num(); i++) {
		if (i > 0) ss << ",";
		ss << children[i]->to_newick();
	}
	ss << ")";
	return ss.str();
}
