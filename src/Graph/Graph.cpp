#include "Graph.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <execution>
int nextId = 0;

void Graph::addEdge(const std::string & source, const std::string & target)
{
  if (m_adjacency_list.find(source) == m_adjacency_list.end()) {
      m_adjacency_list[source] = std::vector<std::string>();
  }
  m_adjacency_list[source].push_back(target);
}

bool Graph::exists(const std::string& element)
{
  for (auto & iter: m_adjacency_list) {
    if (std::find(std::execution::par_unseq, iter.second.begin(), iter.second.end(), element) != iter.second.end())
      return 1;
  }
  return 0;
}

void Graph::printGraph() const
{
  for (const auto& entry : m_adjacency_list) {
    std::cout << entry.first << ": ";
    for (const auto& neighbor : entry.second) {
      std::cout << neighbor << "; ";
    }
    std::cout << std::endl;
  }
}

void Graph::toDotFile(const std::string & filename) const {
  std::ofstream out(filename);
  out << "digraph {\n";
  for (const auto & entry : m_adjacency_list) {
    const std::string & source = entry.first;
    for (const std::string & target : entry.second) {
      out << "  \"" << source << "\" -> \"" << target << "\"\n";
    }
  }
  out << "}\n";
  out.close();
}


/*
  Poset implementation
*/
Poset::Poset(): m_root(nullptr) 
{
}

PosetNode* Poset::insert(const std::string & subset)
{
  if (m_root == nullptr) {
    m_root = new PosetNode(subset);
    m_root->id = nextId++;
    return m_root;
  }

  PosetNode* current = m_root;
  while (current != nullptr) {
    if (subset_is_contained(current->subset, subset)) {
      if (subset_is_contained(subset, current->subset)) {
        return current;
      }
      if (current->left == nullptr) {
        current->left = new PosetNode(subset);
        current->left->id = nextId++;
        return current->left;
      }
      current = current->left;
    } 
    else {
      if (current->right == nullptr) {
        current->right = new PosetNode(subset);
        current->right->id = nextId++;
        return current->right;
      }
      current = current->right;
    }
  }

  return nullptr;
}

bool Poset::subset_is_contained(const std::string & lhs, const std::string & rhs)
{
  return (lhs == rhs);
}

void Poset::freeMemory(PosetNode* node)
{
  if (node == nullptr) 
    return;
  
  freeMemory(node->left);
  freeMemory(node->right);
  
  delete node;
}

void Poset::destroy()
{
  freeMemory(m_root);
  m_root = nullptr;
}

void Poset::toDotFile(PosetNode* node, std::ofstream& out) {
  if (node == nullptr) 
    return;
  out << node->id << " [label=\"" << node->subset << "\"];\n";
  if (node->left != nullptr) {
    out << node->id << " -> " << node->left->id << ";\n";
    toDotFile(node->left, out);
  }
  if (node->right != nullptr) {
    out << node->id << " -> " << node->right->id << ";\n";
    toDotFile(node->right, out);
  }
}

void Poset::writeDotFile(const std::string& filename) {
  std::ofstream out(filename);
  out << "digraph {\n";
  toDotFile(m_root, out);
  out << "}\n";
  out.close();
}