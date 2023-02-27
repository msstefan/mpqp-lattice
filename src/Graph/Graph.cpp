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