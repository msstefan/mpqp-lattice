/*
 * Graph.h
 *
 *      Authors: Stefan S. Mihai, Florin Stoican
 *      Department of Automatic Control and Systems Engineering
 *          at "Politehnica" University of Bucharest
 */

#pragma once

#include <map>
#include <vector>
#include <string>

class Graph {

private:
    std::map<std::string, std::vector<std::string>> m_adjacency_list;

public:
    /**
     * @brief Construct a new Graph object
     * 
     */
    Graph() = default;

    /**
     * @brief Add a new edge in the graph
     * 
     * @param source Source node
     * @param target Target node
     */
    void addEdge(const std::string & source, const std::string & target);

    /**
     * @brief Print the elements of the graph
     * 
     */
    void printGraph() const;

    /**
     * @brief Crate a .dot file for Graphviz
     * 
     */
    void toDotFile(const std::string & filename) const;

    /**
     * @brief Check if an element exists in the graph.
     * 
     * @param element Element to be searched.
     * @return true if it exists
     * @return false if the element cannot be found
     */
    bool exists(const std::string& element);
};