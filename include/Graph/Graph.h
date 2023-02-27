#pragma once

#include <map>
#include <vector>
#include <string>
// #include <boost/archive/binary_oarchive.hpp>
// #include <boost/serialization/export.hpp>

class Graph {

private:
    /**
     * @brief 
     * 
     */
    std::map<std::string, std::vector<std::string>> m_adjacency_list;

    // friend class boost::serialization::access;
    // template<typename Archive>
    // void serialize(Archive & ar, const unsigned int version)
    // {
    //     ar & (*this);
    // }
public:
    /**
     * @brief Construct a new Graph object
     * 
     */
    Graph() = default;

    /**
     * @brief 
     * 
     * @param source 
     * @param target 
     */
    void addEdge(const std::string & source, const std::string & target);

    /**
     * @brief 
     * 
     */
    void printGraph() const;

    /**
     * @brief 
     * 
     */
    void toDotFile(const std::string & filename) const;

    /**
     * @brief 
     * 
     * @param element 
     * @return true 
     * @return false 
     */
    bool exists(const std::string& element);
};

/**
 * ToDo: refactor this section. Memory deallocation needs to be implemented!
 * 
 */

struct PosetNode {
    int id;
    std::string subset;
    PosetNode* left;
    PosetNode* right;

    PosetNode(const std::string & subset) : subset(subset), left(nullptr), right(nullptr), id(0) {}
};

/**
 * @brief 
 * 
 */
class Poset {
private:
    bool subset_is_contained(const std::string & lhs, const std::string & rhs);
    PosetNode* m_root;
    void freeMemory(PosetNode* node);

public:
    /**
     * @brief Construct a new Poset object
     * 
     */
    Poset();

    /**
     * @brief 
     * 
     * @param subset 
     * @return Node* 
     */
    PosetNode* insert(const std::string & subset);

    /**
     * @brief 
     * 
     */
    void destroy();

    /**
     * @brief 
     * 
     * @param node 
     * @param out 
     */
    void toDotFile(PosetNode* node, std::ofstream& out);

    /**
     * @brief 
     * 
     * @param filename 
     */
    void writeDotFile(const std::string& filename);
};
