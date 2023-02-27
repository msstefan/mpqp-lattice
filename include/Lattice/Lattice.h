/*
 * Lattice.h
 *
 *      Authors: Stefan S. Mihai, Florin Stoican
 *      Department of Automatic Control and Systems Engineering
 *          at "Politehnica" University of Bucharest
 */

#include "Polyhedron.h"
#include "Graph.h"

#define LATTICE_DEBUG

#ifndef EIGEN_CONFIG
#define EIGEN_CONFIG
// We set it 1e-8 to avoid some issues with f-vector computation
#define LU_EPS 1e-8
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;
typedef Eigen::Vector<bool, Eigen::Dynamic> VectorXb;
typedef Eigen::IndexedView<MatrixXb, Eigen::internal::AllRange<-1>, Eigen::VectorXi> BlockMatrixXb;
#endif

class Lattice
{
private:
    std::vector<Eigen::VectorXi> m_faceLattice;
    std::vector<Eigen::VectorXi> m_reducedLattice;
    std::vector<Eigen::VectorXi> m_activeIndices;
    std::vector<Eigen::VectorXi> m_reduced_activeIndices;
    std::vector<std::string> m_graphSource;
    std::vector<std::string> m_graphTarget;
    Eigen::VectorXi m_fVector;
    Eigen::VectorXi m_reduced_fVector;
    Eigen::MatrixXd m_AHrep;
    Eigen::VectorXd m_bHrep;
    Eigen::MatrixXd m_V;
    MatrixXb m_incidenceMx;
    int m_nineq, m_dim;
    int m_nvert;

    void computeIncidenceMatrix();
    Eigen::VectorXi closureMap(const Eigen::VectorXi &S);
    std::vector<Eigen::VectorXi> minimalSet(const Eigen::VectorXi &H);
    int contains(const std::vector<Eigen::VectorXi> &v, const Eigen::VectorXi &item);
    void computeHasseDiagram();
    void compute_kSkeleton(unsigned int dim);
    void computeFVector();
    void computeReducedFVector();

public:

    /**
     * @brief Construct a new Lattice object
     * 
     */
    Lattice() = default;

    /**
     * @brief Construct a new Lattice object from the H-representation Ax <= b and the vertex matrix
     * 
     * @param A A matrix of the H-representation
     * @param b b matrix of the H-representation
     * @param V vertex matrix
     */
    Lattice(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::MatrixXd &V);

    /**
     * @brief Computes the entire Hasse diagram
     * 
     * @return std::vector<Eigen::VectorXi> stack with faces of P
     */
    std::vector<Eigen::VectorXi> hasseDiagram();

    /**
     * @brief Computes the f-vector
     * 
     * @return Eigen::VectorXi the f-vector
     */
    Eigen::VectorXi fVector();

    // Eigen::VectorXi reducedFVector();

    /**
     * @brief 
     * 
     * @param dimension 
     * @return std::vector<Eigen::VectorXi> 
     */
    std::vector<Eigen::VectorXi> kSkeleton(unsigned int dimension);

    /**
     * @brief Get the Active Ind Mx object
     * 
     * @return std::vector<Eigen::VectorXi> 
     */
    std::vector<Eigen::VectorXi> getActiveIndMx();

    /**
     * @brief Execution time.
     * 
     */
    std::chrono::duration<double, std::milli> p_executionTime;
    
    /**
     * @brief Face lattice graph
     * 
     */
    Graph latticeGraph;
};