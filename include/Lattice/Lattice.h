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

/**
 * @brief
 *
 */
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

    /**
     * @brief
     *
     */
    void computeIncidenceMatrix();

    /**
     * @brief
     *
     * @param S
     * @return Eigen::VectorXi
     */
    Eigen::VectorXi closureMap(const Eigen::VectorXi &S);

    /**
     * @brief
     *
     * @param H
     * @return std::vector<Eigen::VectorXi>
     */
    std::vector<Eigen::VectorXi> minimalSet(const Eigen::VectorXi &H);

    /**
     * @brief 
     * 
     * @param v 
     * @param item 
     * @return int 
     */
    int contains(const std::vector<Eigen::VectorXi> &v, const Eigen::VectorXi &item);

    /**
     * @brief 
     * 
     */
    void computeHasseDiagram();

    void compute_kSkeleton(unsigned int dim);

    /**
     * @brief 
     * 
     */
    void computeFVector();

    void computeReducedFVector();

public:
    /**
     * @brief Construct a new Lattice object
     * 
     */
    Lattice() = default;

    /**
     * @brief Construct a new Lattice object
     * 
     * @param A 
     * @param b 
     * @param V 
     */
    Lattice(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::MatrixXd &V);

    /**
     * @brief 
     * 
     * @return std::vector<Eigen::VectorXi> 
     */
    std::vector<Eigen::VectorXi> hasseDiagram();

    /**
     * @brief 
     * 
     * @return Eigen::VectorXi 
     */
    Eigen::VectorXi fVector();

    Eigen::VectorXi reducedFVector();

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
     * @brief 
     * 
     */
    std::chrono::duration<double, std::milli> p_executionTime;

    /**
     * @brief 
     * 
     * @param V 
     * @param A 
     * @param b 
     */
    void vrep2hrep(const Eigen::MatrixXd & V, Eigen::MatrixXd & A, Eigen::VectorXd & b);
    
    /**
     * @brief 
     * 
     */
    Graph latticeGraph;
};