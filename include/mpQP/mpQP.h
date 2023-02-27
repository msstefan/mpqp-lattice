/*
 * mpQP.h
 *
 *      Authors: Stefan S. Mihai, Florin Stoican
 *      Department of Automatic Control and Systems Engineering
 *          at "Politehnica" University of Bucharest
 */

#pragma once
#include "Polyhedron.h"
#include "Graph.h"
#include <set>

#ifndef EIGEN_CONFIG
#define EIGEN_CONFIG
// We set it 1e-8 to avoid some issues with f-vector computation
#define LU_EPS 1e-8
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;
typedef Eigen::Vector<bool, Eigen::Dynamic> VectorXb;
typedef Eigen::IndexedView<MatrixXb, Eigen::internal::AllRange<-1>, Eigen::VectorXi> BlockMatrixXb;
#endif

struct mpQPSolution {
    std::string active = "";
    Eigen::MatrixXd CR_AHRep;
    Eigen::VectorXd CR_bHRep;
};

class mpQP {

private:
    int m_nineq, m_predN, m_dx, m_du, m_colsF;

    MatrixXb m_incidenceMx;
    int m_nvert, m_dim;

    Eigen::MatrixXd m_proA;
    Eigen::VectorXd m_prob;
    Eigen::MatrixXd m_proF;
    Eigen::MatrixXd m_proQ;
    Eigen::MatrixXd m_proHt;
    Eigen::MatrixXd m_Qi; 
    Eigen::MatrixXd m_ldom_AHRep;
    Eigen::MatrixXd m_ldom_VRep;
    Eigen::Polyhedron m_cddWrapper;
    Eigen::MatrixXd m_polar_VRep;
    Eigen::MatrixXd m_polar_AHRep;
    Eigen::VectorXd m_polar_bHRep;

    template <typename T>
    Eigen::MatrixXd null_isempty(T & mx, Eigen::MatrixXd & G,  Eigen::FullPivLU<Eigen::MatrixXd> & lu);

    bool zeroIsContained(Eigen::MatrixXd Z, Eigen::VectorXd & b);
    Eigen::VectorXi closureMap(const Eigen::VectorXi &S);
    std::vector<Eigen::VectorXi> minimalSet(const Eigen::VectorXi &H);
    int contains(const std::vector<Eigen::VectorXi> &v, const Eigen::VectorXi &item);

public:
    /**
     * @brief Construct a new mpQP object
     * 
     */
    mpQP() = default;

    /**
     * @brief Reads the input from a binary file.
     * 
     * @param filename the input file
     */
    void readBinaryInput(const char * filename);

    /**
     * @brief Solves the mpQP problem
     * 
     */
    void solve_mpQP();

    /**
     * @brief Compute a critical region
     * 
     * @param CR_AHRep A matrix of the critical region's H-representation
     * @param CR_bHRep b vector of the critical region's H-representation
     * @param active Set of active indices
     * @return 1 if the region can be created and is nonempty
     */
    int computeCR(Eigen::MatrixXd & CR_AHRep, Eigen::VectorXd & CR_bHRep, Eigen::VectorXi & active);

    /**
     * @brief A set o critical regions.
     * 
     */
    std::set<struct mpQPSolution> criticalRegions;  

    /**
     * @brief The solution poset.
     * 
     */
    Graph postree;
};