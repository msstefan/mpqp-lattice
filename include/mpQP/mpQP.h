/*
 * mpQP.h
 *
 *      Authors: Stefan S. Mihai, Florin Stoican
 *      Department of Automatic Control and Systems Engineering
 *          at "Politehnica" University of Bucharest
 */

#pragma once
// #include "Lattice.h"
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
    // Lattice methods
    MatrixXb m_incidenceMx;
    int m_nvert, m_dim;
    //
    Eigen::MatrixXd m_proA;
    Eigen::VectorXd m_prob;
    Eigen::MatrixXd m_proF;
    Eigen::MatrixXd m_proQ;
    Eigen::MatrixXd m_proHt;
    Eigen::MatrixXd m_Qi; 
    //
    Eigen::MatrixXd m_ldom_AHRep;
    Eigen::MatrixXd m_ldom_VRep;
    /* 
    * Eigen::MatrixXd m_ldom_bHRep; //it's actually m_prob
    */
    Eigen::Polyhedron m_cddWrapper;
    //
    Eigen::MatrixXd m_polar_VRep;
    Eigen::MatrixXd m_polar_AHRep;
    Eigen::VectorXd m_polar_bHRep;

    template <typename T>
    Eigen::MatrixXd null_isempty(T & mx, Eigen::MatrixXd & G,  Eigen::FullPivLU<Eigen::MatrixXd> & lu);

    bool zeroIsContained(Eigen::MatrixXd Z, Eigen::VectorXd & b);

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

public:
    /**
     * @brief Construct a new mpQP object
     * 
     */
    mpQP() = default;

    /**
     * @brief 
     * 
     * @param filename 
     */
    void readBinaryInput(const char * filename);

    /**
     * @brief 
     * 
     */
    void solve_mpQP();

    /**
     * @brief 
     * 
     * @param CR_AHRep 
     * @param CR_bHRep 
     * @param active 
     * @return int 
     */
    int computeCR(Eigen::MatrixXd & CR_AHRep, Eigen::VectorXd & CR_bHRep, Eigen::VectorXi & active);

    /**
     * @brief 
     * 
     */
    std::set<struct mpQPSolution> criticalRegions;  

    /**
     * @brief 
     * 
     */
    Graph postree;
};