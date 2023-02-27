/*
 * Lattice.cpp
 *
 *      Authors: Stefan S. Mihai, Florin Stoican
 *      Department of Automatic Control and Systems Engineering
 *          at "Politehnica" University of Bucharest
 */

#include "Lattice.h"
#include <igl/all.h>
#include <igl/find.h>
#include <igl/setdiff.h>
#include <iostream>
#include <execution>
#include <deque>
#include <Eigen/LU>

#include <sstream>
#include <iterator>
#include <string>

#include <omp.h>

#define DEBUG_LATTICE
#define SHOW_HASSE_TIME

#define VectorAppend(X, elem)               \
    (X).conservativeResize((X).rows() + 1); \
    X((X).rows() - 1) = (elem);

Lattice::Lattice(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::MatrixXd &V)
{
    m_AHrep = A;
    m_bHrep = b;
    m_V = V;
    m_nineq = A.rows();
    m_nvert = V.rows(); // TBD
    m_dim = m_AHrep.cols();
    m_incidenceMx.resize(m_nineq, m_nvert);
}

void Lattice::computeIncidenceMatrix()
{
    // Compute abs(A*V' - b) column-wise
    Eigen::MatrixXd AV = m_AHrep * m_V.transpose();
    Eigen::MatrixXd sub = ((-AV).colwise() + m_bHrep).cwiseAbs();
    // Element-wise comparison to dd_almostzero
    m_incidenceMx = (sub.array() <= dd_almostzero);
}

// ToDo: use igl::slice instead of igl::find
Eigen::VectorXi Lattice::closureMap(const Eigen::VectorXi &S)
{
    VectorXb F, V;
    Eigen::VectorXi FIndices, VIndices;
    BlockMatrixXb blockMx = m_incidenceMx(Eigen::all, S);

    F = blockMx.rowwise().all();

    if (!(F.rows()))
        return VIndices;

    igl::find(F, FIndices);
    V = m_incidenceMx(FIndices, Eigen::all).colwise().all();
    igl::find(V, VIndices);

    return VIndices;
}

std::vector<Eigen::VectorXi> Lattice::minimalSet(const Eigen::VectorXi &H)
{
    std::vector<Eigen::VectorXi> G_cal;
    Eigen::VectorXi index, V_H(m_nvert), V_H_label, h, h_v, Hv, dump;
    VectorXb G(m_nvert * m_nvert);
    int v, i, j, iter, len, G_col, h_len, h_v_len, ind;
    bool flag;

    std::fill(std::execution::par_unseq, G.begin(), G.end(), 0);
    index = Eigen::VectorXi::LinSpaced(m_nvert, 0, m_nvert - 1);
    // Compute V \ H
    igl::setdiff(index, H, V_H, dump);
    
    // Create V_H_label to indentify minimal sets
    len = V_H.rows();
    V_H_label.resize(len);
    std::fill(std::execution::par_unseq, V_H_label.begin(), V_H_label.end(), -1);

    iter = 0;
    G_col = 0;
    Hv.resize(H.rows());
    h_v.resize(len);
    for (i = 0; i < len; i++)
    {
        v = V_H(i);
        Hv = H;
        VectorAppend(Hv, v)

            h = closureMap(Hv);

        if (!h.rows())
        {
            V_H_label(i) = 1;
            continue;
        }

        std::sort(std::execution::par_unseq, Hv.begin(), Hv.end());
        
        // Compute h \ Hv
        igl::setdiff(h, Hv, h_v, dump);

        h_len = h.rows();
        h_v_len = h_v.rows();
        if (!h_v_len)
        {
            V_H_label(i) = 0;
            for (j = 0; j < h_len; j++)
                G(h(j) * m_nvert + G_col) = 1;
            G_col++;
        }
        else
        {
            flag = 0;
            for (j = 0; j < h_v_len; j++)
            {
                ind = 0;
                while (V_H(ind++) != h_v(j))
                    ;
                if (V_H_label(ind - 1) <= 0)
                {
                    V_H_label(i) = 1;
                    flag = 1;
                    break;
                }
            }
            if (!flag)
            {
                V_H_label(i) = 0;
                for (j = 0; j < h_len; j++)
                    G(h(j) * m_nvert + G_col) = 1;
                G_col++;
            }
            // h_v.clear();
        }
    }

    for (i = 0; i < G_col; i++)
    {
        Eigen::VectorXi col(m_nvert);
        v = 0;
        for (j = 0; j < m_nvert; j++)
            if (G(j * m_nvert + i))
            {
                col(v++) = j;
            }
        col.conservativeResize(v);
        G_cal.emplace_back(col);
    }

    return G_cal;
}

int Lattice::contains(const std::vector<Eigen::VectorXi> &v, const Eigen::VectorXi &item)
{
    int v_len = v.size();
    int i, j, item_len;
    int found = 0;

    for (i = 0; i < v_len; i++)
    {
        item_len = item.rows();
        if (v[i].size() != item_len)
            continue;

        for (j = 0; j < item_len; j++)
            if (v[i](j) != item(j))
                break;

        if (j == item_len)
            return 1;

        /*if (v[i] == item)
            return 1;*/
    }

    return 0;
}

void Lattice::computeHasseDiagram()
{
    auto t_start = std::chrono::high_resolution_clock::now();
    std::deque<Eigen::VectorXi> Q;
    std::vector<Eigen::VectorXi> L;
    std::vector<Eigen::VectorXi> G;
    Eigen::VectorXi closure;
    int lenG, i, k, sz, found;

    sz = m_nvert * m_nvert;

    m_graphSource.reserve(sz);
    m_graphTarget.reserve(sz);
    L.reserve(sz);
    Q.emplace_back(Eigen::VectorXi());
    computeIncidenceMatrix();

    while (Q.size())
    {
        std::ostringstream queueElemToString;
        // std::copy(Q.front().begin(), Q.front().end(), std::ostream_iterator<int>(queueElemToString, " "));
        queueElemToString << Q.front().transpose();
        G = minimalSet(Q.front());
        Q.pop_front();
        lenG = G.size();

        for (i = 0; i < lenG; i++)
        {
            closure = closureMap(G[i]);

            // ToDo: investigate why the entire polytope appears
            if (closure.rows() == m_nvert)
                continue;

            found = contains(L, closure);

            if (!found)
            {
                Q.emplace_front(closure);
                L.emplace_back(closure);

                // Only the tree
                std::ostringstream closureToString;
                // std::copy(closure.begin(), closure.end(), std::ostream_iterator<int>(closureToString, " "));
                closureToString << closure.transpose();
                m_graphSource.emplace_back(queueElemToString.str());
                m_graphTarget.emplace_back(closureToString.str());
            }
            // Entire graph
            // std::ostringstream closureToString;
            // std::copy(closure.begin(), closure.end(), std::ostream_iterator<int>(closureToString, " "));
            // m_graphSource.emplace_back(queueElemToString.str());
            // m_graphTarget.emplace_back(closureToString.str());
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
#ifdef SHOW_HASSE_TIME
    std::cout << "\n\nHasse-diagram computation time: " << std::chrono::duration<double, std::milli>(t_end - t_start).count() << " ms" << std::endl;
#endif
    m_faceLattice = L;
}

void Lattice::compute_kSkeleton(unsigned int dim)
{
    auto t_start = std::chrono::high_resolution_clock::now();
    std::deque<Eigen::VectorXi> Q;
    std::vector<Eigen::VectorXi> L;
    std::vector<Eigen::VectorXi> G;
    Eigen::VectorXi closure;
    Eigen::VectorXi VIndices;
    VectorXb V;
    
    int lenG, i, k, sz, found;

    dim--;
    sz = m_nvert * m_nvert;

    m_graphSource.reserve(sz);
    m_graphTarget.reserve(sz);
    L.reserve(sz);
    Q.emplace_back(Eigen::VectorXi());
    computeIncidenceMatrix();

    while (Q.size())
    {
        // Source node
        std::ostringstream queueElemToString;
        queueElemToString << Q.front().transpose();

        G = minimalSet(Q.front());
        Q.pop_front();
        lenG = G.size();

        for (i = 0; i < lenG; i++)
        {
            closure = closureMap(G[i]);

            // Test rank condition for the k-skeleton
            V = m_incidenceMx(Eigen::all, closure).rowwise().all();
            igl::find(V, VIndices);
            Eigen::FullPivLU<Eigen::MatrixXd> LU(m_AHrep(VIndices, Eigen::all));
            LU.setThreshold(LU_EPS);

            if (LU.kernel().cols() > dim) {
                continue;
            }

            // ToDo: investigate why the entire polytope appears
            if (closure.rows() == m_nvert)
                continue;

            found = contains(L, closure);

            if (!found)
            {
                Q.emplace_front(closure);
                L.emplace_back(closure);

                // // Only the tree
                // std::ostringstream closureToString;
                // closureToString << closure.transpose(); // Target node
                // latticeGraph.addEdge(queueElemToString.str(), closureToString.str());
            }
            // Only the tree
            std::ostringstream closureToString;
            closureToString << closure.transpose(); // Target node
            latticeGraph.addEdge(queueElemToString.str(), closureToString.str());
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    p_executionTime = std::chrono::duration<double, std::milli>(t_end - t_start);
#ifdef SHOW_HASSE_TIME
    std::cout << "\n[DEBUG] k-Skeleton computation time: " << std::chrono::duration<double, std::milli>(t_end - t_start).count() << " ms" << std::endl;
#endif

    m_faceLattice = L;
}

std::vector<Eigen::VectorXi> Lattice::hasseDiagram()
{
    if (!m_faceLattice.size())
    {
        computeHasseDiagram();

        // std::cout << "\nSource\n";
        // for (int i = 0; i < m_graphSource.size(); i++)
        //     std::cout << "\"" << m_graphSource[i] << "\", ";

        // std::cout << "\nTarget\n";
        // for (int i = 0; i < m_graphTarget.size(); i++)
        //     std::cout << "\"" << m_graphTarget[i] << "\", ";
        // std::cout << "\n";
    }

    return m_faceLattice;
}

std::vector<Eigen::VectorXi> Lattice::kSkeleton(unsigned int dimension)
{
    if (!m_faceLattice.size()) {
        compute_kSkeleton(dimension);

        // std::cout << "\nSource\n";
        // for (int i = 0; i < m_graphSource.size(); i++)
        //     std::cout << "\"" << m_graphSource[i] << "\", ";

        // std::cout << "\nTarget\n";
        // for (int i = 0; i < m_graphTarget.size(); i++)
        //     std::cout << "\"" << m_graphTarget[i] << "\", ";
        // std::cout << "\n";
    }

    return m_faceLattice;
}

void Lattice::computeFVector()
{
    Eigen::VectorXi VIndices;
    VectorXb V;
    int len = m_faceLattice.size();
    
    if (!len) {
        computeHasseDiagram();
        len = m_faceLattice.size();
    }

    m_fVector.resize(m_dim);
    std::fill(std::execution::par_unseq, m_fVector.begin(), m_fVector.end(), 0);
    for (int i = 0; i < len; i++) {
        V = m_incidenceMx(Eigen::all, m_faceLattice[i]).rowwise().all();
        igl::find(V, VIndices);
        m_activeIndices.emplace_back(VIndices);

        // f-dimension rank condition
        Eigen::FullPivLU<Eigen::MatrixXd> LU(m_AHrep(VIndices, Eigen::all));
        LU.setThreshold(LU_EPS);
        m_fVector(m_dim - LU.rank())++;
    }
}

Eigen::VectorXi Lattice::fVector()
{
    if (!m_fVector.rows())
    {
        computeFVector();
    }
    return m_fVector;
}

std::vector<Eigen::VectorXi> Lattice::getActiveIndMx()
{
    return m_activeIndices;
}

void Lattice::vrep2hrep(const Eigen::MatrixXd & V, Eigen::MatrixXd & A, Eigen::VectorXd & b)
{
    // // Find the mean of all vertices
    // Eigen::VectorXd mean = V.colwise().mean();

    // // Subtract the mean from all vertices
    // Eigen::MatrixXd V_centered = V.rowwise() - mean.transpose();

    // // Find the covariance matrix
    // Eigen::MatrixXd cov = (V_Center.adjoint() * V_Center) / double(V.rows() - 1);

    // // Find the eigenvector corresponding to the smallest eigenvalue
    // Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver(cov);
    // Eigen::VectorXd normal = eigenSolver.eigenvectors().col(0);

    // // Normalize the normal vector
    // normal.normalize();

    // // Define the hyperplane coefficients
    // VectorXd A(normal.rows());
    // double b = normal.dot(mean);
    // A = normal;
}