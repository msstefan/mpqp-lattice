/*
 * mpQP.cpp
 *
 *      Authors: Stefan S. Mihai, Florin Stoican
 *      Department of Automatic Control and Systems Engineering
 *          at "Politehnica" University of Bucharest
 */

#include "mpQP.h"
#include <iostream>
#include <fstream>
#include <igl/cat.h>
#include <igl/find.h>
#include <igl/setdiff.h>
#include <execution>
#include <deque>
#include <sstream>
#include <iterator>
#include <string>

#include <Eigen/QR>

#define SHOW_CR_TIME
// #define USE_TRUFFET
#define SHOW_HASSE_TIME
// #define TEST_AQAmx_RANK

#define VectorAppend(X, elem)               \
    (X).conservativeResize((X).rows() + 1); \
    X((X).rows() - 1) = (elem);

// #define TEST_KKTmx_RANK

bool operator<(const struct mpQPSolution & lhs, const struct mpQPSolution & rhs)
{
    return (lhs.active < rhs.active);
}

void mpQP::readBinaryInput(const char * filename)
{
    FILE* fp;
    if ((fp = fopen(filename, "rb")) == NULL) {
        std::cout << "Cannot open the binary file! Execution stopped.\n";
        return;
    }

    double *A, *b, *F, *Ht, *Q;
    size_t bytesRead;

    bytesRead = fread(&m_nineq, sizeof(int), 1, fp);
    bytesRead = fread(&m_predN, sizeof(int), 1, fp);
    bytesRead = fread(&m_dx, sizeof(int), 1, fp);
    bytesRead = fread(&m_du, sizeof(int), 1, fp);
    bytesRead = fread(&m_colsF, sizeof(int), 1, fp);

    A = new double[m_nineq * m_du * m_predN];
    b = new double[m_nineq];
    F = new double [m_nineq * m_colsF];
    Ht = new double[m_du * m_predN * m_colsF];
    Q = new double[m_du * m_predN * m_du * m_predN];
    
    bytesRead = fread(A, m_nineq * m_du * m_predN * sizeof(double), 1, fp);
    bytesRead = fread(b, m_nineq * sizeof(double), 1, fp);
    bytesRead = fread(F, m_nineq * m_colsF * sizeof(double), 1, fp);
    bytesRead = fread(Ht, m_du * m_predN * m_colsF * sizeof(double), 1, fp);
    bytesRead = fread(Q, m_du * m_predN * m_du * m_predN * sizeof(double), 1, fp);

    m_proA  = Eigen::Map<Eigen::Matrix<double,-1,-1,Eigen::RowMajor>> (A, m_nineq, m_du * m_predN);
    m_prob  = Eigen::Map<Eigen::VectorXd> (b, m_nineq);
    m_proF  = Eigen::Map<Eigen::Matrix<double,-1,-1,Eigen::RowMajor>> (F, m_nineq, m_colsF);
    m_proHt = Eigen::Map<Eigen::Matrix<double,-1,-1,Eigen::RowMajor>> (Ht, m_du * m_predN, m_colsF);
    m_proQ  = Eigen::Map<Eigen::Matrix<double,-1,-1,Eigen::RowMajor>> (Q, m_du * m_predN, m_du * m_predN);

    // std::cout << m_proA << std::endl;
    // std::cout << m_prob << std::endl;
    // std::cout << m_proF << std::endl;
    // std::cout << m_proHt << std::endl;
    // std::cout << m_proQ << std::endl;

    // Compute the inverse of m_proQ
    // Or: m_proQ.llt().solve(I);
    m_Qi = m_proQ.inverse();

    // std::cout << m_Qi << std::endl;

    delete[] A;
    delete[] b;
    delete[] F;
    delete[] Ht;
    delete[] Q;

    fclose(fp);

    std::cout << "[INFO] multi-parametric QP description\n\tno. inequalities: " << m_nineq 
              << "\n\tno. states: " << m_dx
              << "\n\tno. inputs: " << m_du
              << "\n\tprediction horizon: " << m_predN
              <<"\n";
}

void mpQP::solve_mpQP()
{
    std::deque<Eigen::VectorXi> Q;
    Graph lattice;
    // std::vector<Eigen::VectorXi> L;
    std::vector<Eigen::VectorXi> G;
    Eigen::VectorXi closure;
    // Critical region matrices
    Eigen::MatrixXd CR_AHRep; 
    Eigen::VectorXd CR_bHRep;
    int lenG, i, k, sz, found, check;

    m_ldom_AHRep.conservativeResize(m_proA.rows(), m_proA.cols() + m_proF.cols());
    m_ldom_AHRep << m_proA, (-1)*m_proF;

    int ldom_cols = m_ldom_AHRep.cols();
    int ldom_rows = m_ldom_AHRep.rows();
    int counter = 0;
    
    // Polar polytope vertices computation
    m_polar_VRep.conservativeResize(ldom_rows, ldom_cols); // Transposed
    /*  
        Tried using pragma omp
            omp_set_num_threads(ldom_cols);
            #pragma omp parallel for
        it's not suitable for small for loops
    */
    for (int i = 0; i < ldom_cols; i++) {
        m_polar_VRep(Eigen::all, i) = m_ldom_AHRep.col(i).array() / m_prob.array();
    }

    std::cout << "[INFO] Starting polar poly computation..\n";
    auto t_start_vertex = std::chrono::high_resolution_clock::now();
    m_cddWrapper.setVrep(m_polar_VRep, Eigen::VectorXd::Ones(ldom_rows));
    m_polar_AHRep = m_cddWrapper.hrep().first;
    m_polar_bHRep = m_cddWrapper.hrep().second;
    auto t_stop_vertex = std::chrono::high_resolution_clock::now();

    std::cout << "[INFO] Polar polytope size descripton:\n";
    std::cout << "\tvertex matrix size: " << m_polar_VRep.rows() << " x " << m_polar_VRep.cols() << "\n";
    std::cout << "\tinequalities: " << m_polar_AHRep.rows() << "\n";
    
    int polar_AHRep_cols = m_polar_AHRep.cols();

    std::ofstream solutionFile("active.txt", std::ios::out);

    auto t_start = std::chrono::high_resolution_clock::now();

    m_nvert = m_polar_VRep.rows();
    sz = m_nvert * m_nvert;

    // L.reserve(sz);
    Q.emplace_back(Eigen::VectorXi());
    
    /*
        Compute the incidence matrix
    */ 
    // Compute abs(A*V' - b) column-wise
    Eigen::MatrixXd AV = m_polar_AHRep * m_polar_VRep.transpose();
    Eigen::MatrixXd sub = ((-AV).colwise() + m_polar_bHRep).cwiseAbs();
    // Element-wise comparison to dd_almostzero
    m_incidenceMx.conservativeResize(m_nineq, m_nvert);
    m_incidenceMx = (sub.array() <= dd_almostzero);

    // Poset postree;
    std::cout << "[INFO] Starting solution search..\n";

    // Check the particular solution
    if (computeCR(CR_AHRep, CR_bHRep, closure)) {
        solutionFile << "\nActive indices of solution: \\emptyset \n";
        // std::cout << "\nActive indices of solution: \\emptyset \n";
        counter++;
    }

    // Cretaing the lattice
    while (Q.size())
    {
        std::ostringstream queueElemToString;
        queueElemToString << Q.front().transpose();

        G = this->minimalSet(Q.front());
        Q.pop_front();
        lenG = G.size();

        for (i = 0; i < lenG; i++) {
            closure = this->closureMap(G[i]);

            // ToDo: investigate why the entire polytope appears
            if (closure.rows() == m_nvert)
                continue;

            std::ostringstream clo;
            clo << closure.transpose(); 
            found = lattice.exists(clo.str());

            if (!found) {
                if (!computeCR(CR_AHRep, CR_bHRep, closure)) {
                    continue;
                }

                Q.emplace_front(closure);
                // L.emplace_back(closure);
                lattice.addEdge(queueElemToString.str(), clo.str());

                solutionFile << "Active indices of solution: " << closure.transpose() << "\n";
                // std::cout << "Active indices of solution: " << closure.transpose() << "\n";
                counter++;

                // struct mpQPSolution sol;
                // sol.active = ai_string.str();
                // sol.CR_AHRep = CR_AHRep;
                // sol.CR_bHRep = CR_bHRep;
                // criticalRegions.insert(sol);
            }
            // Entire graph
            // Add element to the poset
            // std::ostringstream ai_string;
            // ai_string << closure.transpose(); 
            // postree.addEdge(queueElemToString.str(), ai_string.str());
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();

    solutionFile.close();
#ifdef SHOW_HASSE_TIME
    std::cout << "\n\nSolution computation time: " << std::chrono::duration<double, std::milli>(t_end - t_start).count() / 1000.0 << " sec";
    std::cout << "\nFound " << counter << " solutions.";
    std::cout << "\ntotal computation time: " << std::chrono::duration<double, std::milli>(t_end - t_start).count()/ 1000.0 + std::chrono::duration<double, std::milli>(t_stop_vertex - t_start_vertex).count()/ 1000.0 << " ms" << std::endl;
#endif

    // Create dot file
    postree.toDotFile("postree.dot");

    // // Open a binary file for writing
    // std::ofstream ofile("postree.bin", std::ios::binary);
    // // Create a binary archive
    // boost::archive::binary_oarchive oa(ofile);
    // // Serialize the adjacency list
    // oa << postree;
}

int mpQP::computeCR(
    Eigen::MatrixXd & CR_AHRep, 
    Eigen::VectorXd & CR_bHRep, 
    Eigen::VectorXi & active)
{
    if (active.rows() == 0) {
        // Construct the critical regions
        CR_AHRep.conservativeResize(m_proF.rows(), m_proF.cols());
        CR_bHRep.conservativeResize(m_prob.rows(), 1);
        CR_AHRep = -0.5 * m_proA * m_Qi * m_proHt - m_proF;
        CR_bHRep = m_prob;
    }
    else {
        int n = m_proA.rows();
        int active_sz = active.rows();

        Eigen::MatrixXd KKTmx, KKTmx_row_u, KKTmx_row_l;
        Eigen::VectorXi inactive, dump;
        Eigen::VectorXi span = Eigen::VectorXi::LinSpaced(n, 0, n);

        igl::setdiff(span, active, inactive, dump);

        #ifdef TEST_KKTmx_RANK
            // Create the KKT matrix
            KKTmx_row_u.conservativeResize(m_proQ.rows(), m_proQ.cols() + active_sz);
            KKTmx_row_l.conservativeResize(active_sz, m_proA.cols() + active_sz);
            KKTmx.conservativeResize(KKTmx_row_u.rows() + KKTmx_row_l.rows(), KKTmx_row_u.cols());

            KKTmx_row_u << m_proQ, m_proA(active,Eigen::all).transpose();
            KKTmx_row_l << m_proA(active,Eigen::all), Eigen::MatrixXd::Zero(active_sz, active_sz);
            KKTmx << KKTmx_row_u,
                    KKTmx_row_l;

            // std::cout << "[DEBUG]\n" << KKTmx <<"\nActive: " << active.transpose() << "\n\n";

            /* 
                Test the rank of the KKT matrix.
                *******************************************************************************************
                Ideally, FullPivLU class uses the LAPACK library for performing the LU decomposition, 
                which is optimized for performance. So it is the fastest way to determine the rank of
                a matrix in Eigen.
                
                However,  if you don't need the full decomposition, you can use Eigen::ColPivHouseholderQR
                or Eigen::HouseholderQR class, which is faster than FullPivLU when you only need the rank.
                ToDo: see how to use Eigen::HouseholderQR to compute the rank.
                *******************************************************************************************
            */
            Eigen::FullPivLU<Eigen::MatrixXd> lu(KKTmx);
            lu.setThreshold(LU_EPS);

            if (lu.rank() != KKTmx.cols()) {
                return 0;
            }
        #else
            /*
                What if we compute the rank of AQA = m_proA(active,Eigen::all) * m_Qi * m_proA(active,Eigen::all).transpose() 
                or m_proA(active,Eigen::all).transpose() instead.. The KKTmx is KKTmx = [Q A(active,:)'; A(active,:) 0]. Would
                this matrix be rank-deficient only when A(active,:) is not full rank, since Q is the same.
            */
            #ifdef TEST_AQAmx_RANK
                Eigen::MatrixXd mxpro = m_proA(active,Eigen::all) * m_Qi * m_proA(active,Eigen::all).transpose(); // Should be multiplied with (-1), but it won't affect the rank.
            #else
                Eigen::MatrixXd mxpro = m_proA(active,Eigen::all).transpose();
            #endif

            Eigen::FullPivLU<Eigen::MatrixXd> lu(mxpro);
            lu.setThreshold(LU_EPS);
            if (lu.rank() != mxpro.cols()) {
                return 0;
            }

        #endif

        Eigen::MatrixXd Ai = (m_proA(active,Eigen::all) * m_Qi * m_proA(active,Eigen::all).transpose()).inverse();
        Eigen::MatrixXd SA = -1.0 * Ai * (2.0 * m_proF(active,Eigen::all) + m_proA(active,Eigen::all) * m_Qi * m_proHt);
        Eigen::MatrixXd sA = -1.0 * Ai * (2.0 * m_prob(active,Eigen::all)); // m_proc is zero
        Eigen::MatrixXd ATAi = m_proA(active,Eigen::all).transpose() * Ai;
        Eigen::MatrixXd LA = 0.5 * m_Qi * (ATAi * (2.0 * m_proF(active,Eigen::all) + m_proA(active,Eigen::all) * m_Qi * m_proHt) - m_proHt);
        Eigen::MatrixXd lA = 0.5 * m_Qi * (ATAi * (2.0 * m_prob(active,Eigen::all)));

        // Construct the critical regions
        CR_AHRep.conservativeResize(inactive.rows() + SA.rows(), SA.cols());
        CR_bHRep.conservativeResize(inactive.rows() + sA.rows(), 1);
        CR_AHRep << (m_proA(inactive, Eigen::all) * LA - m_proF(inactive,Eigen::all)),
                    -1.0 * SA;
        CR_bHRep << (m_prob(inactive) - m_proA(inactive, Eigen::all) * lA),
                     sA;
        
        #ifdef USE_TRUFFET
            /*
                © "Some Ideas to Test if a Polyhedron is Empty" by Laurent Truffet 
                    IMT-A Dpt. Automatique-Productique-Informatique Nantes, France 
                    email: laurent.truffet@imt-atlantique.fr 
                    April 28, 2020
            */

            int mm, nn;
            bool empty = 0;

            mm = CR_AHRep.rows();
            nn = CR_AHRep.cols();

            Eigen::MatrixXd A1 = CR_AHRep(Eigen::seq(0, mm - nn - 1), Eigen::all);
            Eigen::MatrixXd A2 = CR_AHRep(Eigen::seq(mm - nn, mm - 1), Eigen::all);

            lu.compute(A2);
            if (lu.rank() != nn) {
                #ifdef TRUEF_TRY_SHUFFLE
                Eigen::MatrixXd AA = CR_AHRep(Eigen::lastN(CR_AHRep.rows()).reverse(), Eigen::all);
                Eigen::MatrixXd bb = CR_bHRep(Eigen::lastN(CR_bHRep.rows()).reverse(), Eigen::all);

                CR_AHRep = AA;
                CR_bHRep = bb;
                A1 = CR_AHRep(Eigen::seq(0, mm - nn - 1), Eigen::all);
                A2 = CR_AHRep(Eigen::seq(mm - nn, mm - 1), Eigen::all);
                #else
                lu.compute(CR_AHRep);
                if (lu.rank() != nn)
                     return 0;
                #endif
            }

            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(mm - nn, mm - nn);

            Eigen::MatrixXd R(A1.rows(), A2.cols());
            R = A1 * A2.inverse();
  
            Eigen::VectorXd b1 = CR_bHRep.head(mm - nn);
            Eigen::VectorXd b2 = CR_bHRep.tail(nn); 
            Eigen::MatrixXd G(I.rows(), I.cols() + R.cols());
            G << I, (-1.0) * R;

            if (this->zeroIsContained(G, CR_bHRep) == false) {
                return 0;
            }

            Eigen::MatrixXd nullspace_b1 = this->null_isempty(b1, G, lu);
            Eigen::MatrixXd nullspace_Rb2 = this->null_isempty(R * b2, G, lu);
            Eigen::MatrixXd nullspace_R = this->null_isempty(R, G, lu);

            Eigen::MatrixXd nullspace;
            if (nullspace_b1.rows() > 0) {
                nullspace.conservativeResize(nullspace_b1.rows(), nullspace_b1.cols());
                nullspace << nullspace_b1;
            }
            if (nullspace_Rb2.rows() > 0) {
                nullspace.conservativeResize(nullspace_Rb2.rows(), nullspace_Rb2.cols());
                nullspace << nullspace_Rb2;
            }
            if (nullspace_R.rows() > 0) {
                nullspace.conservativeResize(nullspace_R.rows(), nullspace_R.cols());
                nullspace << nullspace_R;
            }

            if (nullspace.rows() > 0) {
                if (this->zeroIsContained(nullspace.transpose() * G, CR_bHRep) == false) {
                    return 0;
                }
            }

            for (int j = 0; j < nn; j++)
                for (int i = 0; i < mm - nn - 1; i++)
                    for (int ip = i + 1; ip < mm - nn; ip++) {
                        Eigen::VectorXd kp = -1 * R(ip, j) * I(Eigen::all, i) + R(i, j) * I(Eigen::all, ip);
                        if (this->zeroIsContained(kp.transpose() * G, CR_bHRep) == false) {
                            return 0;
                        }
                    }
        /*
                End of Trufet's algorithm.
        */
        #else
            /*
                Use cddlib to create a polyhedron from the H-representation.
                Test if the vertex matrix is nonempty.
            */
           m_cddWrapper.setHrep(CR_AHRep, CR_bHRep);
           if (m_cddWrapper.vrep().first.rows() < 1)
            return 0;
        #endif
    }

    return 1;
}

template <typename T>
Eigen::MatrixXd mpQP::null_isempty(T & mx, Eigen::MatrixXd & G, Eigen::FullPivLU<Eigen::MatrixXd> & lu)
{
    VectorXb index;
    Eigen::MatrixXd nullspace;
    Eigen::VectorXi vec_null_cols;
    Eigen::MatrixXd null_vec_plus;
    Eigen::MatrixXd null_transpose;

    Eigen::MatrixXd null_vec = lu.compute(mx.transpose()).kernel();
    if (null_vec.rows() > 0) {
        null_vec_plus.conservativeResize(null_vec.rows(), 2*null_vec.cols());
        null_vec_plus << null_vec, (-1.0) * null_vec;
        null_transpose = (null_vec_plus.transpose() * G).transpose();
        index = (null_transpose.array() >= 0).rowwise().all();
        igl::find(index, vec_null_cols);
        if (vec_null_cols.rows() > 0) {
            nullspace.conservativeResize(null_vec_plus.rows(), vec_null_cols.rows());
            nullspace = null_vec_plus(Eigen::all, vec_null_cols);
        }
    }

    return nullspace;
}

bool mpQP::zeroIsContained(Eigen::MatrixXd Z, Eigen::VectorXd & b)
{
    if ((b.array() >= 0).all()) {
        return 1;
    }

    int len = Z.rows();
    for (int i = 0; i < len; i++) {

        if ((Z.row(i).array() == 0).all()) {
            continue;
        }

        if ( ( ( !(
                    ((Z.row(i).array() >= 0).all()) && ((Z.row(i) * b) <= 0)
                )
            ) && ( !(
                    ((Z.row(i).array() <= 0).all()) && ((Z.row(i) * b) >= 0)
                )
            ) ) == false ) {
            return false;
        }

    }

    return true;
}

// ToDo: use igl::slice instead of igl::find
Eigen::VectorXi mpQP::closureMap(const Eigen::VectorXi &S)
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

std::vector<Eigen::VectorXi> mpQP::minimalSet(const Eigen::VectorXi &H)
{
    std::vector<Eigen::VectorXi> G_cal;
    Eigen::VectorXi index, V_H(m_nvert), V_H_label, h, h_v, Hv, dump;
    VectorXb G(m_nvert * m_nvert);
    int v, i, j, iter, len, G_col, h_len, h_v_len, ind;
    bool flag;

    std::fill(std::execution::seq, G.begin(), G.end(), 0);
    index = Eigen::VectorXi::LinSpaced(m_nvert, 0, m_nvert - 1);
    // Compute V \ H
    igl::setdiff(index, H, V_H, dump);
    
    // Create V_H_label to indentify minimal sets
    len = V_H.rows();
    V_H_label.resize(len);
    std::fill(std::execution::seq, V_H_label.begin(), V_H_label.end(), -1);

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

        std::sort(std::execution::par, Hv.begin(), Hv.end());
        
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

int mpQP::contains(const std::vector<Eigen::VectorXi> &v, const Eigen::VectorXi &item)
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