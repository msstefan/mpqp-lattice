/*
 * main.cpp
 *
 *      Authors: Stefan S. Mihai, Florin Stoican
 *      Department of Automatic Control and Systems Engineering
 *          at "Politehnica" University of Bucharest
 */

#include "main.h"

int main()
{
    mpQP problem;

    problem.readBinaryInput("../mpQP_test_1.bin");
    problem.solve_mpQP();
  
    return 0;
}