#include <iostream>
#include <gtest/gtest.h>
#include "hw1.h"
int main(int argc, char **argv)
{
    if (false) // make false to run unit-tests
    {
        Matrix mat11{
            {1, 11, 1, 1, 1, 1},
            {0, 0, 1, 0, 1, 0},
            {0, 0, 0, 0, 10, 1},
            {0, 1, 1, 0, 0, 1},
            {0, 0, 0, 1, 10, 1},
            {0, 0, 0, 0, 0, 1}};     
        Matrix mat12{
            {1, 11, 1, 1},
            {0, 0, 1, 0},
            {0, 1, 0, 0},
            {1, 0, 0, 0}};    
        Matrix mat13{
            {1, 12, 55},
            {78, 5, 92},
            {2, 3, 7}};
        algebra::show(mat11);
        algebra::show(algebra::upper_triangular(mat11));
    }
    else
    {
        ::testing::InitGoogleTest(&argc, argv);
        std::cout << "RUNNING TESTS ..." << std::endl;
        int ret{RUN_ALL_TESTS()};
        if (!ret)
            std::cout << "<<<SUCCESS>>>" << std::endl;
        else
            std::cout << "FAILED" << std::endl;
    }
    return 0;
}