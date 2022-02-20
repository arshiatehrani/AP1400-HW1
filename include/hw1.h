#ifndef AP_HW1_H
#define AP_HW1_H

#include<vector>
#include<iostream>
#include<random>
#include<iomanip>
#include<cmath>

using Matrix = std::vector<std::vector<double>>;
namespace algebra
{
    Matrix zeros(size_t, size_t);
    Matrix ones(size_t, size_t);
    Matrix random(size_t, size_t, double, double);
    void show(const Matrix&);
    Matrix multiply(const Matrix&, double);
    Matrix multiply(const Matrix&, const Matrix&);
    Matrix sum(const Matrix&, double);
    Matrix sum(const Matrix&, const Matrix&);
    Matrix transpose(const Matrix&);
    Matrix minor(const Matrix&, size_t, size_t);
    Matrix delcol(const Matrix&, size_t);
    Matrix delrow(const Matrix&, size_t);
    double determinant(const Matrix&);
    Matrix cofactor(const Matrix&);
    Matrix Adjugate(const Matrix&);
    Matrix inverse(const Matrix&);
    Matrix concatenate(const Matrix&, const Matrix&, int);
    Matrix ero_swap(const Matrix&, size_t, size_t);
    Matrix ero_multiply(const Matrix&, size_t, double);
    Matrix ero_sum(const Matrix&, size_t, double, size_t);
    Matrix upper_triangular(const Matrix&);
}
#endif //AP_HW1_H

