#include "hw1.h"
namespace algebra
{
    Matrix zeros(size_t n, size_t m)
    {
        Matrix matrix(n, std::vector<double>(m,0));
        return matrix;
    }

    Matrix ones(size_t n, size_t m)
    {
        Matrix matrix(n, std::vector<double>(m,1));
        return matrix;
    }

    Matrix random(size_t n, size_t m, double min, double max)
    {
        if(min >= max)
        {
            throw std::logic_error("Wrong input arguments!");
        }
        // using random library
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(min, max);

        // using push_back method
        Matrix matrix{};
        for(size_t i{}; i<n; i++)
        {
            std::vector<double> v;
            for(size_t j{}; j<m; j++)
            {
                v.push_back(dist(mt));
            }
            matrix.push_back(v);
        }
        return matrix;
    }

    void show(const Matrix& matrix)
    {
        for (size_t i{}; i < matrix.size(); i++)
        {
            for (size_t j{}; j < matrix[i].size(); j++)
            {
                std::cout << std::fixed << std::setprecision(3) << std::right << std::setw(10) << matrix[i][j] ;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    Matrix multiply(const Matrix& matrix, double c)
    {
        Matrix mat;
        for (size_t i{}; i < matrix.size(); i++)
        {
            std::vector<double> temp;
            for (size_t j{}; j < matrix[i].size(); j++)
            {
                temp.push_back(matrix[i][j] * c);
            }
            mat.push_back(temp);
        }
        return mat;
    }

    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2)
    {
        // n1 -> matrix1.size()
        // m1 -> matrix1[0].size()
        // n2 -> matrix2.size()
        // m2 -> matrix2[0].size()
        // ans -> n1 * m2
        if(matrix1.size() == 0 || matrix2.size() == 0 )
        {
            return matrix1;
        }
        if(matrix1[0].size() != matrix2.size())
        {
            throw std::logic_error("Matrices cannot be multiplied!");
        }
        Matrix mat;
        for(size_t k{}; k < matrix1.size() ; k++)
        {
            std::vector<double> Temp;
            for(size_t i{}; i < matrix2[0].size(); i++)
            {
                double temp{};
                for(size_t j{}; j < matrix1[0].size(); j++)
                {
                    temp+=(matrix1[k][j] * matrix2[j][i]);
                }
                Temp.push_back(temp);
            }
            mat.push_back(Temp);
        }
        return mat;
    }

    Matrix sum(const Matrix& matrix, double c)
    {
        Matrix mat;
        for (size_t i{}; i < matrix.size(); i++)
        {
            std::vector<double> temp;
            for (size_t j{}; j < matrix[i].size(); j++)
            {
                temp.push_back(matrix[i][j] + c);
            }
            mat.push_back(temp);
        }
        return mat;        
    }

    Matrix sum(const Matrix& matrix1, const Matrix& matrix2)
    {
        if(matrix1.size() == 0 && matrix2.size() == 0)
        {
            return matrix1;           
        }
        if((matrix1.size() != matrix2.size()) || (matrix1[0].size() != matrix2[0].size()))
        {
            throw std::logic_error("Matrices cannot be summed!");
        }        
        Matrix mat;
        for (size_t i{}; i < matrix1.size(); i++)
        {
            std::vector<double> temp;
            for (size_t j{}; j < matrix1[i].size(); j++)
            {
                temp.push_back(matrix1[i][j] + matrix2[i][j]);
            }
            mat.push_back(temp);
        }
        return mat;            
    }

    Matrix transpose(const Matrix& matrix)
    {
        if(matrix.size() == 0)
        {
            return matrix;
        }
        Matrix mat;
        for(size_t i{}; i < matrix[0].size(); i++)
        {
            std::vector<double> temp;
            for(size_t j{}; j < matrix.size(); j++)
            {
                temp.push_back(matrix[j][i]);
            }
            mat.push_back(temp);
        }
        return mat;
    }

    Matrix delcol(const Matrix& matrix, size_t c)
    {
        Matrix mat{matrix};
        if(matrix[0].size() > c && c >= 0)
        {
            for(size_t i{}; i < matrix.size(); i++)
            {
                mat[i].erase(mat[i].begin() + c);
            }
            return mat;
        }
        else
        {
            throw std::logic_error("Column cannot be deleted!");
        }
    }

    Matrix delrow(const Matrix& matrix, size_t r)
    {
        Matrix mat{matrix};
        if(matrix.size() > r && r >= 0)
        {
            mat.erase(mat.begin() + r);
            return mat;
        }
        else
        {
            throw std::logic_error("Row cannot be deleted!");
        }
    }

    Matrix minor(const Matrix& matrix, size_t n, size_t m)
    {
        Matrix r_del{delrow(matrix, n)};
        if(matrix.size() != r_del.size())
        {
            Matrix c_del{delcol(r_del, m)};
            if(c_del[0].size() != matrix[0].size())
            {
                return c_del;
            }
            else
            {
                return matrix;
            }
        }
        else
        {
            throw std::logic_error("Minor cannot be achieved!");
        }
    }

    double determinant(const Matrix& matrix)
    {
        if(matrix.size() == 0)
        {
            return 1.0;
        }
        if(matrix.size() != matrix[0].size())
        {
            throw std::logic_error("Determinant cannot be calculated!");
        }
        else
        {
            // 1*1
            if(matrix.size() == 1 )
                return matrix[0][0];

            // 2*2
            if(matrix.size() == 2)
                return ( matrix[0][0] * matrix[1][1] ) - ( matrix[0][1] * matrix[1][0] ); 

            // 3*3
            if(matrix.size() == 3)
            {
                return (matrix[0][0] * determinant(minor(matrix, 0, 0)) - matrix[0][1] * determinant(minor(matrix, 0, 1)) + matrix[0][2] * determinant(minor(matrix, 0, 2))); 
            }

            // n*n 
            size_t i{};
            double det{};
            for (size_t j{}; j < matrix.size(); j++)
            {
                det += pow(-1, i + j) * matrix[i][j] * determinant(minor(matrix, i, j));
            }
            return det;            
        }
   }

    Matrix cofactor(const Matrix& matrix)
    {
        if(matrix.size() != matrix[0].size())
        {
            throw std::logic_error("Cofactor Matrix cannot be calculated!");
        }     

        Matrix c;
        for(size_t i{}; i < matrix.size(); i++)
        {
            std::vector<double> temp;
            for(size_t j{}; j < matrix[0].size(); j++)
            {
                temp.push_back(pow(-1, i + j) * determinant(minor(matrix, i, j)));
            }
            c.push_back(temp);
        }
        return c;
    }

    Matrix Adjugate(const Matrix& matrix)
    {
        return transpose(cofactor(matrix));
    }

    Matrix inverse(const Matrix& matrix)
    {
        if(matrix.size() == 0)
        {
            return matrix;
        }
        if(determinant(matrix) == 0)
        {
            throw std::logic_error("Singular matrix/ has no inverse");
        }
        return multiply(Adjugate(matrix), 1 / determinant(matrix)); 
    }

    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis=0)
    {
        if(axis == 0)
        {
            if(matrix1[0].size() != matrix2[0].size())
            {
                throw std::logic_error("Cannot concatenate two Matrices!");
            }
            else
            {
                Matrix c;
                for(size_t i{}; i < matrix1.size(); i++)
                {
                    std::vector<double> temp1;
                    for(size_t j{}; j < matrix1[0].size(); j++)
                    {
                        temp1.push_back(matrix1[i][j]);
                    }
                    c.push_back(temp1);
                }
                for(size_t i{}; i < matrix2.size(); i++)
                {
                    std::vector<double> temp2;
                    for(size_t j{}; j < matrix2[0].size(); j++)
                    {
                        temp2.push_back(matrix2[i][j]);
                    }
                    c.push_back(temp2);
                }
                return c;
            }
        }
        if(axis == 1)
        {
            if(matrix1.size() != matrix2.size())
            {
                throw std::logic_error("Cannot concatenate two Matrices!");
            }
            else
            {
                Matrix c;
                for(size_t i{}; i < matrix1.size(); i++)
                {
                    std::vector<double> temp;
                    for(size_t j{}; j < matrix1[0].size(); j++)
                    {
                        temp.push_back(matrix1[i][j]);
                    }
                    for(size_t k{}; k < matrix2[0].size(); k++)
                    {
                        temp.push_back(matrix2[i][k]);
                    }                    
                    c.push_back(temp);
                }
                return c;                
            }            
        }
        if(axis != 0 && axis != 1)
        {
            return zeros(matrix1.size(), matrix1[0].size());
        }
        return zeros(matrix1.size(), matrix1[0].size());
    }

    Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2)
    {
        if(matrix.size() > r1 && matrix.size() > r2)
        {
            Matrix mat{matrix};
            std::vector<double> temp;
            for(size_t i{}; i < matrix[0].size(); i++)
            {
                temp.push_back(matrix[r1][i]);
                mat[r1].erase(mat[r1].begin() + i);
                mat[r1].insert(mat[r1].begin() + i, mat[r2][i]);
                mat[r2].erase(mat[r2].begin() + i);
                mat[r2].insert(mat[r2].begin() + i, temp[i]);
            }
            return mat;
        }
        else
        {
            throw std::logic_error("r1 or r2 inputs are out of range");          
        }
    }

    Matrix ero_multiply(const Matrix& matrix, size_t r, double c)
    {
        if(matrix.size() > r)
        {
            Matrix mat;
            for(size_t i{}; i < matrix.size(); i++)
            {
                std::vector<double> temp;
                if(i == r)
                {
                    for(size_t j{}; j < matrix[0].size(); j++)
                    {
                        temp.push_back(matrix[i][j] * c);
                    }
                    mat.push_back(temp);                    
                }
                else
                {
                    for(size_t j{}; j < matrix[0].size(); j++)
                    {
                        temp.push_back(matrix[i][j]);
                    }
                    mat.push_back(temp);
                }
            }            
            return mat;
        }
        else
        {
            throw std::logic_error("Row cannot be Multiplied!");                     
        }
    }

    Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2)
    {
        Matrix mat{matrix};
        for(size_t i{}; i < mat[0].size(); i++)
        {
            mat[r2][i] += mat[r1][i] * c;
        }
        return mat;
    }

    Matrix upper_triangular(const Matrix& matrix)
    {
        if(matrix.size() == 0)
        {
            return matrix;
        }
        if(matrix.size() != matrix[0].size())
        {
            throw std::logic_error("Non-square matrix"); 
        }
        else
        {
            Matrix mat{matrix};
            for(size_t i{}; i < mat.size() - 1; i++)
            {
                // Detecting the zero pivot
                if(mat[i][i] == 0)
                {
                    size_t ind{i + 1};
                    size_t log_err{};
                    double cmp{mat[i + 1][i]};
                    for(size_t k{i + 1}; k < mat.size(); k++)
                    {
                        // Detecting the row to be exchanged
                        if( mat[k][i] != 0 )
                        {
                            cmp = mat[k][i];
                            ind = k;
                            break;                            
                        }
                        // Detecting Singular Matrix Error
                        if (mat[k][i] == 0)
                        {
                            log_err++;
                        }
                        if(log_err == mat.size() - i - 1)
                        {
                            throw std::logic_error("Singular Matrix");
                        }
                    }
                    mat = ero_swap(mat, i, ind);
                }
                    for(size_t j{}; j < mat[0].size() - i - 1; j++)
                    {
                        mat = ero_sum(mat, i, -mat[i + j + 1][i] / mat[i][i], i + j + 1);
                    }
            }
            return mat;
        }
    }
}