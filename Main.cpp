#include <iostream>
#include "Matrix.h"  // Include the header for the Matrix class
#define EPS 1e-15
// Проверка на единичную матрицу
bool isIdentityMatrix(const Matrix& mat) {
    int n = mat.N;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                // Проверка, что на диагонали стоят единички
                if (std::abs(mat(i, j) - 1.0) > EPS) {
                    return false;
                }
            } else {
                // Проверка, что вне диагонали стоят нули
                if (std::abs(mat(i, j)) > EPS) {
                    return false;
                }
            }
        }
    }
    return true;
}

// Функция для нахождения отклонения от единичной матрицы
double findDeviation(const Matrix& mat) {
    int n = mat.N;
    double deviation = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                // Считаем отклонение диагональных элементов от 1
                deviation += std::abs(mat(i, j) - 1.0);
            } else {
                // Считаем отклонение недиагональных элементов от 0
                deviation += std::abs(mat(i, j));
            }
        }
    }
    return deviation;
}

int main() {
    // Ask the user to input the matrix size
    int size;
    std::cout << "Enter the size of the matrix (N x N): ";
    std::cin >> size;

    // Create a matrix of the given size
    Matrix matrix(size);

    // Input the elements of the matrix
    std::cout << "Enter the elements of the matrix (row by row):" << std::endl;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << "Element [" << i << "][" << j << "]: ";
            std::cin >> matrix(i, j);
        }
    }

    // Output the original matrix
    std::cout << "Original Matrix:" << std::endl;
    std::cout << matrix << std::endl;

    // Transpose the matrix
    Matrix transposedMatrix = matrix.transpose();
    std::cout << "Transposed Matrix:" << std::endl;
    std::cout << transposedMatrix << std::endl;

    // Calculate and display the determinant
    double determinant = matrix.determinant();
    std::cout << "Determinant of the matrix: " << determinant << std::endl;

    // Check if the matrix is invertible and print the inverse matrix if it exists
        Matrix inverseMatrix = matrix.inverse();  // Compute the inverse
        std::cout << "Inverse Matrix:" << std::endl;
        std::cout << inverseMatrix << std::endl;


    // Проверяем произведение matrix * inverseMatrix^-1 и inverseMatrix^-1 * matrix на единичность
        Matrix identity1 = matrix * inverseMatrix;
        Matrix identity2 = inverseMatrix * matrix;

        std::cout << "matrix * inverseMatrix:" << std::endl;
        std::cout << identity1 << std::endl;

        std::cout << "inverseMatrix * matrix:" << std::endl;
        std::cout << identity2 << std::endl;

        // Проверка на то, что это единичные матрицы
        if (isIdentityMatrix(identity1) && isIdentityMatrix(identity2)) {
            std::cout << "Matrix matrix and its inverse satisfy matrix * inverseMatrix = inverseMatrix * matrix = I (Identity matrix)." << std::endl;
        } else {
            std::cout << "Matrix matrix and its inverse do NOT satisfy matrix * inverseMatrix = inverseMatrix * matrix = I." << std::endl;

            // Находим отклонение
            double deviation1 = findDeviation(identity1);
            double deviation2 = findDeviation(identity2);

            std::cout << "Deviation of A * A^-1 from Identity Matrix: " << deviation1 << std::endl;
            std::cout << "Deviation of A^-1 * A from Identity Matrix: " << deviation2 << std::endl;
        }

    return 0;
}
