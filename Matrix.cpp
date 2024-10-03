#include "Matrix.h"
#include <iostream>
#include <stdexcept>
#include <iomanip>   // Для форматирования вывода
#include <cmath>     // Для вычислений
#define EPS 1e-12


// Конструктор с размером
Matrix::Matrix(int size) : N(size), data(new double[size * size]()) {}

// Конструктор копирования
Matrix::Matrix(const Matrix& other) : N(other.N), data(new double[other.N * other.N]) {
    for (int i = 0; i < N * N; ++i) {
        data[i] = other.data[i];
    }
}

// Деструктор
Matrix::~Matrix() {
    delete[] data;
}

// Оператор присваивания
Matrix& Matrix::operator=(const Matrix& other) {
    if (this != &other) {
        delete[] data;
        N = other.N;
        data = new double[N * N];
        for (int i = 0; i < N * N; ++i) {
            data[i] = other.data[i];
        }
    }
    return *this;
}

// Операция сложения
Matrix Matrix::operator+(const Matrix& other) const {
    if (N != other.N) throw std::invalid_argument("Размеры матриц должны совпадать");
    Matrix result(N);
    for (int i = 0; i < N * N; ++i) {
        result.data[i] = data[i] + other.data[i];
    }
    return result;
}

// Операция вычитания
Matrix Matrix::operator-(const Matrix& other) const {
    if (N != other.N) throw std::invalid_argument("Размеры матриц должны совпадать");
    Matrix result(N);
    for (int i = 0; i < N * N; ++i) {
        result.data[i] = data[i] - other.data[i];
    }
    return result;
}

// Операция умножения
Matrix Matrix::operator*(const Matrix& other) const {
    if (N != other.N) throw std::invalid_argument("Размеры матриц должны совпадать");
    Matrix result(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            result(i, j) = 0;
            for (int k = 0; k < N; ++k) {
                result(i, j) += (*this)(i, k) * other(k, j);
            }
        }
    }
    return result;
}

// Операция +=
Matrix& Matrix::operator+=(const Matrix& other) {
    if (N != other.N) throw std::invalid_argument("Размеры матриц должны совпадать");
    for (int i = 0; i < N * N; ++i) {
        data[i] += other.data[i];
    }
    return *this;
}

// Операция -=
Matrix& Matrix::operator-=(const Matrix& other) {
    if (N != other.N) throw std::invalid_argument("Размеры матриц должны совпадать");
    for (int i = 0; i < N * N; ++i) {
        data[i] -= other.data[i];
    }
    return *this;
}

// Операция *=
Matrix& Matrix::operator*=(const Matrix& other) {
    *this = *this * other;
    return *this;
}

// Операторы сравнения
bool Matrix::operator==(const Matrix& other) const {
    if (N != other.N) return false;
    for (int i = 0; i < N * N; ++i) {
        if (data[i] != other.data[i]) return false;
    }
    return true;
}

bool Matrix::operator!=(const Matrix& other) const {
    return !(*this == other);
}

bool Matrix::operator<(const Matrix& other) const {
    if (N != other.N) throw std::invalid_argument("Размеры матриц должны совпадать");
    for (int i = 0; i < N * N; ++i) {
        if (data[i] >= other.data[i]) return false;
    }
    return true;
}

bool Matrix::operator<=(const Matrix& other) const {
    return *this < other || *this == other;
}

bool Matrix::operator>(const Matrix& other) const {
    if (N != other.N) throw std::invalid_argument("Размеры матриц должны совпадать");
    for (int i = 0; i < N * N; ++i) {
        if (data[i] <= other.data[i]) return false;
    }
    return true;
}

bool Matrix::operator>=(const Matrix& other) const {
    return *this > other || *this == other;
}

// Оператор доступа к элементам (индексация)
double& Matrix::operator()(int row, int col) {
    if (row < 0 || row >= N || col < 0 || col >= N) {
        throw std::out_of_range("Matrix index out of bounds");
    }
    return data[row * N + col];
}

const double& Matrix::operator()(int row, int col) const {
    if (row < 0 || row >= N || col < 0 || col >= N) {
        throw std::out_of_range("Matrix index out of bounds");
    }
    return data[row * N + col];
}

// Операторы ввода/вывода
std::istream& operator>>(std::istream& in, Matrix& matrix) {
    for (int i = 0; i < matrix.N * matrix.N; ++i) {
        in >> matrix.data[i];
    }
    return in;
}

std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
    for (int i = 0; i < matrix.N; ++i) {
        for (int j = 0; j < matrix.N; ++j) {
            out << std::setw(1) << matrix(i, j) << " ";
        }
        out << std::endl;
    }
    return out;
}

// Транспонирование
Matrix Matrix::transpose() const {
    Matrix result(N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            result(j, i) = (*this)(i, j);
    return result;
}

// Вычисление определителя (метод Гаусса)
double Matrix::determinant() const {
    Matrix temp(*this);
    double det = 1;
    for (int i = 0; i < N; ++i) {
        if (temp(i, i) == 0) {
            bool found = false;
            for (int j = i + 1; j < N; ++j) {
                if (temp(j, i) != 0) {
                    temp.swapRows(i, j);
                    det *= -1;
                    found = true;
                    break;
                }
            }
            if (!found) return 0;  // Если столбец полностью нулевой
        }
        det *= temp(i, i);
        for (int j = i + 1; j < N; ++j) {
            double ratio = temp(j, i) / temp(i, i);
            for (int k = i; k < N; ++k) {
                temp(j, k) -= ratio * temp(i, k);
            }
        }
    }
    return det;
}

// Обратная матрица
Matrix Matrix::inverse() const {
    if ((determinant() < 0+EPS) && (determinant() > 0-EPS)) throw std::invalid_argument("Determinant = 0");
    Matrix result(N);
    Matrix temp(*this);

    // Создание единичной матрицы
    for (int i = 0; i < N; ++i)
        result(i, i) = 1;

    // Приведение к диагональному виду с параллельным преобразованием единичной матрицы
    for (int i = 0; i < N; ++i) {
        double pivot = temp(i, i);
        for (int j = 0; j < N; ++j) {
            temp(i, j) /= pivot;
            result(i, j) /= pivot;
        }
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                double ratio = temp(j, i);
                for (int k = 0; k < N; ++k) {
                    temp(j, k) -= ratio * temp(i, k);
                    result(j, k) -= ratio * result(i, k);
                }
            }
        }
    }
    return result;
}

// Элементарные преобразования

// Обмен строк
void Matrix::swapRows(int i, int j) {
    for (int k = 0; k < N; ++k) {
        std::swap(data[i * N + k], data[j * N + k]);
    }
}

// Обмен столбцов
void Matrix::swapCols(int i, int j) {
    for (int k = 0; k < N; ++k) {
        std::swap(data[k * N + i], data[k * N + j]);
    }
}

// Умножение строки на число
void Matrix::multiplyRow(int i, double scalar) {
    for (int k = 0; k < N; ++k) {
        data[i * N + k] *= scalar;
    }
}

// Умножение столбца на число
void Matrix::multiplyCol(int i, double scalar) {
    for (int k = 0; k < N; ++k) {
        data[k * N + i] *= scalar;
    }
}

// Добавление к строке другой строки, умноженной на число
void Matrix::addRowMultiple(int i, int j, double scalar) {
    for (int k = 0; k < N; ++k) {
        data[i * N + k] += scalar * data[j * N + k];
    }
}

// Добавление к столбцу другого столбца, умноженного на число
void Matrix::addColMultiple(int i, int j, double scalar) {
    for (int k = 0; k < N; ++k) {
        data[k * N + i] += scalar * data[k * N + j];
    }
}

Matrix Matrix::gaussianElimination() const {
    Matrix result = *this;  // Копия исходной матрицы

    int n = result.N;  // Размер матрицы

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; ++i) {
        // Поиск максимального элемента в текущем столбце (для выбора главного элемента)
        double maxElement = result(i, i);
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(result(k, i)) > std::abs(maxElement)) {
                maxElement = result(k, i);
                maxRow = k;
            }
        }

        // Поменять местами текущую строку и строку с максимальным элементом
        result.swapRows(i, maxRow);

        // Обнуляем все строки ниже текущей в текущем столбце
        for (int k = i + 1; k < n; ++k) {
            double coefficient = -result(k, i) / result(i, i);
            result.addRowMultiple(k, i, coefficient);
        }
    }

    // Обратный ход метода Гаусса
    for (int i = n - 1; i >= 0; --i) {
        // Нормализуем строку: делим на ведущий элемент (если он не равен 1)
        if (result(i, i) != 1.0) {
            result.multiplyRow(i, 1.0 / result(i, i));
        }

        // Обнуляем все элементы выше текущего в текущем столбце
        for (int k = i - 1; k >= 0; --k) {
            double coefficient = -result(k, i);
            result.addRowMultiple(k, i, coefficient);
        }
    }

    return result;
}
