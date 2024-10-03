#include <iostream>

class Matrix {
public:
    int N;             // Размер матрицы N x N
    double* data;      // Одномерный массив для хранения элементов матрицы
    // Конструкторы и деструктор
    Matrix(int size);                           // Конструктор с размером
    Matrix(const Matrix& other);                // Конструктор копирования
    ~Matrix();                                  // Деструктор

    // Оператор присваивания
    Matrix& operator=(const Matrix& other);

    // Арифметические операции
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix& operator*=(const Matrix& other);
    Matrix& operator+=(const Matrix& other);
    Matrix& operator-=(const Matrix& other);

    // Операции сравнения
    bool operator==(const Matrix& other) const;
    bool operator!=(const Matrix& other) const;
    bool operator<(const Matrix& other) const;
    bool operator<=(const Matrix& other) const;
    bool operator>(const Matrix& other) const;
    bool operator>=(const Matrix& other) const;

    // Операторы доступа к элементам
    double& operator()(int row, int col);
    const double& operator()(int row, int col) const;

    // Операторы ввода/вывода
    friend std::istream& operator>>(std::istream& in, Matrix& matrix);
    friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix);

    // Другие методы
    Matrix transpose() const;                   // Транспонирование
    double determinant() const;                 // Вычисление определителя
    Matrix inverse() const;                     // Обратная матрица (если существует)

    // Элементарные преобразования над строками и столбцами
    void swapRows(int i, int j);                // Обмен строк
    void swapCols(int i, int j);                // Обмен столбцов
    void multiplyRow(int i, double scalar);     // Умножение строки на число
    void multiplyCol(int i, double scalar);     // Умножение столбца на число
    void addRowMultiple(int i, int j, double scalar);  // Добавление к строке другой строки, умноженной на число
    void addColMultiple(int i, int j, double scalar);  // Добавление к столбцу другого столбца, умноженного на число
    Matrix gaussianElimination() const;  // Gaussian elimination to row echelon form
};
