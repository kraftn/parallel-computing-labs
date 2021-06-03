#ifndef MATRIX_H
#define MATRIX_H

#include <cstdint>
#include <cstddef>
#include <iostream>

#include "vector.h"

// Класс вещественной квадратной матрицы размера n x n
class Matrix {
 public:
  // Конструктор с возможностью задания значения n и выбора
  // между последовательной и параллельной версией
  Matrix(const size_t n, const bool is_parallel);
  // Конструктор копирования
  Matrix(const Matrix& obj);
  // Оператор присваивания
  Matrix& operator=(const Matrix& obj);
  ~Matrix();
  size_t GetN() const;
  // Метод для доступа к элементам матрицы
  double& Value(const ptrdiff_t row, const ptrdiff_t column);
  // Константный метод для доступа к элементам матрицы
  const double& Value(const ptrdiff_t row, const ptrdiff_t column) const;
  // Операция транспонирования
  Matrix& Transpose();
  // Присваивающий оператор сложения двух матриц
  Matrix& operator+=(const Matrix& rhs);
  // Присваивающий оператор вычитания двух матриц
  Matrix& operator-=(const Matrix& rhs);
  // Присваивающий оператор умножения матрицы на число
  Matrix& operator*=(const double rhs);
  // Присваивающий оператор умножения двух матриц
  Matrix& operator*=(const Matrix& rhs);
  // Умножение матрицы на вектор
  Vector operator*(const Vector& rhs) const;
  // Вычисление нормы Фробениуса
  double FrobeniusNorm() const;
  std::ostream& WriteTo(std::ostream& ostrm) const;

 private:
  // Очистить память
  void Clean();
  size_t n_{ 0 };
  bool is_parallel_{ false };
  double** data_{ nullptr };
};

// Сложение двух матриц
Matrix operator+(const Matrix& lhs, const Matrix& rhs);
// Вычитание двух матриц
Matrix operator-(const Matrix& lhs, const Matrix& rhs);
// Умножение матрицы на число
Matrix operator*(const Matrix& lhs, const double rhs);
// Умножение матрицы на число
Matrix operator*(const double lhs, const Matrix& rhs);
// Умножение двух матриц
Matrix operator*(const Matrix& lhs, const Matrix& rhs);
inline std::ostream& operator<<(std::ostream& ostrm, const Matrix& rhs) {
  return rhs.WriteTo(ostrm);
}

#endif
