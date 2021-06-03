#ifndef VECTOR_H
#define VECTOR_H

#include <cstdint>
#include <cstddef>
#include <iostream>

// Класс n-мерного вещественного вектора
class Vector {
 public:
  // Конструктор с возможностью задания длины n и выбора
  // между последовательной и параллельной версией
  Vector(const size_t n, const bool is_parallel);
  // Конструктор копирования
  Vector(const Vector& obj);
  // Оператор присваивания
  Vector& operator=(const Vector& obj);
  ~Vector();
  size_t GetN() const;
  // Оператор индексации
  double& operator[](const ptrdiff_t number);
  // Константная версия оператора индексации
  const double& operator[](const ptrdiff_t number) const;
  // Присваивающий оператор сложения двух векторов
  Vector& operator+=(const Vector& rhs);
  // Присваивающий оператор вычитания двух векторов
  Vector& operator-=(const Vector& rhs);
  // Присваивающий оператор умножения вектора на число
  Vector& operator*=(const double rhs);
  // Скалярное произведение двух векторов
  double DotProduct(const Vector& rhs) const;
  // Вычисление длины вектора
  double Length() const;
  std::ostream& WriteTo(std::ostream& ostrm) const;

 private:
  size_t n_{ 0 };
  bool is_parallel_{ false };
  double* data_{ nullptr };
};

// Сложение двух векторов
Vector operator+(const Vector& lhs, const Vector& rhs);
// Вычитание двух векторов
Vector operator-(const Vector& lhs, const Vector& rhs);
// Умножение вектора на число
Vector operator*(const double lhs, const Vector& rhs);
// Умножение вектора на число
Vector operator*(const Vector& lhs, const double rhs);
inline std::ostream& operator<<(std::ostream& ostrm, Vector& rhs) {
  return rhs.WriteTo(ostrm);
}

#endif
