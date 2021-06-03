#include <stdexcept>
#include <omp.h>
#include <cmath>

#include "vector.h"

Vector::Vector(const size_t n, const bool is_parallel)
  : n_(n), is_parallel_(is_parallel) {
  if (n <= 0) {
    throw std::invalid_argument("Wrong size");
  }
  data_ = new double[n_] { 0.0 };
}



Vector::Vector(const Vector& obj)
  : n_(obj.n_), is_parallel_(obj.is_parallel_), data_(new double[n_] {0.0}) {
  for (ptrdiff_t i_data(0); i_data < obj.n_; i_data += 1) {
    data_[i_data] = obj.data_[i_data];
  }
}



Vector& Vector::operator=(const Vector& obj) {
  // Выполняем присваивание, если адреса двух объектов различны
  if (this != &obj) {
    if (n_ != obj.n_) {
      throw std::invalid_argument("Lengths of vectors aren't equal");
    }
    for (ptrdiff_t i_data(0); i_data < obj.n_; i_data++) {
      data_[i_data] = obj.data_[i_data];
    }
  }
  return *this;
}



Vector::~Vector() {
  delete[] data_;
}



size_t Vector::GetN() const {
  return n_;
}



double& Vector::operator[](const ptrdiff_t number) {
  if (number < 0 || number >= n_) {
    throw std::out_of_range("Index is out of range of Vector");
  }
  return data_[number];
}



const double& Vector::operator[](const ptrdiff_t number) const {
  if (number < 0 || number >= n_) {
    throw std::out_of_range("Index is out of range of Vector");
  }
  return data_[number];
}



Vector& Vector::operator+=(const Vector& rhs) {
  if (n_ != rhs.n_) {
    throw std::invalid_argument("Lengths of vectors aren't equal");
  }
  #pragma omp parallel for if (is_parallel_)
  for (ptrdiff_t i_data(0); i_data < n_; i_data += 1) {
    data_[i_data] += rhs.data_[i_data];
  }
  return *this;
}



Vector& Vector::operator-=(const Vector& rhs) {
  if (n_ != rhs.n_) {
    throw std::invalid_argument("Lengths of vectors aren't equal");
  }
  #pragma omp parallel for if (is_parallel_)
  for (ptrdiff_t i_data(0); i_data < n_; i_data += 1) {
    data_[i_data] -= rhs.data_[i_data];
  }
  return *this;
}



Vector& Vector::operator*=(const double rhs) {
  #pragma omp parallel for if (is_parallel_)
  for (ptrdiff_t i_data(0); i_data < n_; i_data += 1) {
    data_[i_data] *= rhs;
  }
  return *this;
}



double Vector::DotProduct(const Vector& rhs) const {
  if (n_ != rhs.n_) {
    throw std::invalid_argument("Lengths of vectors aren't equal");
  }
  
  double sum(0.0);
  // Найдём в каждой нити значения частичных сумм, после чего сложим их
  #pragma omp parallel for reduction(+:sum) if (is_parallel_)
  for (ptrdiff_t i(0); i < n_; i += 1) {
    sum += data_[i] * rhs.data_[i];
  }

  return sum;
}



double Vector::Length() const {
  // ||x|| = sqrt(<x, x>)
  return sqrt(DotProduct(*this));
}



std::ostream& Vector::WriteTo(std::ostream& ostrm) const {
  for (ptrdiff_t i_data(0); i_data < n_; i_data += 1) {
    ostrm << data_[i_data] << '\n';
  }
  return ostrm;
}



Vector operator+(const Vector& lhs, const Vector& rhs) {
  Vector res(lhs);
  res += rhs;
  return res;
}



Vector operator-(const Vector& lhs, const Vector& rhs) {
  Vector res(lhs);
  res -= rhs;
  return res;
}



Vector operator*(const double lhs, const Vector& rhs) {
  Vector res(rhs);
  res *= lhs;
  return res;
}



Vector operator*(const Vector& lhs, const double rhs) {
  Vector res(lhs);
  res *= rhs;
  return res;
}
