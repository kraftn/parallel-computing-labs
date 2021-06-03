#include <stdexcept>

#include "matrix.h"

Matrix::Matrix(const size_t n, const bool is_parallel)
  : n_(n), is_parallel_(is_parallel) {
  if (n <= 0) {
    throw std::invalid_argument("Wrong size");
  }
  data_ = new double*[n];
  for (ptrdiff_t i_rows(0); i_rows < n; i_rows += 1) {
    data_[i_rows] = new double[n] { 0.0 };
  }
}



Matrix::Matrix(const Matrix& obj)
  : n_(obj.n_), is_parallel_(obj.is_parallel_), data_(new double* [obj.n_]) {
  for (ptrdiff_t i_rows(0); i_rows < obj.n_; i_rows += 1) {
    data_[i_rows] = new double[obj.n_] { 0.0 };
  }

  for (ptrdiff_t i_rows(0); i_rows < obj.n_; i_rows += 1) {
    for (ptrdiff_t i_columns(0); i_columns < obj.n_; i_columns += 1) {
      data_[i_rows][i_columns] = obj.data_[i_rows][i_columns];
    }
  }
}



Matrix& Matrix::operator=(const Matrix& obj) {
  // Выполняем присваивание, если адреса двух объектов различны
  if (this != &obj) {
    if (n_ != obj.n_) {
      throw std::invalid_argument("Sizes of Matrix aren't equal");
    }
    for (ptrdiff_t i_rows(0); i_rows < n_; i_rows += 1) {
      for (ptrdiff_t i_columns(0); i_columns < n_; i_columns += 1) {
        data_[i_rows][i_columns] = obj.data_[i_rows][i_columns];
      }
    }
  }
  return *this;
}



Matrix::~Matrix() {
  Clean();
}



size_t Matrix::GetN() const {
  return n_;
}



double& Matrix::Value(const ptrdiff_t row, const ptrdiff_t column) {
  if (row < 0 || row >= n_ || column < 0 || column >= n_) {
    throw std::out_of_range("Index is out of range of Matrix");
  }
  return data_[row][column];
}



const double& Matrix::Value(const ptrdiff_t row, const ptrdiff_t column) const {
  if (row < 0 || row >= n_ || column < 0 || column >= n_) {
    throw std::out_of_range("Index is out of range of Matrix");
  }
  return data_[row][column];
}



Matrix& Matrix::Transpose() {
  double** new_data = new double*[n_];

  #pragma omp parallel for if (is_parallel_)
  for (ptrdiff_t i_rows(0); i_rows < n_; i_rows += 1) {
    // Создадим новую строку матрицы
    new_data[i_rows] = new double[n_] {0.0};
    for (ptrdiff_t i_columns(0); i_columns < n_; i_columns += 1) {
      new_data[i_rows][i_columns] = data_[i_columns][i_rows];
    }
  }

  // Удалим старые данные
  Clean();
  data_ = new_data;

  return *this;
}



Matrix& Matrix::operator+=(const Matrix& rhs) {
  if (n_ != rhs.n_) {
    throw std::invalid_argument("Sizes of Matrix aren't equal");
  }
  #pragma omp parallel for if (is_parallel_)
  for (ptrdiff_t i_rows(0); i_rows < n_; i_rows += 1) {
    for (ptrdiff_t i_columns(0); i_columns < n_; i_columns += 1) {
      data_[i_rows][i_columns] += rhs.data_[i_rows][i_columns];
    }
  }
  return *this;
}



Matrix& Matrix::operator-=(const Matrix& rhs) {
  if (n_ != rhs.n_) {
    throw std::invalid_argument("Sizes of Matrix aren't equal");
  }
  #pragma omp parallel for if (is_parallel_) 
  for (ptrdiff_t i_rows(0); i_rows < n_; i_rows += 1) {
    for (ptrdiff_t i_columns(0); i_columns < n_; i_columns += 1) {
      data_[i_rows][i_columns] -= rhs.data_[i_rows][i_columns];
    }
  }
  return *this;
}



Matrix& Matrix::operator*=(const double rhs) {
  #pragma omp parallel for if (is_parallel_)
  for (ptrdiff_t i_rows(0); i_rows < n_; i_rows += 1) {
    for (ptrdiff_t i_columns(0); i_columns < n_; i_columns += 1) {
      data_[i_rows][i_columns] *= rhs;
    }
  }
  return *this;
}



Matrix& Matrix::operator*=(const Matrix& rhs) {
  if (n_ != rhs.n_) {
    throw std::invalid_argument("The number of columns of the first matrix isn't equal to the number of rows of the second matrix");
  }

  #pragma omp parallel for if (is_parallel_)
  for (ptrdiff_t i_rows(0); i_rows < n_; i_rows += 1) {
    // Новая строка матрицы
    double* new_data = new double[n_] {0.0};
    for (ptrdiff_t i_columns(0); i_columns < n_; i_columns += 1) {
      for (ptrdiff_t k(0); k < n_; k += 1) {
        new_data[i_columns] += data_[i_rows][k] * rhs.data_[k][i_columns];
      }
    }
    // Удалим старые данные
    delete[] data_[i_rows];
    data_[i_rows] = new_data;
  }

  return *this;
}



Vector Matrix::operator*(const Vector& rhs) const {
  if (n_ != rhs.GetN()) {
    throw std::invalid_argument("Sizes of Matrix and Vector aren't equal");
  }

  Vector res(n_, is_parallel_);
  #pragma omp parallel for if (is_parallel_)
  for (ptrdiff_t i_rows(0); i_rows < n_; i_rows += 1) {
    for (ptrdiff_t k(0); k < n_; k += 1) {
      res[i_rows] += data_[i_rows][k] * rhs[k];
    }
  }

  return res;
}



double Matrix::FrobeniusNorm() const {
  double sum(0.0);
  #pragma omp parallel for reduction(+:sum) if (is_parallel_)
  for (ptrdiff_t i_rows(0); i_rows < n_; i_rows += 1) {
    for (ptrdiff_t i_columns(0); i_columns < n_; i_columns += 1) {
      sum += pow(data_[i_rows][i_columns], 2);
    }
  }
  return sqrt(sum);
}



std::ostream& Matrix::WriteTo(std::ostream& ostrm) const {
  for (ptrdiff_t i_rows(0); i_rows < n_; i_rows += 1) {
    for (ptrdiff_t i_columns(0); i_columns < n_; i_columns += 1) {
      ostrm << data_[i_rows][i_columns] << ' ';
    }
    ostrm << '\n';
  }
  return ostrm;
}



void Matrix::Clean() {
  for (ptrdiff_t i_rows(0); i_rows < n_; i_rows += 1) {
    delete[] data_[i_rows];
  }
  delete[] data_;
}



Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
  Matrix sum(lhs);
  sum += rhs;
  return sum;
}



Matrix operator-(const Matrix& lhs, const Matrix& rhs) {
  Matrix dif(lhs);
  dif -= rhs;
  return dif;
}



Matrix operator*(const Matrix& lhs, const double rhs) {
  Matrix res(lhs);
  res *= rhs;
  return res;
}



Matrix operator*(const double lhs, const Matrix& rhs) {
  Matrix res(rhs);
  res *= lhs;
  return res;
}


Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
  Matrix res(lhs);
  res *= rhs;
  return res;
}
