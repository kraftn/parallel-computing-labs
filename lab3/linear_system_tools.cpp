#include <cstdint>
#include <cstddef>
#include <iostream>
#include <omp.h>

#include "linear_system_tools.h"

std::vector<std::vector<double>> LinearSystemTools::GenerateA(const size_t n) {
  using namespace std;
  // Генератор для выбора случайной строки
  std::uniform_int_distribution<> indexes(0, n - 1);
  // Генератор для выбора значений диагональных элементов матрицы
  std::uniform_real_distribution<> determinant(8.0, 10.0);
  // Генератор для выбора значения, на которое домножаются элементы строки
  std::uniform_real_distribution<> c_gen(-0.1, 0.1);

  vector<vector<double>> a(n);
  // Имеет ли место в матрице диагональное преобладание
  bool cond_diagonal(false);
  while (!cond_diagonal) {
    for (ptrdiff_t i(0); i < n; i += 1) {
      a[i].resize(n, 0.0);
    }

    // Расположить на главной диагонали ненулевые элементы
    for (ptrdiff_t i(0); i < n; i += 1) {
      a[i][i] = determinant(gen_);
    }

    // Прибавить n раз к случайно выбранной строке другую случайно выбранную строку,
    // домноженную на случайное число
    for (ptrdiff_t i_shuffle(0); i_shuffle < n; i_shuffle += 1) {
      int32_t x(indexes(gen_));
      int32_t y(indexes(gen_));
      if (y == x) {
        y = (y + 1) % n;
      }
      for (ptrdiff_t i_column(0); i_column < n; i_column += 1) {
        a[x][i_column] += c_gen(gen_) * a[y][i_column];
      }
    }

    double epsilon(0.000001);
    ptrdiff_t i(0);
    cond_diagonal = true;
    while (cond_diagonal && (i < n)) {
      // На главной диагонали матрицы не должно быть нулевых элементов
      if (abs(a[i][i]) < epsilon) {
        ptrdiff_t i_add(0);
        while ((i_add < n) && (abs(a[i_add][i]) < epsilon)) {
          i_add += 1;
        }
        for (ptrdiff_t i_column(0); i_column < n; i_column += 1) {
          a[i][i_column] += a[i_add][i_column];
        }
      }

      // Проверка на диагональное преобладание
      double sum(0.0);
      for (ptrdiff_t i_column(0); i_column < n; i_column += 1) {
        if (i_column != i) {
          sum += abs(a[i][i_column]);
        }
      }
      cond_diagonal = abs(a[i][i]) > sum;
      i += 1;
    }
  }

  return a;
}



std::vector<double> LinearSystemTools::GenerateB(const double n) {
  using namespace std;

  // Генератор для случайного заполнения вектора правой части СЛАУ
  std::uniform_real_distribution<> elements(-10, 10);
  vector<double> b(n, 0);

  for (ptrdiff_t i(0); i < n; i += 1) {
    b[i] = elements(gen_);
  }

  return b;
}



std::vector<double> LinearSystemTools::GaussianElimination(
    std::vector<std::vector<double>> a, std::vector<double> b,
    const bool is_parallel) const {

  using namespace std;

  size_t n(a.size());
  // Прямой ход метода Гаусса
  for (ptrdiff_t k(0); k < n - 1; k += 1) {
    #pragma omp parallel for if (is_parallel)
    for (ptrdiff_t i(k + 1); i < n; i += 1) {
      double t(a[i][k] / a[k][k]);
      b[i] -= t * b[k];
      for (ptrdiff_t j(k + 1); j < n; j += 1) {
        a[i][j] -= t * a[k][j];
      }
    }
  }

  vector<double> x(n, 0);
  // Обратный ход
  x[n - 1] = b[n - 1] / a[n - 1][n - 1];
  for (ptrdiff_t k(n - 2); k > -1; k -= 1) {
    double sum(0.0);
    #pragma omp parallel for reduction(+:sum) if (is_parallel)
    for (ptrdiff_t j(n - 1); j > k; j -= 1) {
      sum += a[k][j] * x[j];
    }

    x[k] = (b[k] - sum) / a[k][k];
  }

  return x;
}



std::vector<double> LinearSystemTools::JacobiAlgorithm(
    std::vector<std::vector<double>> a, std::vector<double> b,
    const double accuracy, const bool is_parallel) const {

  using namespace std;

  size_t n(a.size());
  vector<double> x(n, 0);
  // Выбор начального приближения решения
  #pragma omp parallel for if (is_parallel)
  for (ptrdiff_t i(0); i < n; i += 1) {
    x[i] = b[i] / a[i][i];
  }

  vector<double> transitional_x(n, 0);
  // Норма разницы между векторами
  double diff(accuracy + 1.0);

  // Вычисление нормы матрицы B
  vector<vector<double>> a_modified(n);
  #pragma omp parallel for if (is_parallel)
  for (ptrdiff_t i(0); i < n; i += 1) {
    a_modified[i].resize(n, 0);
    for (ptrdiff_t j(0); j < n; j += 1) {
      if (i == j) {
        a_modified[i][j] = 0;
      }
      else {
        a_modified[i][j] = -a[i][j] / a[i][i];
      }
    }
  }
  double q(Norm(a_modified, is_parallel));

  while (diff > accuracy) {
    transitional_x = x;

    #pragma omp parallel for if (is_parallel)
    for (ptrdiff_t i(0); i < n; i += 1) {
      double sum(0.0);
      for (ptrdiff_t j(0); j < n; j += 1) {
        if (i != j) {
          sum += a[i][j] * transitional_x[j];
        }
      }
      x[i] = -(sum - b[i]) / a[i][i];
    }

    // Вычисление разницы векторов
    #pragma omp parallel for if (is_parallel)
    for (ptrdiff_t i(0); i < n; i += 1) {
      transitional_x[i] -= x[i];
    }

    diff = Norm(transitional_x, is_parallel) * q / (1 - q);
  }

  return x;
}



std::vector<double> LinearSystemTools::SeidelMethod(
    std::vector<std::vector<double>> a, std::vector<double> b,
    const double accuracy, const double omega, const bool is_parallel) const {

  using namespace std;

  size_t n(a.size());
  vector<double> x(n, 0);
  // Выбор начального приближения решения
  #pragma omp parallel for if (is_parallel)
  for (ptrdiff_t i(0); i < n; i += 1) {
    x[i] = b[i] / a[i][i];
  }

  double diff(accuracy + 1.0);
  vector<double> transitional_x;
  // Определение числа нитей
  int n_threads(omp_get_num_procs());
  while (diff > accuracy) {
    transitional_x = x;

    // У каждой нити свой вектор с текущим приближением решения
    vector<vector<double>> x_threads(n_threads, x);

    #pragma omp parallel for num_threads(n_threads) if (is_parallel)
    for (ptrdiff_t i(0); i < n; i += 1) {
      int i_thread(omp_get_thread_num());

      double sum(0.0);
      for (ptrdiff_t j(0); j < n; j += 1) {
        if (j != i) {
          sum += a[i][j] * x_threads[i_thread][j];
        }
      }

      // Обновим значение в текущей нити
      x_threads[i_thread][i] = (1 - omega) * x_threads[i_thread][i] + omega * (b[i] - sum) / a[i][i];
      // Обновим значение в глобальном векторе
      x[i] = x_threads[i_thread][i];
    }

    // Вычисление разницы векторов
    #pragma omp parallel for if (is_parallel)
    for (ptrdiff_t i(0); i < n; i += 1) {
      transitional_x[i] -= x[i];
    }
    diff = Norm(transitional_x, is_parallel);
  }

  return x;
}



void LinearSystemTools::Write(const std::vector<double>& vec) const {
  size_t n(vec.size());
  for (ptrdiff_t i(0); i < n; i += 1) {
    std::cout << vec[i] << " ";
  }
  std::cout << std::endl;
}



void LinearSystemTools::Write(const std::vector<std::vector<double>>& matrix) const {
  size_t n(matrix.size());
  for (ptrdiff_t i(0); i < n; i += 1) {
    for (ptrdiff_t j(0); j < n; j += 1) {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << std::endl;
  }
}



uint32_t LinearSystemTools::GetSeed() const {
  return seed_;
}



double LinearSystemTools::Norm(const std::vector<double>& vec, const bool is_parallel) const {
  // Реализация нормы-максимум
  size_t n(vec.size());
  int n_threads(omp_get_num_procs());
  std::vector<double> maxes(n_threads, 0.0);

  #pragma omp parallel for num_threads(n_threads) if (is_parallel)
  for (ptrdiff_t i(0); i < n; i += 1) {
    double abs_v(abs(vec[i]));
    int i_thread(omp_get_thread_num());
    if (abs_v > maxes[i_thread]) {
      maxes[i_thread] = abs_v;
    }
  }

  double max(maxes[0]);
  for (ptrdiff_t i(1); i < n_threads; i += 1) {
    if (maxes[i] > max) {
      max = maxes[i];
    }
  }

  return max;
}



double LinearSystemTools::Norm(const std::vector<std::vector<double>>& matrix, const bool is_parallel) const {
  size_t n(matrix.size());
  int n_threads(omp_get_num_procs());
  std::vector<double> maxes(n_threads, 0.0);

  #pragma omp parallel for num_threads(n_threads) if (is_parallel)
  for (ptrdiff_t i(0); i < n; i += 1) {
    double sum(0.0);
    for (ptrdiff_t j(0); j < n; j += 1) {
      sum += abs(matrix[i][j]);
    }
    int i_thread(omp_get_thread_num());
    if (sum > maxes[i_thread]) {
      maxes[i_thread] = sum;
    }
  }

  double max(maxes[0]);
  for (ptrdiff_t i(1); i < n_threads; i += 1) {
    if (maxes[i] > max) {
      max = maxes[i];
    }
  }

  return max;
}
