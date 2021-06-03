#include <cstdint>
#include <cstddef>
#include <cmath>
#include <omp.h>

#include "ode.h"

using namespace std;

vector<vector<double>> Solve(
    const vector<function<double(const double, const vector<double>&)>>& f,
    const vector<double>& alpha,
    const double t0,
    const double tm,
    const double step,
    const bool is_parallel) {

  // Количество потоков
  const int32_t n_threads(is_parallel ? omp_get_num_procs() : 1);

  // Значения переменной t на отрезке [t0; tm] c шагом step
  vector<double> t;
  double curr_t(t0);
  while (curr_t < tm) {
    t.push_back(curr_t);
    curr_t += step;
  }
  t.push_back(tm);

  // Количество точек
  const int32_t total_n_points(t.size());
  // Количество точек в каждой нити, кроме последней
  const int32_t n_points(total_n_points / n_threads);
  // Размерность решаемой системы
  const int32_t n_u(f.size());

  // Значения определённых интегралов
  vector<vector<double>> integrals(total_n_points - 1);
  // Значения искомых функций (u[i][j] - значение функции u[j] в точке t[i])
  vector<vector<double>> u(total_n_points, alpha);
  // Матрица для хоанения истории значений u
  vector<vector<double>> prev_u(total_n_points);

  #pragma omp parallel for if (is_parallel)
  for (ptrdiff_t i_point(0); i_point < total_n_points - 1; i_point += 1) {
    integrals[i_point].resize(n_u, 0.0);
  }

  // Если разница между текущим и прошлым значением вектора u[i] меньше epsilon,
  // то прекращаем вычисления
  const double epsilon(0.000001);
  // Условие продолжения выполнения вычислений
  bool condition(true);

  // Значения сумм интегралов в каждой нити для каждой искомой функции u
  vector<vector<double>> sum(n_threads);
  #pragma omp parallel for if (is_parallel)
  for (ptrdiff_t i_thread(0); i_thread < n_threads; i_thread += 1) {
    sum[i_thread].resize(n_u, 0.0);
  }

  #pragma omp parallel num_threads(n_threads)
  {
    while (condition) {
      // Сохранить найденные значения u
      #pragma omp single nowait
      {
        prev_u = u;
      }

      // Номер нити
      int32_t i_thread(omp_get_thread_num());
      // Номер первого значения перменной t в нити
      int32_t start(i_thread * n_points);
      // Номер последнего значения перменной t в нити
      int32_t end(start + n_points - 1);
      // В последней нити нужно учесть total_n_points % n_threads
      if (i_thread == (n_threads - 1)) {
        end = t.size() - 2;
      }
      // Обнулим значения сумм
      for (ptrdiff_t i_u(0); i_u < n_u; i_u += 1) {
        sum[i_thread][i_u] = 0.0;
      }

      // Вычисление значений интегралов 
      // (Старченко А.В., Берцун В.Н. "Методы параллельных вычислений", с.170)
      for (ptrdiff_t m(start); m <= end; m += 1) {
        double delta_t(t[m + 1] - t[m]);
        for (ptrdiff_t i_u(0); i_u < n_u; i_u += 1) {
          integrals[m][i_u] = delta_t * (f[i_u](t[m], u[m]) + f[i_u](t[m + 1], u[m + 1])) / 2;
          sum[i_thread][i_u] += integrals[m][i_u];
        }
      }

      // Изменение номера последнего значения t
      if (i_thread < (n_threads - 1)) {
        end = start + n_points - 2;
      }

      // Ожидание завершения вычисления значений интегралов во всех нитях
      #pragma omp barrier

      // Вычисление начальных значений u[start]
      for (ptrdiff_t i_u(0); i_u < n_u; i_u += 1) {
        u[start][i_u] = u[0][i_u];
        for (ptrdiff_t i_sum(0); i_sum < i_thread; i_sum += 1) {
          u[start][i_u] += sum[i_sum][i_u];
        }
      }
      // Обновление значений искомых функций при всех значениях аргумента t 
      for (ptrdiff_t m(start); m <= end; m += 1) {
        for (ptrdiff_t i_u(0); i_u < n_u; i_u += 1) {
          u[m + 1][i_u] = u[m][i_u] + integrals[m][i_u];
        }
      }

      #pragma omp barrier

      // Одна из нитей обновляет условие выхода из цикла
      #pragma omp single
      {
        ptrdiff_t m(0);
        double diff(EuclideanMetric(u[m], prev_u[m], is_parallel));
        while (((diff < epsilon) && (m + 1) < total_n_points)) {
          m += 1;
          diff = EuclideanMetric(u[m], prev_u[m], is_parallel);
        }
        condition = diff >= epsilon;
      }

      // Ожидание завершения вычислений во всех нитях
      #pragma omp barrier
    }
  }

  return u;
}



double EuclideanMetric(const vector<double>& x, const vector<double>& y, const bool is_parallel) {
  double res(0.0);
  size_t n(x.size());
  // Найдём в параллельном или последовательном режиме значение суммы
  #pragma omp parallel for reduction(+:res) if (is_parallel)
  for (ptrdiff_t i(0); i < n; i += 1) {
    res += pow(x[i] - y[i], 2);
  }
  return sqrt(res);
}



vector<double> Calculate(
    const function<double(const double)>& u,
    const double t0,
    const double tm,
    const double step) {
  vector<double> ans;
  double curr_t(t0);
  while (curr_t < tm) {
    ans.push_back(u(curr_t));
    curr_t += step;
  }
  ans.push_back(u(tm));
  return ans;
}



vector<vector<double>> Transpose(const vector<vector<double>>& matrix) {
  int32_t n(matrix.size());
  int32_t m(matrix[0].size());
  vector<vector<double>> res(m);
  for (ptrdiff_t i(0); i < m; i += 1) {
    res[i].resize(n, 0.0);
    for (ptrdiff_t j(0); j < n; j += 1) {
      res[i][j] = matrix[j][i];
    }
  }
  return res;
}



double MaximumMetric(const vector<double>& x, const vector<double>& y) {
  double max(abs(x[0] - y[0]));
  double v(0.0);
  size_t n(x.size());

  for (ptrdiff_t i(1); i < n; i += 1) {
    v = abs(x[i] - y[i]);
    if (v > max) {
      max = v;
    }
  }

  return max;
}
