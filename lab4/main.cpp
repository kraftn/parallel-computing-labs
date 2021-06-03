#include <iostream>
#include <omp.h>

#include "ode.h"


// Функции ui для тестового примера

double u1(const double t) {
  return t * t + 2 * t + 1;
}



double u2(const double t) {
  return t + 1;
}



double u3(const double t) {
  return log(t) + t * t;
}



double u4(const double t) {
  return sin(t);
}



double u5(const double t) {
  return cos(t);
}

// Функции fi для тестового примера

double f1(const double t, const std::vector<double>& u) {
  return 2 * u[1];
}



double f2(const double t, const std::vector<double>& u) {
  return u[3] * u[3] + u[4] * u[4];
}



double f3(const double t, const std::vector<double>& u) {
  return 1 / (u[1] - 1) + 2 * t;
}



double f4(const double t, const std::vector<double>& u) {
  return u[4];
}



double f5(const double t, const std::vector<double>& u) {
  return -u[3];
}



int main() {
  using namespace std;

  vector<function<double(const double, const vector<double>&)>> f{ f1, f2, f3, f4, f5 };
  double t0(0.1);
  double tm(10.1);
  vector<double> alpha{ u1(t0), u2(t0), u3(t0), u4(t0), u5(t0) };
  double step(0.001);

  // Найдём численное решение системы ОДУ в последовательном и параллельном режиме
  vector<vector<double>> solved1(Solve(f, alpha, t0, tm, step, false));
  vector<vector<double>> transposed_solved1(Transpose(solved1));
  vector<vector<double>> solved2(Solve(f, alpha, t0, tm, step, true));
  vector<vector<double>> transposed_solved2(Transpose(solved2));

  vector<function<double(const double)>> u{ u1, u2, u3, u4, u5 };
  vector<double> calculated;
  for (ptrdiff_t i_func(0); i_func < 5; i_func += 1) {
    calculated = Calculate(u[i_func], t0, tm, step);
    // Максимальное значение модуля разницы между точными значениями 
    // и значениями численного решения при последовательном подходе
    cout << "Consequential\nu" << i_func + 1 << ": " << MaximumMetric(calculated, transposed_solved1[i_func]) << endl;
    // Максимальное значение модуля разницы между точными значениями 
    // и значениями численного решения при параллельном подходе 
    cout << "Parallel\nu" << i_func + 1 << ": " << MaximumMetric(calculated, transposed_solved2[i_func]) << endl << endl;
  }
  cout << endl;

  // Найдём значение коэффициента ускорения
  double time_sequental(0.0);
  double time_parallel(0.0);
  double start(0.0);
  for (ptrdiff_t i(0); i < 5; i += 1) {
    start = omp_get_wtime();
    Solve(f, alpha, t0, tm, step, false);
    time_sequental += omp_get_wtime() - start;

    start = omp_get_wtime();
    Solve(f, alpha, t0, tm, step, true);
    time_parallel += omp_get_wtime() - start;
  }
  cout << "Coefficient: " << time_sequental / time_parallel << endl;

  char stop('0');
  cin >> stop;
  return 0;
}
