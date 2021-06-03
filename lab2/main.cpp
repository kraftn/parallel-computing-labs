#include <iostream>
#include <omp.h>

#include "vector.h"
#include "matrix.h"

int main() {
  using namespace std;

  // Протестируем основные операции над векторами и матрицами
  Matrix A_test(3, true), B_test(3, true);
  Vector x_test(3, true), y_test(3, true);
  for (ptrdiff_t i(0); i < A_test.GetN(); i += 1) {
    for (ptrdiff_t j(0); j < A_test.GetN(); j += 1) {
      A_test.Value(i, j) = i % A_test.GetN();
    }
  }
  for (ptrdiff_t i(0); i < B_test.GetN(); i += 1) {
    for (ptrdiff_t j(0); j < B_test.GetN(); j += 1) {
      B_test.Value(i, j) = j % B_test.GetN();
    }
  }
  for (ptrdiff_t i(0); i < x_test.GetN(); i += 1) {
    x_test[i] = i;
  }
  for (ptrdiff_t i(0); i < y_test.GetN(); i += 1) {
    y_test[i] = 2 * i;
  }
  cout << "x:\n" << x_test << "y:\n" << y_test << endl;
  cout << "x + y:\n" << x_test + y_test << endl;
  cout << "x - y:\n" << x_test - y_test << endl;
  cout << "2 * x:\n" << 2 * x_test << endl;
  cout << "||x|| = " << x_test.Length() << endl << endl;
  cout << "<x, y> = " << x_test.DotProduct(y_test) << endl << endl;

  cout << "A:\n" << A_test << "B:\n" << B_test << endl;
  A_test.Transpose();
  cout << "A after transposition:\n" << A_test << endl;
  cout << "A + B:\n" << A_test + B_test << endl;
  cout << "A - B:\n" << A_test - B_test << endl;
  cout << "2 * A:\n" << 2 * A_test << endl;
  cout << "A * B:\n" << A_test * B_test << endl;
  cout << "Ax:\n" << A_test * x_test << endl;
  cout << "||A|| = " << A_test.FrobeniusNorm() << endl << endl << endl;

  // Найдём коэффициенты ускорения
  double start(0.0), end(0.0), time1(0.0), time2(0.0);
  size_t length(1000000);
  Vector x(length, true), y(length, false);
  size_t n_experiments(1000);

  // Сравним последовательные операции над векторами и их параллельные аналоги
  cout << "Vector" << endl;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    x += y;
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    y += x;
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "x += y: " << time2 / time1 << endl;

  time1 = 0.0;
  time2 = 0.0;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    x -= y;
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    y -= x;
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "x -= y: " << time2 / time1 << endl;

  time1 = 0.0;
  time2 = 0.0;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    x *= 2;
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    y *= 2;
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "x *= lambda: " << time2 / time1 << endl;

  time1 = 0.0;
  time2 = 0.0;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    x.Length();
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    y.Length();
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "x.Length(): " << time2 / time1 << endl;

  time1 = 0.0;
  time2 = 0.0;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    x.DotProduct(y);
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    y.DotProduct(x);
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "x.DotProduct(y): " << time2 / time1 << endl << endl << endl;

  // Сравним последовательные операции над матрицами и их параллельные аналоги
  n_experiments = 10000;
  length = 100;
  Matrix A(length, true), B(length, false);
  Vector z(length, false);
  time1 = 0.0;
  time2 = 0.0;

  cout << "Matrix" << endl;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    A.Transpose();
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    B.Transpose();
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "A.Transpose(): " << time2 / time1 << endl;

  time1 = 0.0;
  time2 = 0.0;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    A += B;
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    B += A;
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "A += B: " << time2 / time1 << endl;

  time1 = 0.0;
  time2 = 0.0;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    A -= B;
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    B -= A;
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "A -= B: " << time2 / time1 << endl;

  time1 = 0.0;
  time2 = 0.0;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    A *= 2;
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    B *= 2;
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "A *= lambda: " << time2 / time1 << endl;

  time1 = 0.0;
  time2 = 0.0;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    A *= B;
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    B *= A;
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "A *= B: " << time2 / time1 << endl;

  time1 = 0.0;
  time2 = 0.0;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    A * z;
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    B * z;
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "A * x: " << time2 / time1 << endl;

  time1 = 0.0;
  time2 = 0.0;
  for (ptrdiff_t i_experiment(0); i_experiment < n_experiments;
       i_experiment += 1) {
    start = omp_get_wtime();
    A.FrobeniusNorm();
    end = omp_get_wtime();
    time1 += end - start;

    start = omp_get_wtime();
    B.FrobeniusNorm();
    end = omp_get_wtime();
    time2 += end - start;
  }
  cout << "A.FrobeniusNorm(): " << time2 / time1 << endl;

  char stop(' ');
  cin >> stop;
  return 0;
}
