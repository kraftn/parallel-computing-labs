#include <iostream>
#include <omp.h>

#include "linear_system_tools.h"

int main() {
  using namespace std;

  LinearSystemTools tools;
  cout << "Seed: " << tools.GetSeed() << endl << endl;

  size_t n_calc(3);
  vector<vector<double>> a_calc(tools.GenerateA(n_calc));
  vector<double> b_calc(tools.GenerateB(n_calc));
  double accuracy_calc(0.000001);
  double omega_calc(1.5);

  // Демонстрация работы последовательных и параллельных алгоритмов
  cout << "A:" << endl;
  tools.Write(a_calc);
  cout << "b:" << endl;
  tools.Write(b_calc);
  cout << endl;

  cout << "Serial Gaussian elimination:" << endl;
  tools.Write(tools.GaussianElimination(a_calc, b_calc, false));
  cout << "Parallel Gaussian elimination:" << endl;
  tools.Write(tools.GaussianElimination(a_calc, b_calc, true));
  cout << endl;

  cout << "Serial Jacobi algorithm: " << endl;
  tools.Write(tools.JacobiAlgorithm(a_calc, b_calc, accuracy_calc, false));
  cout << "Parallel Jacobi algorithm: " << endl;
  tools.Write(tools.JacobiAlgorithm(a_calc, b_calc, accuracy_calc, true));
  cout << endl;

  cout << "Serial Seidel method: " << endl;
  tools.Write(tools.SeidelMethod(a_calc, b_calc, accuracy_calc, omega_calc, false));
  cout << "Parallel Seidel method: " << endl;
  tools.Write(tools.SeidelMethod(a_calc, b_calc, accuracy_calc, omega_calc, true));
  cout << endl;

  // Вычисление коэффициентов ускорения
  size_t n(50);
  vector<vector<double>> a;
  vector<double> b;
  double time_serial[3]{ 0.0 }, time_parallel[3]{ 0.0 };
  
  double start(0.0);
  double accuracy(0.000001);
  double omega(1.5);
  for (ptrdiff_t i_experiment(0); i_experiment < 1000; i_experiment += 1) {
    a = tools.GenerateA(n);
    b = tools.GenerateB(n);

    start = omp_get_wtime();
    tools.GaussianElimination(a, b, false);
    time_serial[0] += omp_get_wtime() - start;

    start = omp_get_wtime();
    tools.GaussianElimination(a, b, true);
    time_parallel[0] += omp_get_wtime() - start;

    start = omp_get_wtime();
    tools.JacobiAlgorithm(a, b, accuracy, false);
    time_serial[1] += omp_get_wtime() - start;

    start = omp_get_wtime();
    tools.JacobiAlgorithm(a, b, accuracy, true);
    time_parallel[1] += omp_get_wtime() - start;

    start = omp_get_wtime();
    tools.SeidelMethod(a, b, accuracy, omega, false);
    time_serial[2] += omp_get_wtime() - start;

    start = omp_get_wtime();
    tools.SeidelMethod(a, b, accuracy, omega, true);
    time_parallel[2] += omp_get_wtime() - start;
  }
  cout << "Gaussian elimination: " << time_serial[0] / time_parallel[0] << endl;
  cout << "Jacobi algorithm: " << time_serial[1] / time_parallel[1] << endl;
  cout << "Seidel method: " << time_serial[2] / time_parallel[2] << endl;

  char stop(0);
  cin >> stop;
  return 0;
}
