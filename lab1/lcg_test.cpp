#include "lcg.h"
#include "evaluator.h"

#include <iostream>
#include <cstdint>
#include <cmath>
#include <omp.h>


int main() {
  using namespace std;

  Evaluator evaluator;
  int64_t d(10);
  int64_t t(2);
  // Длина тестируемой последовательности
  int64_t n(evaluator.GetPeriod());
  // Количество физических процессоров
  size_t num_procs(static_cast<size_t>(omp_get_num_procs()));
  // Количество повторений для каждого теста
  int64_t n_experiments(3);

  double start(0.0), end(0.0), time1(0.0), time2(0.0), res1(0.0), res2(0.0);
  // Массивы для хранения времени выполнения тестов
  vector<double> uniformity, uniformity_p, series, series_p, interval, interval_p,
    poker, poker_p;

  cout << "Number of threads: " << num_procs << endl;
  cout << "Test\t\t\tSequental (s)\tParallel (s)\tCoefficient" << endl;

  for (int64_t i_experiment(0); i_experiment < n_experiments; i_experiment += 1) {
    // Критерий частот
    start = omp_get_wtime();
    res1 = evaluator.UniformityTest(n, d);
    end = omp_get_wtime();
    time1 = end - start;
    uniformity.push_back(time1);

    // Параллельный вариант критерия частот
    start = omp_get_wtime();
    res2 = evaluator.UniformityTest(n, d, num_procs);
    end = omp_get_wtime();
    time2 = end - start;
    uniformity_p.push_back(time2);
    // Если результаты последовательного и параллельного вариантов не совпадают,
    // запустить ошибку
    if (abs(res1 - res2) > 0.000001) {
      throw logic_error("Error in LCG::UniformityTest");
    }
    cout << "UniformityTest\t\t" << time1 << "\t\t" << time2 << endl;

    // Критерий серий
    start = omp_get_wtime();
    res1 = evaluator.SeriesTest(n / 2, d);
    end = omp_get_wtime();
    time1 = end - start;
    series.push_back(time1);

    // Параллельный вариант критерия серий
    start = omp_get_wtime();
    res2 = evaluator.SeriesTest(n / 2, d, num_procs);
    end = omp_get_wtime();
    time2 = end - start;
    series_p.push_back(time2);
    if (abs(res1 - res2) > 0.000001) {
      throw logic_error("Error in LCG::SeriesTest");
    }
    cout << "SeriesTest\t\t" << time1 << "\t\t" << time2 << endl;

    // Критерий интервалов
    start = omp_get_wtime();
    res1 = evaluator.IntervalTest(n, d, t);
    end = omp_get_wtime();
    time1 = end - start;
    interval.push_back(time1);

    // Параллельный вариант критерия интервалов
    start = omp_get_wtime();
    res2 = evaluator.IntervalTest(n, d, t, num_procs);
    end = omp_get_wtime();
    time2 = end - start;
    interval_p.push_back(time2);
    if (abs(res1 - res2) > 0.000001) {
      throw logic_error("Error in LCG::IntervalTest");
    }
    cout << "IntervalTest\t\t" << time1 << "\t\t" << time2 << endl;

    // Покер-критерий
    start = omp_get_wtime();
    res1 = evaluator.PokerTest(n / 5, d);
    end = omp_get_wtime();
    time1 = end - start;
    poker.push_back(time1);

    // Параллельный вариант покер-критерия
    start = omp_get_wtime();
    res2 = evaluator.PokerTest(n / 5, d, num_procs);
    end = omp_get_wtime();
    time2 = end - start;
    poker_p.push_back(time2);
    if (abs(res1 - res2) > 0.000001) {
      throw logic_error("Error in LCG::PokerTest");
    }
    cout << "PokerTest\t\t" << time1 << "\t\t" << time2 << endl << endl;
  }

  cout << endl << "Mean:" << endl;

  // Найдём среднее время выполнения и коэффициенты ускорения для каждого теста
  double seq(evaluator.Mean(uniformity)), par(evaluator.Mean(uniformity_p));
  cout << "UniformityTest\t\t" << seq << "\t\t" << par << "\t\t" << (seq / par) << endl;

  seq = evaluator.Mean(series);
  par = evaluator.Mean(series_p);
  cout << "SeriesTest\t\t" << seq << "\t\t" << par << "\t\t" << (seq / par) << endl;

  seq = evaluator.Mean(interval);
  par = evaluator.Mean(interval_p);
  cout << "IntervalTest\t\t" << seq << "\t\t" << par << "\t\t" << (seq / par) << endl;

  seq = evaluator.Mean(poker);
  par = evaluator.Mean(poker_p);
  cout << "PokerTest\t\t" << seq << "\t\t" << par << "\t\t" << (seq / par) << endl;

  char stop(' ');
  cin >> stop;
  return 0;
}
