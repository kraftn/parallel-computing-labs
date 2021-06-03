#include "evaluator.h"

#include <cmath>
#include <cstddef>

Evaluator::Evaluator()
  : generator_(new LCG()) {
}



Evaluator::~Evaluator() {
  delete generator_;
  generator_ = nullptr;
}



double Evaluator::UniformityTest(const int64_t n, const int64_t d) {
  // Перед каждым тестом сбросить состояние генератора к начальному
  generator_->Reset();

  // Вектор для хранения количества наблюдений для каждой категории (0, ..., d - 1)
  std::vector<int64_t> counter(d, 0);
  for (int64_t i(0); i < n; i += 1) {
    counter[static_cast<int64_t>(d * generator_->Random())] += 1;
  }

  // Вероятность каждой категории равна 1/d
  std::vector<double> p(d, 1.0 / d);

  return ChiSquaredTest(counter, p, n);
}



double Evaluator::UniformityTest(const int64_t n, const int64_t d, const size_t n_threads) const {
  // Количество элементов в каждой нити, кроме последней
  int64_t step(n / n_threads);
  // Создадим по 1 генератору на каждую нить
  std::vector<LCG*> generators(GetVectorOfGenerators(n_threads, step));

  // Каждая нить имеет свой вектор для хранения количества наблюдений для
  // каждой категории (0, ..., d - 1)
  std::vector<std::vector<int64_t>> counters(n_threads);

  #pragma omp parallel for
  for (ptrdiff_t i_thread(0); i_thread < n_threads; i_thread += 1) {
    counters[i_thread].resize(d, 0);

    // К количеству элементов последней нити добавим n % n_threads 
    int64_t n_stop(step);
    if (i_thread == (n_threads - 1)) {
      n_stop += (n % n_threads);
    }

    for (int64_t i_num(0); i_num < n_stop; i_num += 1) {
      counters[i_thread][static_cast<int64_t>(d * generators[i_thread]->Random())] += 1;
    }

    // Освободим память
    delete generators[i_thread];
    generators[i_thread] = nullptr;
  }

  // Соберём измерения из разных нитей в один массив
  std::vector<int64_t> counter(Merge(counters));
  // Вероятность каждой категории равна 1/d
  std::vector<double> p(d, 1.0 / d);

  return ChiSquaredTest(counter, p, n);
}



double Evaluator::SeriesTest(const int64_t n, const int64_t d) {
  generator_->Reset();

  // Вектор для хранения количества наблюдений для каждой категории
  std::vector<int64_t> counter(d * d, 0);
  for (int64_t i(0); i < n; i += 1) {
    int64_t first_v(static_cast<int64_t>(d * generator_->Random()));
    int64_t second_v(static_cast<int64_t>(d * generator_->Random()));

    // Рассмотрим пару чисел first_v, second_v как запись числа в системе счисления
    // с основанием d
    counter[d * first_v + second_v] += 1;
  }

  // Все категории равновероятны (1/(d*d))
  std::vector<double> p(d * d, 1.0 / (d * d));

  return ChiSquaredTest(counter, p, n);
}



double Evaluator::SeriesTest(const int64_t n, const int64_t d, const size_t n_threads) const {
  // Количество пар в каждой нити, кроме последней
  int64_t step(n / n_threads);
  // Создадим по 1 генератору на каждую нить
  std::vector<LCG*> generators(GetVectorOfGenerators(n_threads, 2 * step));

  // Каждая нить имеет свой вектор для хранения количества наблюдений для
  // каждой категории
  std::vector<std::vector<int64_t>> counters(n_threads);

  #pragma omp parallel for
  for (ptrdiff_t i_thread(0); i_thread < n_threads; i_thread += 1) {
    counters[i_thread].resize(d * d, 0);

    // Учтём последние n % n_threads пар
    int64_t n_stop(step);
    if (i_thread == (n_threads - 1)) {
      n_stop += (n % n_threads);
    }

    for (int64_t i_pair(0); i_pair < n_stop; i_pair += 1) {
      int64_t first_v(static_cast<int64_t>(d * generators[i_thread]->Random()));
      int64_t second_v(static_cast<int64_t>(d * generators[i_thread]->Random()));

      // Рассмотрим пару чисел first_v, second_v как запись числа в системе счисления
      // с основанием d
      counters[i_thread][d * first_v + second_v] += 1;
    }

    // Освободим память
    delete generators[i_thread];
    generators[i_thread] = nullptr;
  }
  

  // Соберём измерения из разных нитей в один массив
  std::vector<int64_t> counter(Merge(counters));
  // Все категории равновероятны (1/(d*d))
  std::vector<double> p(d * d, 1.0 / (d * d));

  return ChiSquaredTest(counter, p, n);
}



double Evaluator::IntervalTest(const int64_t n, const int64_t d, const int64_t t) {
  generator_->Reset();

  /*
  Вектор для хранения количества наблюдений для каждой категории
  counter[i] - количество последовательностей длины (i + 1), 0 <= i < t
  counter[t] - количество последовательностей длиной > t
  */
  std::vector<int64_t> counter(t + 1, 0);
  int64_t previous(0), current(static_cast<int64_t>(d * generator_->Random())), length(1);
  // Переменная для подсчёта количества последовательностей
  int64_t total_number(0);

  // Посчитаем количество последовательностей, состоящих из одинаковых элементов
  for (int64_t i_experiment(1); i_experiment < n; i_experiment += 1) {
    previous = current;
    current = static_cast<int64_t>(d * generator_->Random());
    if (current == previous) {
      length += 1;
    }
    else {
      if (length <= t) {
        counter[length - 1] += 1;
      }
      else {
        counter[t] += 1;
      }
      total_number += 1;
      length = 1;
    }
  }
  if (length <= t) {
    counter[length - 1] += 1;
  }
  else {
    counter[t] += 1;
  }
  total_number += 1;

  // Найдём вероятность каждой категории
  std::vector<double> p(t + 1, 0);
  double sum(0.0);
  for (int64_t i_p(0); i_p < p.size() - 1; i_p += 1) {
    int64_t length(i_p + 1);
    /*
    Пусть x - длина последовательности из одинаковых элементов
    P(x>=m) = 1/(d^(m-1)), так как (m-1) элемент должны быть равны первому элементу
    Найдём P(x=m): P(x=m) = P(x>=m) - P(x>=m+1) = 1/(d^(m-1)) - 1/(d^m) =
    = (d-1)/(d^m), 0 < m < n
    Если m=n, то P(x=n) = 1/(d^(n-1))
    */
    p[i_p] = (d - 1) / pow(d, length);
    sum += p[i_p];
  }
  p[t] = 1 - sum;

  return ChiSquaredTest(counter, p, total_number);
}



double Evaluator::IntervalTest(const int64_t n, const int64_t d,
                               const int64_t t, const size_t n_threads) const {
  // Начальное количество элементов в каждой нити
  int64_t step(n / n_threads);
  // Создадим по 1 генератору на каждую нить
  std::vector<LCG*> generators(n_threads);
  // Первый генератор используется для настройки остальных
  generators[0] = new LCG();
  // Итоговое количество элементов в каждой нити
  std::vector<int64_t> lengths(n_threads, step);

  for (ptrdiff_t i_setting(1); i_setting < n_threads; i_setting += 1) {
    int64_t last(generators[0]->TransitionProcedure(i_setting * step - 1));
    int64_t initial_value(generators[0]->TransitionProcedure(i_setting * step));
    int64_t shift(0);
    // Необходимо, чтобы последовательность из одинаковых элементов не разделялась
    // разными нитями
    while (initial_value == last) {
      shift += 1;
      last = initial_value;
      initial_value = generators[0]->TransitionProcedure(i_setting * step + shift);
    }
    // Сдвинем вправо начальное значение настраиваемого генератора
    lengths[i_setting] -= shift;
    lengths[i_setting - 1] += shift;
    generators[i_setting] = new LCG(initial_value);
  }

  /* 
  Вектора для хранения количества наблюдений для каждой категории
  counters[j][i] - количество последовательностей длины (i + 1) в нити j, 0 <= i < t
  counters[j][t] - количество последовательностей длиной > t в нити j
  */
  std::vector<std::vector<int64_t>> counters(n_threads);
  // У каждой нити свой элемент массива для подсчёта количества последовательностей
  std::vector<int64_t> total_numbers(n_threads, 0);

  #pragma omp parallel for
  for (ptrdiff_t i_thread(0); i_thread < n_threads; i_thread += 1) {
    counters[i_thread].resize(t + 1, 0);

    // Учтём последние n % n_threads элементов
    int64_t n_stop(lengths[i_thread]);
    if (i_thread == (n_threads - 1)) {
      n_stop += (n % n_threads);
    }

    int64_t previous(0), current(static_cast<int64_t>(d * generators[i_thread]->Random())), length(1);

    // Посчитаем количество последовательностей, состоящих из одинаковых элементов
    for (int64_t i_experiment(1); i_experiment < n_stop; i_experiment += 1) {
      previous = current;
      current = static_cast<int64_t>(d * generators[i_thread]->Random());
      if (current == previous) {
        length += 1;
      }
      else {
        if (length <= t) {
          counters[i_thread][length - 1] += 1;
        }
        else {
          counters[i_thread][t] += 1;
        }
        total_numbers[i_thread] += 1;
        length = 1;
      }
    }
    if (length <= t) {
      counters[i_thread][length - 1] += 1;
    }
    else {
      counters[i_thread][t] += 1;
    }
    total_numbers[i_thread] += 1;

    // Освободим память
    delete generators[i_thread];
    generators[i_thread] = nullptr;
  }
  
  std::vector<int64_t> counter(Merge(counters));
  // Найдём общее число последовательностей
  int64_t total_number(0);
  for (ptrdiff_t i_thread(0); i_thread < n_threads; i_thread += 1) {
    total_number += total_numbers[i_thread];
  }

  // Найдём вероятность каждой категории
  std::vector<double> p(t + 1, 0);
  double sum(0.0);
  for (int64_t i_p(0); i_p < p.size() - 1; i_p += 1) {
    int64_t length(i_p + 1);
    /*
    Пусть x - длина последовательности из одинаковых элементов
    P(x>=m) = 1/(d^(m-1)), так как (m-1) элемент должны быть равны первому элементу
    Найдём P(x=m): P(x=m) = P(x>=m) - P(x>=m+1) = 1/(d^(m-1)) - 1/(d^m) =
    = (d-1)/(d^m), 0 < m < n
    Если m=n, то P(x=n) = 1/(d^(n-1))
    */
    p[i_p] = (d - 1) / pow(d, length);
    sum += p[i_p];
  }
  p[t] = 1 - sum;

  return ChiSquaredTest(counter, p, total_number);
}



double Evaluator::PokerTest(const int64_t n, const int64_t d) {
  generator_->Reset();

  // Вектор для хранения количества пятёрок с 1, 2, 3, 4 и 5 различными
  // значениями
  std::vector<int64_t> counter(5, 0);
  int64_t num(0);
  // Количество различных элементов в очередном множестве из пяти элементов
  int64_t total_different(0);

  for (int64_t i_experiment(0); i_experiment < n; i_experiment += 1) {
    total_different = 0;
    // Вектор для хранения того, какие элементы встретились во множестве
    std::vector<bool> mark_num(d, false);

    for (ptrdiff_t i_value(0); i_value < 5; i_value += 1) {
      num = static_cast<int64_t>(d * generator_->Random());
      if (!mark_num[num]) {
        mark_num[num] = true;
        total_different += 1;
      }
    }

    counter[total_different - 1] += 1;
  }

  // Найдём вероятность для каждого вида множеств
  std::vector<double> p(5, 0);
  int64_t n_subsets(0);
  for (ptrdiff_t i_p(0); i_p < p.size(); i_p += 1) {
    n_subsets = i_p + 1;
    // Дональд Э. Кнут "Искусство программирования", том 2 "Получисленные алгоритмы", с.86
    p[i_p] = (StirlingNumber(5, n_subsets) / pow(d, 5)) * Multiplication(d - n_subsets + 1, d);
  }

  return ChiSquaredTest(counter, p, n);
}



double Evaluator::PokerTest(const int64_t n, const int64_t d, const size_t n_threads) const {
  // Количество пятёрок элементов в каждой нити, кроме последней
  int64_t step(n / n_threads);
  // Настроим генераторы
  std::vector<LCG*> generators(GetVectorOfGenerators(n_threads, 5 * step));

  // Массивы для хранения количества пятёрок с 1, 2, 3, 4 и 5 различными
  // значениями
  std::vector<std::vector<int64_t>> counters(n_threads);

  #pragma omp parallel for
  for (ptrdiff_t i_thread(0); i_thread < n_threads; i_thread += 1) {
    counters[i_thread].resize(5, 0);

    // Изменим количество пятёрок в последней нити
    int64_t n_stop(step);
    if (i_thread == (n_threads - 1)) {
      n_stop += (n % n_threads);
    }

    int64_t num(0);
    // Количество различных элементов в очередном множестве из пяти элементов
    int64_t total_different(0);
    for (int64_t i_experiment(0); i_experiment < n_stop; i_experiment += 1) {
      total_different = 0;
      // Вектор для хранения того, какие элементы встретились во множестве
      std::vector<bool> mark_num(d, false);

      for (ptrdiff_t i_value(0); i_value < 5; i_value += 1) {
        num = static_cast<int64_t>(d * generators[i_thread]->Random());
        if (!mark_num[num]) {
          mark_num[num] = true;
          total_different += 1;
        }
      }

      counters[i_thread][total_different - 1] += 1;
    }

    delete generators[i_thread];
    generators[i_thread] = nullptr;
  }

  std::vector<int64_t> counter(Merge(counters));

  // Найдём вероятность для каждого вида множеств
  std::vector<double> p(5, 0);
  int64_t n_subsets(0);
  for (ptrdiff_t i_p(0); i_p < p.size(); i_p += 1) {
    n_subsets = i_p + 1;
    // Дональд Э. Кнут "Искусство программирования", том 2 "Получисленные алгоритмы", с.86
    p[i_p] = (StirlingNumber(5, n_subsets) / pow(d, 5)) * Multiplication(d - n_subsets + 1, d);
  }

  return ChiSquaredTest(counter, p, n);
}


double Evaluator::GetLevel(const double result_test,
                           const int64_t n_freedom) const {
  double step(0.001);

  size_t n_quantiles(static_cast<size_t>(round(1 / step)) - 1);
  std::vector<double> values(n_quantiles, 0);
  // Найдём значения квантилей распределения хи-квадрат во всех точках
  // отрезка [step; 1-step] с шагом step
  for (ptrdiff_t i(0); i < n_quantiles; i += 1) {
    values[i] = ApproximationGoldshtein((i + 1) * step, n_freedom);
  }

  // Определим, к квантили какого уровня наиболее близко значение статистичекого теста
  // result_test
  double min_diff(abs(result_test - values[0]));
  ptrdiff_t i_min(0);
  for (ptrdiff_t i(1); i < n_quantiles; i += 1) {
    double diff(abs(values[i] - result_test));
    if (diff <= min_diff) {
      min_diff = diff;
      i_min = i;
    }
  }

  return (i_min + 1) * step;
}



bool Evaluator::InterpretTestResult(const double result_test,
                                    const int64_t n_freedom) const {
  // Тест пройден, если его результат меньше либо равен 0,1-квантили
  // распределения хи-квадрат
  return result_test <= ApproximationGoldshtein(0.1, n_freedom);
}



double Evaluator::Mean(const std::vector<double>& arr) const {
  double sum(0.0);

  for (ptrdiff_t i(0); i < arr.size(); i += 1) {
    sum += arr[i];
  }

  return sum / arr.size();
}



double Evaluator::ChiSquaredTest(const std::vector<int64_t>& v,
                                 const std::vector<double>& p, const int64_t n) const {
  double ans(0.0);

  for (int64_t i(0); i < v.size(); i += 1) {
    ans += pow(v[i] - n * p[i], 2) / (n * p[i]);
  }

  return ans;
}



double Evaluator::ApproximationGoldshtein(const double level,
                                          const int64_t n_freedom) const {
  // Кобзарь А.И. "Прикладная математическая статистика. Для инженеров и научных работников.", с. 46
  const double a_array[7] { 1.0000886, 0.4713941, 0.0001348028, -0.008553069, 0.00312558, -0.0008426812, 0.00009780499 };
  const double b_array[7] { -0.2237368, 0.02607083, 0.01128186, -0.01153761, 0.005169654, 0.00253001, -0.001450117 };
  const double c_array[7] { -0.01513904, -0.008986007, 0.02277679, -0.01323293, -0.006950356, 0.001060438, 0.001565326 };
  double d(0.0);
  if (level >= 0.5) {
    d = 2.0637 * pow((log(1 / (1 - level)) - 0.16), 0.4274) - 1.5774;
  } else {
    d = -2.0637 * pow((log(1 / level) - 0.16), 0.4274) + 1.5774;
  }

  double ans(0.0);
  for (ptrdiff_t i(0); i < 7; i += 1) {
    ans += pow(n_freedom, -i / 2.0) * pow(d, i) * (a_array[i] + b_array[i] / n_freedom + c_array[i] / (n_freedom * n_freedom));
  }

  return pow(ans, 3) * n_freedom;
}



int64_t Evaluator::Multiplication(const int64_t a, const int64_t b) const {
  int64_t ans(1);
  for (int64_t i(a); i <= b; i += 1) {
    ans *= i;
  }
  return ans;
}



int64_t Evaluator::Factorial(const int64_t n) const {
  return Multiplication(2, n);
}



int64_t Evaluator::StirlingNumber(const int64_t n, const int64_t k) const {
  int64_t ans(0);

  for (int64_t j(0); j <= k; j += 1) {
    if (((j + k) % 2) == 0) {
      ans += Multiplication(j + 1, k) * static_cast<int64_t>(pow(j, n)) / Factorial(k - j);
    }
    else {
      ans += -Multiplication(j + 1, k) * static_cast<int64_t>(pow(j, n)) / Factorial(k - j);
    }
  }

  return ans / Factorial(k);
}



std::vector<LCG*> Evaluator::GetVectorOfGenerators(const size_t n_threads, const int64_t shift) const {
  std::vector<LCG*> generators(n_threads);
  generators[0] = new LCG();
  for (ptrdiff_t i(1); i < n_threads; i += 1) {
    // По процедуре перехода вычислим начальное значение для генератора
    int64_t initial_value(generators[0]->TransitionProcedure(i * shift));
    generators[i] = new LCG(initial_value);
  }

  return generators;
}



std::vector<int64_t> Evaluator::Merge(const std::vector<std::vector<int64_t>>& res_arrays) const {
  size_t n_threads(res_arrays.size());
  int64_t n_categories(res_arrays[0].size());

  std::vector<int64_t> ans(n_categories, 0);
  // В параллельном режиме найдём сумму массивов res_arrays
  #pragma omp parallel for
  for (int64_t i_category(0); i_category < n_categories; i_category += 1) {
    for (ptrdiff_t i_thread(0); i_thread < n_threads; i_thread += 1) {
      ans[i_category] += res_arrays[i_thread][i_category];
    }
  }

  return ans;
}



bool Evaluator::CheckTheoremConditions() const {
  return generator_->CheckTheoremConditions();
}



int64_t Evaluator::GetPeriod() const {
  return generator_->GetPeriod();
}
