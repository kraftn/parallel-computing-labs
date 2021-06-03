#ifndef EVALUATOR_H
#define EVALUATOR_H

#include "lcg.h"

#include <cstdint>
#include <vector>

// Класс для вычисления значений статистических тестов
class Evaluator {
 public:
  Evaluator();
  ~Evaluator();
  // Критерий частот
  double UniformityTest(const int64_t n, const int64_t d);
  // Параллельный вариант критерия частот
  double UniformityTest(const int64_t n, const int64_t d, const size_t n_threads) const;
  // Критерий серий
  double SeriesTest(const int64_t n, const int64_t d);
  // Параллельный вариант критерия серий
  double SeriesTest(const int64_t n, const int64_t d, const size_t n_threads) const;
  // Критерий интервалов
  double IntervalTest(const int64_t n, const int64_t d, const int64_t t);
  // Параллельный вариант критерия интервалов
  double IntervalTest(const int64_t n, const int64_t d, const int64_t t, const size_t n_threads) const;
  // Покер-критерий
  double PokerTest(const int64_t n, const int64_t d);
  // Параллельный вариант покер-критерия
  double PokerTest(const int64_t n, const int64_t d, const size_t n_threads) const;
  // Приближённо определить, квантилью какого уровня распределения хи-квадрат
  // является полученное значение теста (для n_freedom степеней свободы)
  double GetLevel(const double result_test, const int64_t n_freedom) const;
  // Определить, является ли допустимым полученное значение статистического теста
  bool InterpretTestResult(const double result_test,
                           const int64_t n_freedom) const;
  // Метод для вычисления мат. ожидания
  double Mean(const std::vector<double>& arr) const;
  // Проверка, достигается ли максимальный период
  bool CheckTheoremConditions() const;
  // Посчитать период генератора
  int64_t GetPeriod() const;

 private:
  // Вычисление значения хи-квадрат по известным эмпирическим наблюдениям v,
  // количеству измерений n и вероятностям категорий p
  double ChiSquaredTest(const std::vector<int64_t>& v,
                        const std::vector<double>& p, const int64_t n) const;
  // Реализация аппроксимации Голдштейна распределения хи-квадрат
  double ApproximationGoldshtein(const double level,
                                 const int64_t n_freedom) const;
  // Вычисление факториала n!
  int64_t Factorial(const int64_t n) const;
  // Вычисление произведения всех целых чисел, принадлежащих отрезку [a, b] 
  int64_t Multiplication(const int64_t a, const int64_t b) const;
  // Реализация явной формулы для вычисления чисел Стирлинга второго рода
  // (количество неупорядоченных разбиений n-элементного множества на k
  // непустых подмножеств)
  int64_t StirlingNumber(const int64_t n, const int64_t k) const;
  // Метод, возвращающий генераторы, настроенные на работу в n_threads нитях
  std::vector<LCG*> GetVectorOfGenerators(const size_t n_threads, const int64_t shift) const;
  // Метод, возвращающий вектор, являющийся результатом сложения векторов в res_arrays
  std::vector<int64_t> Merge(const std::vector<std::vector<int64_t>>& res_arrays) const;

  // Генератор, тестируемый в последовательном режиме
  LCG* generator_{ nullptr };
};

#endif
