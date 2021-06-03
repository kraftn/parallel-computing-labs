#ifndef LINEAR_SYSTEM_TOOLS_H
#define LINEAR_SYSTEM_TOOLS_H

#include <vector>
#include <random>

class LinearSystemTools {
 public:
  // Метод для генерации матрицы с коэффициентами СЛАУ
  std::vector<std::vector<double>> GenerateA(const size_t n);
  // Метод для генерации вектора с коэффициентами правой части СЛАУ
  std::vector<double> GenerateB(const double n);
  // Прямой (классический) метод Гаусса
  std::vector<double> GaussianElimination(std::vector<std::vector<double>> a,
                                          std::vector<double> b, const bool is_parallel) const;
  // Итерационный метод Якоби-Гаусса
  std::vector<double> JacobiAlgorithm(std::vector<std::vector<double>> a,
                                      std::vector<double> b, const double accuracy, const bool is_parallel) const;
  // Итерационный метод Зейделя с релаксацией
  std::vector<double> SeidelMethod(std::vector<std::vector<double>> a,
                                   std::vector<double> b, const double accuracy, const double omega,
                                   const bool is_parallel) const;
  // Вывод вектора на экран
  void Write(const std::vector<double>& veс) const;
  // Вывод матрицы на экран
  void Write(const std::vector<std::vector<double>>& matrix) const;
  // Получить начальное значение генератора псевдослучайных чисел
  uint32_t GetSeed() const;

 private:
  // Генератор начального значения
  std::random_device rd_;
  // Начальное значение генератора псевдослучайных чисел
  int64_t seed_{ rd_() };
  // Генератор псевдослучайных чисел
  std::default_random_engine gen_{ seed_ };
  // Вычисление нормы вектора
  double Norm(const std::vector<double>& vec, const bool is_parallel) const;
  // Вычисление нормы матрицы
  double Norm(const std::vector<std::vector<double>>& matrix, const bool is_parallel) const;
};


#endif
