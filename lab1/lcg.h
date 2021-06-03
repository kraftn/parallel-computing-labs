#ifndef LCG_H
#define LCG_H

#include <cstdint>
#include <cmath>
#include <map>

// Линейный конгруэнтный генератор HoaglinLCG(2^31 - 1, a=397204094, 0, 1)
class LCG {
 public:
  // Конструктор
  LCG();
  // Конструктор для создания генератора с заданным начальным значением
  LCG(const int64_t initial_value);
  // Деструктор
  ~LCG();
  // Метод для получения очередного действительного случайного числа,
  // принадлежащего [0; 1)
  double Random();
  // Сброс состояния генератора к начальному
  void Reset();
  // Процедура перехода
  int64_t TransitionProcedure(const int64_t k) const;
  // Проверка, достигается ли максимальный период
  bool CheckTheoremConditions() const;
  // Посчитать период генератора
  int64_t GetPeriod() const;

 private:
  // Возведение числа x в степень y по модулю z
  int64_t PowMod(const int64_t x, const int64_t y, const int64_t z) const;
  // Поиск НОД(a, b)
  int64_t GCD(const int64_t a, const int64_t b) const;
  // Проверка, является ли a первообразным элементом по модулю m
  bool PrimitiveElement(const int64_t a, const int64_t m) const;
  // Факторизация числа n
  std::map<int64_t, int64_t> Factorization(const int64_t n) const;

  // Параметры генератора псевдослучайных чисел HoaglinLCG
  const int64_t m_{ static_cast<int64_t>(pow(2, 31)) - 1 }, a_{ 397204094 },
        seed_{ 1 };
  // Текущее случайное значение
  int64_t current_x_{ 0 };
};

#endif
