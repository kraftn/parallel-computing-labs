#ifndef TIMER_H
#define TIMER_H

#include <chrono>

// Класс для измерения времени выполнения участка кода
class Timer {
public:
  // Начать отсчёт времени
  void Start();
  // Найти время выполнения
  double Finish() const;

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> beginning_time_;
};


#endif
