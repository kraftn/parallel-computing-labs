#include "timer.h"


void Timer::Start() {
  // Сохранить текущий момент времени
  beginning_time_ = std::chrono::high_resolution_clock::now();
}



double Timer::Finish() const {
  using namespace std::chrono;
  return duration_cast<duration<double, std::ratio<1>>>(high_resolution_clock::now() - beginning_time_).count();
}
