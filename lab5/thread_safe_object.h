#ifndef THREAD_SAFE_OBJECT_H
#define THREAD_SAFE_OBJECT_H

#include <vector>
#include <mutex>

// Базовый класс потокобезопасного объекта
class ThreadSafeObject {
public:
  virtual ~ThreadSafeObject() = default;
  // Найти значение среднего времени использования ресурса
  double GetAverageCaptureTime();
  // Найти значение среднего времени ожидания ресурса
  double GetAverageWaitingTime();
  // Найти максимальное время использования ресурса
  double GetMaxCaptureTime();
  // Найти максимальное время ожидания
  double GetMaxWaitingTime();
  // Найти минимальное время использования ресурса
  double GetMinCaptureTime();
  // Найти минимальное время ожидания
  double GetMinWaitingTime();
  // Получить общее время использования ресурса
  double GetCaptureTime();

protected:
  // Время использования разделяемого ресурса
  std::vector<double> capture_time_;
  // Время ожидания ресурса
  std::vector<double> waiting_time_;
  // Мьютекс для обеспечения потоковой безопасности
  std::mutex mutex_;
};

#endif
