#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <cstddef>

#include "thread_safe_object.h"

// Класс потокобезопасного вектора
class Vector : public ThreadSafeObject {
public:
  Vector() = default;
  explicit Vector(const size_t size);
  Vector(const Vector& obj);
  Vector& operator=(const Vector& obj);
  ~Vector() = default;
  // Получить количество элементов
  size_t Size();
  // Получить вместимость
  size_t Capacity();
  // Оператор индексации
  int& operator[](const ptrdiff_t number);
  // Изменить размер
  void Resize(size_t size);

private:
  std::vector<int> data_;
};

#endif
