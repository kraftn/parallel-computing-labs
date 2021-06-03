#ifndef QUEUE_H
#define QUEUE_H

#include <cstddef>
#include <queue>

#include "thread_safe_object.h"

// Класс потокобезопасной очереди
class Queue : public ThreadSafeObject {
public:
  Queue() = default;
  Queue(const Queue& obj);
  ~Queue() = default;
  Queue& operator=(const Queue&);
  // Поместить новый элемент
  void Push(const int value);
  // Удалить первый элемент
  void Pop();
  // Получить первый элемент
  int Top();
  // Проверить на пустоту
  bool IsEmpty();
  // Получить количество элементов
  size_t Size();
  // Удалить все элементы
  void Clear();
  
private:
  std::queue<int> data_;
};

#endif
