#ifndef STACK_H
#define STACK_H

#include <cstddef>
#include <stack>

#include "thread_safe_object.h"

// Класс потокобезопасного стека
class Stack : public ThreadSafeObject {
public:
  Stack() = default;
  Stack(const Stack& obj);
  ~Stack() = default;
  Stack& operator=(const Stack& obj);
  // Поместить новый элемент
  void Push(const int value);
  // Проверить на пустоту
  bool IsEmpty();
  // Получить количество элементов
  size_t Size();
  // Удалить элемент
  void Pop();
  // Получить элемент
  int Top();
  // Удалить все элементы
  void Clear();

private:
  std::stack<int> data_;
};

#endif
