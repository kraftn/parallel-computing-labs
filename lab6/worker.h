#ifndef WORKER_H
#define WORKER_H

#include <QObject>

// Класс для нахождения n самых встречающихся слов
class Worker : public QObject {
  Q_OBJECT

 public:
  Worker() = default;
  ~Worker() = default;

 public slots:
  // Начать подсчёт слов
  void DoWork(const QString& data, const int n);

 signals:
  // Вернуть результат
  void SendResult(const QString& result);
};

// Сортировка слов по их частотности
std::vector<std::pair<QString, int>> Sort(const std::map<QString, int>& data);
// Функция для сравнения слов при сортировке
bool Compare(std::pair<QString, int> left, std::pair<QString, int> right);

#endif
