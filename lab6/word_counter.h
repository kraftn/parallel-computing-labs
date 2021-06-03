#ifndef WORD_COUNTER_H
#define WORD_COUNTER_H

#include <map>
#include <vector>

#include <QtWidgets/QMainWindow>
#include <QProcess>
#include <QString>
#include <QThread>

#include "ui_word_counter.h"

// Класс основного окна приложения
class WordCounter : public QMainWindow {
  Q_OBJECT

 public:
  WordCounter(QWidget* parent = 0);
  ~WordCounter();

 public slots:
  // Найти файл для анализа
  void Open();
  // Начать подсчёт слов
  void Search();
  // Получить данные из утилиты mystem
  void GetDataFromMyStem(int exitCode, QProcess::ExitStatus exitStatus);
  // Получить результат подсчёта слов
  void GetResult(const QString& result);
  // Указать местоположение утилиты mystem
  void FindMystem();

signals:
  // Начать подсчёт количества слов
  void LaunchWorker(const QString& process, const int n);

 private:
  Ui::MainWindow ui;
  // Путь до файла
  QString file_path_;
  // Путь до утилиты mystem
  QString mystem_path_;
  // Процесс для вызова программы mystem
  QProcess mystem_;
  // Поток для подсчёта количества слов и их сортировки
  QThread thread_;
};

#endif
