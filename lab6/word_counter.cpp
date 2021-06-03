#include <QFileDialog>
#include <QStringList>
#include <QMessageBox>

#include "word_counter.h"
#include "worker.h"

WordCounter::WordCounter(QWidget* parent)
  : QMainWindow(parent) {
  ui.setupUi(this);

  // Связать нажатия на кнопки с обработчиками
  connect(ui.pushButton_open, SIGNAL(clicked()), this, SLOT(Open()));
  connect(ui.pushButton_search, SIGNAL(clicked()), this, SLOT(Search()));
  connect(ui.pushButton_select_mystem, SIGNAL(clicked()), this, SLOT(FindMystem()));

  connect(&mystem_, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(GetDataFromMyStem(int, QProcess::ExitStatus)));
}



WordCounter::~WordCounter() {
  thread_.quit();
  thread_.wait();
}



void WordCounter::Open() {
  file_path_ = QFileDialog::getOpenFileName(this, "Открыть файл", QString(), "*.txt");
  ui.label_file_path->setText(file_path_);
}



void WordCounter::Search() {
  // Если указан путь до файла и программы
  if (!file_path_.isNull() && !mystem_path_.isNull()) {
    // Если выбрано количество слов
    if (ui.spinBox_number_words->value() != 0) {
      ui.pushButton_search->setEnabled(false);
      ui.spinBox_number_words->setEnabled(false);
      ui.statusbar->showMessage("Обработка документа...");
      ui.textBrowser_output->clear();

      // Запустить новый процесс
      QStringList arguments;
      arguments << "-ldn" << file_path_;
      mystem_.start(mystem_path_, arguments);
    } else {
      QMessageBox::information(this, "Ошибка", "Выберите количество слов.");
    }
  }
  else {
    QMessageBox::information(this, "Ошибка", "Выберите файл и укажите местоположение mystem.");
  }
}



void WordCounter::GetDataFromMyStem(int exitCode, QProcess::ExitStatus exitStatus) {
  // Если mystem корректно завершил работу, запустить подсчёт слов и сортировку
  if (QProcess::NormalExit == exitStatus) {
    Worker* worker = new Worker;

    worker->moveToThread(&thread_);
    connect(this, &WordCounter::LaunchWorker, worker, &Worker::DoWork);
    connect(worker, &Worker::SendResult, this, &WordCounter::GetResult);
    connect(&thread_, &QThread::finished, worker, &QObject::deleteLater);

    thread_.start();
    emit LaunchWorker(QString(mystem_.readAllStandardOutput()), ui.spinBox_number_words->value());
  }
  else {
    QMessageBox::information(this, "Ошибка", "Во время работы mystem произошла ошибка.");
  }
}



void WordCounter::GetResult(const QString& result) {
  // Вывести результат на экран
  ui.textBrowser_output->setText(result);

  ui.pushButton_search->setEnabled(true);
  ui.spinBox_number_words->setEnabled(true);

  ui.statusbar->clearMessage();
}



void WordCounter::FindMystem() {
  mystem_path_ = QFileDialog::getOpenFileName(this, "Найти mystem", QString(), "*.exe");
  ui.label_mystem_path->setText(mystem_path_);
}
