#include <cstddef>

#include "worker.h"

void Worker::DoWork(const QString& data, const int n) {
  // Разделить строку по символу перевода строки
  QStringList words(data.split('\n'));
  // Удалить перевод строки в конце вывода
  words.removeLast();
  for (ptrdiff_t i_word(0); i_word < words.size(); i_word += 1) {
    words[i_word].replace(QString("?"), QString(""));
  }

  // Подсчёт количества слов
  std::map<QString, int> map;
  for (ptrdiff_t i_word(0); i_word < words.size(); i_word += 1) {
    map[words[i_word]] += 1;
  }

  // Отсортировать слова по их частотности, начания с самых употребляемых
  std::vector<std::pair<QString, int>> sorted(Sort(map));

  // Записать результата в строку
  int i_top(0);
  QString result;
  while ((i_top < n) && (i_top < sorted.size())) {
    result += QString::number(sorted[i_top].second) + QString(": ") + sorted[i_top].first + '\n';
    i_top += 1;
  }

  // Отправить результат
  emit SendResult(result);
}



std::vector<std::pair<QString, int>> Sort(const std::map<QString, int>& data) {
  using namespace std;

  // Сформировать вектор пар
  vector<pair<QString, int>> res;
  for (auto i = data.begin(); i != data.end(); i++) {
    res.push_back(make_pair(i->first, i->second));
  }
  // Отсортировать вектор
  sort(res.begin(), res.end(), Compare);

  return res;
}



bool Compare(std::pair<QString, int> left, std::pair<QString, int> right) {
  // Компаратор для сортировки по убыванию
  return left.second > right.second;
}
