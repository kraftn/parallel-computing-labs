#include "stack.h"
#include "timer.h"

#include <exception>


Stack::Stack(const Stack& obj) {
  data_ = obj.data_;
}



Stack& Stack::operator=(const Stack& obj) {
  if (this != &obj) {
    Timer timer;
    timer.Start();
    std::unique_lock<std::mutex> lock(mutex_);
    waiting_time_.push_back(timer.Finish());

    timer.Start();
    data_ = obj.data_;
    capture_time_.push_back(timer.Finish());
  }

  return *this;
}



void Stack::Push(const int value) {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  data_.push(value);
  capture_time_.push_back(timer.Finish());
}



bool Stack::IsEmpty() {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  bool res(data_.empty());
  capture_time_.push_back(timer.Finish());

  return res;
}



size_t Stack::Size() {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  size_t res(data_.size());
  capture_time_.push_back(timer.Finish());

  return res;
}



void Stack::Pop() {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  data_.pop();
  capture_time_.push_back(timer.Finish());
}



int Stack::Top() {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  int res(data_.top());
  capture_time_.push_back(timer.Finish());

  return res;
}



void Stack::Clear() {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  while (!data_.empty()) {
    data_.pop();
  }
  capture_time_.push_back(timer.Finish());
}
