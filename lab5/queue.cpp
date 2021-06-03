#include "queue.h"
#include "timer.h"

#include <exception>


Queue::Queue(const Queue& obj) {
  data_ = obj.data_;
}



Queue& Queue::operator=(const Queue& obj) {
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



void Queue::Push(const int value) {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  data_.push(value);
  capture_time_.push_back(timer.Finish());
}



bool Queue::IsEmpty() {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  bool res(data_.empty());
  capture_time_.push_back(timer.Finish());

  return res;
}



size_t Queue::Size() {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  size_t res(data_.size());
  capture_time_.push_back(timer.Finish());

  return res;
}



void Queue::Pop() {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  if (!data_.empty()) {
    data_.pop();
  }
  capture_time_.push_back(timer.Finish());
}



int Queue::Top() {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  int res(data_.front());
  capture_time_.push_back(timer.Finish());

  return res;
}



void Queue::Clear() {
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
