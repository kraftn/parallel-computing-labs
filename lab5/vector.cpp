#include "vector.h"
#include "timer.h"

#include <exception>


Vector::Vector(const size_t size) {
  data_.resize(size);
}



Vector::Vector(const Vector& obj) {
  data_ = obj.data_;
}



Vector& Vector::operator=(const Vector& obj) {
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



size_t Vector::Size() {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  size_t res(data_.size());
  capture_time_.push_back(timer.Finish());

  return res;
}



size_t Vector::Capacity() {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  size_t res(data_.capacity());
  capture_time_.push_back(timer.Finish());

  return res;
}



int& Vector::operator[](const ptrdiff_t number) {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  int& res(data_.at(number));
  capture_time_.push_back(timer.Finish());

  return res;
}



void Vector::Resize(size_t size) {
  Timer timer;
  timer.Start();
  std::unique_lock<std::mutex> lock(mutex_);
  waiting_time_.push_back(timer.Finish());

  timer.Start();
  data_.resize(size);
  capture_time_.push_back(timer.Finish());
}
