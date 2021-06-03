#include <algorithm>

#include "thread_safe_object.h"


double ThreadSafeObject::GetAverageCaptureTime() {
  return GetCaptureTime() / capture_time_.size();
}



double ThreadSafeObject::GetAverageWaitingTime() {
  std::unique_lock<std::mutex> lock(mutex_);
  double sum(0.0);
  for (ptrdiff_t i(0); i < waiting_time_.size(); i += 1) {
    sum += waiting_time_[i];
  }
  return sum / waiting_time_.size();
}



double ThreadSafeObject::GetMaxCaptureTime() {
  std::unique_lock<std::mutex> lock(mutex_);
  double res(*std::max_element(capture_time_.begin(), capture_time_.end()));
  return res;
}



double ThreadSafeObject::GetMaxWaitingTime() {
  std::unique_lock<std::mutex> lock(mutex_);
  double res(*std::max_element(waiting_time_.begin(), waiting_time_.end()));
  return res;
}



double ThreadSafeObject::GetMinCaptureTime() {
  std::unique_lock<std::mutex> lock(mutex_);
  double res(*std::min_element(capture_time_.begin(), capture_time_.end()));
  return res;
}



double ThreadSafeObject::GetMinWaitingTime() {
  std::unique_lock<std::mutex> lock(mutex_);
  double res(*std::min_element(waiting_time_.begin(), waiting_time_.end()));
  return res;
}



double ThreadSafeObject::GetCaptureTime() {
  std::unique_lock<std::mutex> lock(mutex_);
  double sum(0.0);
  for (ptrdiff_t i(0); i < capture_time_.size(); i += 1) {
    sum += capture_time_[i];
  }
  return sum;
}

