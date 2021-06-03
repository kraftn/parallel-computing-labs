#include <thread>
#include <iostream>
#include <atomic>

#include "stack.h"
#include "queue.h"
#include "vector.h"
#include "timer.h"


const size_t n_threads(100);
Stack thread_safe_stack;
Queue thread_safe_queue;
Vector thread_safe_vector(1);
Timer timer;
double total_time(0.0);
std::atomic<int> counter(0);


void stack_func() {
  thread_safe_stack.Push(1);

  int value(thread_safe_stack.Top());
  bool is_empty(thread_safe_stack.IsEmpty());
  size_t size(thread_safe_stack.Size());
  thread_safe_stack.Pop();

  thread_safe_stack.Push(1);

  if (thread_safe_stack.Size() == n_threads) {
    total_time = timer.Finish();
  }
}



void queue_func() {
  thread_safe_queue.Push(1);

  int value(thread_safe_queue.Top());
  bool is_empty(thread_safe_queue.IsEmpty());
  size_t size(thread_safe_queue.Size());
  thread_safe_queue.Pop();

  thread_safe_queue.Push(1);

  if (thread_safe_queue.Size() == n_threads) {
    total_time = timer.Finish();
  }
}



void vector_func() {
  int value(thread_safe_vector[counter]);
  counter += 1;
  thread_safe_vector.Size();
  thread_safe_vector.Capacity();
  thread_safe_vector.Resize(counter + 1);
  thread_safe_vector[counter] = value + 1;

  if (thread_safe_vector.Size() == (n_threads + 1)) {
    total_time = timer.Finish();
  }
}



int main() {
  using namespace std;
  vector<thread> threads;

  // Тестирование стека
  timer.Start();
  for (ptrdiff_t i(0); i < n_threads; i += 1) {
    threads.push_back(thread(stack_func));
  }
  
  for (ptrdiff_t i(0); i < n_threads; i += 1) {
    threads[i].join();
  }

  cout << "Stack" << endl << endl;
  cout << "Coefficient: " << thread_safe_stack.GetCaptureTime() / total_time << endl << endl;

  cout << "Capture time" << endl;
  cout << "average: " << thread_safe_stack.GetAverageCaptureTime() << endl;
  cout << "min: " << thread_safe_stack.GetMinCaptureTime() << endl;
  cout << "max: " << thread_safe_stack.GetMaxCaptureTime() << endl << endl;

  cout << "Waiting time" << endl;
  cout << "average: " << thread_safe_stack.GetAverageWaitingTime() << endl;
  cout << "min: " << thread_safe_stack.GetMinWaitingTime() << endl;
  cout << "max: " << thread_safe_stack.GetMaxWaitingTime() << endl << endl << endl;


  // Тестирование очереди
  threads.clear();
  timer.Start();
  for (ptrdiff_t i(0); i < n_threads; i += 1) {
    threads.push_back(thread(queue_func));
  }

  for (ptrdiff_t i(0); i < n_threads; i += 1) {
    threads[i].join();
  }

  cout << "Queue" << endl << endl;
  cout << "Coefficient: " << thread_safe_queue.GetCaptureTime() / total_time << endl << endl;

  cout << "Capture time" << endl;
  cout << "average: " << thread_safe_queue.GetAverageCaptureTime() << endl;
  cout << "min: " << thread_safe_queue.GetMinCaptureTime() << endl;
  cout << "max: " << thread_safe_queue.GetMaxCaptureTime() << endl << endl;

  cout << "Waiting time" << endl;
  cout << "average: " << thread_safe_queue.GetAverageWaitingTime() << endl;
  cout << "min: " << thread_safe_queue.GetMinWaitingTime() << endl;
  cout << "max: " << thread_safe_queue.GetMaxWaitingTime() << endl << endl << endl;


  // Тестирование вектора
  threads.clear();
  timer.Start();
  for (ptrdiff_t i(0); i < n_threads; i += 1) {
    threads.push_back(thread(vector_func));
  }

  for (ptrdiff_t i(0); i < n_threads; i += 1) {
    threads[i].join();
  }

  cout << "Vector" << endl << endl;
  cout << "Coefficient: " << thread_safe_vector.GetCaptureTime() / total_time << endl << endl;

  cout << "Capture time" << endl;
  cout << "average: " << thread_safe_vector.GetAverageCaptureTime() << endl;
  cout << "min: " << thread_safe_vector.GetMinCaptureTime() << endl;
  cout << "max: " << thread_safe_vector.GetMaxCaptureTime() << endl << endl;

  cout << "Waiting time" << endl;
  cout << "average: " << thread_safe_vector.GetAverageWaitingTime() << endl;
  cout << "min: " << thread_safe_vector.GetMinWaitingTime() << endl;
  cout << "max: " << thread_safe_vector.GetMaxWaitingTime() << endl;

  char stop(' ');
  std::cin >> stop;
  return 0;
}
