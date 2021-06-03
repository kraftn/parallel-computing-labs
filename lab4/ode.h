#ifndef ODE_H
#define ODE_H

#include <vector>
#include <functional>

// Реализация метода Пикара для решения системы ОДУ
// на отрезке [t0; tm] с функциями правой части f
// и начальными условиями alpha
std::vector<std::vector<double>> Solve(
    const std::vector<std::function<double(const double, const std::vector<double>&)>>& f,
    const std::vector<double>& alpha,
    const double t0,
    const double tm,
    const double step,
    const bool is_parallel);

// Евклидова метрика
double EuclideanMetric(const std::vector<double>& x, const std::vector<double>& y, const bool is_parallel);

// Функция для получения значений функции u на отрезке [t0; tm]
std::vector<double> Calculate(
    const std::function<double(const double)>& u,
    const double t0,
    const double tm,
    const double step);

// Реализация операции транспонирования
std::vector<std::vector<double>> Transpose(const std::vector<std::vector<double>>& matrix);

// Реализация метрики, основанной на макимум-норме
double MaximumMetric(const std::vector<double>& x, const std::vector<double>& y);

#endif
