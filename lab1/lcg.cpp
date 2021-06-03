#include "lcg.h"


LCG::LCG()
  : current_x_(seed_) {
}



LCG::LCG(const int64_t initial_value)
  : current_x_(initial_value) {
}



LCG::~LCG() {
}



double LCG::Random() {
  double double_x(current_x_ / static_cast<double>(m_));
  // В данном генераторе с = 0
  current_x_ = (a_ * current_x_) % m_;
  return double_x;
}



void LCG::Reset() {
  current_x_ = seed_;
}



int64_t LCG::PowMod(const int64_t x, const int64_t y, const int64_t z) const {
  int64_t s(1), c(x), v(y);
  while (v != 0) {
    if (v % 2 == 1) {
      s = (s * c) % z;
    }

    c = (c * c) % z;
    v /= 2;
  }
  return s;
}



int64_t LCG::TransitionProcedure(const int64_t k) const {
  return (current_x_ * PowMod(a_, k, m_)) % m_;
}



bool LCG::CheckTheoremConditions() const {
  // Генератор имеет максимальный период lambda(m), если выполнены
  // следующие условия
  return (GCD(seed_, m_) == 1) && PrimitiveElement(a_, m_);
}



int64_t LCG::GCD(const int64_t a, const int64_t b) const {
  int64_t x(a);
  int64_t y(b);
  if (x > y) {
    x = b;
    y = a;
  }

  while (x != 0) {
    int64_t z(x);
    x = y % x;
    y = z;
  }

  return y;
}



bool LCG::PrimitiveElement(const int64_t a, const int64_t m) const {
  std::map<int64_t, int64_t> factorization(Factorization(m));
  if (factorization.size() > 1) {
    throw std::logic_error("Wrong m in PrimitiveElement");
  }
  int64_t p(factorization.begin()->first);
  int64_t e(factorization.begin()->second);

  // Дональд Э. Кнут "Искусство программирования", том 2 "Получисленные алгоритмы", с.40
  if (p == 2) {
    if ((e == 1) && (a % 2 == 1)) {
      return true;
    }
    if ((e == 2) && (a % 4 == 3)) {
      return true;
    }
    int32_t a_mod_8(a % 8);
    if ((e == 3) && ((a_mod_8 == 3) || (a_mod_8 == 5) || (a_mod_8 == 7))) {
      return true;
    }
    if ((e >= 4) && ((a_mod_8 == 3) || (a_mod_8 == 5))) {
      return true;
    }
  } else {
    factorization = Factorization(p - 1);
    bool condition(true);
    std::map<int64_t, int64_t>::iterator i_divisor(factorization.begin());
    while (condition && (i_divisor != factorization.end())) {
      int64_t q(i_divisor->first);
      condition = PowMod(a, (p - 1) / q, p) != 1;
      ++i_divisor;
    }
    if ((p % 2 == 1) && (e == 1) && (a % p != 0) && condition) {
      return true;
    }
    if ((p % 2 == 1) && (e > 1) && (a % p != 0) && condition && (PowMod(a, p - 1, p * p) != 1)) {
      return true;
    }
  }

  return false;
}



std::map<int64_t, int64_t> LCG::Factorization(const int64_t n) const {
  std::map<int64_t, int64_t> res;

  int64_t curr_n(n), i(2);
  while (i <= sqrt(curr_n)) {

    while ((curr_n % i) == 0) {
      curr_n /= i;
      res[i] += 1;
    }

    i += 1;
  }

  if (curr_n != 1) {
    res[curr_n] += 1;
  }

  return res;
}



int64_t LCG::GetPeriod() const {
  std::map<int64_t, int64_t> factorization(Factorization(m_));
  if (factorization.size() > 1) {
    throw std::logic_error("Wrong m_");
  }
  int64_t p(factorization.begin()->first);
  int64_t e(factorization.begin()->second);

  // Дональд Э. Кнут "Искусство программирования", том 2 "Получисленные алгоритмы", с.40
  if (p == 2) {
    if (e == 1) {
      return 1;
    }
    if (e == 2) {
      return 2;
    }

    return static_cast<int64_t>(pow(2, e - 2));
  }
  return static_cast<int64_t>(pow(p, e - 1)) * (p - 1);
}
