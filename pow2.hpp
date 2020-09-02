#include "common.hpp"



#ifndef POW2_H
#define POW2_H

// Represents the cardinality of a pow2 sized set. Allows div/mod arithmetic operations on indexes.
template <typename T>
struct Pow2 {
  Pow2(uint_fast8_t bits)
      : _bits(bits) {
    assume(bits < 8 * sizeof(T), "Pow2(%u > %u)", unsigned(bits), unsigned(8 * sizeof(T)));
  }

  uint_fast8_t bits() const { return _bits; }
  T value() const { return T(1) << _bits; }
  explicit operator T() const { return value(); }
  T max() const { return value() - T(1); }
  Pow2& operator=(const Pow2&) = default;

  Pow2()
      : _bits(0) {}
  Pow2(const Pow2&) = default;

  friend T operator*(const T& x, const Pow2& y) { return x << y._bits; }
  friend T& operator*=(T& x, const Pow2& y) { return x <<= y._bits; }
  friend T operator/(const T& x, const Pow2& y) { return x >> y._bits; }
  friend T& operator/=(T& x, const Pow2& y) { return x >>= y._bits; }
  friend T operator%(const T& x, const Pow2& y) { return x & y.max(); }
  friend T& operator%=(T& x, const Pow2& y) { return x &= y.max(); }
  Pow2& operator>>=(uint_fast8_t d) {
    _bits -= d;
    return *this;
  }
  Pow2& operator<<=(uint_fast8_t d) {
    _bits += d;
    return *this;
  }
  friend bool operator<(const T& x, const Pow2& y) { return x < y.value(); }
  friend bool operator<=(const T& x, const Pow2& y) { return x < y.value(); }
  friend T operator+(const T& x, const Pow2& y) { return x + y.value(); }
  friend T& operator+=(T& x, const Pow2& y) { return x += y.value(); }
  friend T operator-(const T& x, const Pow2& y) { return x - y.value(); }
  friend T& operator-=(T& x, const Pow2& y) { return x -= y.value(); }

  private:
  uint_fast8_t _bits;
};

#endif
