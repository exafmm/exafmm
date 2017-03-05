#ifndef vec_h
#define vec_h
#include "namespace.h"
#include <ostream>

namespace EXAFMM_NAMESPACE {
  template<int N, typename T>
  class vec {
  private:
    T data[N];
  public:
    vec(){}                                                     // Default constructor
    vec(const T &v) {                                           // Copy constructor (scalar)
      for (int i=0; i<N; i++) data[i] = v;
    }
    vec(const vec &v) {                                         // Copy constructor (vector)
      for (int i=0; i<N; i++) data[i] = v[i];
    }
    ~vec(){}                                                    // Destructor
    const vec &operator=(const T v) {                           // Scalar assignment
      for (int i=0; i<N; i++) data[i] = v;
      return *this;
    }
    const vec &operator+=(const T v) {                          // Scalar compound assignment (add)
      for (int i=0; i<N; i++) data[i] += v;
      return *this;
    }
    const vec &operator-=(const T v) {                          // Scalar compound assignment (subtract)
      for (int i=0; i<N; i++) data[i] -= v;
      return *this;
    }
    const vec &operator*=(const T v) {                          // Scalar compound assignment (multiply)
      for (int i=0; i<N; i++) data[i] *= v;
      return *this;
    }
    const vec &operator/=(const T v) {                          // Scalar compound assignment (divide)
      for (int i=0; i<N; i++) data[i] /= v;
      return *this;
    }
    const vec &operator>=(const T v) {                          // Scalar compound assignment (greater than)
      for (int i=0; i<N; i++) data[i] >= v;
      return *this;
    }
    const vec &operator<=(const T v) {                          // Scalar compound assignment (less than)
      for (int i=0; i<N; i++) data[i] <= v;
      return *this;
    }
    const vec &operator&=(const T v) {                          // Scalar compound assignment (bitwise and)
      for (int i=0; i<N; i++) data[i] &= v;
      return *this;
    }
    const vec &operator|=(const T v) {                          // Scalar compound assignment (bitwise or)
      for (int i=0; i<N; i++) data[i] |= v;
      return *this;
    }
    const vec &operator=(const vec & v) {                       // Vector assignment
      for (int i=0; i<N; i++) data[i] = v[i];
      return *this;
    }
    const vec &operator+=(const vec & v) {                      // Vector compound assignment (add)
      for (int i=0; i<N; i++) data[i] += v[i];
      return *this;
    }
    const vec &operator-=(const vec & v) {                      // Vector compound assignment (subtract)
      for (int i=0; i<N; i++) data[i] -= v[i];
      return *this;
    }
    const vec &operator*=(const vec & v) {                      // Vector compound assignment (multiply)
      for (int i=0; i<N; i++) data[i] *= v[i];
      return *this;
    }
    const vec &operator/=(const vec & v) {                      // Vector compound assignment (divide)
      for (int i=0; i<N; i++) data[i] /= v[i];
      return *this;
    }
    const vec &operator>=(const vec & v) {                      // Vector compound assignment (greater than)
      for (int i=0; i<N; i++) data[i] >= v[i];
      return *this;
    }
    const vec &operator<=(const vec & v) {                      // Vector compound assignment (less than)
      for (int i=0; i<N; i++) data[i] <= v[i];
      return *this;
    }
    const vec &operator&=(const vec & v) {                      // Vector compound assignment (bitwise and)
      for (int i=0; i<N; i++) data[i] &= v[i];
      return *this;
    }
    const vec &operator|=(const vec & v) {                      // Vector compound assignment (bitwise or)
      for (int i=0; i<N; i++) data[i] |= v[i];
      return *this;
    }
    vec operator+(const T v) const {                            // Scalar arithmetic (add)
      return vec(*this) += v;
    }
    vec operator-(const T v) const {                            // Scalar arithmetic (subtract)
      return vec(*this) -= v;
    }
    vec operator*(const T v) const {                            // Scalar arithmetic (multiply)
      return vec(*this) *= v;
    }
    vec operator/(const T v) const {                            // Scalar arithmetic (divide)
      return vec(*this) /= v;
    }
    vec operator>(const T v) const {                            // Scalar arithmetic (greater than)
      return vec(*this) >= v;
    }
    vec operator<(const T v) const {                            // Scalar arithmetic (less than)
      return vec(*this) <= v;
    }
    vec operator&(const T v) const {                            // Scalar arithmetic (bitwise and)
      return vec(*this) &= v;
    }
    vec operator|(const T v) const {                            // Scalar arithmetic (bitwise or)
      return vec(*this) |= v;
    }
    vec operator+(const vec & v) const {                        // Vector arithmetic (add)
      return vec(*this) += v;
    }
    vec operator-(const vec & v) const {                        // Vector arithmetic (subtract)
      return vec(*this) -= v;
    }
    vec operator*(const vec & v) const {                        // Vector arithmetic (multiply)
      return vec(*this) *= v;
    }
    vec operator/(const vec & v) const {                        // Vector arithmetic (divide)
      return vec(*this) /= v;
    }
    vec operator>(const vec & v) const {                        // Vector arithmetic (greater than)
      return vec(*this) >= v;
    }
    vec operator<(const vec & v) const {                        // Vector arithmetic (less than)
      return vec(*this) <= v;
    }
    vec operator&(const vec & v) const {                        // Vector arithmetic (bitwise and)
      return vec(*this) &= v;
    }
    vec operator|(const vec & v) const {                        // Vector arithmetic (bitwise or)
      return vec(*this) |= v;
    }
    vec operator-() const {                                     // Vector arithmetic (negation)
      vec temp;
      for (int i=0; i<N; i++) temp[i] = -data[i];
      return temp;
    }
    T &operator[](int i) {                                      // Indexing (lvalue)
      return data[i];
    }
    const T &operator[](int i) const {                          // Indexing (rvalue)
      return data[i];
    }
    operator       T* ()       {return data;}                   // Type-casting (lvalue)
    operator const T* () const {return data;}                   // Type-casting (rvalue)
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {// Component-wise output stream
      for (int i=0; i<N; i++) s << v[i] << ' ';
      return s;
    }
    friend T sum(const vec & v) {                               // Sum vector
      T temp = 0;
      for (int i=0; i<N; i++) temp += v[i];
      return temp;
    }
    friend T norm(const vec & v) {                              // L2 norm squared
      T temp = 0;
      for (int i=0; i<N; i++) temp += v[i] * v[i];
      return temp;
    }
    friend vec min(const vec & v, const vec & w) {              // Element-wise minimum
      vec temp;
      for (int i=0; i<N; i++) temp[i] = v[i] < w[i] ? v[i] : w[i];
      return temp;
    }
    friend vec max(const vec & v, const vec & w) {              // Element-wise maximum
      vec temp;
      for (int i=0; i<N; i++) temp[i] = v[i] > w[i] ? v[i] : w[i];
      return temp;
    }
    friend T min(const vec & v) {                               // Reduce minimum
      T temp = v[0];
      for (int i=1; i<N; i++) temp = temp < v[i] ? temp : v[i];
      return temp;
    }
    friend T max(const vec & v) {                               // Reduce maximum
      T temp = v[0];
      for (int i=1; i<N; i++) temp = temp > v[i] ? temp : v[i];
      return temp;
    }
    friend vec sin(const vec & v) {                             // Sine function
      vec temp;
      for (int i=0; i<N; i++) temp[i] = sin(v[i]);
      return temp;
    }
    friend vec cos(const vec & v) {                             // Cosine function
      vec temp;
      for (int i=0; i<N; i++) temp[i] = cos(v[i]);
      return temp;
    }
    friend void sincos(vec & s, vec & c, const vec & v) {       // Sine & cosine function
      for (int i=0; i<N; i++) {
	s[i] = sin(v[i]);
	c[i] = cos(v[i]);
      }
    }
    friend vec exp(const vec & v) {                             // Exponential function
      vec temp;
      for (int i=0; i<N; i++) temp[i] = exp(v[i]);
      return temp;
    }
    friend int wrap(vec & v, const vec & w) {                     // Wrap around periodic boundary
      int iw = 0;
      for (int i=0; i<N; i++) {
	if(v[i] < -w[i] / 2) {
	  v[i] += w[i];
	  iw |= 1 << i;
	}
	if(v[i] >  w[i] / 2) {
	  v[i] -= w[i];
	  iw |= 1 << i;
	}
      }
      return iw;
    }
    friend void unwrap(vec & v, const vec & w, const int & iw) {  // Undo wrap around periodic boundary
      for (int i=0; i<N; i++) {
	if((iw >> i) & 1) v[i] += (v[i] > 0 ? -w[i] : w[i]);
      }
    }
  };
}
#endif
