#ifndef vec_h
#define vec_h
#include <ostream>
namespace exafmm {
  template<int N, typename T>
  class vec {
  private:
    T data[N];
  public:
    vec(){}
    vec(const T &v) {
      for (int i=0; i<N; i++) data[i] = v;
    }
    vec(const vec &v) {
      for (int i=0; i<N; i++) data[i] = v[i];
    }
    ~vec(){}
    const vec &operator=(const T v) {
      for (int i=0; i<N; i++) data[i] = v;
      return *this;
    }
    const vec &operator+=(const T v) {
      for (int i=0; i<N; i++) data[i] += v;
      return *this;
    }
    const vec &operator-=(const T v) {
      for (int i=0; i<N; i++) data[i] -= v;
      return *this;
    }
    const vec &operator*=(const T v) {
      for (int i=0; i<N; i++) data[i] *= v;
      return *this;
    }
    const vec &operator/=(const T v) {
      for (int i=0; i<N; i++) data[i] /= v;
      return *this;
    }
    const vec &operator>=(const T v) {
      for (int i=0; i<N; i++) data[i] >= v;
      return *this;
    }
    const vec &operator<=(const T v) {
      for (int i=0; i<N; i++) data[i] <= v;
      return *this;
    }
    const vec &operator&=(const T v) {
      for (int i=0; i<N; i++) data[i] &= v;
      return *this;
    }
    const vec &operator|=(const T v) {
      for (int i=0; i<N; i++) data[i] |= v;
      return *this;
    }
    const vec &operator=(const vec & v) {
      for (int i=0; i<N; i++) data[i] = v[i];
      return *this;
    }
    const vec &operator+=(const vec & v) {
      for (int i=0; i<N; i++) data[i] += v[i];
      return *this;
    }
    const vec &operator-=(const vec & v) {
      for (int i=0; i<N; i++) data[i] -= v[i];
      return *this;
    }
    const vec &operator*=(const vec & v) {
      for (int i=0; i<N; i++) data[i] *= v[i];
      return *this;
    }
    const vec &operator/=(const vec & v) {
      for (int i=0; i<N; i++) data[i] /= v[i];
      return *this;
    }
    const vec &operator>=(const vec & v) {
      for (int i=0; i<N; i++) data[i] >= v[i];
      return *this;
    }
    const vec &operator<=(const vec & v) {
      for (int i=0; i<N; i++) data[i] <= v[i];
      return *this;
    }
    const vec &operator&=(const vec & v) {
      for (int i=0; i<N; i++) data[i] &= v[i];
      return *this;
    }
    const vec &operator|=(const vec & v) {
      for (int i=0; i<N; i++) data[i] |= v[i];
      return *this;
    }
    vec operator+(const T v) const {
      return vec(*this) += v;
    }
    vec operator-(const T v) const {
      return vec(*this) -= v;
    }
    vec operator*(const T v) const {
      return vec(*this) *= v;
    }
    vec operator/(const T v) const {
      return vec(*this) /= v;
    }
    vec operator>(const T v) const {
      return vec(*this) >= v;
    }
    vec operator<(const T v) const {
      return vec(*this) <= v;
    }
    vec operator&(const T v) const {
      return vec(*this) &= v;
    }
    vec operator|(const T v) const {
      return vec(*this) |= v;
    }
    vec operator+(const vec & v) const {
      return vec(*this) += v;
    }
    vec operator-(const vec & v) const {
      return vec(*this) -= v;
    }
    vec operator*(const vec & v) const {
      return vec(*this) *= v;
    }
    vec operator/(const vec & v) const {
      return vec(*this) /= v;
    }
    vec operator>(const vec & v) const {
      return vec(*this) >= v;
    }
    vec operator<(const vec & v) const {
      return vec(*this) <= v;
    }
    vec operator&(const vec & v) const {
      return vec(*this) &= v;
    }
    vec operator|(const vec & v) const {
      return vec(*this) |= v;
    }
    vec operator-() const {
      vec temp;
      for (int i=0; i<N; i++) temp[i] = -data[i];
      return temp;
    }
    T &operator[](int i) {
      return data[i];
    }
    const T &operator[](int i) const {
      return data[i];
    }
    operator       T* ()       {return data;}
    operator const T* () const {return data;}
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {
      for (int i=0; i<N; i++) s << v[i] << ' ';
      return s;
    }
    friend T sum(const vec & v) {
      T temp = 0;
      for (int i=0; i<N; i++) temp += v[i];
      return temp;
    }
    friend T norm(const vec & v) {
      T temp = 0;
      for (int i=0; i<N; i++) temp += v[i] * v[i];
      return temp;
    }
    friend vec min(const vec & v, const vec & w) {
      vec temp;
      for (int i=0; i<N; i++) temp[i] = v[i] < w[i] ? v[i] : w[i];
      return temp;
    }
    friend vec max(const vec & v, const vec & w) {
      vec temp;
      for (int i=0; i<N; i++) temp[i] = v[i] > w[i] ? v[i] : w[i];
      return temp;
    }
    friend T min(const vec & v) {
      T temp = v[0];
      for (int i=1; i<N; i++) temp = temp < v[i] ? temp : v[i];
      return temp;
    }
    friend T max(const vec & v) {
      T temp = v[0];
      for (int i=1; i<N; i++) temp = temp > v[i] ? temp : v[i];
      return temp;
    }
    friend vec sin(const vec & v) {
      vec temp;
      for (int i=0; i<N; i++) temp[i] = sin(v[i]);
      return temp;
    }
    friend vec cos(const vec & v) {
      vec temp;
      for (int i=0; i<N; i++) temp[i] = cos(v[i]);
      return temp;
    }
    friend void sincos(vec & s, vec & c, const vec & v) {
      for (int i=0; i<N; i++) {
        s[i] = sin(v[i]);
        c[i] = cos(v[i]);
      }
    }
    friend vec exp(const vec & v) {
      vec temp;
      for (int i=0; i<N; i++) temp[i] = exp(v[i]);
      return temp;
    }
    friend int wrap(vec & v, const vec & w) {
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
    friend void unwrap(vec & v, const vec & w, const int & iw) {
      for (int i=0; i<N; i++) {
        if((iw >> i) & 1) v[i] += (v[i] > 0 ? -w[i] : w[i]);
      }
    }
  };
}
#endif
