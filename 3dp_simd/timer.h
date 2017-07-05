#ifndef timer_h
#define timer_h
#include <map>
#include <sys/time.h>

namespace exafmm {
  static timeval t;
  static std::map<std::string,timeval> timer;

  void start(std::string event) {
    gettimeofday(&t, NULL);
    timer[event] = t;
  }

  double stop(std::string event) {
    gettimeofday(&t, NULL);
    double eventTime = t.tv_sec - timer[event].tv_sec +
      (t.tv_usec - timer[event].tv_usec) * 1e-6;
    printf("%-20s : %f s\n", event.c_str(), eventTime);
    return eventTime;
  }
}
#endif
