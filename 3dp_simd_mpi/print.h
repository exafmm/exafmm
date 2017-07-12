#ifndef print_h
#define print_h
#include <string>
#include <iomanip>
#include <iostream>

namespace exafmm {
  bool VERBOSE = true;                          //!< Print to screen
  static const int stringLength = 20;           //!< Length of formatted string
  static const int decimal = 7;                 //!< Decimal precision

  void print(std::string s) {
    if (VERBOSE) {
      s += " ";
      std::cout << "--- " << std::setw(stringLength) << std::left
                << std::setfill('-') << s << std::setw(decimal+1) << "-"
                << std::setfill(' ') << std::endl;
    }
  }

  template<typename T>
  void print(std::string s, T v, bool fixed=true) {
    if (VERBOSE) {
      std::cout << std::setw(stringLength) << std::left << s << " : ";
      if(fixed)
        std::cout << std::setprecision(decimal) << std::fixed;
      else
        std::cout << std::setprecision(1) << std::scientific;
      std::cout << v << std::endl;
    }
  }
}
#endif
