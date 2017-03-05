#ifndef logger_h
#define logger_h
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include "namespace.h"
#include <pthread.h>
#include <queue>
#include <stdint.h>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <vector>

#if EXAFMM_USE_PAPI
#include <cstring>
#include <papi.h>
#endif

namespace EXAFMM_NAMESPACE {
   //! Timer and Tracer logger
  namespace logger {
    typedef std::map<std::string,double> Timer;                 //!< Map of timer event name to timed value
    typedef Timer::iterator              T_iter;                //!< Iterator of timer event name map

    Timer           beginTimer;                                 //!< Timer base value
    Timer           timer;                                      //!< Timings of all events
    int stringLength = 20;                                      //!< Max length of event name
    int decimal = 7;                                            //!< Decimal precision
    bool verbose = false;                                       //!< Print to screen
    const char * path = "./";                                   //!< Path to save files

    //! Timer function
    double get_time() {
      struct timeval tv;                                        // Time value
      gettimeofday(&tv, NULL);                                  // Get time of day in seconds and microseconds
      return double(tv.tv_sec)+double(tv.tv_usec)*1e-6;         // Combine seconds and microseconds and return
    }

    //! Cycle counter
    inline uint64_t get_cycle() {
      uint32_t low = 0, high = 0;                               // Define low and high 32 bits of cycle counter
#if !__FUJITSU
      asm volatile ("rdtsc" : "=a" (low), "=d" (high));         // Call rdtsc
#endif
      return (uint64_t(high) << 32) | uint64_t(low);            // Return 64 bit cycle counter
    }

    //! Cycle counter with thread ID
    inline uint64_t get_cycle(uint32_t * id) {
      uint32_t low = 0, high = 0;                               // Define low and high 32 bits of cycle counter
      if (!id) return 0;                                        // Count only for valid thread ID
#if !__FUJITSU
      asm volatile ("rdtscp" : "=a" (low), "=d" (high), "=c" (*id));// Call rdtscp
#endif
      return (uint64_t(high) << 32) | uint64_t(low);            // Return 64 bit cycle counter
    }

    //! Print message to standard output
    inline void printTitle(std::string title) {
      if (verbose) {                                            // If verbose flag is true
	title += " ";                                           //  Append space to end of title
	std::cout << "--- " << std::setw(stringLength)          //  Align string length
		  << std::left                                  //  Left shift
		  << std::setfill('-')                          //  Set to fill with '-'
		  << title << std::setw(10) << "-"              //  Fill until end of line
		  << std::setfill(' ') << std::endl;            //  Set back to fill with ' '
      }                                                         // End if for verbose flag
    }

    //! Start timer for given event
    inline void startTimer(std::string event) {
      beginTimer[event] = get_time();                           // Get time of day and store in beginTimer
    }

    //! Print timings of a specific event
    inline void printTime(std::string event) {
      if (verbose) {                                            // If verbose flag is true
	std::cout << std::setw(stringLength) << std::left       //  Set format
		  << event << " : " << std::setprecision(decimal) << std::fixed
		  << timer[event] << " s" << std::endl;         //  Print event and timer
      }                                                         // End if for verbose flag
    }

    //! Stop timer for given event
    double stopTimer(std::string event, int print=1) {
      double endTimer = get_time();                             // Get time of day and store in endTimer
      timer[event] += endTimer - beginTimer[event];             // Accumulate event time to timer
      if (verbose && print) printTime(event);                   // Print event and timer to screen
      return endTimer - beginTimer[event];                      // Return the event time
    }

    //! Write timings of all events
    inline void writeTime(int mpirank=0) {
      std::stringstream name;                                   // File name
      name << path << "time" << std::setfill('0') << std::setw(6) // Set format
	   << mpirank << ".dat";                                // Create file name for timer
      std::ofstream timerFile(name.str().c_str());              // Open timer log file
      for (T_iter E=timer.begin(); E!=timer.end(); E++) {       // Loop over all events
	timerFile << std::setw(stringLength) << std::left       //  Set format
		  << E->first << " " << E->second << std::endl; //  Print event and timer
      }                                                         // End loop over all events
      timerFile.close();                                        // Close timer log file
      timer.clear();                                            // Clear timer
    }

    //! Erase single event in timer
    inline void resetTimer(std::string event) {
      timer.erase(event);                                       // Erase event from timer
    }

    //! Erase all events in timer
    inline void resetTimer() {
      timer.clear();                                            // Clear timer
    }
  }
}
#endif
